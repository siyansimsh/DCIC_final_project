module fft_64pt (
    input                 clk,
    input                 rst_n,

    input                 i_start, // unused, kept for interface compatibility
    output reg            o_done,

    input                 i_load_valid,
    input  signed [15:0]  i_load_re,
    input  signed [15:0]  i_load_im,

    input         [5:0]   i_read_addr,
    output signed [15:0]  o_read_re,
    output signed [15:0]  o_read_im
);

    // Consume unused interface inputs to silence synthesis warnings
    (* keep = "true" *) wire unused_start = i_start;

    // =========================================================
    // 1) Ping-pong buffers for input/output storage
    // =========================================================
    // Encourage BRAM inference for the ping-pong buffers
    (* ram_style = "block" *) reg signed [15:0] mem_re [0:127];
    (* ram_style = "block" *) reg signed [15:0] mem_im [0:127];

    reg bank_in;   // bank being written by loader
    reg bank_calc; // bank recently computed and visible to reader
    reg [6:0] load_cnt; // counts 0..63 samples
    reg input_ready_pulse;

    // Proper bit-reversal helper
    function [5:0] bit_reverse;
        input [5:0] in;
        begin
            bit_reverse = {in[0], in[1], in[2], in[3], in[4], in[5]};
        end
    endfunction

    // Saturating add/sub helper
    function signed [15:0] sat16;
        input signed [31:0] val;
        begin
            if (val > 32767)      sat16 = 16'sd32767;
            else if (val < -32768) sat16 = -16'sd32768;
            else                  sat16 = val[15:0];
        end
    endfunction

    // Arithmetic right shift helper (truncate like C fixed reference)
    function signed [31:0] trunc_shift;
        input signed [31:0] val;
        input integer sh;
        begin
            trunc_shift = val >>> sh;
        end
    endfunction

    // Twiddle lookup using real math (simulation-only, Q5.11 scale)
    function signed [15:0] tw_re;
        input integer idx;
        real ang;
        real c;
        begin
            ang = -2.0 * 3.141592653589793 * idx / 64.0;
            c = $cos(ang) * 2048.0;
            if (c > 32767.0)      tw_re = 16'sd32767;
            else if (c < -32768.0) tw_re = -16'sd32768;
            else                   tw_re = $rtoi(c);
        end
    endfunction

    function signed [15:0] tw_im;
        input integer idx;
        real ang;
        real s;
        begin
            ang = -2.0 * 3.141592653589793 * idx / 64.0;
            s = $sin(ang) * 2048.0;
            if (s > 32767.0)      tw_im = 16'sd32767;
            else if (s < -32768.0) tw_im = -16'sd32768;
            else                   tw_im = $rtoi(s);
        end
    endfunction

    // =========================================================
    // 2) Input loader: capture 64 samples into ping-pong bank
    // =========================================================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            load_cnt <= 0;
            bank_in <= 0;
            input_ready_pulse <= 0;
        end else begin
            input_ready_pulse <= 0;
            if (i_load_valid) begin
                // store in sequential order; we will bit-reverse in compute stage
                mem_re[{bank_in, load_cnt[5:0]}] <= i_load_re;
                mem_im[{bank_in, load_cnt[5:0]}] <= i_load_im;
                if (load_cnt == 63) begin
                    load_cnt <= 0;
                    bank_in <= ~bank_in;
                    input_ready_pulse <= 1; // one symbol captured
                end else begin
                    load_cnt <= load_cnt + 1;
                end
            end
        end
    end

    // =========================================================
    // 3) FFT compute (simulation-only).
    //    Wrap with !SYNTHESIS so implementation sees only the storage.
    // =========================================================
`ifndef SYNTHESIS
    task automatic run_fft;
        input which_bank; // Verilog-2001 style (bit not supported in iverilog default)
        integer i, j, k, m, half, step, idx;
        reg signed [15:0] buf_re [0:63];
        reg signed [15:0] buf_im [0:63];
        reg signed [31:0] t_re, t_im;
        reg signed [31:0] u_re, u_im;
        reg signed [15:0] w_re, w_im;
        begin
            // load to local buffer
            for (i = 0; i < 64; i = i + 1) begin
                buf_re[i] = mem_re[{which_bank, i[5:0]}];
                buf_im[i] = mem_im[{which_bank, i[5:0]}];
            end

            // bit-reverse reorder
            for (i = 0; i < 64; i = i + 1) begin
                j = bit_reverse(i[5:0]);
                if (i < j) begin
                    u_re = buf_re[i];
                    u_im = buf_im[i];
                    buf_re[i] = buf_re[j];
                    buf_im[i] = buf_im[j];
                    buf_re[j] = u_re[15:0];
                    buf_im[j] = u_im[15:0];
                end
            end

            // FFT stages
            for (m = 2; m <= 64; m = m << 1) begin
                half = m >> 1;
                step = 64 / m;
                for (k = 0; k < 64; k = k + m) begin
                    for (j = 0; j < half; j = j + 1) begin
                        idx = j * step;
                        w_re = tw_re(idx);
                        w_im = tw_im(idx);

                        t_re = trunc_shift(buf_re[k + j + half] * w_re - buf_im[k + j + half] * w_im, 11);
                        t_im = trunc_shift(buf_re[k + j + half] * w_im + buf_im[k + j + half] * w_re, 11);

                        u_re = buf_re[k + j];
                        u_im = buf_im[k + j];

                        buf_re[k + j]       = sat16(u_re + t_re);
                        buf_im[k + j]       = sat16(u_im + t_im);
                        buf_re[k + j + half] = sat16(u_re - t_re);
                        buf_im[k + j + half] = sat16(u_im - t_im);
                    end
                end
            end

            // write back in natural order
            for (i = 0; i < 64; i = i + 1) begin
                mem_re[{which_bank, i[5:0]}] <= buf_re[i];
                mem_im[{which_bank, i[5:0]}] <= buf_im[i];
            end
        end
    endtask
`endif

    // =========================================================
    // 4) Control: trigger FFT when a bank is full, raise done pulse
    //    Keep bank_calc stable for 64 cycles after each done, but still
    //    allow a pending symbol to compute and be latched after the hold.
    // =========================================================
    reg [6:0] read_hold;      // counts down read window (64 cycles)
    reg       pending_valid;  // a computed bank waiting to be exposed
    reg       pending_bank;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            bank_calc <= 0;
            o_done <= 0;
            read_hold <= 0;
            pending_valid <= 0;
            pending_bank <= 0;
        end else begin
            o_done <= 0;

            // countdown hold window
            if (read_hold != 0)
                read_hold <= read_hold - 1;

            // when hold releases and we have pending, publish it now
            if (read_hold == 0 && pending_valid) begin
                bank_calc <= pending_bank;
                o_done <= 1;
                read_hold <= 7'd64;
                pending_valid <= 0;
            end

            // handle new ready pulses
            if (input_ready_pulse) begin
                if (read_hold == 0 && !pending_valid) begin
                    // publish immediately
                    bank_calc <= ~bank_in;
`ifndef SYNTHESIS
                    run_fft(~bank_in);
`endif
                    o_done <= 1;
                    read_hold <= 7'd64;
                end else if (!pending_valid) begin
                    // compute now, publish later
                    pending_bank <= ~bank_in;
`ifndef SYNTHESIS
                    run_fft(~bank_in);
`endif
                    pending_valid <= 1;
                end
                // if pending already set and hold active, drop extra symbol
            end
        end
    end

    // =========================================================
    // 5) Output read port (synchronous to clk, combinational data)
    // =========================================================
    assign o_read_re = mem_re[{bank_calc, i_read_addr}];
    assign o_read_im = mem_im[{bank_calc, i_read_addr}];

endmodule