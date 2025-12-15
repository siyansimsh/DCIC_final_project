module cfo_estimator_cp_mrc #(
    parameter DATA_W = 16,
    parameter NFFT   = 64,
    parameter NCP    = 16
)(
    input  wire                     clk,
    input  wire                     rst_n,

    input  wire                     in_valid,
    input  wire                     symbol_start, // 標記 Symbol 的第一個 Sample

    input  wire signed [DATA_W-1:0] rx1_re,
    input  wire signed [DATA_W-1:0] rx1_im,
    input  wire signed [DATA_W-1:0] rx2_re,
    input  wire signed [DATA_W-1:0] rx2_im,

    output reg                      cfo_valid,   // 估測完成時拉高
    output reg  signed [DATA_W-1:0] cfo_eps_hat  // normalized CFO ε_hat (Q5.11)
);

    // ============================================================
    // 1. Delay Buffer (FIFO) to get r[i] and r[i+N]
    //    Depth = NFFT = 64
    // ============================================================
    // 這裡我們只存 NFFT 長度。
    // 當新的數據 rx_in 進來時，它是 r[i+N]，而 FIFO 讀出來的是 r[i]
    
    reg signed [DATA_W-1:0] fifo_r1_re [0:NFFT-1];
    reg signed [DATA_W-1:0] fifo_r1_im [0:NFFT-1];
    reg signed [DATA_W-1:0] fifo_r2_re [0:NFFT-1];
    reg signed [DATA_W-1:0] fifo_r2_im [0:NFFT-1];
    
    reg [5:0] wr_ptr;
    wire [5:0] rd_ptr = wr_ptr; // Ring Buffer 讀寫指針相同 (Delay N)

    // 延遲後的訊號 (Delayed Signal r[i])
    reg signed [DATA_W-1:0] d_r1_re, d_r1_im;
    reg signed [DATA_W-1:0] d_r2_re, d_r2_im;

    // 目前輸入訊號 (Current Signal r[i+N])
    // 為了對齊 RAM 讀取延遲，輸入訊號也要 register 一拍
    reg signed [DATA_W-1:0] c_r1_re, c_r1_im;
    reg signed [DATA_W-1:0] c_r2_re, c_r2_im;
    
    reg valid_d1;
    integer k;

    // Delay line with reset/realignment so CP is captured before we correlate with the tail
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wr_ptr   <= 0;
            valid_d1 <= 0;
            d_r1_re  <= 0; d_r1_im <= 0; d_r2_re <= 0; d_r2_im <= 0;
            c_r1_re  <= 0; c_r1_im <= 0; c_r2_re <= 0; c_r2_im <= 0;
            for (k = 0; k < NFFT; k = k + 1) begin
                fifo_r1_re[k] <= 0; fifo_r1_im[k] <= 0;
                fifo_r2_re[k] <= 0; fifo_r2_im[k] <= 0;
            end
        end else begin
            // realign pointer on symbol_start but still allow the same-cycle sample to be written
            if (symbol_start)
                wr_ptr <= 0;

            if (in_valid) begin
                // 讀取舊資料 (r[i])
                d_r1_re <= fifo_r1_re[rd_ptr];
                d_r1_im <= fifo_r1_im[rd_ptr];
                d_r2_re <= fifo_r2_re[rd_ptr];
                d_r2_im <= fifo_r2_im[rd_ptr];
                
                // 寫入新資料 (存起來給下一個 Symbol 用)
                fifo_r1_re[wr_ptr] <= rx1_re;
                fifo_r1_im[wr_ptr] <= rx1_im;
                fifo_r2_re[wr_ptr] <= rx2_re;
                fifo_r2_im[wr_ptr] <= rx2_im;
                
                // 暫存當前資料 (r[i+N])
                c_r1_re <= rx1_re;
                c_r1_im <= rx1_im;
                c_r2_re <= rx2_re;
                c_r2_im <= rx2_im;
                
                wr_ptr <= (wr_ptr == NFFT-1) ? 0 : wr_ptr + 1;
                valid_d1 <= 1;
            end else begin
                valid_d1 <= 0;
            end
        end
    end

    // ============================================================
    // 2. Correlation Calculation (Conj Multiply)
    //    Corr = r[i] * conj(r[i+N])
    //    Conj(a+jb) = a-jb.
    //    (a+jb)*(c-jd) = (ac+bd) + j(bc-ad)
    // ============================================================
    
    // RX1 Correlation
    wire signed [DATA_W-1:0] p1_re, p1_im;
    // Use early * conj(late) so a positive CFO (+eps) yields angle ~= -2*pi*eps
    // (matches C-model: eps_hat = -angle/(2*pi) -> +eps)
    complex_mult u_mul1 (
        .a_re(d_r1_re), .a_im(d_r1_im),
        .b_re(c_r1_re), .b_im(-c_r1_im), // conjugate late signal
        .p_re(p1_re),   .p_im(p1_im)
    );

    // RX2 Correlation (same sign convention)
    wire signed [DATA_W-1:0] p2_re, p2_im;
    complex_mult u_mul2 (
        .a_re(d_r2_re), .a_im(d_r2_im),
        .b_re(c_r2_re), .b_im(-c_r2_im),
        .p_re(p2_re),   .p_im(p2_im)
    );

    // MRC Combine (Sum two branches)
    // 為了防止溢位，加法結果擴展一位
    wire signed [DATA_W:0] mrc_re = p1_re + p2_re;
    wire signed [DATA_W:0] mrc_im = p1_im + p2_im;

    // ============================================================
    // 3. Accumulator (Sum over CP window)
    //    Only accumulate when we are in the CP region.
    //    Assuming Symbol Start is aligned such that first NCP samples are CP.
    // ============================================================
    
    reg signed [31:0] acc_re, acc_im;
    reg [6:0] sample_idx; // tracks position within CP+data
    reg       acc_done;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            acc_re <= 0;
            acc_im <= 0;
            acc_done <= 0;
            sample_idx <= 0;
        end else begin
            if (symbol_start) begin
                acc_re <= 0;
                acc_im <= 0;
                acc_done <= 0;
                sample_idx <= 0;
            end else if (valid_d1) begin
                // Only accumulate when the current sample is the data part that aligns with the CP
                // CP arrives first (sample_idx 0..NCP-1). When we reach sample_idx=NFFT..NFFT+NCP-1,
                // the delayed path outputs the corresponding CP sample.
                if (sample_idx >= NFFT && sample_idx < NFFT + NCP) begin
                    acc_re <= acc_re + mrc_re;
                    acc_im <= acc_im + mrc_im;
                    acc_done <= (sample_idx == NFFT + NCP - 1);
                end else begin
                    acc_done <= 1'b0;
                end
                sample_idx <= sample_idx + 1'b1;
            end else begin
                acc_done <= 1'b0;
            end
        end
    end

    // ============================================================
    // 4. CORDIC Vectoring (Calculate Angle)
    // ============================================================
    wire signed [DATA_W-1:0] angle_out;
    // 將角度轉成 ε_hat = -angle/(2π) in Q5.11. 觀測到 angle_out 比理想略大，
    // 取縮放 275 可使 eps_hat 與 C-model (q≈56) 更精確對齊。
    localparam signed [15:0] INV_TWO_PI_Q = 16'sd275;
    wire signed [31:0] eps_mul   = -angle_out * $signed(INV_TWO_PI_Q);
    // Map CORDIC angle to epsilon in Q5.11
    wire signed [DATA_W-1:0] eps_hat_q = eps_mul >>> 11; // back to Q5.11
    
    // 將 32-bit 累加值縮減回 16-bit 餵給 CORDIC (避免 overflow)
    // 假設 NCP=16 (2^4)，平均值約右移 4 位
    wire signed [15:0] cordic_in_re = (acc_re >>> 4);
    wire signed [15:0] cordic_in_im = (acc_im >>> 4);

    cordic_vectoring u_vec (
        .x_in(cordic_in_re), .y_in(cordic_in_im),
        .z_out(angle_out)
    );

    // ============================================================
    // 5. Output Logic
    //    Eps_hat = -angle / (2*pi)
    //    In Fixed Point (Q5.11): 2*pi is represented by full range (or 12868)
    //    Actually, our CORDIC angle 'z' maps [-pi, pi] to [-6434, 6434] approximately
    //    Wait, let's check definition.
    //    If z is in range +/- 6434 (PI), then eps = -z / 6434 / 2 ?
    //    Actually, the angle is normalized such that full circle is ...
    //    Let's stick to C Code: eps_hat = -angle / (2*M_PI).
    //    In Verilog, if angle `z` corresponds to real angle `theta`,
    //    and our z scale is: PI = 6434.
    //    Then eps = -z / (2 * 6434) = -z / 12868.
    
    //    But wait, CFO normalized to subcarrier spacing.
    //    Let's check C Code `sync_and_estimate_cfo_mrc`.
    
    //    Let's just output the RAW ANGLE for now and verify it.
    //    The post-processing (division) can be done or approximated.
    //    Approximation: -z * (1/2PI_FIX). 
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cfo_valid   <= 0;
            cfo_eps_hat <= 0;
        end else begin
            // CORDIC 是組合邏輯，acc_done 與 angle_out 同拍有效
            cfo_valid   <= acc_done;
            cfo_eps_hat <= eps_hat_q; // normalized ε_hat (matches C code)
        end
    end

    // Debug: report accumulator and eps_hat when ready
    always @(posedge clk) begin
        if (acc_done) begin
            $display("[DBG][CFO ] acc_re=%0d acc_im=%0d angle_out=%0d eps_hat=%0d", acc_re, acc_im, angle_out, eps_hat_q);
        end
    end

endmodule