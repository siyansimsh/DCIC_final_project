module mimo_ofdm_rx_top #(
    parameter DATA_W = 16,
    parameter DEBUG_BYPASS_ROT = 0 // enable rotator (set 1 to force phase=0)
)(
    input  wire                 clk,
    input  wire                 rst_n,

    // ADC input
    input  wire                 in_valid,
    input  wire signed [DATA_W-1:0] adc1_re, adc1_im,
    input  wire signed [DATA_W-1:0] adc2_re, adc2_im,

    // [新增] H_inv 輸入 (從 Testbench 讀入 Golden Channel)
    // 雖然這樣做不符合實際晶片應用(實際要有通道估測)，但這是驗證流程的標準做法
    input  wire signed [DATA_W-1:0] h00_re, h00_im,
    input  wire signed [DATA_W-1:0] h01_re, h01_im,
    input  wire signed [DATA_W-1:0] h10_re, h10_im,
    input  wire signed [DATA_W-1:0] h11_re, h11_im,

    // Output
    output wire                 out_valid,
    output wire [3:0]           rx_bits,
    // Debug taps (for testbench visibility)
    output wire                 zf_valid_tap,
    output wire                 cpe_phase_valid_tap,
    output wire signed [DATA_W-1:0] zf_x1_re_tap,
    output wire signed [DATA_W-1:0] zf_x1_im_tap,
    output wire signed [DATA_W-1:0] zf_x2_re_tap,
    output wire signed [DATA_W-1:0] zf_x2_im_tap,
    output wire signed [DATA_W-1:0] cpe_phase_err_tap,
    output wire signed [DATA_W-1:0] nco_acc_tap
);

`ifndef SYNTHESIS
    // Debug counters (throttle prints to first symbol)
    reg [15:0] dbg_in_cnt, dbg_fft_in_cnt, dbg_fft_out_cnt, dbg_zf_cnt, dbg_cpe_cnt;
`endif

    // ====================================================
    // 0. CP-based Synchronization (auto-detect CP start)
    // ====================================================
    wire                 sync_valid;
    wire                 sync_symbol_start;
    wire signed [DATA_W-1:0] sync1_re, sync1_im;
    wire signed [DATA_W-1:0] sync2_re, sync2_im;

    sync_block #( .DATA_W(DATA_W), .NFFT(64), .NCP(16) ) u_sync (
        .clk(clk), .rst_n(rst_n),
        .in_valid(in_valid),
        .in1_re(adc1_re), .in1_im(adc1_im),
        .in2_re(adc2_re), .in2_im(adc2_im),
        .out_valid(sync_valid),
        .symbol_start(sync_symbol_start),
        .out1_re(sync1_re), .out1_im(sync1_im),
        .out2_re(sync2_re), .out2_im(sync2_im),
        .dbg_best_idx()
    );

    // ====================================================
    // 1. CFO Estimation
    // ====================================================
    wire                 cfo_est_valid;
    wire signed [DATA_W-1:0] cfo_eps_hat;

    cfo_estimator_cp_mrc #( .DATA_W(DATA_W) ) u_cfo_est (
        .clk(clk), .rst_n(rst_n),
        .in_valid(sync_valid), .symbol_start(sync_symbol_start),
        .rx1_re(sync1_re), .rx1_im(sync1_im),
        .rx2_re(sync2_re), .rx2_im(sync2_im),
        .cfo_valid(cfo_est_valid),
        .cfo_eps_hat(cfo_eps_hat)
    );

    // ====================================================
    // 1.5 Symbol Buffer for CFO latency hiding
    // ====================================================
    // We need CFO estimate before rotating the symbol, so buffer CP+data then replay once eps_hat is ready.
    reg [6:0] buf_wr_idx;      // 0..79 (16 CP + 64 data)
    reg [6:0] buf_rd_idx;
    reg [6:0] captured_count;  // number of samples captured
    reg [6:0] play_len;        // how many samples to replay
    reg       collecting;      // capturing sync outputs
    reg       sync_valid_d;    // track falling edge of sync_valid
    reg       collected_done;  // asserted when we already latched a full symbol
    reg       playing;         // replaying to rotator/FFT
    reg       playing_d;       // edge detect for debug
    reg       played_once;     // prevent multiple replays of same buffer
    reg       cfo_ready;       // latched when estimator finishes

    // Encourage BRAM inference for the symbol buffer
    (* ram_style = "block" *) reg signed [DATA_W-1:0] buf1_re [0:79];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] buf1_im [0:79];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] buf2_re [0:79];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] buf2_im [0:79];

    // Capture raw (post-sync) samples
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            buf_wr_idx <= 0;
            captured_count <= 0;
            play_len <= 0;
            collecting <= 0;
            collected_done <= 0;
            cfo_ready <= 0;
            sync_valid_d <= 0;
        end else begin
            sync_valid_d <= sync_valid;

            if (sync_symbol_start) begin
                // start a fresh capture window
                collecting <= 1;
                buf_wr_idx <= 0;
                captured_count <= 0;
                play_len <= 0;
                collected_done <= 0;
                cfo_ready <= 0;

                if (sync_valid) begin
                    buf1_re[0] <= sync1_re;
                    buf1_im[0] <= sync1_im;
                    buf2_re[0] <= sync2_re;
                    buf2_im[0] <= sync2_im;
                    buf_wr_idx <= 7'd1;
                    captured_count <= 7'd1;
                    play_len <= 7'd1;
                end

            end else if (collecting && sync_valid) begin
                // normal capture path
                buf1_re[buf_wr_idx] <= sync1_re;
                buf1_im[buf_wr_idx] <= sync1_im;
                buf2_re[buf_wr_idx] <= sync2_re;
                buf2_im[buf_wr_idx] <= sync2_im;
                captured_count <= buf_wr_idx + 1'b1;

                if (buf_wr_idx == 7'd79) begin
                    collecting <= 0;
                    collected_done <= 1'b1;
                    play_len <= 7'd80; // full CP+data captured
                    $display("[DBG][CAP ] T=%t captured full symbol: count=%0d", $time, buf_wr_idx + 1'b1);
                end else begin
                    buf_wr_idx <= buf_wr_idx + 1'b1;
                    if (buf_wr_idx + 1'b1 >= 7'd79) begin
                        collected_done <= 1'b1;
                        play_len <= buf_wr_idx + 1'b1;
                    end
                end

            end else if (collecting && sync_valid_d && !sync_valid) begin
                // Fallback: sync_valid dropped early; latch what we have
                collecting <= 0;
                play_len <= captured_count;
                if (captured_count != 0)
                    collected_done <= 1'b1;
                $display("[DBG][CAP ] T=%t sync_valid fell, captured_count=%0d", $time, captured_count);
            end

            if (cfo_est_valid) begin
                cfo_ready <= 1'b1;
                // do not force play_len beyond what we actually captured
                $display("[DBG][CFO_RDY] T=%t cfo_ready=1 captured_count=%0d collected_done=%0d play_len=%0d", $time, captured_count, collected_done, play_len);
            end
        end
    end

    // Playback controller
    reg rot_in_valid;
    reg signed [DATA_W-1:0] rot_in1_re, rot_in1_im;
    reg signed [DATA_W-1:0] rot_in2_re, rot_in2_im;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            playing <= 0;
            buf_rd_idx <= 0;
            rot_in_valid <= 0;
            rot_in1_re <= 0; rot_in1_im <= 0;
            rot_in2_re <= 0; rot_in2_im <= 0;
            played_once <= 0;
        end else begin
            rot_in_valid <= 0;
            playing_d <= playing;

            if (sync_symbol_start) begin
                playing <= 0;
                played_once <= 0;
                buf_rd_idx <= 0;
            end

            // Start playback when CFO ready and full symbol captured (80 samples)
            if (!playing && !played_once && cfo_ready && collected_done && play_len != 0) begin
                playing <= 1;
                buf_rd_idx <= 0;
            end

            // Debug: report when playback starts (expect CP[0])
            if (!playing_d && playing) begin
                $display("[DBG][PLAY] T=%t start play_len=%0d captured_count=%0d wr=%0d rd=%0d", $time, play_len, captured_count, buf_wr_idx, buf_rd_idx);
            end

            if (playing) begin
                rot_in_valid <= 1;
                rot_in1_re <= buf1_re[buf_rd_idx];
                rot_in1_im <= buf1_im[buf_rd_idx];
                rot_in2_re <= buf2_re[buf_rd_idx];
                rot_in2_im <= buf2_im[buf_rd_idx];

                if (buf_rd_idx == play_len - 1'b1) begin
                    playing <= 0;
                    played_once <= 1'b1;
                end else begin
                    buf_rd_idx <= buf_rd_idx + 1'b1;
                end
            end
        end
    end

    // ====================================================
    // 2. NCO Loop & Rotation
    // ====================================================
    wire signed [DATA_W-1:0] nco_phase;
    wire signed [DATA_W-1:0] cpe_phase_err; // 來自後端 CPE 的回授
    wire              cpe_phase_valid; // 來自後端 CPE 的 valid

    // NCO Controller (累積相位)
    nco_phase_ctrl #( .PHASE_W(DATA_W) ) u_nco_ctrl (
        .clk(clk), .rst_n(rst_n),
        .cfo_valid(cfo_est_valid),
        .cfo_eps_hat(cfo_eps_hat),
        .phase_update_en(rot_in_valid), // use buffered stream for phase stepping
        .phase_err_valid(cpe_phase_valid),
        .phase_err(cpe_phase_err),
        .phase_out(nco_phase),
        .acc_tap(nco_acc_tap)
    );

    // Rotator (修正時域訊號)
    wire rot_valid;
    wire signed [DATA_W-1:0] rot1_re, rot1_im;
    wire signed [DATA_W-1:0] rot2_re, rot2_im;

    wire signed [DATA_W-1:0] rot_phase = DEBUG_BYPASS_ROT ? {DATA_W{1'b0}} : nco_phase;

    // Apply -phase to cancel positive CFO ramp (NCO accumulates +phase_inc)
    cfo_rotator #( .DATA_W(DATA_W), .PHASE_W(DATA_W), .NEGATE_PHASE(1) ) u_rot (
        .clk(clk), .rst_n(rst_n),
        .in_valid(rot_in_valid),
        .in1_re(rot_in1_re), .in1_im(rot_in1_im),
        .in2_re(rot_in2_re), .in2_im(rot_in2_im),
        .phase(rot_phase),
        .out_valid(rot_valid),
        .out1_re(rot1_re), .out1_im(rot1_im),
        .out2_re(rot2_re), .out2_im(rot2_im)
    );

`ifndef SYNTHESIS
    // Debug: dump the 64 samples actually loaded into the FFT (post-rotator, CP removed)
    integer dbg_fft_in_file;
    integer dbg_fft_out_file;
    integer dbg_zf_out_file;
    reg [6:0] dbg_fft_in_idx;
    reg [6:0] dbg_fft_out_idx;
    reg [6:0] dbg_zf_out_idx;
    initial begin
        dbg_fft_in_file = $fopen("rtl_fft_in_dump.hex", "w");
        dbg_fft_out_file = $fopen("rtl_fft_out_dump.hex", "w");
        dbg_zf_out_file  = $fopen("rtl_zf_out_dump.hex", "w");
        dbg_fft_in_idx = 0;
        dbg_fft_out_idx = 0;
        dbg_zf_out_idx  = 0;
    end
    always @(posedge clk) begin
        if (!rst_n) begin
            dbg_fft_in_idx <= 0;
            dbg_fft_out_idx <= 0;
            dbg_zf_out_idx  <= 0;
        end else if (sync_symbol_start) begin
            dbg_fft_in_idx  <= 0; // restart per symbol
            dbg_fft_out_idx <= 0;
            dbg_zf_out_idx  <= 0;
        end else if (fft_load_valid && dbg_fft_in_idx < 64) begin
            // Format: re1 im1 re2 im2 (each 16-bit) per line
            $fdisplay(dbg_fft_in_file, "%04h%04h%04h%04h", rot1_re[15:0], rot1_im[15:0], rot2_re[15:0], rot2_im[15:0]);
            dbg_fft_in_idx <= dbg_fft_in_idx + 1'b1;
        end

        // Dump FFT outputs as seen by ZF input ordering
        if (fft_data_valid_q && dbg_fft_out_idx < 64) begin
            $fdisplay(dbg_fft_out_file, "%04h%04h%04h%04h", fft1_re_q[15:0], fft1_im_q[15:0], fft2_re_q[15:0], fft2_im_q[15:0]);
            dbg_fft_out_idx <= dbg_fft_out_idx + 1'b1;
        end

        // Dump ZF outputs before CPE
        if (zf_valid_out && dbg_zf_out_idx < 64) begin
            $fdisplay(dbg_zf_out_file, "%04h%04h%04h%04h", X1_re[15:0], X1_im[15:0], X2_re[15:0], X2_im[15:0]);
            dbg_zf_out_idx <= dbg_zf_out_idx + 1'b1;
        end
    end
`endif

    // ====================================================
    // 3. FFT (Instantiate x2)
    // ====================================================
    wire fft1_done, fft2_done;
    wire signed [DATA_W-1:0] fft1_out_re, fft1_out_im;
    wire signed [DATA_W-1:0] fft2_out_re, fft2_out_im;
    
    // FFT Read Controller Signals
    reg [5:0] fft_read_addr;
    reg       fft_read_en;

    // Drop CP before feeding FFT
    reg [6:0] sym_sample_idx;
    wire      fft_load_valid = rot_valid && (sym_sample_idx > 15) && (sym_sample_idx < 80);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sym_sample_idx <= 0;
        end else if (sync_symbol_start) begin
            sym_sample_idx <= 0;
        end else if (rot_valid) begin
            sym_sample_idx <= sym_sample_idx + 1'b1;
        end
    end
    
    // RX1 FFT
    fft_64pt u_fft1 (
        .clk(clk), .rst_n(rst_n),
        .i_start(1'b0), // Auto-start by load
        .o_done(fft1_done),
        .i_load_valid(fft_load_valid),
        .i_load_re(rot1_re), .i_load_im(rot1_im),
        .i_read_addr(fft_read_addr),
        .o_read_re(fft1_out_re), .o_read_im(fft1_out_im)
    );

    // RX2 FFT
    fft_64pt u_fft2 (
        .clk(clk), .rst_n(rst_n),
        .i_start(1'b0),
        .o_done(fft2_done),
        .i_load_valid(fft_load_valid),
        .i_load_re(rot2_re), .i_load_im(rot2_im),
        .i_read_addr(fft_read_addr),
        .o_read_re(fft2_out_re), .o_read_im(fft2_out_im)
    );

    // ====================================================
    // 3.5 FFT Readout Logic (Bridge to ZF)
    // ====================================================
    
    reg [2:0] state_read;
    localparam S_WAIT = 0, S_READ = 1;
    reg fft_data_valid_d; // asserted when registered FFT data are valid (raw)
    reg fft_data_valid_q; // 1-cycle delayed valid aligned with pipelined FFT data
    reg fft_done_pending; // latched FFT done when we are still reading previous frame

    // 將 FFT 讀出資料對齊有效位元，避免位址/valid 交錯
    reg signed [DATA_W-1:0] fft1_re_r, fft1_im_r;
    reg signed [DATA_W-1:0] fft2_re_r, fft2_im_r;
    reg signed [DATA_W-1:0] fft1_re_q, fft1_im_q;
    reg signed [DATA_W-1:0] fft2_re_q, fft2_im_q;
    reg [5:0] fft_read_addr_q;

    wire start_fft_read = (state_read == S_WAIT) && (fft_done_pending || (fft1_done && fft2_done));

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state_read <= S_WAIT;
            fft_read_addr <= 0;
            fft_read_en <= 0; // 這是控制讀取地址的 enable
            fft_data_valid_d <= 0;
            fft_done_pending <= 0;
            fft1_re_r <= 0; fft1_im_r <= 0;
            fft2_re_r <= 0; fft2_im_r <= 0;
        end else begin
            case (state_read)
                S_WAIT: begin
                    fft_read_addr <= 0;
                    fft_read_en <= 0;
                    fft_data_valid_d <= 0;
                    // 當兩顆 FFT 都做完時 (o_done 拉高)，進入讀取狀態
                    if (start_fft_read) begin 
                        state_read <= S_READ;
                        fft_read_en <= 1; // 開始讀取
                        fft_done_pending <= 0; // 消耗掉 pending
                    end
                end
                S_READ: begin
                    // 先擷取目前位址的資料，再在循環尾端更新位址
                    fft1_re_r <= fft1_out_re;
                    fft1_im_r <= fft1_out_im;
                    fft2_re_r <= fft2_out_re;
                    fft2_im_r <= fft2_out_im;
                    fft_data_valid_d <= 1'b1;

                    // Address Counter (0 ~ 63)
                    if (fft_read_addr == 63) begin
                        state_read <= S_WAIT;
                        fft_read_en <= 0; // 讀完最後一筆，停止
                        fft_read_addr <= 0;
                    end else begin
                        fft_read_addr <= fft_read_addr + 1;
                    end

                    // 若在讀取期間又收到新的 FFT done，先暫存，待會讀完再觸發下一輪
                    if (fft1_done && fft2_done)
                        fft_done_pending <= 1;
                end
            endcase
        end
    end

    // Align valid/address with registered FFT outputs for ZF and debug/dumps
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            fft_data_valid_q <= 1'b0;
            fft1_re_q <= 0; fft1_im_q <= 0;
            fft2_re_q <= 0; fft2_im_q <= 0;
            fft_read_addr_q  <= 6'd0;
        end else begin
            fft_data_valid_q <= fft_data_valid_d; // align with pipelined FFT data
            fft1_re_q <= fft1_re_r; fft1_im_q <= fft1_im_r;
            fft2_re_q <= fft2_re_r; fft2_im_q <= fft2_im_r;
            fft_read_addr_q  <= fft_read_addr;
        end
    end

    // ====================================================
    // 4. MIMO ZF
    // ====================================================
    wire zf_valid_out;
    wire signed [DATA_W-1:0] X1_re, X1_im, X2_re, X2_im;

    // Debug taps
    assign zf_valid_tap = zf_valid_out;
    assign zf_x1_re_tap = X1_re;
    assign zf_x1_im_tap = X1_im;
    assign zf_x2_re_tap = X2_re;
    assign zf_x2_im_tap = X2_im;

    mimo_zf_2x2 #( .DATA_W(DATA_W) ) u_zf (
        .clk(clk), .rst_n(rst_n),
        .in_valid(fft_data_valid_q), // 接 FFT Read En (aligned with registered data)
        .Y1_re(fft1_re_q), .Y1_im(fft1_im_q), // 已對齊 valid 的 FFT 輸出
        .Y2_re(fft2_re_q), .Y2_im(fft2_im_q),
        // H_inv 來自 Top Level Input
        .h00_re(h00_re), .h00_im(h00_im),
        .h01_re(h01_re), .h01_im(h01_im),
        .h10_re(h10_re), .h10_im(h10_im),
        .h11_re(h11_re), .h11_im(h11_im),
        
        .out_valid(zf_valid_out),
        .X1_re(X1_re), .X1_im(X1_im),
        .X2_re(X2_re), .X2_im(X2_im)
    );

    // ====================================================
    // 5. CPE Tracker
    // ====================================================
    wire cpe_out_valid;
    wire signed [DATA_W-1:0] X1c_re, X1c_im, X2c_re, X2c_im;

    cpe_tracker #( .DATA_W(DATA_W), .NFFT(64) ) u_cpe (
        .clk(clk), .rst_n(rst_n),
        .in_valid(zf_valid_out),
        .X1_re(X1_re), .X1_im(X1_im),
        .X2_re(X2_re), .X2_im(X2_im),
        
        .out_valid(cpe_out_valid),
        .X1c_re(X1c_re), .X1c_im(X1c_im),
        .X2c_re(X2c_re), .X2c_im(X2c_im),
        
        .phase_err_valid(cpe_phase_valid),
        .phase_err(cpe_phase_err) // 回授給 NCO
    );

    // Debug taps
    assign cpe_phase_valid_tap = cpe_phase_valid;
    assign cpe_phase_err_tap = cpe_phase_err;

    // ====================================================
    // 6. Demap & Output
    // ====================================================
    qpsk_demap #( .DATA_W(DATA_W) ) u_demap (
        .clk(clk), .rst_n(rst_n),
        .in_valid(cpe_out_valid),
        .X1_re(X1c_re), .X1_im(X1c_im),
        .X2_re(X2c_re), .X2_im(X2c_im),
        .out_valid(out_valid),
        .bits_out(rx_bits)
    );

    // =========================================================
    // Debug Probes: 追蹤 Valid 訊號流向 (修正版)
    // =========================================================
`ifndef SYNTHESIS
    reg fft_done_latched;
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            fft_done_latched <= 0;
            dbg_in_cnt <= 0;
            dbg_fft_in_cnt <= 0;
            dbg_fft_out_cnt <= 0;
            dbg_zf_cnt <= 0;
            dbg_cpe_cnt <= 0;
        end else begin
            if (in_valid)
                $display("[DBG][ADC ] T=%t idx=%0d in1=(%0d,%0d) in2=(%0d,%0d)", $time, dbg_in_cnt, adc1_re, adc1_im, adc2_re, adc2_im);
            if (in_valid)
                dbg_in_cnt <= dbg_in_cnt + 1'b1;

            if (sync_symbol_start)
                $display("[DBG][SYNC] T=%t best_symbol_start asserted (sym_sample_idx=%0d)", $time, sym_sample_idx);

            if (cfo_est_valid)
                $display("[DBG][CFO ] T=%t eps_hat=%0d", $time, cfo_eps_hat);

            if (rot_valid && fft_load_valid && dbg_fft_in_cnt < 80) begin
                $display("[DBG][FFT_IN] T=%t idx=%0d load_re1=%0d load_im1=%0d load_re2=%0d load_im2=%0d", 
                    $time, sym_sample_idx, rot1_re, rot1_im, rot2_re, rot2_im);
                dbg_fft_in_cnt <= dbg_fft_in_cnt + 1'b1;
            end

            if (fft1_done) begin
                $display("[DBG][FFT ] T=%t: FFT1 Done Pulse!", $time);
                fft_done_latched <= 1;
            end
            if (fft_data_valid_q && dbg_fft_out_cnt < 64) begin
                $display("[DBG][FFT_OUT] T=%t addr=%0d y1=(%0d,%0d) y2=(%0d,%0d)",
                    $time, fft_read_addr_q, fft1_re_q, fft1_im_q, fft2_re_q, fft2_im_q);
                dbg_fft_out_cnt <= dbg_fft_out_cnt + 1'b1;
            end

            if (zf_valid_out && dbg_zf_cnt < 64) begin
                $display("[DBG][ZF  ] T=%t k=%0d X1=(%0d,%0d) X2=(%0d,%0d)",
                    $time, dbg_zf_cnt, X1_re, X1_im, X2_re, X2_im);
                dbg_zf_cnt <= dbg_zf_cnt + 1'b1;
            end

            if (cpe_out_valid && dbg_cpe_cnt < 64) begin
                $display("[DBG][CPE ] T=%t k=%0d X1c=(%0d,%0d) X2c=(%0d,%0d) phase_err_valid=%b phase_err=%0d", 
                    $time, dbg_cpe_cnt, X1c_re, X1c_im, X2c_re, X2c_im, cpe_phase_valid, cpe_phase_err);
                dbg_cpe_cnt <= dbg_cpe_cnt + 1'b1;
            end
        end
    end
`endif

endmodule