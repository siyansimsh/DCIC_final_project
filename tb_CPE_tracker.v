`timescale 1ns/1ps

module tb_cpe_tracker;

    // 參數設定
    parameter DATA_W = 16;
    parameter PHASE_W = 16;
    parameter CLK_PERIOD = 10;
    parameter NFFT = 64;

    // 訊號宣告
    reg clk, rst_n, in_valid;
    reg signed [DATA_W-1:0] X1_re, X1_im, X2_re, X2_im;

    wire out_valid, phase_err_valid;
    wire signed [DATA_W-1:0] X1c_re, X1c_im, X2c_re, X2c_im;
    wire signed [PHASE_W-1:0] phase_err;

    // 記憶體 (Input: ZF Out, Expected: Final Out)
    reg [63:0] mem_zf_in   [0:NFFT-1]; // 讀取 golden_zf_out.hex
    reg [63:0] mem_final   [0:NFFT-1]; // 讀取 golden_final_out.hex

    integer i, err_cnt;

    // DUT 實例化
    cpe_tracker #(
        .DATA_W(DATA_W), .PHASE_W(PHASE_W), .NFFT(NFFT)
    ) u_dut (
        .clk(clk), .rst_n(rst_n), .in_valid(in_valid),
        .X1_re(X1_re), .X1_im(X1_im), .X2_re(X2_re), .X2_im(X2_im),
        .out_valid(out_valid),
        .X1c_re(X1c_re), .X1c_im(X1c_im), .X2c_re(X2c_re), .X2c_im(X2c_im),
        .phase_err_valid(phase_err_valid), .phase_err(phase_err)
    );

    always #(CLK_PERIOD/2) clk = ~clk;

    // 驗證流程
    initial begin
        // 1. 讀取檔案
        // 注意：golden_zf_out 是 CPE 的輸入
        $readmemh("golden_zf_out.hex", mem_zf_in); 
        $readmemh("golden_final_out.hex", mem_final);

        // 2. 初始化
        clk = 0; rst_n = 0; in_valid = 0; err_cnt = 0;
        X1_re=0; X1_im=0; X2_re=0; X2_im=0;
        
        #(CLK_PERIOD*5); rst_n = 1; #(CLK_PERIOD*2);

        $display("-------------------------------------------");
        $display("Starting CPE Tracker Verification");
        $display("NOTE: CPE has Ping-Pong buffer delay (64 cycles)");
        $display("-------------------------------------------");

        // 3. 灌入資料 (第一回合：真實資料)
        // 這會存入 Bank 0，算出誤差
        for (i = 0; i < NFFT; i = i + 1) begin
            @(posedge clk);
            #1; // hold time
            in_valid = 1;
            {X1_re, X1_im, X2_re, X2_im} = mem_zf_in[i];
        end

        // 4. 灌入資料 (第二回合：推擠資料)
        // 因為是 Ping-Pong，我們需要繼續餵入 valid 訊號，
        // 這樣 Bank 0 的資料才會被讀出來做修正。
        // 這裡我們餵 0 也可以，反正我們只比對第一組輸出。
        $display("Feeding Dummy Data to flush Ping-Pong buffer...");
        for (i = 0; i < NFFT; i = i + 1) begin
            @(posedge clk);
            #1;
            in_valid = 1;
            X1_re=0; X1_im=0; X2_re=0; X2_im=0; // Dummy
        end

        // 停止輸入
        @(posedge clk); in_valid = 0;

        #(CLK_PERIOD*20);
        
        $display("-------------------------------------------");
        if (err_cnt == 0) $display("TEST PASSED! CPE Output matches Golden.");
        else              $display("TEST FAILED! Errors: %d", err_cnt);
        $display("-------------------------------------------");
        $stop;
    end

    // 自動比對
    integer out_idx = 0;
    reg signed [DATA_W-1:0] exp_re1, exp_im1, exp_re2, exp_im2;
    integer TOL = 2; // 容許誤差

    // 取絕對值函數
    function integer abs(input integer v);
        abs = (v < 0) ? -v : v;
    endfunction

    always @(posedge clk) begin
        if (out_valid) begin
            // 由於我們灌了兩輪資料，out_valid 會拉高兩次 (共128 cycles)
            // 我們只比對前 64 筆 (對應真實資料)
            if (out_idx < NFFT) begin
                {exp_re1, exp_im1, exp_re2, exp_im2} = mem_final[out_idx];

                if (abs(X1c_re - exp_re1) > TOL || abs(X1c_im - exp_im1) > TOL ||
                    abs(X2c_re - exp_re2) > TOL || abs(X2c_im - exp_im2) > TOL) begin
                    $display("[ERR] Idx %d: Got(%d,%d) Exp(%d,%d)", 
                             out_idx, X1c_re, X1c_im, exp_re1, exp_im1);
                    err_cnt = err_cnt + 1;
                end
            end
            out_idx = out_idx + 1;
        end
    end

endmodule