`timescale 1ns/1ps

module tb_mimo_zf;

    // ============================================================
    // 1. 參數設定
    // ============================================================
    parameter DATA_W = 16;
    parameter CLK_PERIOD = 10; // 100MHz
    parameter NFFT = 64;

    // ============================================================
    // 2. 訊號宣告
    // ============================================================
    reg clk;
    reg rst_n;
    reg in_valid;

    // Y Input
    reg signed [DATA_W-1:0] Y1_re, Y1_im;
    reg signed [DATA_W-1:0] Y2_re, Y2_im;

    // H_inv Input
    reg signed [DATA_W-1:0] h00_re, h00_im;
    reg signed [DATA_W-1:0] h01_re, h01_im;
    reg signed [DATA_W-1:0] h10_re, h10_im;
    reg signed [DATA_W-1:0] h11_re, h11_im;

    // Outputs
    wire out_valid;
    wire signed [DATA_W-1:0] X1_re, X1_im;
    wire signed [DATA_W-1:0] X2_re, X2_im;

    // ============================================================
    // 3. 記憶體宣告 (儲存 Golden Vectors)
    // ============================================================
    // FFT Out (Y): 格式 RE1 IM1 RE2 IM2 (共 4*16 = 64 bits)
    reg [63:0] mem_fft [0:NFFT-1];
    
    // ZF Out (Expected X): 格式 RE1 IM1 RE2 IM2
    reg [63:0] mem_zf_gold [0:NFFT-1];

    // H Inv: 格式 RE IM (共 32 bits)，檔案有 4 行
    reg [31:0] mem_h [0:3]; 

    integer i, err_cnt;

    // ============================================================
    // 4. DUT 實例化
    // ============================================================
    mimo_zf_2x2 #(
        .DATA_W(DATA_W)
    ) u_dut (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(in_valid),
        
        .Y1_re(Y1_re), .Y1_im(Y1_im),
        .Y2_re(Y2_re), .Y2_im(Y2_im),

        .h00_re(h00_re), .h00_im(h00_im),
        .h01_re(h01_re), .h01_im(h01_im),
        .h10_re(h10_re), .h10_im(h10_im),
        .h11_re(h11_re), .h11_im(h11_im),

        .out_valid(out_valid),
        .X1_re(X1_re), .X1_im(X1_im),
        .X2_re(X2_re), .X2_im(X2_im)
    );

    // ============================================================
    // 5. 時脈產生
    // ============================================================
    always #(CLK_PERIOD/2) clk = ~clk;

    // ============================================================
    // 6. 測試流程 (Stimulus)
    // ============================================================
    initial begin
        // --- 6.1 讀取檔案 ---
        // 請確保這些 hex 檔跟 tb_mimo_zf.v 在同一個目錄，或使用絕對路徑
        $readmemh("golden_fft_out.hex", mem_fft);
        $readmemh("golden_h_inv.hex", mem_h);
        $readmemh("golden_zf_out.hex", mem_zf_gold);

        // --- 6.2 初始化 ---
        clk = 0;
        rst_n = 0;
        in_valid = 0;
        err_cnt = 0;
        Y1_re=0; Y1_im=0; Y2_re=0; Y2_im=0;
        
        // --- 6.3 載入 H 矩陣 (固定值) ---
        // 根據 C Code 寫入順序：H00, H01, H10, H11
        {h00_re, h00_im} = mem_h[0];
        {h01_re, h01_im} = mem_h[1];
        {h10_re, h10_im} = mem_h[2];
        {h11_re, h11_im} = mem_h[3];

        $display("-------------------------------------------");
        $display("Loading H Matrix:");
        $display("H00 = %d + j%d", h00_re, h00_im);
        $display("H01 = %d + j%d", h01_re, h01_im);
        $display("-------------------------------------------");

        // --- 6.4 Reset Release ---
        #(CLK_PERIOD * 5);
        rst_n = 1;
        #(CLK_PERIOD * 2);

        // --- 6.5 開始灌入資料 ---
        $display("Starting Input Data Feed...");
        $display("Idx |   Y1_re   Y1_im   |   Y2_re   Y2_im   | Expected H*Y approx?"); // 表頭

        for (i = 0; i < NFFT; i = i + 1) begin
            @(posedge clk);
            #1;
            in_valid = 1;
            {Y1_re, Y1_im, Y2_re, Y2_im} = mem_fft[i];

            // [新增] 印出目前送入 DUT 的數值
            $display("%3d | %7d %7d | %7d %7d", i, Y1_re, Y1_im, Y2_re, Y2_im);
        end

        // 結束輸入
        @(posedge clk);
        #1;
        in_valid = 0;
        Y1_re=0; Y1_im=0; Y2_re=0; Y2_im=0;

        // 等待 pipeline 輸出完成
        #(CLK_PERIOD * 20);
        
        $display("-------------------------------------------");
        if (err_cnt == 0) 
            $display("TEST PASSED! All outputs match Golden Vectors.");
        else 
            $display("TEST FAILED! Total Errors: %d", err_cnt);
        $display("-------------------------------------------");
        
        $stop;
    end

    // ============================================================
    // 7. 自動比對 (Scoreboard)
    // ============================================================
    integer out_idx = 0;
    reg signed [DATA_W-1:0] exp_x1_re, exp_x1_im, exp_x2_re, exp_x2_im;
    
    // 定義容許誤差 (Tolerance)
    // 定點數運算在 C 跟 Verilog 的截位方式 (Floor vs Round) 
    // 可能導致 1~2 LSB 的差異，這是正常的。
    integer TOL = 2; 

    always @(posedge clk) begin
        if (out_valid) begin
            // 讀取對應的 Golden Output
            {exp_x1_re, exp_x1_im, exp_x2_re, exp_x2_im} = mem_zf_gold[out_idx];

            // 比對 X1 Real
            if ( abs(X1_re - exp_x1_re) > TOL ) begin
                $display("[ERROR] Idx %0d: X1_re Exp=%d, Got=%d", out_idx, exp_x1_re, X1_re);
                err_cnt = err_cnt + 1;
            end
            // 比對 X1 Imag
            if ( abs(X1_im - exp_x1_im) > TOL ) begin
                $display("[ERROR] Idx %0d: X1_im Exp=%d, Got=%d", out_idx, exp_x1_im, X1_im);
                err_cnt = err_cnt + 1;
            end
            // 比對 X2 Real
            if ( abs(X2_re - exp_x2_re) > TOL ) begin
                $display("[ERROR] Idx %0d: X2_re Exp=%d, Got=%d", out_idx, exp_x2_re, X2_re);
                err_cnt = err_cnt + 1;
            end
            // 比對 X2 Imag
            if ( abs(X2_im - exp_x2_im) > TOL ) begin
                $display("[ERROR] Idx %0d: X2_im Exp=%d, Got=%d", out_idx, exp_x2_im, X2_im);
                err_cnt = err_cnt + 1;
            end

            out_idx = out_idx + 1;
        end
    end

    // 輔助函數：取絕對值
    function integer abs;
        input integer val;
        begin
            abs = (val < 0) ? -val : val;
        end
    endfunction

endmodule