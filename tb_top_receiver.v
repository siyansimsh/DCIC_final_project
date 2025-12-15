`timescale 1ns/1ps

module tb_top_receiver;

    // 1. 參數設定
    parameter DATA_W = 16;
    parameter CLK_PERIOD = 10; // 100MHz
    parameter NFFT = 64;
    parameter NCP  = 16;
    // 根據 C Code，Stream Length = NFFT + 3*NCP = 64 + 48 = 112
    // 但我們只看一個 Symbol 的有效輸出 (64點)
    
    // 2. 訊號宣告
    reg clk;
    reg rst_n;
    reg in_valid;
    
    // ADC Inputs (Time Domain)
    reg signed [DATA_W-1:0] adc1_re, adc1_im;
    reg signed [DATA_W-1:0] adc2_re, adc2_im;

    // H_inv Inputs (Frequency Domain)
    reg signed [DATA_W-1:0] h00_re, h00_im;
    reg signed [DATA_W-1:0] h01_re, h01_im;
    reg signed [DATA_W-1:0] h10_re, h10_im;
    reg signed [DATA_W-1:0] h11_re, h11_im;

    // Outputs
    wire out_valid;
    wire [3:0] rx_bits;
    // Debug taps
    wire zf_valid_dbg;
    wire cpe_phase_valid_dbg;
    wire signed [DATA_W-1:0] zf_x1_re_dbg, zf_x1_im_dbg, zf_x2_re_dbg, zf_x2_im_dbg;
    wire signed [DATA_W-1:0] cpe_phase_err_dbg;
    wire signed [DATA_W-1:0] nco_acc_dbg;

    // 3. 記憶體宣告 (Golden Vectors)
    // Raw ADC with CP+data (golden_adc_time.hex) for sync input
    localparam MEM_ADC_LEN = NFFT + NCP;
    reg [63:0] mem_adc_time   [0:MEM_ADC_LEN-1];
    reg [63:0] mem_final_out  [0:NFFT-1]; // CPE 修正後的黃金輸出 (X1c, X2c)
    reg [63:0] mem_zf_out     [0:NFFT-1]; // ZF 輸出黃金 (X1, X2)
    reg [31:0] mem_h_inv      [0:3];   // 讀取 golden_h_inv.hex
    // 注意：我們沒有直接比對 Demap 後的 Bits (因為 C Code 沒產生 bit hex)
    // 但我們可以觀察波形，或者自己加 logic 比對
    // 為了簡單起見，我們先跑通波形，確認 valid 有拉起，且沒有紅字

    integer i;

    // 4. DUT 實例化 (Top Level)
    mimo_ofdm_rx_top #( .DATA_W(DATA_W) ) u_top (
        .clk(clk),
        .rst_n(rst_n),
        .in_valid(in_valid),
        .adc1_re(adc1_re), .adc1_im(adc1_im),
        .adc2_re(adc2_re), .adc2_im(adc2_im),
        
        // H_inv 連接
        .h00_re(h00_re), .h00_im(h00_im),
        .h01_re(h01_re), .h01_im(h01_im),
        .h10_re(h10_re), .h10_im(h10_im),
        .h11_re(h11_re), .h11_im(h11_im),
        
        .out_valid(out_valid),
        .rx_bits(rx_bits),
        .zf_valid_tap(zf_valid_dbg),
        .cpe_phase_valid_tap(cpe_phase_valid_dbg),
        .zf_x1_re_tap(zf_x1_re_dbg),
        .zf_x1_im_tap(zf_x1_im_dbg),
        .zf_x2_re_tap(zf_x2_re_dbg),
        .zf_x2_im_tap(zf_x2_im_dbg),
        .cpe_phase_err_tap(cpe_phase_err_dbg),
        .nco_acc_tap(nco_acc_dbg)
    );

    // 5. 時脈產生
    always #(CLK_PERIOD/2) clk = ~clk;

    // 6. 主要測試流程
    initial begin
        // --- 讀取檔案 ---
        $readmemh("golden_adc_time.hex", mem_adc_time);
        $readmemh("golden_h_inv.hex", mem_h_inv);
        $readmemh("golden_final_out.hex", mem_final_out);
        $readmemh("golden_zf_out.hex", mem_zf_out);

        // --- 初始化 ---
        clk = 0;
        rst_n = 0;
        in_valid = 0;
        adc1_re=0; adc1_im=0; adc2_re=0; adc2_im=0;
        
        // --- 載入 H 矩陣 ---
        {h00_re, h00_im} = mem_h_inv[0];
        {h01_re, h01_im} = mem_h_inv[1];
        {h10_re, h10_im} = mem_h_inv[2];
        {h11_re, h11_im} = mem_h_inv[3];
        
        $display("-------------------------------------------");
        $display("   TOP LEVEL VERIFICATION START            ");
        $display("-------------------------------------------");
        $display("Loading H Matrix:");
        $display("H00 = %d + j%d", h00_re, h00_im);

        // --- Reset Release ---
        #(CLK_PERIOD * 10);
        rst_n = 1;
        #(CLK_PERIOD * 5);

        // --- 開始灌入時域資料 (ADC Data) ---
        $display("Feeding ADC Data...");
        
        // 模擬一個 Symbol 的輸入 (CP + NFFT)
        // 直接送 golden_adc_time (包含 CP+data)，sync_block 會自行找最佳對位
        for (i = 0; i < NCP + NFFT; i = i + 1) begin
            @(posedge clk);
            #1;
            in_valid = 1;
            {adc1_re, adc1_im, adc2_re, adc2_im} = mem_adc_time[i];
        end

        // ============================================================
        // [新增] 灌入 Dummy Data 以推擠 CPE Ping-Pong Buffer
        // ============================================================
        $display("Feeding Dummy Data to flush pipeline...");
        for (i = 0; i < NFFT + 20; i = i + 1) begin // 多灌一點比較保險
            @(posedge clk);
            #1;
            in_valid = 1;
            adc1_re=0; adc1_im=0; adc2_re=0; adc2_im=0; // 輸入 0
        end
        // ============================================================
        
        // 停止輸入，等待 Pipeline 處理
        @(posedge clk);
        in_valid = 0;
        adc1_re=0; adc1_im=0; adc2_re=0; adc2_im=0;

        $display("Data Feed Done. Waiting for Output...");
        
        // 等待足夠長的時間 (Latency: Rotator + FFT(64+) + ZF(1) + CPE(64+))
        // 估計至少需要 200~300 cycles
        #(CLK_PERIOD * 500);
        
        $display("-------------------------------------------");
        $display("   SIMULATION FINISHED                     ");
        $display("-------------------------------------------");
        $finish;
    end
    
    // 7. 比對 Output：把 golden_final_out 重新 demap，比對 rx_bits
    integer out_idx;
    integer mismatch_cnt;
    integer matched_cnt;
    integer zf_idx;
    integer zf_mismatch_cnt;
    integer zf_matched_cnt;
    reg signed [15:0] g_x1_re, g_x1_im, g_x2_re, g_x2_im;
    reg [3:0] expected_bits;
    integer adc_feed_idx;

    initial begin
        out_idx = 0;
        mismatch_cnt = 0;
        matched_cnt = 0;
        zf_idx = 0;
        zf_mismatch_cnt = 0;
        zf_matched_cnt = 0;
        adc_feed_idx = 0;
    end

    // ZF golden compare: check complex values before CPE
    reg signed [15:0] gz_x1_re, gz_x1_im, gz_x2_re, gz_x2_im;
    // 追加：統計 ZF 與 Golden 的誤差 (sum / max)
    integer abs_err_x1_re, abs_err_x1_im, abs_err_x2_re, abs_err_x2_im;
    integer max_err_x1_re, max_err_x1_im, max_err_x2_re, max_err_x2_im;

    function integer abs16;
        input signed [15:0] v;
        begin
            abs16 = (v < 0) ? -v : v;
        end
    endfunction

    initial begin
        abs_err_x1_re = 0; abs_err_x1_im = 0; abs_err_x2_re = 0; abs_err_x2_im = 0;
        max_err_x1_re = 0; max_err_x1_im = 0; max_err_x2_re = 0; max_err_x2_im = 0;
    end
    always @(posedge clk) begin
        if (zf_valid_dbg && zf_idx < NFFT) begin
            {gz_x1_re, gz_x1_im, gz_x2_re, gz_x2_im} = mem_zf_out[zf_idx];
            $display("[TB][ZF ] t=%t k=%0d DUT=(%0d,%0d)(%0d,%0d) GOLD=(%0d,%0d)(%0d,%0d)", $time, zf_idx,
                zf_x1_re_dbg, zf_x1_im_dbg, zf_x2_re_dbg, zf_x2_im_dbg,
                gz_x1_re, gz_x1_im, gz_x2_re, gz_x2_im);
            // 累積誤差統計
            abs_err_x1_re = abs_err_x1_re + abs16(gz_x1_re - zf_x1_re_dbg);
            abs_err_x1_im = abs_err_x1_im + abs16(gz_x1_im - zf_x1_im_dbg);
            abs_err_x2_re = abs_err_x2_re + abs16(gz_x2_re - zf_x2_re_dbg);
            abs_err_x2_im = abs_err_x2_im + abs16(gz_x2_im - zf_x2_im_dbg);
            if (abs16(gz_x1_re - zf_x1_re_dbg) > max_err_x1_re) max_err_x1_re = abs16(gz_x1_re - zf_x1_re_dbg);
            if (abs16(gz_x1_im - zf_x1_im_dbg) > max_err_x1_im) max_err_x1_im = abs16(gz_x1_im - zf_x1_im_dbg);
            if (abs16(gz_x2_re - zf_x2_re_dbg) > max_err_x2_re) max_err_x2_re = abs16(gz_x2_re - zf_x2_re_dbg);
            if (abs16(gz_x2_im - zf_x2_im_dbg) > max_err_x2_im) max_err_x2_im = abs16(gz_x2_im - zf_x2_im_dbg);

            if ((gz_x1_re !== zf_x1_re_dbg) || (gz_x1_im !== zf_x1_im_dbg) ||
                (gz_x2_re !== zf_x2_re_dbg) || (gz_x2_im !== zf_x2_im_dbg)) begin
                zf_mismatch_cnt = zf_mismatch_cnt + 1;
                $display("[ZFCHK][%0d] T=%t DUT=(%d,%d)(%d,%d) GOLD=(%d,%d)(%d,%d)",
                    zf_idx, $time,
                    zf_x1_re_dbg, zf_x1_im_dbg, zf_x2_re_dbg, zf_x2_im_dbg,
                    gz_x1_re, gz_x1_im, gz_x2_re, gz_x2_im);
            end else begin
                zf_matched_cnt = zf_matched_cnt + 1;
            end
            zf_idx = zf_idx + 1;
        end
    end

    always @(posedge clk) begin
        if (in_valid) begin
            $display("[TB][ADC] t=%t idx=%0d adc1=(%0d,%0d) adc2=(%0d,%0d)", $time, adc_feed_idx, adc1_re, adc1_im, adc2_re, adc2_im);
            adc_feed_idx = adc_feed_idx + 1;
        end
        if (out_valid && out_idx < NFFT) begin
            {g_x1_re, g_x1_im, g_x2_re, g_x2_im} = mem_final_out[out_idx];
            expected_bits[3] = (g_x1_re < 0);
            expected_bits[2] = (g_x1_im < 0);
            expected_bits[1] = (g_x2_re < 0);
            expected_bits[0] = (g_x2_im < 0);

            if (expected_bits !== rx_bits) begin
                mismatch_cnt = mismatch_cnt + 1;
                $display("[CHECK][%0d] Time %t: DUT=%b GOLD=%b", out_idx, $time, rx_bits, expected_bits);
            end else begin
                $display("[CHECK][%0d] Time %t: Match %b", out_idx, $time, rx_bits);
                matched_cnt = matched_cnt + 1;
            end
            out_idx = out_idx + 1;
        end
    end

    // 在模擬結束前報告總結果
    initial begin
        wait(rst_n === 1'b1);
        wait(out_idx == NFFT);
        $display("================================================");
        $display("Golden compare done: %0d mismatches / %0d matches over %0d outputs", mismatch_cnt, matched_cnt, NFFT);
        if (mismatch_cnt == 0)
            $display("STATUS: PASS");
        else
            $display("STATUS: FAIL");
        $display("ZF compare: %0d mismatches / %0d matches over %0d outputs", zf_mismatch_cnt, zf_matched_cnt, NFFT);
        $display("ZF avg abs err X1_re=%0f X1_im=%0f X2_re=%0f X2_im=%0f", 
            abs_err_x1_re / 64.0, abs_err_x1_im / 64.0, abs_err_x2_re / 64.0, abs_err_x2_im / 64.0);
        $display("ZF max abs err X1_re=%0d X1_im=%0d X2_re=%0d X2_im=%0d", 
            max_err_x1_re, max_err_x1_im, max_err_x2_re, max_err_x2_im);
        $display("================================================");
    end

    // 觀察 CPE 迴授給 NCO 的相位修正值以及 NCO 內部累加器
    // 注意：僅為 debug，不改功能
    always @(posedge clk) begin
        if (cpe_phase_valid_dbg) begin
            $display("[DBG][CPE->NCO] T=%t phase_err=%d nco_acc_before=%d", $time,
                cpe_phase_err_dbg, nco_acc_dbg);
        end
    end
    
    // 波形紀錄
    initial begin
        $dumpfile("top_level.vcd");
        $dumpvars(0, tb_top_receiver);
    end

endmodule