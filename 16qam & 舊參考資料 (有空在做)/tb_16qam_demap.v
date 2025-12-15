`timescale 1ns/1ps

module tb_qpsk_demap;

    reg clk;
    reg rst_n;
    
    // Input buffers
    reg signed [15:0] fft_out_re [0:63];
    reg signed [15:0] fft_out_im [0:63];
    reg golden_bits [0:255]; // 64 symbols * 4 bits = 256 bits
    
    // Interface
    reg i_valid;
    reg signed [15:0] i_re, i_im;
    wire o_valid;
    wire o_bit0, o_bit1, o_bit2, o_bit3;
    
    integer i;
    integer err_cnt;

    // DUT Instantiation
    qam16_demap u_demap (  // 改用 16-QAM 模組
        .clk(clk),
        .rst_n(rst_n),
        .i_valid(i_valid),
        .i_re(i_re),
        .i_im(i_im),
        .o_valid(o_valid),
        .o_bit0(o_bit0),
        .o_bit1(o_bit1),
        .o_bit2(o_bit2), // 連接新 port
        .o_bit3(o_bit3)  // 連接新 port
    );

    // Clock Gen
    initial clk = 0;
    always #5 clk = ~clk;

    // =========================================================
    // Feed Data Process
    // =========================================================
    initial begin
        $display("========================================");
        $display("   Step 10: 16qam Demapping Simulation   ");
        $display("========================================");

        // 1. Load Data
        $readmemh("step9_output_re_16qam.txt", fft_out_re);
        $readmemh("step9_output_im_16qam.txt", fft_out_im);
        $readmemb("step10_output_16qam.txt", golden_bits);

        // 2. Reset
        i_valid = 0; i_re = 0; i_im = 0;
        rst_n = 1; 
        #20 rst_n = 0;
        #20 rst_n = 1;
        #20; // Wait for stable

        $display("[TB] Feeding 64 FFT outputs...");

        // 3. Feed Data (Synchronized to Negedge to be safe)
        for (i = 0; i < 64; i = i + 1) begin
            @(negedge clk); // Change inputs on falling edge
            i_valid = 1;
            i_re = fft_out_re[i];
            i_im = fft_out_im[i];

            $display("[TB Feed] i=%d, re=%h, im=%h", i, fft_out_re[i], fft_out_im[i]);
            
            // Debug print for input
            //$display("[TB Feed] i=%d, data=(%h, %h)", i, fft_out_re[i], fft_out_im[i]);
        end
        
        @(negedge clk);
        i_valid = 0;
        i_re = 0;
        i_im = 0;
        
        // Wait for output to flush
    
        if (err_cnt == 0)
            $display("RESULT: PASS (All 256 bits matched)");
        else
            $display("RESULT: FAIL (Total bit errors: %d)", err_cnt);
            
        $finish;
    end

    // =========================================================
    // Check Process (Event-Driven)
    // =========================================================
    integer check_idx = 0;
    reg g_b0, g_b1, g_b2, g_b3;

    // 初始化 error counter
    initial err_cnt = 0;

    always @(posedge clk) begin
        if (rst_n && o_valid) begin
            
            // 1. 取得 Golden
            g_b0 = golden_bits[check_idx*4];
            g_b1 = golden_bits[check_idx*4 + 1];
            g_b2 = golden_bits[check_idx*4 + 2];
            g_b3 = golden_bits[check_idx*4 + 3];

            // 2. Debug Print
            if (check_idx < 5) begin
                 $display("[Check %2d] DUT=(%b%b%b%b) | Golden=(%b%b%b%b)", 
                          check_idx, o_bit0, o_bit1, o_bit2, o_bit3, g_b0, g_b1, g_b2, g_b3);
            end

            // 3. Compare (檢查所有 4 個 bits)
            if (o_bit0 !== g_b0) begin
                $display("  [ERROR] Symbol %d (Bit 0): Out=%b, Golden=%b", check_idx, o_bit0, g_b0);
                err_cnt = err_cnt + 1;
            end
            if (o_bit1 !== g_b1) begin
                $display("  [ERROR] Symbol %d (Bit 1): Out=%b, Golden=%b", check_idx, o_bit1, g_b1);
                err_cnt = err_cnt + 1;
            end
            if (o_bit2 !== g_b2) begin
                $display("  [ERROR] Symbol %d (Bit 2): Out=%b, Golden=%b", check_idx, o_bit2, g_b2);
                err_cnt = err_cnt + 1;
            end
            if (o_bit3 !== g_b3) begin
                $display("  [ERROR] Symbol %d (Bit 3): Out=%b, Golden=%b", check_idx, o_bit3, g_b3);
                err_cnt = err_cnt + 1;
            end
            
            // 4. Advance Index
            check_idx = check_idx + 1;
        end
    end

    // --- 波形輸出 ---
    initial begin
        $dumpfile("qam16_demap.vcd");
        $dumpvars(0, tb_qpsk_demap); 
    end

endmodule