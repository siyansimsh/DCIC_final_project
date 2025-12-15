module mimo_zf_2x2 #(
    parameter DATA_W = 16,
    parameter FRAC   = 11
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,

    // 接收訊號 Y (Frequency Domain)
    input  wire signed [DATA_W-1:0] Y1_re, Y1_im,
    input  wire signed [DATA_W-1:0] Y2_re, Y2_im,

    // H 反矩陣輸入 (H_inv)
    input  wire signed [DATA_W-1:0] h00_re, h00_im, // H_inv[0][0]
    input  wire signed [DATA_W-1:0] h01_re, h01_im, // H_inv[0][1]
    input  wire signed [DATA_W-1:0] h10_re, h10_im, // H_inv[1][0]
    input  wire signed [DATA_W-1:0] h11_re, h11_im, // H_inv[1][1]

    output reg                      out_valid,
    output reg signed [DATA_W-1:0] X1_re, X1_im,
    output reg signed [DATA_W-1:0] X2_re, X2_im
);

    // ============================================================
    // 1. 定義中間乘積訊號 (Wires)
    // ============================================================
    // 公式:
    // X1 = (H00 * Y1) + (H01 * Y2)
    // X2 = (H10 * Y1) + (H11 * Y2)
    
    // Term P00 = H00 * Y1
    wire signed [DATA_W-1:0] p00_re, p00_im;
    // Term P01 = H01 * Y2
    wire signed [DATA_W-1:0] p01_re, p01_im;
    // Term P10 = H10 * Y1
    wire signed [DATA_W-1:0] p10_re, p10_im;
    // Term P11 = H11 * Y2
    wire signed [DATA_W-1:0] p11_re, p11_im;

    // ============================================================
    // 2. 實例化 4 個複數乘法器
    // ============================================================
    // 注意：這裡假設 complex_mult 是純組合邏輯，無 clk
    
    // 運算 X1 的部分
    complex_mult #(.WIDTH(DATA_W), .FRAC(FRAC)) u_mul_00 (
        .a_re(h00_re), .a_im(h00_im), .b_re(Y1_re), .b_im(Y1_im),
        .p_re(p00_re), .p_im(p00_im)
    );

    complex_mult #(.WIDTH(DATA_W), .FRAC(FRAC)) u_mul_01 (
        .a_re(h01_re), .a_im(h01_im), .b_re(Y2_re), .b_im(Y2_im),
        .p_re(p01_re), .p_im(p01_im)
    );

    // 運算 X2 的部分
    complex_mult #(.WIDTH(DATA_W), .FRAC(FRAC)) u_mul_10 (
        .a_re(h10_re), .a_im(h10_im), .b_re(Y1_re), .b_im(Y1_im),
        .p_re(p10_re), .p_im(p10_im)
    );

    complex_mult #(.WIDTH(DATA_W), .FRAC(FRAC)) u_mul_11 (
        .a_re(h11_re), .a_im(h11_im), .b_re(Y2_re), .b_im(Y2_im),
        .p_re(p11_re), .p_im(p11_im)
    );

    // ============================================================
    // 3. 飽和加法邏輯 (Combinational Function)
    // ============================================================
    // 為了符合 C code 的 sat_add，我們必須處理溢位
    // 兩個 16-bit 相加可能變成 17-bit
    function signed [DATA_W-1:0] sat_add;
        input signed [DATA_W-1:0] a;
        input signed [DATA_W-1:0] b;
        reg signed [DATA_W:0] sum; // 17 bits

        // 為了安全，建議定義明確的 signed 邊界變數
        reg signed [DATA_W:0] max_limit;
        reg signed [DATA_W:0] min_limit;

        begin
            sum = a + b;
            // 透過賦值給 signed reg，確保它是 signed
            max_limit = $signed({1'b0, {(DATA_W-1){1'b1}}}); //  32767
            min_limit = $signed({1'b1, {(DATA_W-1){1'b0}}}); // -32768

            // 檢查是否超過 16-bit 範圍 (32767 ~ -32768)
            if (sum > max_limit)
                sat_add = max_limit[DATA_W-1:0];
            else if (sum < min_limit)
                sat_add = min_limit[DATA_W-1:0];
            else
                sat_add = sum[DATA_W-1:0];
        end
    endfunction

    // ============================================================
    // 4. 時序邏輯：加法與輸出暫存 (Latency = 1)
    // ============================================================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;
            X1_re <= 0; X1_im <= 0;
            X2_re <= 0; X2_im <= 0;
        end else begin
            // 1. Valid 訊號延遲 1 個 Cycle (因為運算也花 1 個 Cycle 存入 Reg)
            out_valid <= in_valid;

            if (in_valid) begin
                // 2. 執行加法並存入暫存器 (Pipelining)
                // X1 = P00 + P01
                X1_re <= sat_add(p00_re, p01_re);
                X1_im <= sat_add(p00_im, p01_im);

                // X2 = P10 + P11
                X2_re <= sat_add(p10_re, p11_re);
                X2_im <= sat_add(p10_im, p11_im);
            end else begin
                // 當沒有有效輸入時，保持 0 或最後值皆可，這裡清零較安全
                X1_re <= 0; X1_im <= 0;
                X2_re <= 0; X2_im <= 0;
            end
        end
    end

endmodule