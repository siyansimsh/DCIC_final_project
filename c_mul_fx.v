module complex_mult #(
    parameter WIDTH = 16,
    parameter FRAC  = 11
)(
    input  wire signed [WIDTH-1:0] a_re, a_im,
    input  wire signed [WIDTH-1:0] b_re, b_im,
    output reg  signed [WIDTH-1:0] p_re, p_im
);

    // 1. 中間乘積需要 32 bits (16+16)
    // 這裡使用 wire 讓合成器自動推論出 DSP Block
    wire signed [2*WIDTH-1:0] ac = a_re * b_re;
    wire signed [2*WIDTH-1:0] bd = a_im * b_im;
    wire signed [2*WIDTH-1:0] ad = a_re * b_im;
    wire signed [2*WIDTH-1:0] bc = a_im * b_re;

    // 2. 加減運算 (保留一位進位以防溢位，變成 33 bits)
    wire signed [2*WIDTH:0] re_long = ac - bd;
    wire signed [2*WIDTH:0] im_long = ad + bc;

    // 3. 右移與飽和處理 (Combinational Logic)
    // 回復為截斷：直接算術右移，不做四捨五入
    // Q5.11 * Q5.11 = Q10.22 -> 右移 11 -> Q10.11
    wire signed [2*WIDTH:0] re_shifted = re_long >>> FRAC;
    wire signed [2*WIDTH:0] im_shifted = im_long >>> FRAC;

    // 定義最大最小值 (16-bit signed)
    localparam signed [WIDTH-1:0] MAX_VAL = {1'b0, {(WIDTH-1){1'b1}}}; // 32767
    localparam signed [WIDTH-1:0] MIN_VAL = {1'b1, {(WIDTH-1){1'b0}}}; // -32768

    always @(*) begin
        // 檢查是否超過 16-bit 範圍
        // 如果高位元 (不含最後 16 位) 不是全 0 或全 1，代表溢位
        
        // 實部飽和
        if (re_shifted > MAX_VAL)      p_re = MAX_VAL;
        else if (re_shifted < MIN_VAL) p_re = MIN_VAL;
        else                           p_re = re_shifted[WIDTH-1:0];

        // 虛部飽和
        if (im_shifted > MAX_VAL)      p_im = MAX_VAL;
        else if (im_shifted < MIN_VAL) p_im = MIN_VAL;
        else                           p_im = im_shifted[WIDTH-1:0];
    end

endmodule