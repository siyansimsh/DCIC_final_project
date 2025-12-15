module cordic_vectoring #(
    parameter WIDTH = 16
)(
    // 純組合邏輯 (如果要跑高頻，建議改成 Pipeline)
    input  wire signed [WIDTH-1:0] x_in,
    input  wire signed [WIDTH-1:0] y_in,

    output reg  signed [WIDTH-1:0] z_out // 輸出的角度 (Theta)
);
    // CORDIC Table (Q5.11)
    wire signed [WIDTH-1:0] atan_table [0:11];
    assign atan_table[0]  = 16'd1608;
    assign atan_table[1]  = 16'd949;
    assign atan_table[2]  = 16'd501;
    assign atan_table[3]  = 16'd254;
    assign atan_table[4]  = 16'd127;
    assign atan_table[5]  = 16'd63;
    assign atan_table[6]  = 16'd31;
    assign atan_table[7]  = 16'd15;
    assign atan_table[8]  = 16'd7;
    assign atan_table[9]  = 16'd3;
    assign atan_table[10] = 16'd1;
    assign atan_table[11] = 16'd0;

    reg signed [WIDTH-1:0] x [0:12];
    reg signed [WIDTH-1:0] y [0:12];
    reg signed [WIDTH-1:0] z [0:12];

    integer i;
    always @(*) begin
        // 1. 初始化與預旋轉 (Pre-rotation)
        // 修正 x < 0 (第二、三象限) 的情況
        // 必須對應 C Code: if (xi < 0) { xi = -xi; yi = -yi; zi = FIX_PI; }
        // FIX_PI = 6434 (根據 C header)
        if (x_in < 0) begin
            x[0] = -x_in;
            y[0] = -y_in;
            // 加上 +PI 以保留原本第二、三象限的角度 (原本寫成 -PI 會導致角度反號)
            z[0] = 16'd6434; // +PI (+180度)
        end else begin
            x[0] = x_in;
            y[0] = y_in;
            z[0] = 0;
        end

        // 2. CORDIC 迭代 (邏輯不變)
        for (i = 0; i < 12; i = i + 1) begin
            // Vectoring 模式：目標是把 y 消成 0
            if (y[i] > 0) begin 
                // 順時針轉
                x[i+1] = x[i] + (y[i] >>> i);
                y[i+1] = y[i] - (x[i] >>> i);
                z[i+1] = z[i] + atan_table[i];
            end else begin
                // 逆時針轉
                x[i+1] = x[i] - (y[i] >>> i);
                y[i+1] = y[i] + (x[i] >>> i);
                z[i+1] = z[i] - atan_table[i];
            end
        end
        
        // 輸出角度
        z_out = z[12];
    end

endmodule