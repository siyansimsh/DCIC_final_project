module cordic_rotator #(
    parameter WIDTH = 16
)(
    input  wire                    clk,   // 若要 Pipeline 需加上 clk
    input  wire                    rst_n,
    input  wire signed [WIDTH-1:0] x_in,
    input  wire signed [WIDTH-1:0] y_in,
    input  wire signed [WIDTH-1:0] theta_in, // 角度
    output reg  signed [WIDTH-1:0] x_out,
    output reg  signed [WIDTH-1:0] y_out
);

    // Consume clock/reset to avoid "unused port" warnings in synthesis
    (* keep = "true" *) wire clk_used  = clk;
    (* keep = "true" *) wire rstn_used = rst_n;

    // CORDIC Table (atan(2^-i)) in Q5.11
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

    // 展開 12 級運算 (這裡用組合邏輯示意，實際上建議加 Pipeline Regs)
    reg signed [WIDTH-1:0] x [0:12];
    reg signed [WIDTH-1:0] y [0:12];
    reg signed [WIDTH-1:0] z [0:12];

    integer i;
    always @(*) begin
        // 初始化迭代起點，避免 X 造成整段輸出未知
        x[0] = x_in;
        y[0] = y_in;
        z[0] = theta_in;

        for (i = 0; i < 12; i = i + 1) begin
            // 若 z[i] >= 0 (目標角度是正的)，逆時針；否則順時針
            if (z[i] >= 0) begin
                x[i+1] = x[i] - (y[i] >>> i);
                y[i+1] = y[i] + (x[i] >>> i);
                z[i+1] = z[i] - atan_table[i];
            end else begin
                x[i+1] = x[i] + (y[i] >>> i);
                y[i+1] = y[i] - (x[i] >>> i);
                z[i+1] = z[i] + atan_table[i];
            end
        end
    end

    // 最後乘上 Gain (0.607)
    // 簡單做法： (X * 1243) >> 11
    // 這裡你需要實例化上面的 complex_mult 或者簡單的乘法器
    // 為了簡化範例，這裡直接寫運算
    wire signed [31:0] x_scaled = x[12] * $signed(16'd1243);
    wire signed [31:0] y_scaled = y[12] * $signed(16'd1243);

    always @(*) begin
        x_out = x_scaled >>> 11;
        y_out = y_scaled >>> 11;
    end

endmodule