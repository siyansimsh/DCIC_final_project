module qpsk_demap #(
    parameter DATA_W = 16
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,
    input  wire signed [DATA_W-1:0] X1_re,
    input  wire signed [DATA_W-1:0] X1_im,
    input  wire signed [DATA_W-1:0] X2_re,
    input  wire signed [DATA_W-1:0] X2_im,

    output reg                      out_valid,
    output reg [3:0]                bits_out   // [b0_ch1,b1_ch1,b0_ch2,b1_ch2]
);
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;
            bits_out  <= 4'b0;
        end else begin
            out_valid <= in_valid;
            if (in_valid) begin
                bits_out[3] <= (X1_re < 0); // b0_ch1
                bits_out[2] <= (X1_im < 0); // b1_ch1
                bits_out[1] <= (X2_re < 0); // b0_ch2
                bits_out[0] <= (X2_im < 0); // b1_ch2
            end
        end
    end
endmodule
