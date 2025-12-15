module cfo_rotator #(
    parameter DATA_W  = 16,
    parameter PHASE_W = 16,
    parameter NEGATE_PHASE = 0  // set to 1 to apply -phase for rotation (CFO compensation)
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,

    input  wire signed [DATA_W-1:0] in1_re, in1_im,
    input  wire signed [DATA_W-1:0] in2_re, in2_im,
    input  wire signed [PHASE_W-1:0] phase,   // 來自 NCO 的補償角度

    output reg                      out_valid,
    output reg signed [DATA_W-1:0] out1_re, out1_im,
    output reg signed [DATA_W-1:0] out2_re, out2_im
);

    // 實例化兩個 CORDIC Rotator
    // 注意：這裡假設 phase 是補償角度；若估測出的是 +theta 而 cordic_rotator 需要 -theta，可用 NEGATE_PHASE 切換

    wire signed [PHASE_W-1:0] phase_eff = NEGATE_PHASE ? -phase : phase;

    wire signed [DATA_W-1:0] rot1_re, rot1_im;
    wire signed [DATA_W-1:0] rot2_re, rot2_im;

    cordic_rotator #( .WIDTH(DATA_W) ) u_rot1 (
        .clk(clk), .rst_n(rst_n),
        .x_in(in1_re), .y_in(in1_im), .theta_in(phase_eff),
        .x_out(rot1_re), .y_out(rot1_im)
    );

    cordic_rotator #( .WIDTH(DATA_W) ) u_rot2 (
        .clk(clk), .rst_n(rst_n),
        .x_in(in2_re), .y_in(in2_im), .theta_in(phase_eff),
        .x_out(rot2_re), .y_out(rot2_im)
    );

    // 輸出暫存 (Pipeline Register)
`ifndef SYNTHESIS
    // Debug: print first few rotations to observe effect of phase sign
    integer dbg_rot_cnt;
`endif

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 0;
            out1_re <= 0; out1_im <= 0;
            out2_re <= 0; out2_im <= 0;
`ifndef SYNTHESIS
            dbg_rot_cnt = 0;
`endif
        end else begin
            out_valid <= in_valid;
            out1_re   <= rot1_re;
            out1_im   <= rot1_im;
            out2_re   <= rot2_re;
            out2_im   <= rot2_im;

`ifndef SYNTHESIS
            if (in_valid && dbg_rot_cnt < 8) begin
                $display("[DBG][CFO_ROT] T=%0t step=%0d phase_eff=%0d in1=(%0d,%0d) out1=(%0d,%0d) in2=(%0d,%0d) out2=(%0d,%0d)",
                         $time, dbg_rot_cnt, phase_eff, in1_re, in1_im, rot1_re, rot1_im, in2_re, in2_im, rot2_re, rot2_im);
                dbg_rot_cnt = dbg_rot_cnt + 1;
            end
`endif
        end
    end

endmodule