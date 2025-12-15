module qam16_demap (
    input  wire               clk,
    input  wire               rst_n,
    
    input  wire               i_valid,
    input  wire signed [15:0] i_re,
    input  wire signed [15:0] i_im,
    
    output reg                o_valid,
    output reg                o_bit0, // I - Sign
    output reg                o_bit1, // I - Mag
    output reg                o_bit2, // Q - Sign
    output reg                o_bit3  // Q - Mag
);

    // 16-QAM Thresholds (from C Code)
    // 0.632455532 * 32768 = 20724
    localparam signed [15:0] TH_P2 = 16'd20724;
    localparam signed [15:0] TH_N2 = -16'd20724;
    localparam signed [15:0] MAX_16 = 16'h7FFF;
    localparam signed [15:0] MIN_16 = 16'h8000;

    // Internal signals for scaled values
    reg signed [15:0] I_scaled, Q_scaled;
    
    // Temporary expanded variables for saturation logic
    reg signed [21:0] I_temp, Q_temp; // 16 bits + 6 bits shift

    always @(*) begin
        // Step 1: Scale Up (<< 6)
        I_temp = {{6{i_re[15]}}, i_re} <<< 6; // Re-scaling to compensate FFT right shifts
        Q_temp = {{6{i_im[15]}}, i_im} <<< 6;

        // Step 2: Saturation
        if (I_temp > MAX_16) I_scaled = MAX_16;
        else if (I_temp < MIN_16) I_scaled = MIN_16;
        else I_scaled = I_temp[15:0];

        if (Q_temp > MAX_16) Q_scaled = MAX_16;
        else if (Q_temp < MIN_16) Q_scaled = MIN_16;
        else Q_scaled = Q_temp[15:0];
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            o_valid <= 1'b0;
            o_bit0 <= 0; o_bit1 <= 0;
            o_bit2 <= 0; o_bit3 <= 0;
        end else begin
            o_valid <= i_valid;
            
            if (i_valid) begin
                // I Component (b0, b1)
                // b0: Sign (0 for positive, 1 for negative) -> Inverse of normal sign bit logic in C? 
                // C Code: (I > 0) ? 0 : 1.  Verilog: I[15] is 0 for pos, 1 for neg. So b0 = I[15].
                o_bit0 <= I_scaled[15]; 
                $display("[DUT Recv] i_re=%h, i_im=%h", i_re, i_im);
                
                // b1: Magnitude (Outer const = 1, Inner = 0)
                if (I_scaled > TH_P2 || I_scaled < TH_N2)
                    o_bit1 <= 1'b1;
                else
                    o_bit1 <= 1'b0;

                // Q Component (b2, b3)
                // b2: Sign
                o_bit2 <= Q_scaled[15];
                
                // b3: Magnitude
                if (Q_scaled > TH_P2 || Q_scaled < TH_N2)
                    o_bit3 <= 1'b1;
                else
                    o_bit3 <= 1'b0;
            end
        end
    end

endmodule