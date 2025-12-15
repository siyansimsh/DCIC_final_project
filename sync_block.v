// CP-based synchronization block for 2x2 MIMO.
// Buffers one OFDM symbol (CP + data), runs CP autocorrelation
// across offsets 0..NCP-1 (MRC across two RX chains), picks the
// best start index, then streams out the aligned CP+data window
// with a one-cycle symbol_start pulse.
module sync_block #(
    parameter DATA_W = 16,
    parameter NFFT   = 64,
    parameter NCP    = 16
)(
    input  wire                     clk,
    input  wire                     rst_n,
    input  wire                     in_valid,
    input  wire signed [DATA_W-1:0] in1_re,
    input  wire signed [DATA_W-1:0] in1_im,
    input  wire signed [DATA_W-1:0] in2_re,
    input  wire signed [DATA_W-1:0] in2_im,

    output reg                      out_valid,
    output reg                      symbol_start,
    output reg  signed [DATA_W-1:0] out1_re,
    output reg  signed [DATA_W-1:0] out1_im,
    output reg  signed [DATA_W-1:0] out2_re,
    output reg  signed [DATA_W-1:0] out2_im,
    output wire [4:0]               dbg_best_idx
);
    localparam BUF_LEN = NFFT + 2*NCP; // 96: guard + CP + data
    localparam SYM_LEN = NFFT + NCP;  // length to output (CP + data)

    // Buffer first CP+data span (prepadded with NCP guard to allow positive offset detection)
    reg signed [DATA_W-1:0] buf1_re [0:BUF_LEN-1];
    reg signed [DATA_W-1:0] buf1_im [0:BUF_LEN-1];
    reg signed [DATA_W-1:0] buf2_re [0:BUF_LEN-1];
    reg signed [DATA_W-1:0] buf2_im [0:BUF_LEN-1];

    reg [6:0] fill_cnt;

    // Correlation search
    // Search window starts anywhere up to the last index that keeps a full CP and NFFT in range.
    localparam SEARCH_MAX = BUF_LEN - NFFT - NCP; // 16 when NFFT=64, NCP=16
    reg [4:0] m_idx; // 0..SEARCH_MAX
    reg [4:0] i_idx; // 0..NCP-1
    reg signed [47:0] sum_re, sum_im;
    reg [95:0] best_metric;
    reg [4:0] best_idx;
    assign dbg_best_idx = best_idx;

    // Exponential moving average of metrics (alpha=1/16)
    reg [95:0] metric_ema [0:SEARCH_MAX];
    integer j;

    // Streaming aligned symbol
    reg [6:0] out_idx;

    // Simple Verilog state encoding (avoid typedef for Icarus compatibility)
    localparam S_IDLE   = 3'd0;
    localparam S_FILL   = 3'd1;
    localparam S_CORR   = 3'd2;
    localparam S_STREAM = 3'd3;
    reg [2:0] state;

    // Products for current tap (combinational)
    wire signed [DATA_W-1:0] a1_re = buf1_re[m_idx + i_idx];
    wire signed [DATA_W-1:0] a1_im = buf1_im[m_idx + i_idx];
    wire signed [DATA_W-1:0] b1_re = buf1_re[m_idx + i_idx + NFFT];
    wire signed [DATA_W-1:0] b1_im = buf1_im[m_idx + i_idx + NFFT];

    wire signed [DATA_W-1:0] a2_re = buf2_re[m_idx + i_idx];
    wire signed [DATA_W-1:0] a2_im = buf2_im[m_idx + i_idx];
    wire signed [DATA_W-1:0] b2_re = buf2_re[m_idx + i_idx + NFFT];
    wire signed [DATA_W-1:0] b2_im = buf2_im[m_idx + i_idx + NFFT];

    // r * conj(r_shift)
    wire signed [31:0] p1_re = a1_re * b1_re + a1_im * b1_im;
    wire signed [31:0] p1_im = a1_im * b1_re - a1_re * b1_im;
    wire signed [31:0] p2_re = a2_re * b2_re + a2_im * b2_im;
    wire signed [31:0] p2_im = a2_im * b2_re - a2_re * b2_im;

    wire signed [32:0] p_re = p1_re + p2_re; // MRC sum
    wire signed [32:0] p_im = p1_im + p2_im;

    wire signed [47:0] sum_re_next = sum_re + {{15{p_re[32]}}, p_re};
    wire signed [47:0] sum_im_next = sum_im + {{15{p_im[32]}}, p_im};

    wire [95:0] metric_next = (sum_re_next * sum_re_next) + (sum_im_next * sum_im_next);

    // EMA helper wires (avoid block-local declarations)
    wire signed [96:0] ema_diff = {1'b0, metric_next} - {1'b0, metric_ema[m_idx]};
    wire [95:0] ema_candidate   = metric_ema[m_idx] + (ema_diff >>> 4);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state        <= S_IDLE;
            fill_cnt     <= 0;
            m_idx        <= 0;
            i_idx        <= 0;
            sum_re       <= 0;
            sum_im       <= 0;
            best_metric  <= 0;
            best_idx     <= 0;
            out_idx      <= 0;
            out_valid    <= 0;
            symbol_start <= 0;
            out1_re      <= 0; out1_im <= 0; out2_re <= 0; out2_im <= 0;
            for (j = 0; j <= SEARCH_MAX; j = j + 1) begin
                metric_ema[j] <= 0;
            end
            // Clear buffer to avoid Xs in correlation window (prepended guard stays zero)
            for (j = 0; j < BUF_LEN; j = j + 1) begin
                buf1_re[j] <= 0; buf1_im[j] <= 0;
                buf2_re[j] <= 0; buf2_im[j] <= 0;
            end
        end else begin
            // default strobes
            out_valid    <= 0;
            symbol_start <= 0;

            case (state)
                S_IDLE: begin
                    fill_cnt <= 0;
                    if (in_valid) begin
                        // Prepend a zero guard of NCP samples; store incoming stream starting at offset NCP
                        buf1_re[NCP] <= in1_re; buf1_im[NCP] <= in1_im;
                        buf2_re[NCP] <= in2_re; buf2_im[NCP] <= in2_im;
                        fill_cnt     <= 1;
                        state        <= S_FILL;
                    end
                end
                S_FILL: begin
                    if (in_valid) begin
                        // Shift writes by NCP so that CP sits at index 16 and correlation peak lands at best_idx=16
                        buf1_re[fill_cnt + NCP] <= in1_re; buf1_im[fill_cnt + NCP] <= in1_im;
                        buf2_re[fill_cnt + NCP] <= in2_re; buf2_im[fill_cnt + NCP] <= in2_im;
                        if (fill_cnt == SYM_LEN-1) begin // captured NCP+NFFT samples
                            state   <= S_CORR;
                            m_idx   <= 0;
                            i_idx   <= 0;
                            sum_re  <= 0;
                            sum_im  <= 0;
                            best_metric <= 0;
                            best_idx    <= 0;
                        end
                        fill_cnt <= fill_cnt + 1'b1;
                    end
                end
                S_CORR: begin
                    // accumulate one tap per cycle
                    sum_re <= sum_re_next;
                    sum_im <= sum_im_next;

                    if (i_idx == NCP-1) begin
                        // EMA update: ema = ema + (metric_next - ema)/16
                        // Use filtered value for best metric selection
                        metric_ema[m_idx] <= metric_ema[m_idx] + (ema_diff >>> 4);
                        if (ema_candidate > best_metric) begin
                            best_metric <= ema_candidate;
                            best_idx    <= m_idx;
                        end
                        i_idx  <= 0;
                        sum_re <= 0;
                        sum_im <= 0;
                        if (m_idx == SEARCH_MAX) begin
                            state   <= S_STREAM;
                            out_idx <= 0;
                        end else begin
                            m_idx <= m_idx + 1'b1;
                        end
                    end else begin
                        i_idx <= i_idx + 1'b1;
                    end
                end
                S_STREAM: begin
                    if (out_idx < SYM_LEN) begin
                        out_valid    <= 1;
                        if (out_idx == 0)
                            symbol_start <= 1;
                        out1_re <= buf1_re[best_idx + out_idx];
                        out1_im <= buf1_im[best_idx + out_idx];
                        out2_re <= buf2_re[best_idx + out_idx];
                        out2_im <= buf2_im[best_idx + out_idx];
                        out_idx <= out_idx + 1'b1;
                    end else begin
                        // done with this symbol, prepare for next capture
                        out_valid <= 0;
                        state <= S_IDLE;
                        fill_cnt <= 0;
                        out_idx <= 0;
                    end
                end
            endcase
        end
    end

    // Debug: report chosen best_idx when symbol_start pulses
    always @(posedge clk) begin
        if (symbol_start)
            $display("[DBG][SYNC] best_idx=%0d T=%t", best_idx, $time);
    end
endmodule
