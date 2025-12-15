module nco_phase_ctrl #(
    parameter PHASE_W = 16  // Q5.11
)(
    input  wire                     clk,
    input  wire                     rst_n,

    // 來自 CFO Estimator
    input  wire                     cfo_valid,
    input  wire signed [PHASE_W-1:0] cfo_eps_hat, // normalized epsilon_hat

    // 控制訊號
    input  wire                     phase_update_en, // 每個 ADC Sample 更新一次
    
    // 來自 CPE Tracker 的回授
    input  wire                     phase_err_valid, // CPE residual valid pulse (per symbol)
    input  wire signed [PHASE_W-1:0] phase_err,      // CPE Residual Error

    // 輸出
    output reg  signed [PHASE_W-1:0] phase_out,      // 給 Rotator 用
    output wire signed [PHASE_W-1:0] acc_tap         // Debug tap
);

    // ============================================================
    // 1. 計算相位增量 (Phase Increment)
    //    Delta_Phi = 2*PI * Eps / NFFT
    //    In Q5.11: 2*PI = 12868 (approx)
    //    NFFT = 64
    //    K = 12868 / 64 = 201.0625 -> round to 201
    //    phase_inc = Eps * 201
    // ============================================================
    
    // 儲存穩定的 phase_inc
    reg signed [PHASE_W-1:0] phase_inc;
    
    // 乘法運算可能需要較大位寬
    // Eps (16-bit) * 201 (8-bit) -> 24-bit -> Shift back ?? 
    // 注意：Eps 已經是 Q5.11 格式的數值 (例如 0.03 -> 61)
    // 相位也是 Q5.11 (PI -> 6434)
    // 所以直接乘常數即可，不需要再 Shift
    // Delta = 61 * 201 = 12261 (太大了！這幾乎是 2PI)
    
    // [更正邏輯]
    // 我們的 Eps 是 normalized to Subcarrier Spacing.
    // 0.03 * 2PI / 64 = 0.03 * 0.098 = 0.0029 rad
    // 0.0029 rad in Q5.11 (x 2048) = 6
    // 如果 Eps=61 (0.03), 則 phase_inc = 61 * 201 / 2048 = 6.
    // 所以應該是 * 201 之後要右移 11 (Q_FRAC)
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            phase_inc <= 0;
        end else if (cfo_valid) begin
            // phase_inc = (eps * 201) >> 11
            // 為了保留精度，先轉 32-bit
            phase_inc <= (cfo_eps_hat * $signed(16'd201)) >>> 11;
        end
    end

    // ============================================================
    // 2. 相位累加器 (Phase Accumulator)
    //    phi[n] = phi[n-1] + phase_inc
    //    Rotator applies -phi to cancel CFO
    // ============================================================
    
    // 定義 PI 範圍以進行 Wrap-around
    localparam signed [PHASE_W-1:0] PI_POS = 16'd6434;
    localparam signed [PHASE_W-1:0] PI_NEG = -16'd6434;
    localparam signed [PHASE_W-1:0] TWO_PI = 16'd12868;

    reg signed [PHASE_W-1:0] acc;
    reg signed [PHASE_W:0]   next_acc; // 多一位用來偵測溢位

`ifndef SYNTHESIS
    // 單次 debug 限制輸出，避免洗版
    integer dbg_step_cnt;
`endif

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            acc <= 0;
`ifndef SYNTHESIS
            dbg_step_cnt = 0;
`endif
        end else begin
            // 優先處理 CPE 回授 (每個 Symbol 一次)
            // 這裡假設 phase_err 在 phase_update_en 為 low 的空檔更新
            // 或者直接疊加
            
            // 簡化邏輯：
            next_acc = acc;

            // Reset phase accumulator when a new CFO estimate is latched
            if (cfo_valid)
                next_acc = 0;

            // 每個 sample 的 CFO 相位累加
            if (phase_update_en) begin
                // Accumulate phase directly; negative phase_inc yields negative ramp
                next_acc = acc + phase_inc;

`ifndef SYNTHESIS
                // 限制輸出 8 次，觀察前幾個 sample 的相位走向
                if (dbg_step_cnt < 8) begin
                    $display("[DBG][NCO ] T=%0t step=%0d phase_inc=%0d acc_before=%0d acc_after=%0d", $time, dbg_step_cnt, phase_inc, acc, next_acc);
                    dbg_step_cnt = dbg_step_cnt + 1;
                end
`endif
            end

            // 在一個 symbol 結束時套用 CPE 殘餘誤差 (一次性 pulse)
            if (phase_err_valid) begin
                // loop_gain = 1; 若需放緩可右移 1~2 位
                next_acc = next_acc - phase_err;
            end

            // Wrap Around Logic (-PI ~ PI)
            if (next_acc > PI_POS)      acc <= next_acc - TWO_PI;
            else if (next_acc < PI_NEG) acc <= next_acc + TWO_PI;
            else                        acc <= next_acc[PHASE_W-1:0];
        end
    end

    // 輸出
    always @(posedge clk) begin
        phase_out <= acc;
    end

    // Debug tap
    assign acc_tap = acc;

endmodule