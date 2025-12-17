module cpe_tracker #(
    parameter DATA_W  = 16,
    parameter PHASE_W = 16,
    parameter NFFT    = 64
)(
    input  wire                     clk,
    input  wire                     rst_n,

    input  wire                     in_valid,
    // input  wire                  symbol_start, // 簡化：我們用內部計數器即可
    // input  wire                  symbol_end,

    input  wire signed [DATA_W-1:0] X1_re, X1_im,
    input  wire signed [DATA_W-1:0] X2_re, X2_im,

    output reg                      out_valid,
    output reg  signed [DATA_W-1:0] X1c_re, X1c_im,
    output reg  signed [DATA_W-1:0] X2c_re, X2c_im,

    output reg                      phase_err_valid,
    output reg  signed [PHASE_W-1:0] phase_err
);

    // ============================================================
    // 1. Vectoring: 計算輸入訊號的角度 (Stage 1)
    // ============================================================
    wire signed [PHASE_W-1:0] ang1, ang2;
    
    cordic_vectoring #(.WIDTH(DATA_W)) u_vec1 (
        .x_in(X1_re), .y_in(X1_im), .z_out(ang1)
    );
    cordic_vectoring #(.WIDTH(DATA_W)) u_vec2 (
        .x_in(X2_re), .y_in(X2_im), .z_out(ang2)
    );

    // ============================================================
    // 2. 計算相位誤差 (Diff from QPSK angles) (Stage 2)
    // ============================================================
    // QPSK 理想角度: 45 (1608), 135 (4824), -135 (-4824), -45 (-1608)
    // 技巧：先減去 45度 (1608)，然後對 90度 (3217) 取餘數
    localparam signed [PHASE_W-1:0] QPSK_45  = 16'sd1608;
    localparam signed [PHASE_W-1:0] QPSK_90  = 16'sd3217;
    localparam signed [PHASE_W-1:0] QPSK_45N = -16'sd1608;

    // 強健版 calc_diff (all-signed constants to avoid unsigned wrap)
    function signed [PHASE_W-1:0] calc_diff;
        input signed [PHASE_W-1:0] theta;
        reg signed [PHASE_W-1:0] tmp;
        begin
            // 1. 減去 45 度 (1608)
            tmp = theta - QPSK_45;
            
            // 2. 使用迴圈將角度強制拉回 +/- 45 度 (1608) 範圍內
            // 90度 = 3217
            // 我們允許誤差稍微大一點點，但不能超過 90 度太多
            
            // 處理正向溢位
            if (tmp > QPSK_45)  tmp = tmp - QPSK_90;
            if (tmp > QPSK_45)  tmp = tmp - QPSK_90;
            if (tmp > QPSK_45)  tmp = tmp - QPSK_90;
            
            // 處理負向溢位
            if (tmp < QPSK_45N) tmp = tmp + QPSK_90;
            if (tmp < QPSK_45N) tmp = tmp + QPSK_90;
            if (tmp < QPSK_45N) tmp = tmp + QPSK_90;
            
            calc_diff = tmp;
        end
    endfunction

    wire signed [PHASE_W-1:0] diff1 = calc_diff(ang1);
    wire signed [PHASE_W-1:0] diff2 = calc_diff(ang2);
    wire signed [PHASE_W-1:0] current_total_diff = diff1 + diff2;

    // ============================================================
    // 3. 累加器與 Ping-Pong Buffer 控制
    // ============================================================
    reg signed [31:0] sum_err;
    reg signed [31:0] sum_next; // 暫存計算用，避免漏掉最後一筆
    reg [5:0]         cnt;       // 0~63 (write counter)
    reg               bank;      // 0 or 1 (write bank)

    // Read-side controls (decoupled from in_valid so one symbol can be output without a second symbol)
    reg               read_active;
    reg [5:0]         read_cnt;
    reg               read_bank;
    wire              start_read;
    
    // 雙緩衝 RAM (128 depth, 64-bit width)
    // Bank 0: addr 0~63, Bank 1: addr 64~127
    // Encourage BRAM inference for the ping-pong buffers
    (* ram_style = "block" *) reg signed [DATA_W-1:0] ram_x1_re [0:127];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] ram_x1_im [0:127];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] ram_x2_re [0:127];
    (* ram_style = "block" *) reg signed [DATA_W-1:0] ram_x2_im [0:127];

    // Dedicated RAM write port (reset-free for RAM inference)
    always @(posedge clk) begin
        if (rst_n && in_valid) begin
            ram_x1_re[{bank, cnt}] <= X1_re;
            ram_x1_im[{bank, cnt}] <= X1_im;
            ram_x2_re[{bank, cnt}] <= X2_re;
            ram_x2_im[{bank, cnt}] <= X2_im;
        end
    end

    reg signed [PHASE_W-1:0] avg_err_reg; // 算好的平均誤差

    // Write & phase accumulation
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            cnt <= 0;
            bank <= 0;
            sum_err <= 0;
            phase_err_valid <= 0;
            phase_err <= 0;
            avg_err_reg <= 0;
        end else if (in_valid) begin
            // 3.2 累加誤差
            if (cnt == 0) begin
                sum_err <= current_total_diff; // 第一筆直接存
                cnt <= cnt + 1;               // 前進到下一筆
                phase_err_valid <= 0;
            end else if (cnt == 63) begin
                // 3.3 結算 Symbol (包含最後一筆誤差)
                // 先計算加上最後一筆的總和，再用來產生平均值
                // 並保證 sum_err 也更新為完整總和 (雖然下一輪會覆寫)
                sum_next = sum_err + current_total_diff;
                sum_err <= sum_next;

                cnt <= 0;
                // 完成的 Bank = bank；下一輪寫 ~bank
                bank <= ~bank; // 切換 Bank
 
                avg_err_reg <= sum_next >>> 6; // 除以 64 (one symbol = 64 tones)
                phase_err <= sum_next >>> 6;
                phase_err_valid <= 1;

            end else begin
                sum_err <= sum_err + current_total_diff;
                cnt <= cnt + 1;
                phase_err_valid <= 0;
            end
        end else begin
            phase_err_valid <= 0;
        end
    end

    // ============================================================
    // 4. 讀出修正 (Correction Stage)
    // ============================================================
    // 當正在寫入 Bank 1 時，我們讀取 Bank 0，反之亦然
    // 使用另一個 counter 來讀取，或者利用 data pipeline 特性
    
    reg signed [DATA_W-1:0] r_x1_re, r_x1_im;
    reg signed [DATA_W-1:0] r_x2_re, r_x2_im;
    reg                     valid_d1;

    assign start_read = in_valid && (cnt == 6'd63);

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            read_active <= 0;
            read_cnt    <= 0;
            read_bank   <= 0;
            r_x1_re <= 0; r_x1_im <= 0;
            r_x2_re <= 0; r_x2_im <= 0;
            valid_d1 <= 0;
        end else begin
            valid_d1 <= 0; // default low

            if (start_read) begin
                read_active <= 1;
                read_bank   <= bank;
                read_cnt    <= 0;
            end else if (read_active) begin
                r_x1_re <= ram_x1_re[{read_bank, read_cnt}];
                r_x1_im <= ram_x1_im[{read_bank, read_cnt}];
                r_x2_re <= ram_x2_re[{read_bank, read_cnt}];
                r_x2_im <= ram_x2_im[{read_bank, read_cnt}];
                valid_d1 <= 1;

                if (read_cnt == 63) begin
                    read_cnt    <= 0;
                    read_active <= 0;
                end else begin
                    read_cnt <= read_cnt + 1'b1;
                end
            end
        end
    end

    // CORDIC Rotation (修正)
    // 修正角度 = -avg_err
    wire signed [PHASE_W-1:0] rot_angle = -avg_err_reg;

    wire signed [DATA_W-1:0] x1_fixed_re, x1_fixed_im;
    wire signed [DATA_W-1:0] x2_fixed_re, x2_fixed_im;

    cordic_rotator #(.WIDTH(DATA_W)) u_rot1 (
        .clk(clk), .rst_n(rst_n),
        .x_in(r_x1_re), .y_in(r_x1_im), .theta_in(rot_angle),
        .x_out(x1_fixed_re), .y_out(x1_fixed_im)
    );

    cordic_rotator #(.WIDTH(DATA_W)) u_rot2 (
        .clk(clk), .rst_n(rst_n),
        .x_in(r_x2_re), .y_in(r_x2_im), .theta_in(rot_angle),
        .x_out(x2_fixed_re), .y_out(x2_fixed_im)
    );

    // 輸出暫存
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 0;
            X1c_re <= 0; X1c_im <= 0;
            X2c_re <= 0; X2c_im <= 0;
        end else begin
            out_valid <= valid_d1;
            X1c_re <= x1_fixed_re;
            X1c_im <= x1_fixed_im;
            X2c_re <= x2_fixed_re;
            X2c_im <= x2_fixed_im;
        end
    end

    // Debug: 印出算出來的平均誤差
    always @(posedge clk) begin
        if (cnt == 63 && in_valid) begin
            $display("[DEBUG] Symbol Done. Sum_Err = %d, Avg_Err (Shift 6) = %d", sum_err, sum_err >>> 6);
        end
    end
    // 在 cpe_tracker.v 內
    always @(posedge clk) begin
        if (cnt == 0 && in_valid) begin
            $display("[DEBUG] ang1 = %d, diff1 = %d", ang1, diff1);
        end
    end

endmodule