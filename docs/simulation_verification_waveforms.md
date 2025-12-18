# Simulation & Verification Results（RTL 波形驗證 / GTKWave）

本章節提供一個「可直接貼進期末報告」的 RTL 驗證寫法，並說明如何用 GTKWave 產生波形證據。

---

## 1. 波形產生方式（VCD → GTKWave）

### 1.1 使用 iverilog 產生 VCD（Windows / PowerShell）

以下指令假設你已安裝 `iverilog` 並已加入 PATH，且在專案根目錄執行（同一層能看到 `tb_top_receiver.v` 與各種 `golden_*.hex`）。

檢查版本：

```powershell
iverilog -V
```

（A）Top-level 端到端波形：

```powershell
iverilog -g2012 -o sim_tb_top_receiver.exe tb_top_receiver.v Top_level_receiver.v sync_block.v CFO_estimator_cp_mrc.v nco_phase_ctrl.v cfo_rotator.v cordic_rotator.v cordic_vectoring.v fft_64pt.v mimo_zf_2x2.v CPE_tracker.v QPSK_demap_mimo.v c_mul_fx.v
vvp sim_tb_top_receiver.exe
gtkwave tb_top_receiver.vcd
```

（B）ZF 子模組波形：

```powershell
iverilog -g2012 -o sim_tb_mimo_zf.exe tb_mimo_zf.v mimo_zf_2x2.v
vvp sim_tb_mimo_zf.exe
gtkwave tb_mimo_zf.vcd
```

（C）CPE 子模組波形：

```powershell
iverilog -g2012 -o sim_tb_cpe_tracker.exe tb_CPE_tracker.v CPE_tracker.v
vvp sim_tb_cpe_tracker.exe
gtkwave tb_cpe_tracker.vcd
```

> 如果你跑 `vvp` 時出現 `$readmemh` 找不到 `golden_*.hex`，代表工作目錄不對：請確認你是在專案根目錄跑，或把 hex 檔放到同一層，或改成用相對路徑（例如 `vectors/golden_*.hex`）。

本專案 testbench 已加入 `$dumpfile/$dumpvars`，執行模擬後會在專案目錄產生下列波形檔：

- `tb_top_receiver.vcd`：端到端 Top-level（sync → CFO/CPE → FFT → ZF → CPE → demap）
- `tb_mimo_zf.vcd`：2x2 MIMO ZF 子模組
- `tb_cpe_tracker.vcd`：CPE tracker 子模組

使用 GTKWave：

1. 開啟 `gtkwave`
2. `File → Open New Tab...` 選擇 `*.vcd`
3. 在左側階層找到 `tb_*`，把下列訊號加入波形視窗
4. 對複數/定點訊號建議設定顯示：`Data Format → Signed Decimal`（或 Hex）

> 建議：VCD 檔不要 push 上 GitHub（檔案會很大）；用 README/報告放「截圖」即可。

---

## 2. Top-level 驗證（tb_top_receiver）— 建議截圖與解釋

### 2.0（快速版）tb_top_receiver.vcd 最短截圖清單（建議 12 條）

如果你只想「最快截圖」且能覆蓋整條鏈路（Sync→CFO→FFT→ZF→CPE→Demap），打開 `tb_top_receiver.vcd` 後，在 GTKWave 將以下訊號加入波形視窗即可（建議顯示格式：Signed Decimal；valid/bit 類維持 Binary）。

**A. TB 層（stimulus/基本控制）**

- `clk`
- `rst_n`
- `in_valid`

**B. u_top 內（pipeline 里程碑）**

- `sync_symbol_start`
- `sync_valid`
- `cfo_est_valid`
- `cfo_eps_hat`
- `fft_load_valid`
- `fft_data_valid_q`
- `zf_valid_dbg`（或 `zf_valid_tap`）
- `cpe_phase_valid_dbg`（或 `cpe_phase_valid_tap`）
- `out_valid`
- `rx_bits[3:0]`

> 找訊號路徑提示：`tb_top_receiver` → `u_top`（module: `mimo_ofdm_rx_top`）。

**（可選加分，若你想多一張“數值證據”）**

- ZF 複數輸出：`zf_x1_re_dbg`, `zf_x1_im_dbg`, `zf_x2_re_dbg`, `zf_x2_im_dbg`
- CPE/NCO：`cpe_phase_err_dbg`, `nco_acc_dbg`

### 圖 1：Reset 與輸入資料灌入（Stimulus 正確）

**建議加入訊號**

- `clk`, `rst_n`, `in_valid`
- `adc1_re`, `adc1_im`, `adc2_re`, `adc2_im`

**你要在報告說什麼（可直接貼）**

> 本測試在 `rst_n` 解除後開始灌入 `in_valid`，並依序送入一個 OFDM symbol 的時域取樣（包含 CP+data）。波形中可觀察到 `in_valid` 於輸入區間連續為 1，且 ADC I/Q 數值隨測試向量變化，表示 stimulus 餵入流程正確。

### 圖 2：ZF 輸出 valid 與複數輸出（子鏈路正確性）

**建議加入訊號**

- `zf_valid_dbg`（或 `zf_valid_tap`）
- `zf_x1_re_dbg`, `zf_x1_im_dbg`, `zf_x2_re_dbg`, `zf_x2_im_dbg`

**你要在報告說什麼**

> 當 `zf_valid_dbg` 於每個 subcarrier 對應時刻拉高，ZF 模組輸出複數符號 `X1/X2`。此 valid 訊號呈現連續 64 筆（對應 NFFT=64），表示頻域處理鏈已輸出完整一個 OFDM symbol。

### 圖 3：CPE phase error / NCO accumulator（CPE 行為合理）

**建議加入訊號**

- `cpe_phase_valid_dbg`（或 `cpe_phase_valid_tap`）
- `cpe_phase_err_dbg`
- `nco_acc_dbg`

**你要在報告說什麼**

> CPE tracker 會根據解調符號計算 phase error，並透過 NCO accumulator 更新相位補償量。波形中可觀察 `cpe_phase_valid_dbg` 對齊輸出節拍，`cpe_phase_err_dbg` 與 `nco_acc_dbg` 產生隨時間更新的數值，符合「估測誤差 → 更新補償」的閉迴路行為。

### 圖 4：最終輸出 out_valid / rx_bits（端到端輸出證據）

**建議加入訊號**

- `out_valid`
- `rx_bits[3:0]`

**你要在報告說什麼**

> 當 `out_valid` 拉高時，代表接收端完成一筆輸出（2 streams 的 QPSK demap 結果，合併為 4 bits）。波形可觀察到 `out_valid` 於 symbol 輸出區間連續出現 64 筆，表示從輸入到輸出之 pipeline 對齊正常。

###（加分）圖 5：Sync（CP-based）對齊證據（symbol_start / sync_valid）

這張圖可以證明「你真的有做 OFDM 對位」，通常很加分。

**在哪裡找訊號（GTKWave 層級提示）**

- `tb_top_receiver` → `u_top` →（module: `mimo_ofdm_rx_top`）

**建議加入訊號**

- `in_valid`（來自 tb）
- `sync_symbol_start`
- `sync_valid`
- （可選）`sym_sample_idx`、`collecting`、`playing`

**你要在報告說什麼**

> Sync block 會利用 CP correlation 自動偵測一個 OFDM symbol 的起點。波形中可觀察到 `sync_symbol_start` 於正確時刻產生脈衝，接著 `sync_valid` 連續為 1（對應切出 CP+data 或 data 區間），表示 symbol 對位成功，後續 CFO/FFT 鏈路使用此對齊後資料。

###（加分）圖 6：CFO estimation 生效證據（cfo_valid / eps_hat）

這張圖可以把「估測 → 套用補償」的流程講清楚。

**在哪裡找訊號（GTKWave 層級提示）**

- `tb_top_receiver` → `u_top`

**建議加入訊號**

- `cfo_est_valid`
- `cfo_eps_hat`
- （可選）`cfo_ready`、`playing`（表示 buffer replay 開始）

**你要在報告說什麼**

> CFO estimator 在偵測到 symbol start 後進行 CFO 估測，並於 `cfo_est_valid` 拉高時輸出 `cfo_eps_hat`。波形可看到 `cfo_eps_hat` 在 valid 時刻更新並維持穩定，代表 CFO 估測完成；後續資料由 buffer replay 餵入 rotator/FFT 以套用相位補償。

###（加分）圖 7：Rotator → FFT 的資料流對齊（rot_valid / fft_load_valid / fft_done / fft_data_valid_q）

這張圖能證明你掌握 pipeline 對齊與 FFT frame 邊界。

**在哪裡找訊號（GTKWave 層級提示）**

- `tb_top_receiver` → `u_top`

**建議加入訊號**

- `rot_valid`
- `fft_load_valid`（CP removed 後餵 FFT 的 valid）
- `fft1_done`, `fft2_done`
- `fft_data_valid_q`
- `fft_read_addr_q`（或 `fft_read_addr`）
- `fft1_re_q`, `fft1_im_q`, `fft2_re_q`, `fft2_im_q`

**你要在報告說什麼**

> 經 CFO rotator 修正後，系統以 `fft_load_valid` 將去除 CP 的 64 筆資料送入兩顆 FFT（RX1/RX2）。當 `fft1_done` 與 `fft2_done` 皆完成後，讀出控制器開始以 `fft_read_addr(_q)` 依序讀出 64 個 subcarriers，並以 `fft_data_valid_q` 標示頻域輸出有效。此 valid 也同時作為 ZF 的輸入節拍，確保頻域處理鏈路對齊。

---

## 3. ZF 子模組驗證（tb_mimo_zf）— 建議截圖與解釋

**建議加入訊號**

- `in_valid`, `Y1_re/Y1_im`, `Y2_re/Y2_im`
- `out_valid`, `X1_re/X1_im`, `X2_re/X2_im`

**報告文字（可直接貼）**

> 本測試以 `golden_fft_out.hex` 作為 ZF 輸入（頻域 Y），並以 `golden_h_inv.hex` 提供固定通道反矩陣。當 `out_valid` 拉高時，ZF 輸出 `X1/X2` 與 `golden_zf_out.hex` 進行容許誤差（TOL=2 LSB）比對。波形中可見 `in_valid` 與 `out_valid` 的節拍一致，且輸出複數符號於 valid 週期穩定。

---

## 4. CPE 子模組驗證（tb_cpe_tracker）— 建議截圖與解釋

**建議加入訊號**

- `in_valid`, `X1_re/X1_im`, `X2_re/X2_im`
- `phase_err_valid`, `phase_err`
- `out_valid`, `X1c_re/X1c_im`, `X2c_re/X2c_im`

**報告文字（可直接貼）**

> CPE tracker 使用 ping-pong buffer，因此需額外餵入一輪 dummy valid 以推擠資料，讓第一輪資料被讀出並完成修正輸出。波形中可觀察到 `in_valid` 連續兩段（第一段為真實資料、第二段為 dummy），而 `out_valid` 於延遲後輸出第一輪對應的 64 筆修正結果。輸出 `X1c/X2c` 與 `golden_final_out.hex` 以 TOL=2 LSB 進行比對。

---

## 5. 截圖建議（GTKWave）

為了讓老師一眼看懂，建議每張截圖都做到：

- 畫面左上角包含檔名（`tb_*.vcd`）或階層（`tb_*`），證明是你的 RTL 波形
- 用 marker 標示：`rst_n` 解除、`in_valid` 開始、第一個 `out_valid` 出現的位置
- 定點訊號顯示成 Signed Decimal（或 Hex，但要在文字說明位寬/符號）

建議輸出 3～6 張截圖即可：

- Top-level：Stimulus、ZF valid+輸出、CPE error/NCO、最終 out_valid/rx_bits
- 子模組：ZF、CPE 各 1 張（可選，但很加分）
