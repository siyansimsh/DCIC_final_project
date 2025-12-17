# DCIC Final Project — 2x2 MIMO OFDM Receiver (FPGA RTL)

本專案為 DCIC 期末專題的 FPGA/RTL 實作與驗證素材，目標是完成一個 2x2 MIMO OFDM 接收端的主要處理鏈（同步、CFO/CPE 校正、FFT、MIMO ZF、Demap），並可在 Vivado 完整跑完 Synthesis/Implementation 與產生 timing / utilization 報告。

---

## 期末摘要（報告/影片口白）

本專案完成 2x2 MIMO OFDM Receiver 的 RTL 版本，整體流程可在 Vivado 完成 Synthesis 與 Implementation，並輸出設計資源使用與時序結果。由於未取得真實板端 pinout 與外部 I/O 時序參數（tCO/board skew 等），本專案以 placeholder XDC 進行 STA，並在報告中清楚標註時序假設與限制。

### Timing（Post-Implementation）

- Clock：`sys_clk`（period 20.000 ns / 50 MHz；placeholder）
- Setup：WNS **+0.709 ns**、TNS **0.000 ns**、Failing Endpoints **0**
- Hold：WHS **+0.012 ns**、THS **0.000 ns**、Failing Endpoints **0**

> 註：在缺乏真實 I/O 時序時，primary input 的 hold 檢查容易因 placeholder `set_input_delay -min 0` 造成大量不具代表性的 hold violations；因此本專案採用「保留 setup 檢查、但對 primary inputs 關掉 hold 檢查」的方式，讓 STA 結果反映設計內部可控的時序收斂狀態。若要上板/正式 signoff，需改回真實 I/O constraints。

### Utilization（Post-Implementation）

| Resource | Used | Available | Utilization |
|---|---:|---:|---:|
| LUT | 14306 | 230400 | 6.21% |
| LUTRAM | 160 | 101760 | 0.16% |
| FF | 17022 | 460800 | 3.69% |
| BRAM | 2 | 312 | 0.64% |
| DSP | 60 | 1728 | 3.47% |
| IO | 298 | 464 | 64.22% |
| BUFG | 1 | 544 | 0.18% |

### 報告用圖

- 系統方塊圖/報告用：[方塊圖(報告用).jpg](%E6%96%B9%E5%A1%8A%E5%9C%96(%E5%A0%B1%E5%91%8A%E7%94%A8).jpg)
- 硬體模組規劃：[硬體模組規劃.jpg](%E7%A1%AC%E9%AB%94%E6%A8%A1%E7%B5%84%E8%A6%8F%E5%8A%83.jpg)
- Verilog 模組列表：[verilog設計模組.jpg](verilog%E8%A8%AD%E8%A8%88%E6%A8%A1%E7%B5%84.jpg)
- 驗證平台/黃金向量：[驗證平台和黃金向量.jpg](%E9%A9%97%E8%AD%89%E5%B9%B3%E5%8F%B0%E5%92%8C%E9%BB%83%E9%87%91%E5%90%91%E9%87%8F.jpg)

### 建議補上的 Vivado 截圖（放進 docs/images 後連結就會生效）

- Timing Summary（Post-Implementation）：[docs/images/timing_summary.png](docs/images/timing_summary.png)
- Utilization（Post-Implementation）：[docs/images/utilization.png](docs/images/utilization.png)
- Messages/DRC 摘要（可選）：[docs/images/messages_drc.png](docs/images/messages_drc.png)

> 注意：本 repo 目前使用 **placeholder XDC** 進行靜態時序分析（STA）。若要上板/做真實 signoff，必須改成真實的板端 pinout（LOC）、IOSTANDARD 與 I/O delay 假設。

---

## 專案內容

### 主要 RTL 模組（Verilog）

- `Top_level_receiver.v`：接收端頂層串接（sync / CFO / FFT / ZF / CPE / demap 等）
- `sync_block.v`：同步/緩衝
- `cfo_estimator_cp_mrc.v`：CFO 估測
- `nco_phase_ctrl.v`：NCO 相位/累加控制
- `cfo_rotator.v`、`cordic_rotator.v`、`cordic_vectoring.v`：旋轉/向量化（CORDIC）
- `fft_64pt.v`：64-point FFT
- `mimo_zf_2x2.v`：2x2 MIMO Zero-Forcing
- `CPE_tracker.v`：CPE phase tracker / correction
- `QPSK_demap_mimo.v`：QPSK demapper（MIMO）

### Testbench / 波形

- `tb_top_receiver.v`：整體接收端測試平台（如有使用）
- `tb_CPE_tracker.v`、`tb_mimo_zf.v`：子模組測試
- `top_level.vcd`、`wave*/`：波形輸出（建議不推到 GitHub）

### C 參考/模型（用於驗證/比對）

- `ofdm_final_qpsk.c`、`fixed_point_main.c`：參考流程/固定點相關程式（依實際使用情況）

### Vivado 專案

- `DCIC_final_ofdm/`：Vivado project 相關資料（含 XDC 等）
- `2x2_mimo_cfo_ofdm/`：可能的另一組 Vivado 工作資料（視你使用狀況）

---

## 如何跑（Vivado）

1. 開啟 Vivado
2. 開啟專案：`DCIC_final_ofdm/DCIC_final_ofdm.xpr`
3. 確認 constraints 已加入：
   - 若你使用 placeholder constraint：`temp_placeholder.xdc`
   - 或使用專案內的 `time.xdc`（依你的 run 設定）
4. 依序執行：
   - Run Synthesis
   - Run Implementation
5. 產出報告：
   - Timing Summary（Setup/Hold/Pulse Width）
   - Utilization
   - DRC Report

---

## Placeholder XDC 說明（重要）

檔案：`temp_placeholder.xdc`

- `create_clock` 目前使用較寬鬆的 clock period（例如 20ns / 50MHz）以便在不大改 RTL 的前提下完成 implementation。
- 因缺乏真實板端 I/O 時序（外部 tCO、板級延遲、skew），對 primary inputs 採用 placeholder `set_input_delay`。
- 若遇到大量 **input→同步 FF 的 HOLD violations**（常見於 placeholder `-min 0` 的情境），本專案採用：
  - **保留 setup 檢查，但對 primary inputs 關掉 hold 檢查**（避免不具實際意義的爆量 hold fail）

> 若要上板/正式 signoff：請移除 `set_false_path -hold ...` 並改用真實的 `set_input_delay -min/-max`（依 ADC/外部介面規格與板級延遲）

---

## 版本管理（Git）

本 repo 提供 `.gitignore`，用於避免提交 Vivado 大量產物（`.runs/.sim/.cache/.Xil`）、波形檔（`*.vcd`）與本機編譯輸出（`*.exe`）。

若你曾經把波形/產物加入追蹤，需先從 Git index 移除（不刪本機檔）：

```bash
git rm --cached top_level.vcd
git rm --cached -r simv
```

---

## 已知限制

- 未提供真實板端 pinout / I/O timing：目前以 placeholder constraints 產生 STA 報告。
- 若需更高 Fmax 或更嚴格時序：可能需要針對關鍵路徑做 pipeline/架構調整（本期末版本以流程可交付為優先）。

---

## 作者/課程

- DCIC Final Project
- 2x2 MIMO OFDM Receiver (RTL + Vivado)
