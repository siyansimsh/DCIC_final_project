// ------------------------------------------------------------
//  Final Project: 2x2 MIMO-OFDM Receiver Simulation (Float)
//  Features: 
//    1. 2x2 MIMO Channel with Common CFO [cite: 9]
//    2. MRC-based CFO Estimation [cite: 11]
//    3. Time-domain Compensation [cite: 12]
//    4. Zero-Forcing (ZF) Detection 
// ------------------------------------------------------------

// ------------------------------------------------------------
//  Final Project: 2x2 MIMO-OFDM Receiver Simulation (Fixed Version)
// ------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808
#endif

// ====== 系統參數 ======
#define NUM_TX 2
#define NUM_RX 2
#define NFFT 64u
#define NCP  16u       
#define SEARCH_WINDOW_LEN (NCP + NFFT + NCP) // prev-CP + (CP+data)
#define BITS_PER_SYM 2u // QPSK
#define TOTAL_BITS   (1u<<20) 
#define SYNC_SMOOTH_ALPHA 0.0625 

// ====== 複數運算與結構 ======
typedef struct { double re, im; } cpx; // Complex number
static inline cpx C(double r,double i){ return (cpx){r,i}; } // Construct complex function
static inline cpx c_add(cpx a,cpx b){ return C(a.re+b.re, a.im+b.im);}  
static inline cpx c_sub(cpx a,cpx b){ return C(a.re-b.re, a.im-b.im);}  
static inline cpx c_mul(cpx a,cpx b){ return C(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);}
static inline cpx c_div(cpx a,cpx b){ 
    double d = b.re*b.re + b.im*b.im; 
    return C((a.re*b.re + a.im*b.im)/d, (a.im*b.re - a.re*b.im)/d);
}
static inline cpx c_conj(cpx a){ return C(a.re, -a.im);}  
static inline cpx c_scale(cpx a,double s){ return C(a.re*s, a.im*s);} 
static inline double c_abs_sq(cpx a){ return a.re*a.re + a.im*a.im; }

// ====== 輔助函數 (RNG, FFT) ======
static inline uint32_t prng_u32(void){ static uint32_t s=0x12345678; s = 1664525u*s + 1013904223u; return s; }
static inline double rand_uniform(void){ return (prng_u32()>>8) * (1.0/16777216.0); }
static inline double randn(void){ double u1=rand_uniform()+1e-12, u2=rand_uniform(); return sqrt(-2*log(u1)) * cos(2*M_PI*u2); }

// FFT
static unsigned ilog2u(size_t n){ unsigned p=0; while((1ULL<<p)<n) ++p; return p; }
static size_t bit_reverse(size_t x, unsigned bits){
    size_t r=0; for(unsigned i=0;i<bits;++i){ r=(r<<1)|(x&1); x>>=1; } return r;
}
void fft_inplace(cpx *a, size_t N, int inverse){
    unsigned bits = ilog2u(N);
    for(size_t i=0;i<N;++i){ size_t j=bit_reverse(i,bits); if(j>i){ cpx t=a[i]; a[i]=a[j]; a[j]=t; } }
    for(size_t m=2; m<=N; m<<=1){ 
        double theta = (inverse?2.0:-2.0)*M_PI/m; 
        cpx wp = C(cos(theta), sin(theta));
        for(size_t k=0; k<N; k+=m){ 
            cpx w = C(1.0, 0.0);
            for(size_t j=0; j<m/2; ++j){
                cpx t = c_mul(w, a[k+j+m/2]);
                cpx u = a[k+j];
                a[k+j] = c_add(u, t);
                a[k+j+m/2] = c_sub(u, t);
                w = c_mul(w, wp);
            }
        }
    }
    if(inverse){ for(size_t i=0;i<N;++i) a[i]=c_scale(a[i], 1.0/N); }
}

// ====== 模組 1: QPSK ======
static inline cpx qpsk_map(uint8_t b0,uint8_t b1){
    double s = 1.0/M_SQRT2; 
    return C((b0? -s:s), (b1? -s:s)); 
}
static inline void qpsk_demap(cpx s, uint8_t *b0, uint8_t *b1){
    *b0 = (s.re < 0) ? 1 : 0;
    *b1 = (s.im < 0) ? 1 : 0;
}

// ====== 模組 2: Channel + CFO ======
static void apply_mimo_channel_cfo(
    const cpx *tx1, const cpx *tx2, 
    cpx *rx1, cpx *rx2, 
    size_t len, 
    double SNRdB, 
    double cfo_norm,
    cpx H[2][2], 
    size_t time_offset
) {

    // [FIX 1] 訊號功率校正：IFFT 後時域訊號功率約為 1/NFFT
    double signal_power = 1.0 / (double)NFFT; 
    double noise_power = signal_power * pow(10.0, -SNRdB/10.0);
    double sigma = sqrt(noise_power/2.0);

    for(size_t i=0; i<len; ++i){
        cpx x1 = tx1[i]; cpx x2 = tx2[i];
        
        // MIMO Mixing
        cpx y1 = c_add(c_mul(H[0][0], x1), c_mul(H[0][1], x2));
        cpx y2 = c_add(c_mul(H[1][0], x1), c_mul(H[1][1], x2));

        // CFO Rotation
        double total_time_idx = (double)(time_offset + i);
        double phase = 2.0 * M_PI * (cfo_norm / (double)NFFT) * total_time_idx;
        cpx rot = C(cos(phase), sin(phase));
        y1 = c_mul(y1, rot);
        y2 = c_mul(y2, rot);

        // AWGN
        y1.re += sigma*randn(); y1.im += sigma*randn();
        y2.re += sigma*randn(); y2.im += sigma*randn();

        rx1[i] = y1; rx2[i] = y2;
    }
}

// ====== 模組 3: Sync & CFO Estimation (MRC) ======
static size_t sync_symbol_robust_mimo(
    const cpx *r1_win,
    const cpx *r2_win,
    size_t win_len,
    double *profile_ema
) {
    double best_metric = -1.0;
    size_t best_idx = 0;
    size_t m = 0;

    for (; m + NFFT + NCP <= win_len; ++m) {
        cpx corr = C(0,0);
        for (size_t i = 0; i < NCP; ++i) {
            size_t idx_cp   = m + i;
            size_t idx_tail = m + i + NFFT;
            corr = c_add(corr, c_mul(r1_win[idx_cp], c_conj(r1_win[idx_tail])));
            corr = c_add(corr, c_mul(r2_win[idx_cp], c_conj(r2_win[idx_tail])));
        }

        double metric = corr.re*corr.re + corr.im*corr.im;
        double old_avg = profile_ema[m];
        double new_avg = old_avg + (metric - old_avg) * SYNC_SMOOTH_ALPHA;
        profile_ema[m] = new_avg;

        if (new_avg > best_metric) {
            best_metric = new_avg;
            best_idx = m;
        }
    }

    // decay the tail
    for (; m < win_len; ++m) {
        profile_ema[m] = profile_ema[m] * (1.0 - SYNC_SMOOTH_ALPHA);
    }

    return best_idx;
}

static double estimate_cfo_at_idx(
    const cpx *r1_win,
    const cpx *r2_win,
    size_t start_idx
) {
    cpx corr_total = C(0,0);
    for (size_t i = 0; i < NCP; ++i) {
        size_t idx_cp   = start_idx + i;
        size_t idx_tail = start_idx + i + NFFT;
        corr_total = c_add(corr_total, c_mul(r1_win[idx_cp], c_conj(r1_win[idx_tail])));
        corr_total = c_add(corr_total, c_mul(r2_win[idx_cp], c_conj(r2_win[idx_tail])));
    }
    double angle   = atan2(corr_total.im, corr_total.re);
    double eps_hat = -angle / (2.0 * M_PI);
    if (eps_hat >= 0.5) eps_hat -= 1.0;
    if (eps_hat < -0.5) eps_hat += 1.0;
    return eps_hat;
}




// ====== 模組 4: CFO Compensation ======
// [FIX 3] 支援標準化 CFO 輸入與相位連續性
static void apply_cfo_compensation_continuous(
    cpx *r1, cpx *r2, 
    size_t len, 
    double cfo_est_norm, 
    double *phase_track 
){
    // 每 sample 旋轉量: 2 * pi * epsilon / N
    double phase_inc = (2.0 * M_PI * cfo_est_norm) / (double)NFFT;

    for(size_t i=0; i<len; ++i){
        double phi = *phase_track;
        cpx rot = C(cos(phi), sin(phi));
        
        r1[i] = c_mul(r1[i], rot);
        r2[i] = c_mul(r2[i], rot);

        // 接收端頻率偏高 -> 相位增加 -> 補償要減
        *phase_track -= phase_inc;
        
        if(*phase_track > M_PI) *phase_track -= 2*M_PI;
        if(*phase_track < -M_PI) *phase_track += 2*M_PI;
    }
}

// ====== 模組 5: MIMO ZF ======
static void mimo_zf_detect(cpx y1, cpx y2, cpx H[2][2], cpx *x1_est, cpx *x2_est){
    cpx det = c_sub(c_mul(H[0][0], H[1][1]), c_mul(H[0][1], H[1][0]));
    
    cpx h_inv[2][2];
    h_inv[0][0] = c_div(H[1][1], det);
    h_inv[0][1] = c_div(c_scale(H[0][1], -1.0), det);
    h_inv[1][0] = c_div(c_scale(H[1][0], -1.0), det);
    h_inv[1][1] = c_div(H[0][0], det);

    *x1_est = c_add(c_mul(h_inv[0][0], y1), c_mul(h_inv[0][1], y2));
    *x2_est = c_add(c_mul(h_inv[1][0], y1), c_mul(h_inv[1][1], y2));
}

// CPE Correction
// 修改回傳值為 double
static double correct_residual_phase(cpx *x1, cpx *x2, size_t n){
    double total_phase_err = 0;
    for(size_t k=0; k<n; ++k){
         // Stream 1
         double ang = atan2(x1[k].im, x1[k].re); 
         double diff = ang - (M_PI/4.0);
         while(diff < 0) diff += (M_PI/2.0);
         diff = fmod(diff, M_PI/2.0); 
         if(diff > M_PI/4.0) diff -= (M_PI/2.0); 
         total_phase_err += diff;

         // Stream 2
         ang = atan2(x2[k].im, x2[k].re);
         diff = ang - (M_PI/4.0);
         while(diff < 0) diff += (M_PI/2.0);
         diff = fmod(diff, M_PI/2.0);
         if(diff > M_PI/4.0) diff -= (M_PI/2.0);
         total_phase_err += diff;
    }
    // 計算平均相位誤差
    double avg_err = total_phase_err / (2.0 * n);
    
    // 進行旋轉補償
    cpx rot = C(cos(-avg_err), sin(-avg_err));
    for(size_t k=0; k<n; ++k){
        x1[k] = c_mul(x1[k], rot);
        x2[k] = c_mul(x2[k], rot);
    }

    // [關鍵] 回傳剛剛計算出的誤差
    return avg_err; 
}

// =============================================================
//                       MAIN PROGRAM
// =============================================================
int main(){
    printf("=== 2x2 MIMO-OFDM Receiver Simulation ===\n");
    printf("NFFT: %u, CP: %u, Modulation: QPSK\n", NFFT, NCP);

    const size_t sym_len = NFFT + NCP;
    cpx *tx1_time = malloc(sym_len * sizeof(cpx));
    cpx *tx2_time = malloc(sym_len * sizeof(cpx));
    
    cpx *rx1_prev   = calloc(sym_len, sizeof(cpx));
    cpx *rx2_prev   = calloc(sym_len, sizeof(cpx));
    cpx *rx1_curr   = malloc(sym_len * sizeof(cpx));
    cpx *rx2_curr   = malloc(sym_len * sizeof(cpx));
    cpx *rx1_win    = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    cpx *rx2_win    = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    double *sync_profile = calloc(SEARCH_WINDOW_LEN, sizeof(double));

    // ====== [插入 1] 監測變數宣告 ======
    double max_abs_adc = 0.0;   // 監測 ADC 輸入 (RX Stream)
    double max_abs_fft = 0.0;   // 監測 FFT 輸出
    double max_abs_zf = 0.0;    // 監測 ZF 解調後 (Equalized)
    double max_abs_det = 0.0;   // 監測 ZF 行列式值 (決定除法位寬)

    uint8_t *bits_tx = malloc(TOTAL_BITS);

    for(size_t i=0; i<TOTAL_BITS; ++i) bits_tx[i] = (prng_u32() >> 16) & 1;

    double cfo_actual = 0.03; // Normalized CFO (Δf / Δf_subcarrier), 可設為0 
    double snr_list[] = {0, 5, 10, 15, 20, 25};
    
    printf("Testing with CFO = %.3f (Normalized)\n", cfo_actual);
    printf("SNR(dB) | BER (MIMO)\n");
    printf("--------------------\n");

    // ====== 生成 Rayleigh Fading Channel ======
        // 放在這裡代表：針對這一個 SNR 測試點，通道是固定的 (Block Fading)
        cpx H_ref[2][2]; 

        double scale = 1.0 / sqrt(2.0); // 能量正規化因子
        H_ref[0][0] = C(randn() * scale, randn() * scale);
        H_ref[0][1] = C(randn() * scale, randn() * scale);
        H_ref[1][0] = C(randn() * scale, randn() * scale);
        H_ref[1][1] = C(randn() * scale, randn() * scale);
    

    for(int s=0; s<6; ++s){
        double snr = snr_list[s];
        double rx_nco_phase = 0.0; 
        size_t errors = 0;
        size_t total_bits_processed = 0;
        size_t global_time_counter = 0; 
        
        memset(rx1_prev, 0, sym_len*sizeof(cpx));
        memset(rx2_prev, 0, sym_len*sizeof(cpx));
        memset(sync_profile, 0, SEARCH_WINDOW_LEN*sizeof(double));

        size_t num_syms = (TOTAL_BITS / (BITS_PER_SYM * NUM_TX * NFFT)); 
        size_t bit_idx = 0;

        // [FIX 4] 記憶體配置移至此處 (Scope Fix)
        cpx *x1_arr = malloc(NFFT*sizeof(cpx)); 
        cpx *x2_arr = malloc(NFFT*sizeof(cpx));


        // 確認通道
        if (s == 0) {
            printf("[System] Generated Random Channel:\n");
            printf("H[0][0]=%.3f+j%.3f, H[0][1]=%.3f+j%.3f\n", H_ref[0][0].re, H_ref[0][0].im, H_ref[0][1].re, H_ref[0][1].im);
            printf("H[1][0]=%.3f+j%.3f, H[1][1]=%.3f+j%.3f\n", H_ref[1][0].re, H_ref[1][0].im, H_ref[1][1].re, H_ref[1][1].im);
        }
        
        // ==============================================

        // --- System Warm-up ---
        {
            cpx w_tx1_freq[NFFT], w_tx2_freq[NFFT];
            cpx w_tx1_time[NFFT+NCP], w_tx2_time[NFFT+NCP];
            cpx w_rx1[NFFT+NCP], w_rx2[NFFT+NCP];
            cpx H_ref[2][2];
            
            for(size_t k=0; k<NFFT; ++k){
                w_tx1_freq[k] = qpsk_map(prng_u32()&1, prng_u32()&1);
                w_tx2_freq[k] = qpsk_map(prng_u32()&1, prng_u32()&1);
            }
            cpx t1[NFFT], t2[NFFT];
            memcpy(t1, w_tx1_freq, sizeof(t1)); fft_inplace(t1, NFFT, 1);
            memcpy(t2, w_tx2_freq, sizeof(t2)); fft_inplace(t2, NFFT, 1);
            memcpy(w_tx1_time, t1+NFFT-NCP, NCP*sizeof(cpx)); memcpy(w_tx1_time+NCP, t1, NFFT*sizeof(cpx));
            memcpy(w_tx2_time, t2+NFFT-NCP, NCP*sizeof(cpx)); memcpy(w_tx2_time+NCP, t2, NFFT*sizeof(cpx));
            apply_mimo_channel_cfo(w_tx1_time, w_tx2_time, w_rx1, w_rx2, sym_len, snr, cfo_actual, H_ref, global_time_counter);
            global_time_counter += sym_len;
            memcpy(rx1_prev, w_rx1, sym_len*sizeof(cpx));
            memcpy(rx2_prev, w_rx2, sym_len*sizeof(cpx));
        }

        // --- Symbol Loop ---
        for(size_t sym=0; sym < num_syms; ++sym){
            // 1. Transmitter
            cpx tx1_freq[NFFT], tx2_freq[NFFT];
            for(size_t k=0; k<NFFT; ++k){
                tx1_freq[k] = qpsk_map(bits_tx[bit_idx], bits_tx[bit_idx+1]);
                tx2_freq[k] = qpsk_map(bits_tx[bit_idx+2], bits_tx[bit_idx+3]);
                bit_idx += 4;
            }
            cpx t1[NFFT], t2[NFFT];
            memcpy(t1, tx1_freq, sizeof(t1)); fft_inplace(t1, NFFT, 1);
            memcpy(t2, tx2_freq, sizeof(t2)); fft_inplace(t2, NFFT, 1);
            memcpy(tx1_time, t1+NFFT-NCP, NCP*sizeof(cpx)); memcpy(tx1_time+NCP, t1, NFFT*sizeof(cpx));
            memcpy(tx2_time, t2+NFFT-NCP, NCP*sizeof(cpx)); memcpy(tx2_time+NCP, t2, NFFT*sizeof(cpx));

            /* ====== [DEBUG CHECK] 檢查 TX 能量 ======
            if (sym == 0) {
                double tx_energy = 0.0;
                for(size_t k=0; k<sym_len; ++k) tx_energy += c_abs_sq(tx1_time[k]);
                printf("\n[DEBUG] Sym 0 TX Energy: %.5f (Should be approx %.2f)\n", 
                       tx_energy, (double)sym_len/NFFT); // 預期約 1.25
                
                // 檢查前幾個值
                printf("[DEBUG] TX Sample[0]: %.4f + j%.4f\n", tx1_time[0].re, tx1_time[0].im);
            }
            // =======================================*/

            // 2. Channel
            apply_mimo_channel_cfo(tx1_time, tx2_time, rx1_curr, rx2_curr, sym_len, snr, cfo_actual, H_ref, global_time_counter);
            global_time_counter += sym_len; 

            // Build search window: prev tail (16) + current (80)
            memcpy(rx1_win, rx1_prev + (sym_len - NCP), NCP*sizeof(cpx));
            memcpy(rx2_win, rx2_prev + (sym_len - NCP), NCP*sizeof(cpx));
            memcpy(rx1_win + NCP, rx1_curr, sym_len*sizeof(cpx));
            memcpy(rx2_win + NCP, rx2_curr, sym_len*sizeof(cpx));

            // 3. Sync (robust EMA) & CFO Est
            size_t best_idx = sync_symbol_robust_mimo(rx1_win, rx2_win, SEARCH_WINDOW_LEN, sync_profile);
            double cfo_est_val = estimate_cfo_at_idx(rx1_win, rx2_win, best_idx);

            static double cfo_est_smooth = 0.0;
            if (sym == 0) cfo_est_smooth = cfo_est_val;
            else          cfo_est_smooth = 0.9 * cfo_est_smooth + 0.1 * cfo_est_val;

            double cfo_used = cfo_est_smooth;

            // 4. Compensation (aligned window)
            cpx rx1_aligned[NFFT], rx2_aligned[NFFT];
            for(size_t k=0; k<NFFT; ++k){
                size_t idx = best_idx + NCP + k;
                if(idx >= SEARCH_WINDOW_LEN) idx = SEARCH_WINDOW_LEN-1; 
                rx1_aligned[k] = rx1_win[idx];
                rx2_aligned[k] = rx2_win[idx];
            }

            // Update prev buffer for next symbol
            memcpy(rx1_prev, rx1_curr, sym_len*sizeof(cpx));
            memcpy(rx2_prev, rx2_curr, sym_len*sizeof(cpx));

            // ====== [插入 2] 監測 ADC/Compensation 後數值 ======
            for(size_t k=0; k<NFFT; ++k){
                double val1 = fmax(fabs(rx1_aligned[k].re), fabs(rx1_aligned[k].im));
                if(val1 > max_abs_adc) max_abs_adc = val1;
            }
            
            // 執行資料區段補償
            apply_cfo_compensation_continuous(rx1_aligned, rx2_aligned, NFFT, cfo_used, &rx_nco_phase);
            
            // 手動推進 CP 區段相位 (維持連續性)
            // Phase change = (2*pi*eps/N) * NCP.  We subtract because compensation reverses channel.
            rx_nco_phase -= (2.0 * M_PI * cfo_used / (double)NFFT) * (double)NCP;

            // 5. FFT
            fft_inplace(rx1_aligned, NFFT, 0);
            fft_inplace(rx2_aligned, NFFT, 0);

            // ====== [插入 3] 監測 FFT 後數值 ======
            for(size_t k=0; k<NFFT; ++k){
                double val1 = fmax(fabs(rx1_aligned[k].re), fabs(rx1_aligned[k].im));
                if(val1 > max_abs_fft) max_abs_fft = val1;
            }

            // 6. MIMO Detection & Demap
            // ====== [插入 4] 監測 ZF 後數值 ======
            for(size_t k=0; k<NFFT; ++k){
                // 執行 ZF 偵測
                mimo_zf_detect(rx1_aligned[k], rx2_aligned[k], H_ref, &x1_arr[k], &x2_arr[k]);
                // 注意：這裡不能有 "}"，迴圈還沒結束！

                // 監測輸出最大值 (這段必須在迴圈內)
                double val_out = fmax(fabs(x1_arr[k].re), fabs(x1_arr[k].im));
                if(val_out > max_abs_zf) max_abs_zf = val_out;
                
                // 額外手算一下 Determinant 大小，這對硬體除法器很重要
                cpx det = c_sub(c_mul(H_ref[0][0], H_ref[1][1]), c_mul(H_ref[0][1], H_ref[1][0]));
                double abs_det = sqrt(det.re*det.re + det.im*det.im);
                if(abs_det > max_abs_det) max_abs_det = abs_det;

            } // <--- 正確的結束位置在這裡！
            

            // Phase Correction
            double residual_err = correct_residual_phase(x1_arr, x2_arr, NFFT);
            // [關鍵修正]：將殘餘誤差「扣」回 NCO 相位
            // 如果 residual_err 是正的，代表訊號轉過頭了，我們要把 NCO 的相位減回來
            // 0.1 是一個 Loop Gain (迴路增益)，避免修正過度震盪，通常 0.1~0.5 都可以
            //rx_nco_phase -= (0.5 * residual_err);
            double loop_gain = (sym < 5) ? 1.0 : 0.1; // 前5個 symbol 全力修正，之後慢慢追蹤
            rx_nco_phase -= (loop_gain * residual_err); 

            // 確保相位在 -PI ~ PI 之間 (雖然不做也行，但為了數值穩定建議做)
            if(rx_nco_phase > M_PI) rx_nco_phase -= 2*M_PI;
            if(rx_nco_phase < -M_PI) rx_nco_phase += 2*M_PI;

            // Demap & Error Count
            size_t rx_bit_idx = (sym * 4 * NFFT);
            for(size_t k=0; k<NFFT; ++k){
                uint8_t b0, b1, b2, b3;
                qpsk_demap(x1_arr[k], &b0, &b1);
                qpsk_demap(x2_arr[k], &b2, &b3);
                
                if (sym >= 20){ // Skip first 20 symbols
                    if(b0 != bits_tx[rx_bit_idx++]) errors++;
                    if(b1 != bits_tx[rx_bit_idx++]) errors++;
                    if(b2 != bits_tx[rx_bit_idx++]) errors++;
                    if(b3 != bits_tx[rx_bit_idx++]) errors++;
                    total_bits_processed += 4;
                }
            }
        } 
        
        free(x1_arr); 
        free(x2_arr);
        
        printf("%5.1f dB | %.5e\n", snr, (double)errors/total_bits_processed);
    }
    
    // ====== [插入 5] 印出分析報告 (請加在 main 函數最後 return 0 之前) ======
    printf("\n=============================================\n");
    printf("       FIXED-POINT DYNAMIC RANGE REPORT      \n");
    printf("=============================================\n");
    // 假設我們想用 Q Format，這裡幫你算出需要的整數 bits (包含 sign bit)
    // 公式: Integer Bits = ceil(log2(Max_Val)) + 1 (for sign)
    int adc_bits = (int)ceil(log2(max_abs_adc + 1.0)) + 1;
    int fft_bits = (int)ceil(log2(max_abs_fft + 1.0)) + 1;
    int zf_bits  = (int)ceil(log2(max_abs_zf + 1.0)) + 1;

    printf("Max ADC Input:   %10.4f | Suggest Int Bits: %d (e.g., Q%d.%d)\n", 
           max_abs_adc, adc_bits, adc_bits, 16-adc_bits);
    printf("Max FFT Output:  %10.4f | Suggest Int Bits: %d (e.g., Q%d.%d)\n", 
           max_abs_fft, fft_bits, fft_bits, 16-fft_bits);
    printf("Max ZF Output:   %10.4f | Suggest Int Bits: %d (e.g., Q%d.%d)\n", 
           max_abs_zf, zf_bits, zf_bits, 16-zf_bits);
    printf("Min Determinant: %10.4f (Check for division stability)\n", 
           max_abs_det); // 這裡其實應該監測 Min 避免除以零，但先看 Max 確定範圍也好
    printf("=============================================\n");
    
    
    
    
    free(tx1_time); free(tx2_time);
    free(rx1_prev); free(rx2_prev);
    free(rx1_curr); free(rx2_curr);
    free(rx1_win); free(rx2_win);
    free(sync_profile);
    free(bits_tx);

    return 0;
}