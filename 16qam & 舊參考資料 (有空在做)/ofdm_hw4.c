// ------------------------------------------------------------
//  COM519000 – Homework #3 (Fixed-point Simulation)
//  Base : Homework #2 (OFDM System Simulation)
//  Target : ISO C11
//
//  - TX (float) + Channel (float) : Steps 1-7
//  - RX (fixed-point) : Steps 8-10
//  - RX (float) : Benchmark
//
//  *** MODIFICATION: Using Robust Sync Algorithm (Long-term EMA) ***
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

/* ====== Complex type & helpers (FLOAT) ====== */
typedef struct { double re, im; } cpx;
static inline cpx C(double r,double i){ return (cpx){r,i}; }
static inline cpx c_add(cpx a,cpx b){ return C(a.re+b.re, a.im+b.im);}  
static inline cpx c_sub(cpx a,cpx b){ return C(a.re-b.re, a.im-b.im);}  
static inline cpx c_mul(cpx a,cpx b){ return C(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);}
static inline cpx c_conj(cpx a){ return C(a.re, -a.im);}  
static inline cpx c_scale(cpx a,double s){ return C(a.re*s, a.im*s);}  

// =============================================================
//                       HW3: FIXED-POINT TOOLKIT
// =============================================================
#define FX_FBITS 15  // Number of fractional bits (12-bit Q3.9 format)
typedef int16_t fx_t; 
typedef int32_t fx_acc_t; 
typedef int64_t fx_metric_t; 
typedef struct { fx_t re, im; } cpx_fx; 
// *** 12-bit (Q1.11) Saturation limits ***
#define FX_MAX ( (fx_t)( (1 << (FX_FBITS)) - 1) )   
#define FX_MIN ( (fx_t)(-(1 << (FX_FBITS)) ) )     

static inline fx_t double_to_fx(double d) { // convert double to fixed-point
    d *= (1 << FX_FBITS);
    if (d > FX_MAX) d = FX_MAX;
    if (d < FX_MIN) d = FX_MIN;
    return (fx_t)d;
}
static inline double fx_to_double(fx_t x) {
    return (double)x / (double)(1 << FX_FBITS); // convert back to double
}
static inline fx_t fx_add(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a + (int32_t)b; // addition with saturation
    if (temp > FX_MAX) temp = FX_MAX;
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline fx_t fx_sub(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a - (int32_t)b; // subtraction with saturation
    if (temp > FX_MAX) temp = FX_MAX;
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline fx_t fx_mul(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a * (int32_t)b; // multiplication
    temp = (temp + (1 << (FX_FBITS - 1))) >> FX_FBITS;
    if (temp > FX_MAX) temp = FX_MAX; 
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline cpx_fx c_mul_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_sub(fx_mul(a.re, b.re), fx_mul(a.im, b.im)), fx_add(fx_mul(a.re, b.im), fx_mul(a.im, b.re)) }; // a * b
}
static inline cpx_fx c_add_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_add(a.re, b.re), fx_add(a.im, b.im) }; // a + b
}
static inline cpx_fx c_sub_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_sub(a.re, b.re), fx_sub(a.im, b.im) }; // a - b
}
static inline cpx_fx c_conj_fx(cpx_fx a) {
    return (cpx_fx){ a.re, fx_sub(0, a.im) }; // -a.im
}
// =============================================================
//                 END OF HW3 FIXED-POINT TOOLKIT
// =============================================================


/* ====== FFT & helpers (FLOAT) ====== */
static unsigned ilog2u(size_t n){
    unsigned p = 0; while((1ULL<<p) < n) ++p; return p;
}
static size_t bit_reverse(size_t x, unsigned bits){
    size_t r = 0;
    for(unsigned i=0;i<bits;++i){ r = (r<<1) | (x & 1); x >>= 1; } 
    return r;
}
static void bit_reverse_permute(cpx *a, size_t N){
    unsigned bits = ilog2u(N);
    for(size_t i=0;i<N;++i){
        size_t j = bit_reverse(i, bits);
        if(j>i){ cpx t=a[i]; a[i]=a[j]; a[j]=t; }
    }
}
void fft_inplace(cpx *a, size_t N, int inverse){
    bit_reverse_permute(a, N);
    for(size_t m=2; m<=N; m<<=1){ 
        double theta = (inverse ? +2.0 : -2.0) * M_PI / (double)m; 
        double wpr = cos(theta), wpi = sin(theta);
        for(size_t k=0; k<N; k+=m){ 
            double wr = 1.0, wi = 0.0; 
            for(size_t j=0; j<m/2; ++j){
                size_t i0 = k + j;         
                size_t i1 = i0 + m/2;      
                double tr = wr * a[i1].re - wi * a[i1].im;
                double ti = wr * a[i1].im + wi * a[i1].re;
                double ur = a[i0].re, ui = a[i0].im;
                a[i0].re = ur + tr; a[i0].im = ui + ti;
                a[i1].re = ur - tr; a[i1].im = ui - ti;
                double nwr = wr*wpr - wi*wpi;
                double nwi = wr*wpi + wi*wpr;
                wr = nwr; wi = nwi;
            }
        }
    }
    if(inverse){
        double invN = 1.0 / (double)N;
        for(size_t i=0;i<N;++i){ a[i].re *= invN; a[i].im *= invN; }
    }
} 
#define FFT(a,N)  fft_inplace(a,N,0)
#define IFFT(a,N) fft_inplace(a,N,1)


// ──────────────────────────────────────────────────────────────
// Global OFDM parameters
// ----------------------------------------------------------------
#define NFFT    64u               
#define NCP     (NFFT/4)          
#define BITS_PER_SYM 2u             
#define TOTAL_BITS   (1u<<25)     
#define TOTAL_SYMS   (TOTAL_BITS / BITS_PER_SYM)          
#define OFDM_SYMBOLS (TOTAL_SYMS / NFFT)                 

// --- *** NEW: Robust Sync Parameters *** ---
// EMA (IIR filter) alpha parameter
// alpha = 1 / (2^SHIFT)
#define SYNC_SMOOTH_SHIFT 4 // alpha = 1/16 = 0.0625
#define SYNC_SMOOTH_ALPHA (1.0 / (double)(1 << SYNC_SMOOTH_SHIFT))
// --- *********************************** ---


// =============================================================
//               HW3: FIXED-POINT RX FUNCTIONS (Steps 8-10)
// =============================================================

/* ====== HW3: CP Helper (FIXED-POINT) ====== */
static void remove_cp_fx(const cpx_fx *in, cpx_fx *out){
    memcpy(out, in+NCP, NFFT*sizeof(cpx_fx));
}

/* ====== HW3: Step 8: Timing sync (FIXED-POINT, NAIVE) ====== */
static size_t sync_symbol_fx_naive(const cpx_fx *rx, size_t rx_len) {
    fx_metric_t best_metric = -1;
    size_t best_idx = 0;
    for (size_t m = 0; m + NFFT + NCP <= rx_len; ++m) {
        fx_acc_t acc_re_fx = 0; fx_acc_t acc_im_fx = 0;
        for (size_t i = 0; i < NCP; ++i) {
            cpx_fx p = c_mul_fx(rx[m + i], c_conj_fx(rx[m + i + NFFT]));
            acc_re_fx += p.re; acc_im_fx += p.im;
        }
        fx_metric_t metric = (fx_metric_t)acc_re_fx * acc_re_fx + (fx_metric_t)acc_im_fx * acc_im_fx;
        if (metric > best_metric) {
            best_metric = metric; best_idx = m;
        }
    }
    return best_idx;
}

/* ====== *** NEW: Step 8: Timing sync (FIXED-POINT, ROBUST) *** ====== */
//
static size_t sync_symbol_robust_fx(const cpx_fx *rx, size_t rx_len, fx_metric_t *sync_profile_avg, int enable_debug_print) {
    
    fx_metric_t best_metric = -1; // metric 初始化為最小值
    size_t best_idx = 0; // 最佳索引初始化為 0
    size_t m; // 用於迴圈的索引變數

    for (m = 0; m + NFFT + NCP <= rx_len; ++m) {
        // 1. Calculate current metric (same as naive)
        fx_acc_t acc_re_fx = 0; fx_acc_t acc_im_fx = 0;
        for (size_t i = 0; i < NCP; ++i) {
            cpx_fx p = c_mul_fx(rx[m + i], c_conj_fx(rx[m + i + NFFT]));
            acc_re_fx += p.re; acc_im_fx += p.im;
            // [修正後的 printf]
            // 只有在 main 函式允許時 (黃金向量 symbol)
            // 並且只在 m=0 時 (為了匹配 Verilog TB )
            if (enable_debug_print && m == 0) {
                printf("[C] m=%zu, i=%zu, p=(%hd, %hd), acc=(%d, %d)\n",
                        m, i, p.re, p.im, acc_re_fx, acc_im_fx);
            }
        }
        fx_metric_t current_metric = (fx_metric_t)acc_re_fx * acc_re_fx + (fx_metric_t)acc_im_fx * acc_im_fx;

        // 2. Update the persistent average profile (EMA: Exponential Moving Average)
        // avg_new = avg_old + (sample_new - avg_old) >> SHIFT
        // (這是 avg = (1-a)*old + a*new 的硬體高效實作)
        fx_metric_t old_avg = sync_profile_avg[m];
        fx_metric_t diff = current_metric - old_avg;
        fx_metric_t new_avg = old_avg + (diff >> SYNC_SMOOTH_SHIFT);
        sync_profile_avg[m] = new_avg;

        // 3. Find the peak in the *updated averaged* profile
        if (new_avg > best_metric) {
            best_metric = new_avg;
            best_idx = m;
        }
    }
    
    // 讓未計算的部分 "衰減"
    for (; m < rx_len; m++) {
        sync_profile_avg[m] = sync_profile_avg[m] - (sync_profile_avg[m] >> SYNC_SMOOTH_SHIFT);
    }

    return best_idx;
}


/* ====== HW3: Step 9: FFT (FIXED-POINT) ====== */
//
static void bit_reverse_permute_fx(cpx_fx *a, size_t N){
    unsigned bits = ilog2u(N);
    for(size_t i=0;i<N;++i){
        size_t j = bit_reverse(i, bits);
        if(j>i){ cpx_fx t=a[i]; a[i]=a[j]; a[j]=t; }
    }
}
void fft_inplace_fx(cpx_fx *a, size_t N, int inverse){ 
    bit_reverse_permute_fx(a, N);
    for(size_t m=2; m<=N; m<<=1){ 
        double theta = (inverse ? +2.0 : -2.0) * M_PI / (double)m; 
        fx_t wpr_fx = double_to_fx(cos(theta));
        fx_t wpi_fx = double_to_fx(sin(theta));
        for(size_t k=0; k<N; k+=m){ 
            fx_t wr_fx = double_to_fx(1.0); fx_t wi_fx = 0;
            for(size_t j=0; j<m/2; ++j){
                size_t i0 = k + j; size_t i1 = i0 + m/2;
                fx_t tr_fx = fx_sub(fx_mul(wr_fx, a[i1].re), fx_mul(wi_fx, a[i1].im));
                fx_t ti_fx = fx_add(fx_mul(wr_fx, a[i1].im), fx_mul(wi_fx, a[i1].re));
                cpx_fx u = a[i0]; cpx_fx t = {tr_fx, ti_fx};
                a[i0] = c_add_fx(u, t);
                a[i1] = c_sub_fx(u, t);
                a[i0].re >>= 1; a[i0].im >>= 1;
                a[i1].re >>= 1; a[i1].im >>= 1;
                fx_t nwr_fx = fx_sub(fx_mul(wr_fx, wpr_fx), fx_mul(wi_fx, wpi_fx));
                fx_t nwi_fx = fx_add(fx_mul(wr_fx, wpi_fx), fx_mul(wi_fx, wpr_fx));
                wr_fx = nwr_fx; wi_fx = nwi_fx;
            }
        }
    }
}
#define FFT_FX(a,N)  fft_inplace_fx(a,N,0)
#define IFFT_FX(a,N) fft_inplace_fx(a,N,1)


/* ====== HW3: Step 10: QAM demapping (FIXED-POINT) ====== */
//
static inline void qpsk_demap_fx(cpx_fx s, uint8_t *b0, uint8_t *b1){
    fx_t I = s.re; fx_t Q = s.im;
    *b0 = (I < 0) ? 1 : 0; *b1 = (Q < 0) ? 1 : 0; 
    if(*b0 == 1) *b1 ^= 1;
}
// =============================================================
//               END OF HW3 FIXED-POINT RX FUNCTIONS
// =============================================================


// ---------- RNG helpers (TX/Channel - float) ------------------
static inline uint32_t prng_u32(void){ static uint32_t s=0x12345678; s = 1664525u*s + 1013904223u; return s; }
static inline double rand_uniform(void){ return (prng_u32()>>8) * (1.0/16777216.0); }
static inline double randn(void){ double u1=rand_uniform()+1e-12, u2=rand_uniform(); return sqrt(-2*log(u1)) * cos(2*M_PI*u2); }

// ---------- QPSK mapper (TX - float) --------------------------
static inline cpx qpsk_map(uint8_t b0,uint8_t b1){
    double norm = 1.0 / M_SQRT2;
    double I = (b0==0)? 1.0 : -1.0;     
    double Q = (b1==b0)? 1.0 : -1.0;    
    return c_scale(C(I,Q), norm);
}

// ---------- QPSK demapper (Benchmark RX - float) --------------
static inline void qpsk_demap(cpx s, uint8_t *b0, uint8_t *b1){
    double I=s.re, Q=s.im; *b0 = (I<0)?1:0; *b1 = (Q<0)?1:0; if(*b0==1) *b1 ^=1; }

// ---------- AWGN channel (Channel - float) --------------------
static void awgn(cpx *buf,size_t len,double SNRdB){
    double Ps=0.0;
    for(size_t i=0;i<len;++i) Ps += buf[i].re*buf[i].re + buf[i].im*buf[i].im;
    Ps /= (double)len;
    double N0 = Ps / pow(10.0,SNRdB/10.0); 
    double sigma = sqrt(N0/2.0);
    for(size_t i=0;i<len;++i){ buf[i].re += sigma*randn(); buf[i].im += sigma*randn(); }
}

// ---------- CP helpers (TX & Benchmark RX - float) ------------
static void add_cp(const cpx *in,cpx *out){
    memcpy(out,             in+NFFT-NCP, NCP*sizeof(cpx));
    memcpy(out+NCP,         in,          NFFT*sizeof(cpx));
}
static void remove_cp(const cpx *in,cpx *out){
    memcpy(out, in+NCP, NFFT*sizeof(cpx));
}

// ---------- Timing sync (Benchmark RX - float, NAIVE) ----------------
// (保留舊的函式，但不再使用它)
static size_t sync_symbol_naive(const cpx *rx, size_t rx_len){
    double best_metric = -1.0; size_t best_idx = 0;
    for (size_t m = 0; m + NFFT + NCP <= rx_len; ++m) {
        double acc_re = 0.0, acc_im = 0.0;
        for (size_t i = 0; i < NCP; ++i) {
            cpx p = c_mul(rx[m + i], c_conj(rx[m + i + NFFT]));
            acc_re += p.re; acc_im += p.im;
        }
        double metric = acc_re * acc_re + acc_im * acc_im;
        if (metric > best_metric) { // 找到更好的 metric
            best_metric = metric; best_idx = m; // 更新最佳索引
        }
    }
    return best_idx;
}

/* ====== *** NEW: Timing sync (Benchmark RX - float, ROBUST) *** ====== */ 
//使用指數移動平均法來計算同步符號的最佳起始索引
static size_t sync_symbol_robust_float(const cpx *rx, size_t rx_len, double *sync_profile_avg) {
    
    double best_metric = -1.0;
    size_t best_idx = 0;
    size_t m;

    for (m = 0; m + NFFT + NCP <= rx_len; ++m) {
        // 1. Calculate current metric
        double acc_re = 0.0, acc_im = 0.0;
        for (size_t i = 0; i < NCP; ++i) { 
            cpx p = c_mul(rx[m + i], c_conj(rx[m + i + NFFT])); // 乘法並取共軛
            acc_re += p.re; acc_im += p.im; // 累加acc_re與acc_im
        }
        double current_metric = acc_re * acc_re + acc_im * acc_im; // 計算當前的metric

        // 2. Update the persistent average profile (EMA)
        // avg = (alpha * new) + ((1-alpha) * old)
        double old_avg = sync_profile_avg[m]; // 取得舊的平均值
        double new_avg = (SYNC_SMOOTH_ALPHA * current_metric) + ((1.0 - SYNC_SMOOTH_ALPHA) * old_avg); // 計算新的平均值
        sync_profile_avg[m] = new_avg; // 更新平均值

        // 3. Find the peak in the *updated averaged* profile
        if (new_avg > best_metric) { // 如果新的平均值更大，更新最佳metric和索引
            best_metric = new_avg; // 更新最佳metric
            best_idx = m; // 更新最佳索引
        }
    }

    // 讓未計算的部分 "衰減"
    for (; m < rx_len; m++) {
        sync_profile_avg[m] = (1.0 - SYNC_SMOOTH_ALPHA) * sync_profile_avg[m]; 
    }

    return best_idx; // 返回最佳索引
}


// ---------- BER counter (float & fixed) -----------------------
static uint64_t ber_count(const uint8_t *a,const uint8_t *b,size_t n){ uint64_t e=0; for(size_t i=0;i<n;++i) if(a[i]!=b[i]) ++e; return e; }

// ---------- CSV dump helpers ----------------------------------
static void dump_csv_complex(const char *fn,const cpx *v,size_t N){
    FILE *f=fopen(fn,"w"); if(!f){perror("fopen"); return;} fprintf(f,"re,im\n");
    for(size_t i=0;i<N;++i) fprintf(f,"%.9g,%.9g\n",v[i].re,v[i].im); fclose(f);
}
static void dump_csv_complex_fx(const char *fn,const cpx_fx *v,size_t N){
    FILE *f=fopen(fn,"w"); if(!f){perror("fopen"); return;} fprintf(f,"re,im\n");
    for(size_t i=0;i<N;++i) fprintf(f,"%.9g,%.9g\n",fx_to_double(v[i].re),fx_to_double(v[i].im)); fclose(f);
}


// =============================================================
//                       MAIN SIMULATION (FIXED)
// =============================================================
int main(int argc,char **argv){
    (void)argc; (void)argv;
    // (*** 1. 搜尋視窗現在是固定的 ***)
    // 我們將搜尋 "上一個 symbol 的 CP" + "這一個 symbol"
    // 總搜尋長度 = NCP + (NFFT + NCP)
    const size_t SEARCH_WINDOW_LEN = NCP + (NFFT + NCP); 

    // --- *** ADD THIS (1/4) *** ---
    // 宣告黃金向量的檔案指標
    FILE *f_step7_in_re = fopen("step7_input_re.txt", "w");
    FILE *f_step7_in_im = fopen("step7_input_im.txt", "w");
    FILE *f_step8_out = fopen("step8_output.txt", "w");
    FILE *f_step9_out_re = fopen("step9_output_re.txt", "w");
    FILE *f_step9_out_im = fopen("step9_output_im.txt", "w");   
    FILE *f_step10_out = fopen("step10_output.txt", "w");
    if (!f_step7_in_re || !f_step7_in_im || !f_step8_out || !f_step9_out_re || !f_step9_out_im || !f_step10_out) {
        perror("Failed to open golden vector files"); return 1;
    }
    int has_dumped_vectors = 0; // 確保只儲存一次的旗標
    // --- ************************ ---
    
    printf("HW3 Fixed-Point Sim (Robust Sync) | N=%u, CP=%u, FX_FBITS=%d, SYNC_ALPHA=1/%d\n",
           NFFT,NCP,FX_FBITS,(1<<SYNC_SMOOTH_SHIFT));

    // --- Allocate FLOAT buffers (TX & RX Benchmark) ---
    cpx *tx_freq  = malloc(NFFT*sizeof(cpx));
    cpx *tx_time  = malloc(NFFT*sizeof(cpx));
    cpx *tx_cp    = malloc((NFFT+NCP)*sizeof(cpx));
    cpx *rx_cp    = malloc((NFFT+NCP)*sizeof(cpx));
    cpx *rx_time  = malloc(NFFT*sizeof(cpx));
    cpx *rx_freq  = malloc(NFFT*sizeof(cpx));
    cpx *tx_waveform_dump = malloc((NFFT+NCP)*3*sizeof(cpx));
    cpx *rx_waveform_dump = malloc((NFFT+NCP)*3*sizeof(cpx));
    cpx *rx_const_dump    = malloc(NFFT*3*sizeof(cpx));
    cpx *rx_cp_prev = malloc((NFFT+NCP)*sizeof(cpx));
    // (*** 2. 搜尋 stream 大小固定 ***)
    cpx *search_stream = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    uint8_t *bits_tx = malloc(TOTAL_BITS);
    uint8_t *bits_rx_float = malloc(TOTAL_BITS); 
    
    // --- HW3: Allocate FIXED-POINT buffers (RX) ---
    cpx_fx *rx_time_fx = malloc(NFFT*sizeof(cpx_fx));
    cpx_fx *rx_freq_fx = malloc(NFFT*sizeof(cpx_fx));
    cpx_fx *rx_const_dump_fx = malloc(NFFT*3*sizeof(cpx_fx));
    // (*** 3. 搜尋 stream 大小固定 ***)
    cpx_fx *search_stream_fx = malloc(SEARCH_WINDOW_LEN * sizeof(cpx_fx));
    uint8_t *bits_rx_fx = malloc(TOTAL_BITS);
    
    // --- *** NEW: Allocate Robust Sync State Buffers *** ---
    // (*** 4. 平均緩衝區大小固定 ***)
    double *sync_profile_avg_float = calloc(SEARCH_WINDOW_LEN, sizeof(double));
    fx_metric_t *sync_profile_avg_fx = calloc(SEARCH_WINDOW_LEN, sizeof(fx_metric_t));

    // (*** 檢查 OOM ***)
    if(!tx_freq||!tx_time||!tx_cp||!rx_cp||!rx_time||!rx_freq||!bits_tx||!bits_rx_float||
       !rx_cp_prev||!search_stream||!rx_time_fx||!rx_freq_fx||
       !rx_const_dump_fx||!search_stream_fx||!bits_rx_fx ||
       !sync_profile_avg_float || !sync_profile_avg_fx ) { 
        fprintf(stderr,"OOM\n"); return 1; 
    }

    // Generate random bitstream
    for(size_t i=0;i<TOTAL_BITS;i+=32){ uint32_t r=prng_u32(); for(int b=0;b<32&&i+b<TOTAL_BITS;++b) bits_tx[i+b]=(r>>b)&1; }

    // SNR sweep 
    const double SNR_dB_list[] = {0,3,6,9,12,15}; 
    const size_t NSNR = sizeof(SNR_dB_list)/sizeof(SNR_dB_list[0]);

    FILE *fp_ber = fopen("ber_vs_snr_fixed_robust.csv", "w"); 
    if (!fp_ber) { perror("fopen"); return 1; } 
    fprintf(fp_ber, "SNR(dB),BER_float,BER_fixed,Incorrect_Idx_Percent\n");

    printf("SNR(dB) | BER (float) | BER (fixed) | Incorrect Idx (%%)\n");
    printf("-----------------------------------------------------------\n");

    // --- Main SNR loop ---//
    for(size_t si=0; si<NSNR; ++si){
        double SNRdB = SNR_dB_list[si];
        size_t bit_idx = 0;
        
        uint64_t error_bits_float = 0;
        uint64_t error_bits_fx = 0;
        uint64_t incorrect_indices_count = 0;
        int saved_symbols_count = 0; 

        // --- *** NEW: Reset Sync State for each SNR *** ---
        memset(sync_profile_avg_float, 0, SEARCH_WINDOW_LEN * sizeof(double));
        memset(sync_profile_avg_fx, 0, SEARCH_WINDOW_LEN * sizeof(fx_metric_t));
        memset(rx_cp_prev, 0, (NFFT+NCP)*sizeof(cpx));

        for(size_t sym=0; sym<OFDM_SYMBOLS; ++sym){
            
            // ==============================================
            // TX 側 (Steps 1-7): 保持浮點數
            // ==============================================
            for(size_t k=0;k<NFFT;++k){ 
                tx_freq[k] = qpsk_map(bits_tx[bit_idx], bits_tx[bit_idx+1]);
                bit_idx += 2; //serial to 64xparallel
            }
            memcpy(tx_time, tx_freq, NFFT*sizeof(cpx));
            IFFT(tx_time,NFFT); 
            add_cp(tx_time, tx_cp); 
            memcpy(rx_cp, tx_cp,(NFFT+NCP)*sizeof(cpx));
            awgn(rx_cp,NFFT+NCP,SNRdB);
            
            // (*** 5. 構建固定長度的 search_stream ***)
            // [上一個 symbol 的 CP (16)] + [當前 symbol (80)]
            size_t correct_idx_in_search = NCP; // 正確的峰值 *總是* 在索引 16
            
            // 隨機偏移 (0-16)
            size_t timing_offset = prng_u32() % (NCP + 1);
            
            // 構建: [上一個 symbol 的尾巴 (16-offset)] + [當前 symbol (80)] + [一點雜訊 (offset)]
            // 為了簡化，我們只模擬固定大小
            memcpy(search_stream, rx_cp_prev + (NFFT+NCP) - NCP, NCP * sizeof(cpx));
            memcpy(search_stream + NCP, rx_cp, (NFFT + NCP) * sizeof(cpx));
            
            memcpy(rx_cp_prev, rx_cp, (NFFT+NCP)*sizeof(cpx));

            // ==============================================
            // "ADC" 轉換點: float -> fixed-point
            // ==============================================
            for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) {
                search_stream_fx[i].re = double_to_fx(search_stream[i].re);
                search_stream_fx[i].im = double_to_fx(search_stream[i].im);
            }

            // ==============================================
            // RX 側 (Floating-Point Benchmark)
            // ==============================================
            {
                // --- *** Call Robust Sync (固定長度) *** ---
                size_t start_float_idx = sync_symbol_robust_float(search_stream, 
                                                                  SEARCH_WINDOW_LEN, 
                                                                  sync_profile_avg_float);
                
                // (*** 6. [BUG FIX] 直接從 search_stream 中移除 CP ***)
                // (*** 刪除所有 'start_float' 校正邏輯 ***)
                remove_cp(search_stream + start_float_idx, rx_time);
                
                memcpy(rx_freq, rx_time, NFFT*sizeof(cpx));
                FFT(rx_freq, NFFT);

                for(size_t k=0;k<NFFT;++k){
                    uint8_t b0,b1; 
                    qpsk_demap(rx_freq[k], &b0, &b1); 
                    bits_rx_float[(sym*NFFT+k)*2]   = b0; 
                    bits_rx_float[(sym*NFFT+k)*2+1] = b1; 
                }
                /*
                // (Save float dumps... 邏輯不變)
                if ((fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) && saved_symbols_count < 3) {
                    if (sym == 0) { 
                        memcpy(tx_waveform_dump + saved_symbols_count*(NFFT+NCP), tx_cp, (NFFT+NCP)*sizeof(cpx));
                    }
                    // (*** 7. [BUG FIX] 儲存從 search_stream 對齊的波形 ***)
                    memcpy(rx_waveform_dump + saved_symbols_count*(NFFT+NCP), search_stream + start_float_idx, (NFFT+NCP)*sizeof(cpx));
                    memcpy(rx_const_dump + saved_symbols_count*NFFT, rx_freq, NFFT*sizeof(cpx));
                }
                */    
            }
            
            // ==============================================
            // RX 側 (Fixed-Point)
            // ==============================================
            {   
                int enable_debug_print = 0;
                // 檢查是否為 3.0dB 且尚未 dump 過
                if ( (fabs(SNRdB-3.0)<1e-6) && (has_dumped_vectors == 0) ) {
                    enable_debug_print = 1; // 只有這一輪
                    printf(">>> [DEBUG PRINT ENABLED FOR SYMBOL %zu]\n", sym);
                }
                // --- ***************************** ---
                // --- *** Call Robust Sync (固定長度) *** ---
                size_t start_fx_idx = sync_symbol_robust_fx(search_stream_fx, 
                                                            SEARCH_WINDOW_LEN, 
                                                            sync_profile_avg_fx,
                                                            enable_debug_print);
                
                // --- HW3 指標計算 ---
                if (sym > 100 && start_fx_idx != correct_idx_in_search) {
                     incorrect_indices_count++;
                }
                
                // (*** 8. [BUG FIX] 直接從 search_stream_fx 中移除 CP ***)
                remove_cp_fx(search_stream_fx + start_fx_idx, rx_time_fx);
                
                memcpy(rx_freq_fx, rx_time_fx, NFFT*sizeof(cpx_fx));
                FFT_FX(rx_freq_fx, NFFT);

                for(size_t k=0;k<NFFT;++k){ 
                    uint8_t b0,b1; 
                    qpsk_demap_fx(rx_freq_fx[k], &b0, &b1); 
                    bits_rx_fx[(sym*NFFT+k)*2]   = b0; 
                    bits_rx_fx[(sym*NFFT+k)*2+1] = b1; 
                }

                // 在 SNR=3.0dB 且尚未儲存過時，儲存黃金向量
                // (fabs(SNRdB-3.0)<1e-6) 確保我們在 3.0 dB 時觸發
                // [修改] 在 dumping 區塊內
                if ( (fabs(SNRdB-3.0)<1e-6) && (has_dumped_vectors == 0) ) {
                    printf(">>> [SNR %.1f] Dumping Golden Vectors (Symbol %zu)...\n", SNRdB, sym);
                    
                    // 1. Step 7 Input (保持不變: 寫入原始的 search stream)
                    for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) {
                        fprintf(f_step7_in_re, "%04hx\n", (uint16_t)search_stream_fx[i].re);
                        fprintf(f_step7_in_im, "%04hx\n", (uint16_t)search_stream_fx[i].im);
                    }

                    // 2. Step 8 Output (保持不變: 寫入 sync profile)
                    for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) {
                        fprintf(f_step8_out, "%016llx\n", (uint64_t)sync_profile_avg_fx[i]);
                    }

                    // --- [關鍵修正開始] ---
                    // 為了驗證 FFT 硬體模組，我們必須排除 Sync Jitter 的影響。
                    // Verilog Testbench 固定從 index 32 (跳過 16 Old CP + 16 Current CP) 開始讀取。
                    // 所以這裡我們創建一個 "Perfect" 的 FFT 參考數據，強制從 index 16 (Current Symbol Start) 開始切。
                    
                    cpx_fx rx_time_perfect_fx[NFFT];
                    cpx_fx rx_freq_perfect_fx[NFFT];
                    
                    // 強制從 search_stream_fx[16] 開始移除 CP
                    // 這等同於讀取 search_stream_fx 的 index [32] 到 [95]
                    // 這樣就跟 Verilog 的 input_data_re[32+i] 完全對齊了！
                    remove_cp_fx(search_stream_fx + 16, rx_time_perfect_fx);
                    
                    // 複製到頻域 buffer
                    memcpy(rx_freq_perfect_fx, rx_time_perfect_fx, NFFT*sizeof(cpx_fx));
                    
                    // 執行 FFT (產生乾淨的 Golden Pattern)
                    FFT_FX(rx_freq_perfect_fx, NFFT);

                    // 3. Step 9 Output (儲存這個完美對齊的 FFT 結果)
                    for (size_t k = 0; k < NFFT; ++k) {
                        fprintf(f_step9_out_re, "%04hx\n", (uint16_t)rx_freq_perfect_fx[k].re);
                        fprintf(f_step9_out_im, "%04hx\n", (uint16_t)rx_freq_perfect_fx[k].im);
                    }
                    // --- [關鍵修正結束] ---

                    // 4. Step 10 Output (Output based on perfect FFT)
                    for (size_t k = 0; k < NFFT; ++k) {
                        uint8_t b0, b1;
                        qpsk_demap_fx(rx_freq_perfect_fx[k], &b0, &b1);
                        fprintf(f_step10_out, "%d\n", b0);
                        fprintf(f_step10_out, "%d\n", b1);
                    }

                    has_dumped_vectors = 1;
                    printf(">>> Golden Vectors Dumped (Forced Alignment for Verification).\n");
                }
                
                /*
                // (Save fixed-point dumps... 邏輯不變)
                if ((fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) && saved_symbols_count < 3) {
                    memcpy(rx_const_dump_fx + saved_symbols_count*NFFT, rx_freq_fx, NFFT*sizeof(cpx_fx));
                }*/
            }
            
            if ((fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) && saved_symbols_count < 3) {
                saved_symbols_count++;
            }
        }
        
        // --- 計算 HW3 指標 ---
        error_bits_float = ber_count(bits_tx, bits_rx_float, TOTAL_BITS);
        error_bits_fx = ber_count(bits_tx, bits_rx_fx, TOTAL_BITS);
        double BER_float = (double)error_bits_float / (double)TOTAL_BITS;
        double BER_fx = (double)error_bits_fx / (double)TOTAL_BITS;
        double incorrect_idx_percent = (double)incorrect_indices_count / (double)(OFDM_SYMBOLS - 101) * 100.0;

        printf("%5.1f dB | %.5e | %.5e | %8.3f%%      \n", 
               SNRdB, BER_float, BER_fx, incorrect_idx_percent);
        fprintf(fp_ber, "%.1f,%.5e,%.5e,%.5f\n", SNRdB, BER_float, BER_fx, incorrect_idx_percent);
        
        /*
        // --- Dump CSVs (float & fixed) ---
        if (fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) {
            char f1[64], f2[64], f3[64], f4[64];
            snprintf(f1,sizeof(f1),"waveform_tx_SNR%.0f_robust.csv", SNRdB);
            snprintf(f2,sizeof(f2),"waveform_rx_float_SNR%.0f_robust.csv", SNRdB);
            snprintf(f3,sizeof(f3),"const_rx_float_SNR%.0f_robust.csv",    SNRdB);
            snprintf(f4,sizeof(f4),"const_rx_fixed_SNR%.0f_robust.csv",    SNRdB);
            
            dump_csv_complex(f1, tx_waveform_dump, (NFFT+NCP)*3);
            dump_csv_complex(f2, rx_waveform_dump, (NFFT+NCP)*3);
            dump_csv_complex(f3, rx_const_dump,    NFFT*3);
            dump_csv_complex_fx(f4, rx_const_dump_fx, NFFT*3);
        }*/
    }

    fclose(fp_ber); 
    printf("-----------------------------------------------------------\n");
    // printf("Robust sync BER results saved to ber_vs_snr_fixed_robust.csv\n");

    // 關閉黃金向量的檔案
    fclose(f_step7_in_re);
    fclose(f_step7_in_im);
    fclose(f_step8_out);
    fclose(f_step9_out_re);
    fclose(f_step9_out_im);
    fclose(f_step10_out);
    printf("Golden vector files closed.\n");

    // --- Free all buffers ---
    free(tx_freq); free(tx_time); free(tx_cp); free(rx_cp); free(rx_time); free(rx_freq);
    free(tx_waveform_dump); free(rx_waveform_dump); free(rx_const_dump);
    free(rx_cp_prev); free(search_stream); free(bits_tx); free(bits_rx_float);
    free(rx_time_fx); free(rx_freq_fx); free(rx_const_dump_fx);
    free(search_stream_fx); free(bits_rx_fx);
    free(sync_profile_avg_float);
    free(sync_profile_avg_fx);
    
    return 0;
}