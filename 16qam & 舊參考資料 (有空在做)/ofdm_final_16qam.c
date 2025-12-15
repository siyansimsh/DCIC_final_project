// ------------------------------------------------------------
//  COM519000 – Homework #3 (Fixed-point Simulation)
//  Base : Homework #2 (OFDM System Simulation)
//  Target : ISO C11
//
//  *** MODIFICATION: 16-QAM Bonus Version ***
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
// *** 可變參數，例如 7 (8-bit), 11 (12-bit), 15 (16-bit) ***
#define FX_FBITS 15
typedef int16_t fx_t; // 基礎容器 (16-bit)
typedef int32_t fx_acc_t; 
typedef int64_t fx_metric_t; 
typedef struct { fx_t re, im; } cpx_fx; 

// *** Word Length 飽和極限 (由 FX_FBITS 決定) ***
#define FX_MAX ( (fx_t)( (1 << (FX_FBITS)) - 1) )
#define FX_MIN ( (fx_t)(-(1 << (FX_FBITS)) ) )

// --- *** NEW 16-QAM Fixed-Point Demapper Constants *** ---
// 決策邊界: 0, +2/sqrt(10), -2/sqrt(10)
// NORM_2 = 2.0 / sqrt(10.0) = 0.63245...
#define QAM16_FX_BOUND_P2 ( (fx_t)( (0.632455532 * (1 << FX_FBITS)) ) )
#define QAM16_FX_BOUND_N2 ( (fx_t)( (-0.632455532 * (1 << FX_FBITS)) ) )
#define QAM16_FX_BOUND_0  ( (fx_t)0 )

// --- 轉換函式 ---
static inline fx_t double_to_fx(double d) {
    d *= (1 << FX_FBITS);
    if (d > FX_MAX) d = FX_MAX; // 飽和 
    if (d < FX_MIN) d = FX_MIN; // 飽和
    return (fx_t)d;
}
static inline double fx_to_double(fx_t x) {
    return (double)x / (double)(1 << FX_FBITS);
}

// --- 定點數運算 (使用飽和宏) ---
static inline fx_t fx_add(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a + (int32_t)b;
    if (temp > FX_MAX) temp = FX_MAX;
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline fx_t fx_sub(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a - (int32_t)b;
    if (temp > FX_MAX) temp = FX_MAX;
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline fx_t fx_mul(fx_t a, fx_t b) {
    int32_t temp = (int32_t)a * (int32_t)b;
    temp = (temp + (1 << (FX_FBITS - 1))) >> FX_FBITS;
    if (temp > FX_MAX) temp = FX_MAX;
    if (temp < FX_MIN) temp = FX_MIN;
    return (fx_t)temp;
}
static inline cpx_fx c_mul_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_sub(fx_mul(a.re, b.re), fx_mul(a.im, b.im)), fx_add(fx_mul(a.re, b.im), fx_mul(a.im, b.re)) };
}
static inline cpx_fx c_add_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_add(a.re, b.re), fx_add(a.im, b.im) };
}
static inline cpx_fx c_sub_fx(cpx_fx a, cpx_fx b) {
    return (cpx_fx){ fx_sub(a.re, b.re), fx_sub(a.im, b.im) };
}
static inline cpx_fx c_conj_fx(cpx_fx a) {
    return (cpx_fx){ a.re, fx_sub(0, a.im) };
}
// =============================================================
//                 END OF HW3 FIXED-POINT TOOLKIT
// =============================================================


/* ====== FFT & helpers (FLOAT) ====== */
// ( ... 浮點數 FFT 函式 ... )
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
// *** MODIFIED FOR 16-QAM ***
#define BITS_PER_SYM 4u           
#define TOTAL_BITS   (1u<<25)     
#define TOTAL_SYMS   (TOTAL_BITS / BITS_PER_SYM)          
#define OFDM_SYMBOLS (TOTAL_SYMS / NFFT)                 

// --- Robust Sync Parameters ---
#define SYNC_SMOOTH_SHIFT 4 
#define SYNC_SMOOTH_ALPHA (1.0 / (double)(1 << SYNC_SMOOTH_SHIFT))
// --- *********************************** ---


// =============================================================
//               HW3: FIXED-POINT RX FUNCTIONS (Steps 8-10)
// =============================================================

/* ====== HW3: CP Helper (FIXED-POINT) ====== */
static void remove_cp_fx(const cpx_fx *in, cpx_fx *out){
    memcpy(out, in+NCP, NFFT*sizeof(cpx_fx));
}

/* ====== HW3: Step 8: Timing sync (FIXED-POINT, ROBUST) ====== */
static size_t sync_symbol_robust_fx(const cpx_fx *rx, size_t rx_len, fx_metric_t *sync_profile_avg) {
    fx_metric_t best_metric = -1;
    size_t best_idx = 0;
    size_t m;
    for (m = 0; m + NFFT + NCP <= rx_len; ++m) {
        fx_acc_t acc_re_fx = 0; fx_acc_t acc_im_fx = 0;
        for (size_t i = 0; i < NCP; ++i) {
            cpx_fx p = c_mul_fx(rx[m + i], c_conj_fx(rx[m + i + NFFT]));
            acc_re_fx += p.re; acc_im_fx += p.im;
        }
        fx_metric_t current_metric = (fx_metric_t)acc_re_fx * acc_re_fx + (fx_metric_t)acc_im_fx * acc_im_fx;
        fx_metric_t old_avg = sync_profile_avg[m];
        fx_metric_t diff = current_metric - old_avg;
        fx_metric_t new_avg = old_avg + (diff >> SYNC_SMOOTH_SHIFT);
        sync_profile_avg[m] = new_avg;
        if (new_avg > best_metric) {
            best_metric = new_avg;
            best_idx = m;
        }
    }
    for (; m < rx_len; m++) {
        sync_profile_avg[m] = sync_profile_avg[m] - (sync_profile_avg[m] >> SYNC_SMOOTH_SHIFT);
    }
    return best_idx;
}


/* ====== HW3: Step 9: FFT (FIXED-POINT) ====== */
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


/* ====== ***  Step 10: 16-QAM demapping (FIXED-POINT) *** ====== */
//
static inline void qam16_demap_fx(cpx_fx s, uint8_t *b0, uint8_t *b1, uint8_t *b2, uint8_t *b3){

    // *** BUG FIX: Compensate for 6-stage scaling (2^6 = 64) in FFT_FX ***
    // 我們必須將訊號 << 6 (乘以 64) 來還原幅度
    // NFFT=64 -> 6 stages -> << 6
    // (使用 32-bit 暫存器來防止 shift 溢位)
    int32_t I_scaled = (int32_t)s.re << 6;
    int32_t Q_scaled = (int32_t)s.im << 6;

    // 將 32-bit 暫存器飽和回我們的 word length (FX_MAX/FX_MIN)
    if (I_scaled > FX_MAX) I_scaled = FX_MAX;
    if (I_scaled < FX_MIN) I_scaled = FX_MIN;
    if (Q_scaled > FX_MAX) Q_scaled = FX_MAX;
    if (Q_scaled < FX_MIN) Q_scaled = FX_MIN;
    
    // 現在 I 和 Q 具有正確的幅度
    fx_t I = (fx_t)I_scaled;
    fx_t Q = (fx_t)Q_scaled;

    // --- I component decision (Boundaries: N2, 0, P2) ---
    // b0 (Sign bit)
    *b0 = (I > QAM16_FX_BOUND_0) ? 0 : 1;
    // b1 (Magnitude bit)
    *b1 = (I > QAM16_FX_BOUND_P2 || I < QAM16_FX_BOUND_N2) ? 1 : 0;
    
    // --- Q component decision (Boundaries: N2, 0, P2) ---
    // b2 (Sign bit)
    *b2 = (Q > QAM16_FX_BOUND_0) ? 0 : 1;
    // b3 (Magnitude bit)
    *b3 = (Q > QAM16_FX_BOUND_P2 || Q < QAM16_FX_BOUND_N2) ? 1 : 0;
}

// =============================================================
//               END OF HW3 FIXED-POINT RX FUNCTIONS
// =============================================================


// ---------- RNG helpers (TX/Channel - float) ------------------
static inline uint32_t prng_u32(void){ static uint32_t s=0x12345678; s = 1664525u*s + 1013904223u; return s; }
static inline double rand_uniform(void){ return (prng_u32()>>8) * (1.0/16777216.0); }
static inline double randn(void){ double u1=rand_uniform()+1e-12, u2=rand_uniform(); return sqrt(-2*log(u1)) * cos(2*M_PI*u2); }

// ---------- *** NEW: 16-QAM mapper (TX - float) *** -----------
//
static inline cpx qam16_map(uint8_t b0, uint8_t b1, uint8_t b2, uint8_t b3) {
    // Gray coding for 4 levels (-3, -1, 1, 3)
    // b0/b1 for I, b2/b3 for Q
    double I = (b0 == 0) ? ((b1 == 0) ? 1.0 : 3.0) : ((b1 == 0) ? -1.0 : -3.0);
    double Q = (b2 == 0) ? ((b3 == 0) ? 1.0 : 3.0) : ((b3 == 0) ? -1.0 : -3.0);
    
    // Normalization factor for average power of 1. Avg power is 10.
    double norm = 1.0 / sqrt(10.0);
    return c_scale(C(I, Q), norm);
}


// ---------- *** NEW: 16-QAM demapper (Benchmark RX - float) *** ---
//
static inline void qam16_demap(cpx s, uint8_t *b0, uint8_t *b1, uint8_t *b2, uint8_t *b3) {
    double norm = 1.0 / sqrt(10.0);
    // De-normalize to make decision boundaries at -2, 0, 2
    double I = s.re / norm;
    double Q = s.im / norm;
    
    // Decision boundaries for I component (at -2, 0, 2)
    *b0 = (I > 0) ? 0 : 1;
    *b1 = (fabs(I) > 2.0) ? 1 : 0;
    
    // Decision boundaries for Q component (at -2, 0, 2)
    *b2 = (Q > 0) ? 0 : 1;
    *b3 = (fabs(Q) > 2.0) ? 1 : 0;
}


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

// ---------- Timing sync (Benchmark RX - float, ROBUST) --------

static size_t sync_symbol_robust_float(const cpx *rx, size_t rx_len, double *sync_profile_avg) {
    double best_metric = -1.0;
    size_t best_idx = 0;
    size_t m;
    for (m = 0; m + NFFT + NCP <= rx_len; ++m) {
        double acc_re = 0.0, acc_im = 0.0;
        for (size_t i = 0; i < NCP; ++i) {
            cpx p = c_mul(rx[m + i], c_conj(rx[m + i + NFFT]));
            acc_re += p.re; acc_im += p.im;
        }
        double current_metric = acc_re * acc_re + acc_im * acc_im;
        double old_avg = sync_profile_avg[m];
        double new_avg = (SYNC_SMOOTH_ALPHA * current_metric) + ((1.0 - SYNC_SMOOTH_ALPHA) * old_avg);
        sync_profile_avg[m] = new_avg;
        if (new_avg > best_metric) {
            best_metric = new_avg;
            best_idx = m;
        }
    }
    for (; m < rx_len; m++) {
        sync_profile_avg[m] = (1.0 - SYNC_SMOOTH_ALPHA) * sync_profile_avg[m];
    }
    return best_idx;
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
//                       MAIN SIMULATION
// =============================================================
int main(int argc,char **argv){
    (void)argc; (void)argv;
    // (*** 1. 搜尋視窗現在是固定的 ***)
    const size_t SEARCH_WINDOW_LEN = NCP + (NFFT + NCP);
    
    FILE *f_step7_in_re = fopen("step7_input_re_16qam.txt", "w");
    FILE *f_step7_in_im = fopen("step7_input_im_16qam.txt", "w");
    FILE *f_step8_out = fopen("step8_output_16qam.txt", "w");
    FILE *f_step9_out_re = fopen("step9_output_re_16qam.txt", "w");
    FILE *f_step9_out_im = fopen("step9_output_im_16qam.txt", "w");
    FILE *f_step10_out = fopen("step10_output_16qam.txt", "w");
    if (!f_step7_in_im ||!f_step7_in_re || !f_step8_out || !f_step9_out_im || !f_step7_in_re || !f_step10_out) {
        perror("Failed to open golden vector files"); return 1;
    }
    int has_dumped_vectors = 0; // 確保只儲存一次的旗標
    
    printf("HW3 Fixed-Point Sim (16-QAM Robust Sync) | N=%u, CP=%u, FX_FBITS=%d, SYNC_ALPHA=1/%d\n",
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
    cpx *search_stream = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    uint8_t *bits_tx = malloc(TOTAL_BITS);
    uint8_t *bits_rx_float = malloc(TOTAL_BITS); 
    
    // --- HW3: Allocate FIXED-POINT buffers (RX) ---
    cpx_fx *rx_time_fx = malloc(NFFT*sizeof(cpx_fx));
    cpx_fx *rx_freq_fx = malloc(NFFT*sizeof(cpx_fx));
    cpx_fx *rx_const_dump_fx = malloc(NFFT*3*sizeof(cpx_fx));
    cpx_fx *search_stream_fx = malloc(SEARCH_WINDOW_LEN * sizeof(cpx_fx));
    uint8_t *bits_rx_fx = malloc(TOTAL_BITS);
    
    // --- Allocate Robust Sync State Buffers ---
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
    // (*** 16-QAM 需要更高的 SNR 來達到低 BER ***)
    const double SNR_dB_list[] = {0, 3, 6, 9, 12, 15, 18, 21}; 
    const size_t NSNR = sizeof(SNR_dB_list)/sizeof(SNR_dB_list[0]);

    // *** NEW FILENAME ***
    FILE *fp_ber = fopen("ber_vs_snr_fixed_robust_16qam.csv", "w");
    if (!fp_ber) { perror("fopen"); return 1; }
    fprintf(fp_ber, "SNR(dB),BER_float,BER_fixed,Incorrect_Idx_Percent\n");

    printf("SNR(dB) | BER (float) | BER (fixed) | Incorrect Idx (%%)\n");
    printf("-----------------------------------------------------------\n");
    
    for(size_t si=0; si<NSNR; ++si){
        double SNRdB = SNR_dB_list[si];
        size_t bit_idx = 0;
        
        uint64_t error_bits_float = 0;
        uint64_t error_bits_fx = 0;
        uint64_t incorrect_indices_count = 0;
        int saved_symbols_count = 0; 

        // --- Reset Sync State for each SNR ---
        memset(sync_profile_avg_float, 0, SEARCH_WINDOW_LEN * sizeof(double));
        memset(sync_profile_avg_fx, 0, SEARCH_WINDOW_LEN * sizeof(fx_metric_t));
        memset(rx_cp_prev, 0, (NFFT+NCP)*sizeof(cpx));

        for(size_t sym=0; sym<OFDM_SYMBOLS; ++sym){
            
            // ==============================================
            // TX 側 (Steps 1-7): 保持浮點數
            // ==============================================
            // *** MODIFIED FOR 16-QAM ***
            for(size_t k=0;k<NFFT;++k){ 
                tx_freq[k] = qam16_map(bits_tx[bit_idx], bits_tx[bit_idx+1], bits_tx[bit_idx+2], bits_tx[bit_idx+3]);
                bit_idx += 4;
            }
            
            memcpy(tx_time, tx_freq, NFFT*sizeof(cpx));
            IFFT(tx_time,NFFT); 
            add_cp(tx_time, tx_cp); 
            memcpy(rx_cp, tx_cp,(NFFT+NCP)*sizeof(cpx));
            awgn(rx_cp,NFFT+NCP,SNRdB);
            
            // (*** 構建 search_stream ***)
            size_t correct_idx_in_search = NCP; 
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
                // (*** Sync & FFT 邏輯不變 ***)
                size_t start_float_idx = sync_symbol_robust_float(search_stream, 
                                                                  SEARCH_WINDOW_LEN, 
                                                                  sync_profile_avg_float);
                remove_cp(search_stream + start_float_idx, rx_time);
                memcpy(rx_freq, rx_time, NFFT*sizeof(cpx));
                FFT(rx_freq, NFFT);

                // *** MODIFIED FOR 16-QAM ***
                uint8_t b0,b1,b2,b3;
                for(size_t k=0;k<NFFT;++k){
                    qam16_demap(rx_freq[k], &b0, &b1, &b2, &b3); 
                    bits_rx_float[(sym*NFFT+k)*4]   = b0; 
                    bits_rx_float[(sym*NFFT+k)*4+1] = b1; 
                    bits_rx_float[(sym*NFFT+k)*4+2] = b2; 
                    bits_rx_float[(sym*NFFT+k)*4+3] = b3; 
                }
                /*
                // (Save float dumps... 邏輯不變)
                if ((fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) && saved_symbols_count < 3) {
                    if (sym == 0) { 
                        memcpy(tx_waveform_dump + saved_symbols_count*(NFFT+NCP), tx_cp, (NFFT+NCP)*sizeof(cpx));
                    }
                    memcpy(rx_waveform_dump + saved_symbols_count*(NFFT+NCP), search_stream + start_float_idx, (NFFT+NCP)*sizeof(cpx));
                    memcpy(rx_const_dump + saved_symbols_count*NFFT, rx_freq, NFFT*sizeof(cpx));
                }*/
            }
            
            // ==============================================
            // RX 側 (Fixed-Point)
            // ==============================================
            {
                // (*** Sync & FFT 邏輯不變 ***)
                size_t start_fx_idx = sync_symbol_robust_fx(search_stream_fx, 
                                                            SEARCH_WINDOW_LEN, 
                                                            sync_profile_avg_fx);
                
                if (sym > 100 && start_fx_idx != correct_idx_in_search) {
                     incorrect_indices_count++;
                }
                
                remove_cp_fx(search_stream_fx + start_fx_idx, rx_time_fx);
                memcpy(rx_freq_fx, rx_time_fx, NFFT*sizeof(cpx_fx));
                FFT_FX(rx_freq_fx, NFFT);

                // *** MODIFIED FOR 16-QAM ***
                uint8_t b0,b1,b2,b3;
                for(size_t k=0;k<NFFT;++k){ 
                    qam16_demap_fx(rx_freq_fx[k], &b0, &b1, &b2, &b3); 
                    bits_rx_fx[(sym*NFFT+k)*4]   = b0; 
                    bits_rx_fx[(sym*NFFT+k)*4+1] = b1; 
                    bits_rx_fx[(sym*NFFT+k)*4+2] = b2; 
                    bits_rx_fx[(sym*NFFT+k)*4+3] = b3; 
                }

                // --- ***「黃金向量儲存」的 IF 區塊 *** ---
                // [修改] 在 dumping 區塊內
                if ( (fabs(SNRdB-3.0)<1e-6) && (has_dumped_vectors == 0) ) {
                    printf(">>> [SNR %.1f] Dumping Golden Vectors (Symbol %zu)...\n", SNRdB, sym);

                    // 0. 先計算完美的 FFT 資料 (移到最前面)
                    cpx_fx rx_time_perfect_fx[NFFT];
                    cpx_fx rx_freq_perfect_fx[NFFT];
                    
                    // 強制從 search_stream_fx[16] 開始移除 CP (完美對齊)
                    remove_cp_fx(search_stream_fx + 16, rx_time_perfect_fx);
                    memcpy(rx_freq_perfect_fx, rx_time_perfect_fx, NFFT * sizeof(cpx_fx));
                    FFT_FX(rx_freq_perfect_fx, NFFT);


                    // 1. Step 7 Input (Verilog RX 的總輸入)
                    for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) {
                        fprintf(f_step7_in_re, "%04hx\n", (uint16_t)search_stream_fx[i].re);
                        fprintf(f_step7_in_im, "%04hx\n", (uint16_t)search_stream_fx[i].im);
                    }

                    // 2. Step 8 Output (Sync 電路的輸出)
                    for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) {
                        fprintf(f_step8_out, "%016llx\n", (uint64_t)sync_profile_avg_fx[i]);
                    }

                    // 3. Step 9 Output (FFT 電路的輸出)
                    // [修正] 改用 rx_freq_perfect_fx 
                    for (size_t k = 0; k < NFFT; ++k) {
                        fprintf(f_step9_out_re, "%04hx\n", (uint16_t)rx_freq_perfect_fx[k].re);
                        fprintf(f_step9_out_im, "%04hx\n", (uint16_t)rx_freq_perfect_fx[k].im);
                    }

                    // 4. Step 10 Output (Demapper 電路的輸出)
                    // [修正] 使用 rx_freq_perfect_fx 計算 Bits
                    for (size_t k = 0; k < NFFT; ++k) {
                        uint8_t b0, b1, b2, b3;
                        qam16_demap_fx(rx_freq_perfect_fx[k], &b0, &b1, &b2, &b3);
                        
                        fprintf(f_step10_out, "%d\n", b0);
                        fprintf(f_step10_out, "%d\n", b1);
                        fprintf(f_step10_out, "%d\n", b2);
                        fprintf(f_step10_out, "%d\n", b3);

                        for (size_t k = 0; k < NFFT; ++k) {
                        // ... demap ...
                        if (k < 5) {
                             printf("C_Check: k=%d, re=%04hx, b0=%d\n", k, (uint16_t)rx_freq_perfect_fx[k].re, b0);
                        }
                        // ...
                    }
                    }

                    has_dumped_vectors = 1; 
                    printf(">>> Golden Vectors Dumped (All synchronized).\n");
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
        // --- Dump CSVs (float & fixed) *** NEW FILENAMES *** ---
        if (fabs(SNRdB-3.0)<1e-6 || fabs(SNRdB-15.0)<1e-6) {
            char f1[64], f2[64], f3[64], f4[64];
            snprintf(f1,sizeof(f1),"waveform_tx_SNR%.0f_16qam.csv", SNRdB);
            snprintf(f2,sizeof(f2),"waveform_rx_float_SNR%.0f_16qam.csv", SNRdB);
            snprintf(f3,sizeof(f3),"const_rx_float_SNR%.0f_16qam.csv",    SNRdB);
            snprintf(f4,sizeof(f4),"const_rx_fixed_SNR%.0f_16qam.csv",    SNRdB);
            
            dump_csv_complex(f1, tx_waveform_dump, (NFFT+NCP)*3);
            dump_csv_complex(f2, rx_waveform_dump, (NFFT+NCP)*3);
            dump_csv_complex(f3, rx_const_dump,    NFFT*3);
            dump_csv_complex_fx(f4, rx_const_dump_fx, NFFT*3);
        }*/
    }
    
    fclose(fp_ber); 
    printf("-----------------------------------------------------------\n");
    // printf("16-QAM Robust sync BER results saved to ber_vs_snr_fixed_robust_16qam.csv\n");

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