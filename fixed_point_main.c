// =============================================================
//  Fixed-Point 2x2 MIMO-OFDM Receiver (Verified BER + Golden Vector)
//  Based on user provided fixed_point_main.c
// =============================================================

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

// =============================================================
//  1. Fixed-Point Library (Q5.11)
// =============================================================
#define Q_FRAC 11
#define ONE_FIX (1 << Q_FRAC)  // 2048
#define INV_TWO_PI   ((fix16)327)
#define MAX_FIX 32767
#define MIN_FIX -32768

typedef int16_t fix16;
typedef int32_t fix32;

static inline fix16 sat_add(fix16 a, fix16 b) {
    int32_t res = (int32_t)a + b;
    if (res > MAX_FIX) return MAX_FIX;
    if (res < MIN_FIX) return MIN_FIX;
    return (fix16)res;
}

static inline fix16 sat_sub(fix16 a, fix16 b) {
    int32_t res = (int32_t)a - b;
    if (res > MAX_FIX) return MAX_FIX;
    if (res < MIN_FIX) return MIN_FIX;
    return (fix16)res;
}

static inline fix16 DOUBLE_TO_FIX(double d) {
    double val = d * ONE_FIX;
    if (val > MAX_FIX) return MAX_FIX;
    if (val < MIN_FIX) return MIN_FIX;
    return (fix16)(val + (val >= 0 ? 0.5 : -0.5));
}

static inline double FIX_TO_DOUBLE(fix16 x) {
    return (double)x / ONE_FIX;
}

typedef struct { fix16 re; fix16 im; } cpx_fix;

static inline cpx_fix C_FIX(fix16 r, fix16 i) {
    cpx_fix c; c.re = r; c.im = i; return c;
}
static inline cpx_fix c_add_fix(cpx_fix a, cpx_fix b) {
    cpx_fix c;
    c.re = sat_add(a.re, b.re);
    c.im = sat_add(a.im, b.im);
    return c;
}
static inline cpx_fix c_sub_fix(cpx_fix a, cpx_fix b) {
    cpx_fix c;
    c.re = sat_sub(a.re, b.re);
    c.im = sat_sub(a.im, b.im);
    return c;
}
static inline cpx_fix c_mul_fix(cpx_fix a, cpx_fix b) {
    fix32 re = ((fix32)a.re * b.re - (fix32)a.im * b.im) >> Q_FRAC;
    fix32 im = ((fix32)a.re * b.im + (fix32)a.im * b.re) >> Q_FRAC;
    return C_FIX((fix16)re, (fix16)im);
}

// =============================================================
//  2. CORDIC Algorithm
// =============================================================
#define FIX_PI 6434
#define FIX_TWO_PI 12868
#define FIX_HALF_PI 3217
#define FIX_QUARTER_PI 1608
#define CORDIC_ITER 12
#define CORDIC_GAIN 1243

const fix16 cordic_angles[] = {
    1608, 949, 501, 254, 127, 63, 31, 15, 7, 3, 1, 0
};

void cordic_rotate(fix16 *x, fix16 *y, fix16 z) {
    fix32 xi = *x; fix32 yi = *y; fix32 zi = z;
    fix32 xt, yt;
    for (int i = 0; i < CORDIC_ITER; i++) {
        if (zi < 0) {
            xt = xi + (yi >> i); yt = yi - (xi >> i); zi += cordic_angles[i];
        } else {
            xt = xi - (yi >> i); yt = yi + (xi >> i); zi -= cordic_angles[i];
        }
        xi = xt; yi = yt;
    }
    *x = (fix16)((xi * CORDIC_GAIN) >> Q_FRAC);
    *y = (fix16)((yi * CORDIC_GAIN) >> Q_FRAC);
}

fix16 cordic_vectoring(fix16 x, fix16 y) {
    fix32 xi = x; fix32 yi = y; fix32 zi = 0;
    fix32 xt, yt;
    if (xi < 0) { xi = -xi; yi = -yi; zi = FIX_PI; }
    
    for (int i = 0; i < CORDIC_ITER; i++) {
        if (yi > 0) {
            xt = xi + (yi >> i); yt = yi - (xi >> i); zi += cordic_angles[i];
        } else {
            xt = xi - (yi >> i); yt = yi + (xi >> i); zi -= cordic_angles[i];
        }
        xi = xt; yi = yt;
    }
    return (fix16)zi;
}

// =============================================================
//  3. Fixed-Point Modules
// =============================================================
#define NFFT 64
#define NCP  16
#define STREAM_LEN (NFFT + 3*NCP)
#define SEARCH_WINDOW_LEN (NCP + NFFT + NCP)

// 3.1 FFT (Fixed)
void fft_fixed(cpx_fix *a, int inverse) {
    int j = 0;
    for (int i = 0; i < NFFT - 1; i++) {
        if (i < j) {
            cpx_fix t = a[i];
            a[i] = a[j];
            a[j] = t;
        }
        int k = NFFT / 2;
        while (k <= j) { j -= k; k /= 2; }
        j += k;
    }
    for (int m = 2; m <= NFFT; m <<= 1) {
        double theta_step = (inverse ? 2.0 : -2.0) * M_PI / m;
        cpx_fix wm = C_FIX((fix16)(cos(theta_step)*ONE_FIX), (fix16)(sin(theta_step)*ONE_FIX));
        for (int k = 0; k < NFFT; k += m) {
            cpx_fix w = C_FIX(ONE_FIX, 0);
            for (int j = 0; j < m / 2; j++) {
                cpx_fix t = c_mul_fix(w, a[k + j + m / 2]);
                cpx_fix u = a[k + j];
                a[k + j] = c_add_fix(u, t);
                a[k + j + m / 2] = c_sub_fix(u, t);
                w = c_mul_fix(w, wm);
            }
        }
    }
}

// 3.2 Matrix Inversion
int calc_zf_matrix_fixed(cpx_fix H[2][2], cpx_fix H_inv[2][2]) {
    cpx_fix t1 = c_mul_fix(H[0][0], H[1][1]);
    cpx_fix t2 = c_mul_fix(H[0][1], H[1][0]);
    cpx_fix det = c_sub_fix(t1, t2);
    
    fix32 det_mag_sq = ((fix32)det.re * det.re + (fix32)det.im * det.im) >> Q_FRAC;
    if (det_mag_sq < 10) return 0; 
    
    fix32 inv_factor = (ONE_FIX * ONE_FIX) / det_mag_sq; 
    cpx_fix det_inv;
    det_inv.re = (fix16)(((fix32)det.re * inv_factor) >> Q_FRAC);
    det_inv.im = (fix16)((-(fix32)det.im * inv_factor) >> Q_FRAC); 
    
    H_inv[0][0] = c_mul_fix(H[1][1], det_inv);
    H_inv[0][1] = c_mul_fix(C_FIX(-H[0][1].re, -H[0][1].im), det_inv); 
    H_inv[1][0] = c_mul_fix(C_FIX(-H[1][0].re, -H[1][0].im), det_inv);
    H_inv[1][1] = c_mul_fix(H[0][0], det_inv);
    return 1;
}

// 3.3 MIMO ZF
void mimo_zf_fixed(cpx_fix y1, cpx_fix y2, cpx_fix H_inv[2][2], cpx_fix *x1, cpx_fix *x2) {
    *x1 = c_add_fix(c_mul_fix(H_inv[0][0], y1), c_mul_fix(H_inv[0][1], y2));
    *x2 = c_add_fix(c_mul_fix(H_inv[1][0], y1), c_mul_fix(H_inv[1][1], y2));
}

// 3.4 Phase Tracking
fix16 correct_residual_phase_fixed(cpx_fix *x1, cpx_fix *x2, size_t n) {
    fix32 total_err = 0;
    for(size_t k=0; k<n; ++k) {
        fix16 ang = cordic_vectoring(x1[k].re, x1[k].im); 
        fix16 diff = ang - FIX_QUARTER_PI; 
        while(diff > FIX_HALF_PI/2)  diff -= FIX_HALF_PI;
        while(diff < -FIX_HALF_PI/2) diff += FIX_HALF_PI;
        total_err += diff;

        ang = cordic_vectoring(x2[k].re, x2[k].im);
        diff = ang - FIX_QUARTER_PI;
        while(diff > FIX_HALF_PI/2)  diff -= FIX_HALF_PI;
        while(diff < -FIX_HALF_PI/2) diff += FIX_HALF_PI;
        total_err += diff;
    }
    //printf("[DEBUG C] Idx 0 Input: (%d, %d)\n", x1[0].re, x1[0].im);
    fix16 avg_err = (fix16)(total_err >> 7);
    //printf("[DEBUG C] Avg_Err = %d\n", avg_err); // Debug Line
    fix16 neg_err = -avg_err;
    for(size_t k=0; k<n; ++k){
        cordic_rotate(&x1[k].re, &x1[k].im, neg_err);
        cordic_rotate(&x2[k].re, &x2[k].im, neg_err);
    }
    //printf("[DEBUG C] Idx 0 Output: (%d, %d)\n", x1[0].re, x1[0].im);
    return avg_err;
}

// 3.5 Demap
void qpsk_demap_fixed(cpx_fix s, uint8_t *b0, uint8_t *b1){
    *b0 = (s.re < 0) ? 1 : 0;
    *b1 = (s.im < 0) ? 1 : 0;
}

// =============================================================
//  4. Floating Point Helpers (For Environment Simulation)
// =============================================================
typedef struct { double re, im; } cpx;
static inline cpx C(double r,double i){ return (cpx){r,i}; }
static inline cpx c_add(cpx a,cpx b){ return C(a.re+b.re, a.im+b.im);}  
static inline cpx c_sub(cpx a,cpx b){ return C(a.re-b.re, a.im-b.im);}  
static inline cpx c_mul(cpx a,cpx b){ return C(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);}
static inline cpx c_scale(cpx a,double s){ return C(a.re*s, a.im*s);} 
static inline cpx c_conj(cpx a){ return C(a.re, -a.im); }
static inline uint32_t prng_u32(void){ static uint32_t s=0x12345678; s = 1664525u*s + 1013904223u; return s; }
static inline double rand_uniform(void){ return (prng_u32()>>8) * (1.0/16777216.0); }
static inline double randn(void){ double u1=rand_uniform()+1e-12, u2=rand_uniform(); return sqrt(-2*log(u1)) * cos(2*M_PI*u2); }
static inline cpx qpsk_map(uint8_t b0,uint8_t b1){
    double s = 1.0/M_SQRT2; return C((b0? -s:s), (b1? -s:s)); 
}

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

static void apply_mimo_channel_cfo(const cpx *tx1, const cpx *tx2, cpx *rx1, cpx *rx2, 
    size_t len, double SNRdB, double cfo_norm, cpx H[2][2], size_t time_offset) {
    double signal_power = 1.0 / (double)NFFT; 
    double noise_power = signal_power * pow(10.0, -SNRdB/10.0);
    double sigma = sqrt(noise_power/2.0);

    for(size_t i=0; i<len; ++i){
        cpx y1 = c_add(c_mul(H[0][0], tx1[i]), c_mul(H[0][1], tx2[i]));
        cpx y2 = c_add(c_mul(H[1][0], tx1[i]), c_mul(H[1][1], tx2[i]));
        double phase = 2.0 * M_PI * (cfo_norm / (double)NFFT) * (double)(time_offset + i);
        cpx rot = C(cos(phase), sin(phase));
        y1 = c_mul(y1, rot); y2 = c_mul(y2, rot);
        y1.re += sigma*randn(); y1.im += sigma*randn();
        y2.re += sigma*randn(); y2.im += sigma*randn();
        rx1[i] = y1; rx2[i] = y2;
    }
}

static size_t sync_symbol_robust_mimo_float(
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
        double new_avg = old_avg + (metric - old_avg) * (1.0 / 16.0);
        profile_ema[m] = new_avg;
        if (new_avg > best_metric) {
            best_metric = new_avg;
            best_idx = m;
        }
    }

    for (; m < win_len; ++m) {
        profile_ema[m] = profile_ema[m] * (1.0 - 1.0/16.0);
    }

    return best_idx;
}

static double estimate_cfo_at_idx_float(
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

static size_t sync_and_estimate_cfo_mrc(
    const cpx *rx1_stream, const cpx *rx2_stream,
    size_t stream_len,
    double *profile,
    double *est_cfo_norm
) {
    const cpx *r1 = rx1_stream + NCP;
    const cpx *r2 = rx2_stream + NCP;
    for (size_t i = 0; i < stream_len; ++i) profile[i] = 0.0;
    size_t m = 0;
    cpx corr_total = C(0,0);
    for (size_t i = 0; i < NCP; ++i) {
        size_t idx_cp   = m + i;
        size_t idx_tail = m + i + NFFT;
        cpx r1_cp = r1[idx_cp]; cpx r1_tail = r1[idx_tail];
        cpx r2_cp = r2[idx_cp]; cpx r2_tail = r2[idx_tail];
        corr_total = c_add(corr_total, c_mul(r1_cp, c_conj(r1_tail)));
        corr_total = c_add(corr_total, c_mul(r2_cp, c_conj(r2_tail)));
    }
    double angle   = atan2(corr_total.im, corr_total.re);
    double eps_hat = -angle / (2.0 * M_PI);
    if (eps_hat >= 0.5) eps_hat -= 1.0;
    if (eps_hat < -0.5) eps_hat += 1.0;
    *est_cfo_norm = eps_hat;
    static int cfo_debug_printed = 0;
    if (!cfo_debug_printed) {
        printf("\n=============================================\n");
        printf("       CFO ESTIMATOR GOLDEN VECTORS          \n");
        printf("=============================================\n");
        // 印出複數相關值的實部與虛部 (這對應 Verilog 的 acc_re, acc_im)
        printf("Golden_Corr_Re: %d\n", (int)corr_total.re);
        printf("Golden_Corr_Im: %d\n", (int)corr_total.im);
        
        // 印出計算出的浮點角度
        printf("Golden_Angle_Rad: %.5f\n", angle);
        
        // 預估 Verilog CORDIC 應該跑出的數值 (Z)
        // Verilog PI = 6434, 所以 Z = Angle * (6434 / PI)
        int expected_z = (int)(angle * (6434.0 / M_PI));
        printf("Expected_Verilog_Angle_Out (Z): %d\n", expected_z);
        printf("=============================================\n");
        
        cfo_debug_printed = 1;
    }
    // -------------------------------------
    return NCP;
}

static void apply_cfo_compensation_continuous(
    cpx *r1, cpx *r2, 
    size_t len, 
    double cfo_est_norm, 
    double *phase_track 
){
    double phase_inc = (2.0 * M_PI * cfo_est_norm) / (double)NFFT;
    for(size_t i=0; i<len; ++i){
        double phi = *phase_track;
        cpx rot = C(cos(phi), sin(phi));
        r1[i] = c_mul(r1[i], rot);
        r2[i] = c_mul(r2[i], rot);
        *phase_track -= phase_inc;
        if(*phase_track > M_PI) *phase_track -= 2*M_PI;
        if(*phase_track < -M_PI) *phase_track += 2*M_PI;
    }
}

// =============================================================
//  Helper to Write Hex (Inserted Utility)
// =============================================================
void write_cpx_hex(FILE *fp, cpx_fix val) {
    if (fp) fprintf(fp, "%04hX%04hX ", (unsigned short)val.re, (unsigned short)val.im);
}

// =============================================================
//  5. Main
// =============================================================
#define TOTAL_BITS (1<<20)
#define BITS_PER_SYM 2
#define NUM_TX 2

int main(){
    printf("=== Fixed-Point 2x2 MIMO-OFDM Receiver (Q5.11) + Golden Vector ===\n");
    
    // 檔案指標
    FILE *f_adc = fopen("golden_adc_time.hex", "w"); // Raw ADC (CP+data) for RTL sync input
    if (!f_adc) { perror("open golden_adc_time.hex"); return 1; }
    else { printf("[INFO] opened golden_adc_time.hex\n"); }
    FILE *f_in  = fopen("golden_input_time.hex", "w"); // Input to FFT (aligned, post-sync)
    FILE *f_fft = fopen("golden_fft_out.hex", "w"); // Output of FFT
    FILE *f_zf  = fopen("golden_zf_out.hex", "w"); // Output of MIMO ZF
    FILE *f_pt  = fopen("golden_final_out.hex", "w"); // Output after Phase Tracking
    FILE *f_h   = fopen("golden_h_inv.hex", "w"); // H_inv Matrix

    // New reference dumps for RTL correlation (single symbol capture)
    FILE *f_ref_in   = fopen("ref_fft_in.hex", "w");   // CP removed + CFO compensated time samples
    FILE *f_ref_fft  = fopen("ref_fft_out.hex", "w");  // FFT output of the above

    if (!f_ref_in)  { perror("open ref_fft_in.hex");  return 1; }
    if (!f_ref_fft) { perror("open ref_fft_out.hex"); return 1; }

    size_t sym_len = NFFT + NCP;
    cpx *rx1_float = malloc(sym_len * sizeof(cpx));
    cpx *rx2_float = malloc(sym_len * sizeof(cpx));
    cpx *rx1_prev = calloc(sym_len, sizeof(cpx));
    cpx *rx2_prev = calloc(sym_len, sizeof(cpx));
    cpx *rx1_win  = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    cpx *rx2_win  = malloc(SEARCH_WINDOW_LEN * sizeof(cpx));
    
    cpx_fix x1_arr[NFFT], x2_arr[NFFT];
    cpx_fix rx1_aligned_fix[NFFT]; 
    cpx_fix rx2_aligned_fix[NFFT];
    
    uint8_t *bits_tx = malloc(TOTAL_BITS);
    for(size_t i=0; i<TOTAL_BITS; ++i) bits_tx[i] = (prng_u32() >> 16) & 1;

    double cfo_actual = 0.03;
    double snr_list[] = {0, 5, 10, 15, 20, 25};
    size_t num_syms = TOTAL_BITS / (BITS_PER_SYM * NUM_TX * NFFT);

    // [Step A] 生成 Rayleigh Channel (Float)
    cpx H_ref[2][2];
    double scale = 1.0 / sqrt(2.0);
    H_ref[0][0] = C(randn()*scale, randn()*scale);
    H_ref[0][1] = C(randn()*scale, randn()*scale);
    H_ref[1][0] = C(randn()*scale, randn()*scale);
    H_ref[1][1] = C(randn()*scale, randn()*scale);
    
    printf("Random Channel: H00=%.2f+j%.2f\n", H_ref[0][0].re, H_ref[0][0].im);

    // [Step B] 將 H 轉換為 Fixed Point
    cpx_fix H_fix[2][2], H_inv_fix[2][2];
    H_fix[0][0] = C_FIX(DOUBLE_TO_FIX(H_ref[0][0].re), DOUBLE_TO_FIX(H_ref[0][0].im));
    H_fix[0][1] = C_FIX(DOUBLE_TO_FIX(H_ref[0][1].re), DOUBLE_TO_FIX(H_ref[0][1].im));
    H_fix[1][0] = C_FIX(DOUBLE_TO_FIX(H_ref[1][0].re), DOUBLE_TO_FIX(H_ref[1][0].im));
    H_fix[1][1] = C_FIX(DOUBLE_TO_FIX(H_ref[1][1].re), DOUBLE_TO_FIX(H_ref[1][1].im));
    
    // [Step C] 預先計算反矩陣
    if (!calc_zf_matrix_fixed(H_fix, H_inv_fix)) {
        printf("Error: Singular Matrix in Fixed Point!\n");
        return 1;
    }

    // --- Write H_inv Matrix (Golden) ---
    // 修改：移除 write_cpx_hex，改用緊密排列 fprintf
    // 格式：ReIm (共32-bit，無空白)
    // fprintf(f_h, "// H_inv in order: H00, H01, H10, H11 (Re Im)\n"); // 註解可以拿掉以免 readmemh 讀到
    fprintf(f_h, "%04hX%04hX\n", H_inv_fix[0][0].re, H_inv_fix[0][0].im);
    fprintf(f_h, "%04hX%04hX\n", H_inv_fix[0][1].re, H_inv_fix[0][1].im);
    fprintf(f_h, "%04hX%04hX\n", H_inv_fix[1][0].re, H_inv_fix[1][0].im);
    fprintf(f_h, "%04hX%04hX\n", H_inv_fix[1][1].re, H_inv_fix[1][1].im);
        printf("SNR(dB) | BER (Fixed)\n");
        printf("---------------------\n");

    int wrote_adc = 0;
    for(int s=0; s<6; ++s){
        double snr = snr_list[s];
        size_t errors = 0, total_bits_processed = 0;
        size_t bit_idx = 0;
        size_t global_time = 0;
        
        double nco_phase = 0;
        double cfo_smooth = 0;
        double sync_profile[SEARCH_WINDOW_LEN];
        for (size_t i = 0; i < SEARCH_WINDOW_LEN; ++i) sync_profile[i] = 0.0;
        
        // 設定觸發條件：只在 SNR=25dB 時輸出前 2 個 Symbol
        int capture_enabled_for_snr = (snr == 25.0); 

        for(size_t sym=0; sym < num_syms; ++sym){
            int do_write = (capture_enabled_for_snr && sym < 1);

            // =========================================================
            // 1. TX & Channel
            // =========================================================
            cpx t1_f[NFFT], t2_f[NFFT];
            for(int k=0; k<NFFT; k++) {
                t1_f[k] = qpsk_map(bits_tx[bit_idx], bits_tx[bit_idx+1]);
                t2_f[k] = qpsk_map(bits_tx[bit_idx+2], bits_tx[bit_idx+3]);
                bit_idx += 4;
            }
            fft_inplace(t1_f, NFFT, 1);
            fft_inplace(t2_f, NFFT, 1);
            
            cpx tx1_with_cp[NFFT+NCP], tx2_with_cp[NFFT+NCP];
            memcpy(tx1_with_cp, t1_f + NFFT - NCP, NCP * sizeof(cpx));
            memcpy(tx1_with_cp + NCP, t1_f, NFFT * sizeof(cpx));
            memcpy(tx2_with_cp, t2_f + NFFT - NCP, NCP * sizeof(cpx));
            memcpy(tx2_with_cp + NCP, t2_f, NFFT * sizeof(cpx));

            apply_mimo_channel_cfo(tx1_with_cp, tx2_with_cp, rx1_float, rx2_float, NFFT+NCP, snr, cfo_actual, H_ref, global_time);
            global_time += (NFFT+NCP);

            // =========================================================
            // 2. Dump raw ADC (CP+data) for RTL (first capture only)
            //    Tie this to the same capture condition (do_write)
            // =========================================================
            if (!wrote_adc && do_write) {
                for (int k = 0; k < NFFT + NCP; ++k) {
                    fprintf(f_adc, "%04hX%04hX%04hX%04hX\n",
                        (unsigned short)DOUBLE_TO_FIX(rx1_float[k].re), (unsigned short)DOUBLE_TO_FIX(rx1_float[k].im),
                        (unsigned short)DOUBLE_TO_FIX(rx2_float[k].re), (unsigned short)DOUBLE_TO_FIX(rx2_float[k].im));
                }
                wrote_adc = 1;
                printf("[INFO] wrote ADC dump (sym %zu, snr %.1f)\n", sym, snr);
            }

            // =========================================================
            // 3. Sync (robust EMA) & CFO Estimation/Compensation
            // =========================================================
            size_t prev_copy = (sym == 0) ? 0 : NCP;
            if (prev_copy) {
                memcpy(rx1_win, rx1_prev + (sym_len - prev_copy), prev_copy * sizeof(cpx));
                memcpy(rx2_win, rx2_prev + (sym_len - prev_copy), prev_copy * sizeof(cpx));
            }
            size_t curr_copy = SEARCH_WINDOW_LEN - prev_copy;
            if (curr_copy > sym_len) curr_copy = sym_len;
            memcpy(rx1_win + prev_copy, rx1_float, curr_copy * sizeof(cpx));
            memcpy(rx2_win + prev_copy, rx2_float, curr_copy * sizeof(cpx));
            size_t filled = prev_copy + curr_copy;
            for (size_t i = filled; i < SEARCH_WINDOW_LEN; ++i) {
                rx1_win[i] = C(0.0, 0.0);
                rx2_win[i] = C(0.0, 0.0);
            }

            size_t best_idx = sync_symbol_robust_mimo_float(rx1_win, rx2_win, SEARCH_WINDOW_LEN, sync_profile);
            double eps_hat = estimate_cfo_at_idx_float(rx1_win, rx2_win, best_idx);

            if (sym == 0 && capture_enabled_for_snr) {
                int eps_q_dbg = (int)llround(eps_hat * ONE_FIX);
                printf("[INFO][CMODEL] snr=%.1f sync_start_index=%zu estimated_cfo=%.9f (q5.11=%d)\n",
                       snr, best_idx, eps_hat, eps_q_dbg);
            }

            if (sym == 0) cfo_smooth = eps_hat;
            else cfo_smooth = 0.9 * cfo_smooth + 0.1 * eps_hat;

            cpx y1_aligned[NFFT], y2_aligned[NFFT];
            for (size_t k = 0; k < NFFT; ++k) {
                size_t idx = best_idx + NCP + k;
                if (idx >= SEARCH_WINDOW_LEN) idx = SEARCH_WINDOW_LEN - 1;
                y1_aligned[k] = rx1_win[idx];
                y2_aligned[k] = rx2_win[idx];
            }

            apply_cfo_compensation_continuous(y1_aligned, y2_aligned, NFFT, cfo_smooth, &nco_phase);
            nco_phase -= (2.0 * M_PI * cfo_smooth / (double)NFFT) * (double)NCP;
            if (nco_phase > M_PI) nco_phase -= 2.0 * M_PI;
            if (nco_phase < -M_PI) nco_phase += 2.0 * M_PI;

            for (int k = 0; k < NFFT; ++k) {
                rx1_aligned_fix[k].re = DOUBLE_TO_FIX(y1_aligned[k].re);
                rx1_aligned_fix[k].im = DOUBLE_TO_FIX(y1_aligned[k].im);
                rx2_aligned_fix[k].re = DOUBLE_TO_FIX(y2_aligned[k].re);
                rx2_aligned_fix[k].im = DOUBLE_TO_FIX(y2_aligned[k].im);

                if (do_write) {
                    fprintf(f_in, "%04hX%04hX%04hX%04hX\n",
                        (unsigned short)rx1_aligned_fix[k].re, (unsigned short)rx1_aligned_fix[k].im,
                        (unsigned short)rx2_aligned_fix[k].re, (unsigned short)rx2_aligned_fix[k].im);
                    fprintf(f_ref_in, "%04hX%04hX%04hX%04hX\n",
                        (unsigned short)rx1_aligned_fix[k].re, (unsigned short)rx1_aligned_fix[k].im,
                        (unsigned short)rx2_aligned_fix[k].re, (unsigned short)rx2_aligned_fix[k].im);
                }
            }

            if (do_write) {
                int eps_q = (int)llround(eps_hat * ONE_FIX);
                printf("[INFO][CMODEL] sync_start_index=%zu estimated_cfo=%.9f (q5.11=%d)\n",
                       best_idx, eps_hat, eps_q);
            }

            memcpy(rx1_prev, rx1_float, sym_len * sizeof(cpx));
            memcpy(rx2_prev, rx2_float, sym_len * sizeof(cpx));

            // =========================================================
            // 4. FFT & GOLDEN DUMP
            // =========================================================
            fft_fixed(rx1_aligned_fix, 0);
            fft_fixed(rx2_aligned_fix, 0);

            // [Golden] Dump FFT Output
            if(do_write) {
            // 格式: RX1_RE RX1_IM RX2_RE RX2_IM (64-bit, 無空白)
                for(int k = 0; k < NFFT; ++k) {
                    fprintf(f_fft, "%04hX%04hX%04hX%04hX\n", 
                        (unsigned short)rx1_aligned_fix[k].re, (unsigned short)rx1_aligned_fix[k].im, 
                        (unsigned short)rx2_aligned_fix[k].re, (unsigned short)rx2_aligned_fix[k].im);
                    fprintf(f_ref_fft, "%04hX%04hX%04hX%04hX\n", 
                        (unsigned short)rx1_aligned_fix[k].re, (unsigned short)rx1_aligned_fix[k].im, 
                        (unsigned short)rx2_aligned_fix[k].re, (unsigned short)rx2_aligned_fix[k].im);
                }
            }

            // =========================================================
            // 5. MIMO ZF & GOLDEN DUMP
            // =========================================================
            for (int k = 0; k < NFFT; ++k) {
                mimo_zf_fixed(rx1_aligned_fix[k], rx2_aligned_fix[k], H_inv_fix, &x1_arr[k], &x2_arr[k]);
                
                // [Golden] Dump ZF Output
                if(do_write) {
                // 格式: RX1_RE RX1_IM RX2_RE RX2_IM (64-bit, 無空白)
                fprintf(f_zf, "%04hX%04hX%04hX%04hX\n", 
                    (unsigned short)x1_arr[k].re, (unsigned short)x1_arr[k].im, 
                    (unsigned short)x2_arr[k].re, (unsigned short)x2_arr[k].im);
                }
            }

            // =========================================================
            // 6. Residual Phase Tracking (CRITICAL LOGIC KEPT)
            // =========================================================
            fix16 phase_err_fix = correct_residual_phase_fixed(x1_arr, x2_arr, NFFT);
            
            // [Golden] Dump Final Output (Corrected)
            if(do_write) {
                for (int k = 0; k < NFFT; ++k) {
                    fprintf(f_pt, "%04hX%04hX%04hX%04hX\n", 
                        (unsigned short)x1_arr[k].re, (unsigned short)x1_arr[k].im, 
                        (unsigned short)x2_arr[k].re, (unsigned short)x2_arr[k].im);
                }
            }

            // --- 重要：這裡必須保留 NCO 更新，不然 BER 會爛掉 ---
            double phase_err = (double)phase_err_fix / ONE_FIX;
            double loop_gain = (sym < 5) ? 1.0 : 0.1;
            nco_phase -= loop_gain * phase_err;
            if (nco_phase > M_PI) nco_phase -= 2.0 * M_PI;
            if (nco_phase < -M_PI) nco_phase += 2.0 * M_PI;

            // =========================================================
            // 7. Demap & Error Count
            // =========================================================
            size_t rx_bit_idx = (sym * 4 * NFFT); 
            for(size_t k=0; k<NFFT; ++k){
                uint8_t b0, b1, b2, b3;
                qpsk_demap_fixed(x1_arr[k], &b0, &b1);
                qpsk_demap_fixed(x2_arr[k], &b2, &b3);
                
                if (sym >= 20){ // Skip Warm-up
                    if(b0 != bits_tx[rx_bit_idx++]) errors++;
                    if(b1 != bits_tx[rx_bit_idx++]) errors++;
                    if(b2 != bits_tx[rx_bit_idx++]) errors++;
                    if(b3 != bits_tx[rx_bit_idx++]) errors++;
                    total_bits_processed += 4;
                }
            }
        }
        printf("%5.1f dB | %.5e   | %s\n", snr, (double)errors/total_bits_processed, 
               capture_enabled_for_snr ? "Golden Vectors Generated!" : "");
    }

    fclose(f_adc); fclose(f_in); fclose(f_fft); fclose(f_zf); fclose(f_pt); fclose(f_h);
    
    free(rx1_float); free(rx2_float);
    free(rx1_prev); free(rx2_prev);
    free(rx1_win);  free(rx2_win);
    free(bits_tx);
    return 0;
}