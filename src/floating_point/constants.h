/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "defines.h"
#include "structs.h"

/* DCT */
extern const Complex dct2_16[16];

/* Ari coder */
extern const LC3_INT16 ari_tns_order_cf[2][9];
extern const LC3_INT16 ari_tns_freq_cf[8][18];
extern const LC3_INT ari_spec_lookup_fl[4096];
extern const LC3_INT16 ari_spec_cumfreq_fl[64][18];
extern const LC3_INT ari_spec_bits_fl[64][17];

/* SNS */
extern const LC3_FLOAT sns_W[6];
extern const LC3_FLOAT *sns_preemph_all[6];
extern const LC3_FLOAT sns_LFCB[8][32];
extern const LC3_FLOAT sns_HFCB[8][32];
extern const LC3_INT   pvq_enc_A[16][11];
extern const LC3_FLOAT idct_lookup[M][M];

/* 12.8 kHz resampler */
extern const LC3_FLOAT lp_scale_factors[6];

extern const LC3_INT32 resamp_params[][4];
extern const LC3_FLOAT *lp_filter[6];
extern const LC3_FLOAT    highpass50_filt_b[3];
extern const LC3_FLOAT    highpass50_filt_a[3];
extern const LC3_INT   up_fac[6];

/* TNS */
extern const LC3_FLOAT quants_pts_tns[17];
extern const LC3_INT   huff_bits_tns[8][17];
extern const LC3_INT   order1_tns[8];
extern const LC3_INT   order2_tns[8];
extern const LC3_FLOAT lagw_tns[9];
extern const LC3_FLOAT quants_pts_tns[17];
extern const LC3_FLOAT quants_thr_tns[18];

/* SNS */
extern const LC3_FLOAT sns_vq_far_adj_gains_fl[8];
extern const LC3_FLOAT sns_vq_near_adj_gains_fl[4];
extern const LC3_FLOAT sns_vq_reg_lf_adj_gains_fl[4];
extern const LC3_FLOAT q_g_sns[6];
extern const LC3_FLOAT sns_vq_reg_adj_gains_fl[2];
extern const LC3_FLOAT sns_dec_gains[4][8];

/* Global Gain */
extern const LC3_INT   gg_p1[6];
extern const LC3_INT   gg_p2[6];
extern const LC3_INT   gg_p3[6];
extern const LC3_FLOAT gg_c[6];
extern const LC3_FLOAT gg_d[6];

/* Olpa */
extern const LC3_FLOAT olpa_down2[5];
extern const LC3_FLOAT olpa_acw[98];

/* LTPF */
extern const LC3_FLOAT conf_inter_filter_48[4][12];
extern const LC3_FLOAT conf_inter_filter_32[4][8];
extern const LC3_FLOAT conf_inter_filter_24[4][6];
extern const LC3_FLOAT conf_inter_filter_16[4][4];
extern const LC3_FLOAT conf_tilt_filter_48[4][11];
extern const LC3_FLOAT conf_tilt_filter_32[4][7];
extern const LC3_FLOAT conf_tilt_filter_24[4][5];
extern const LC3_FLOAT conf_tilt_filter_16[4][3];
extern const LC3_FLOAT inter4_1[33];
extern const LC3_FLOAT enc_inter_filter[4][4];

/* Bandwidth Detector */
extern const LC3_INT threshold_quiet[4];
extern const LC3_INT threshold_brickwall[4];
extern const LC3_INT  brickwall_dist[4];
extern const LC3_INT  BW_warp_idx_start_16k[4];
extern const LC3_INT  BW_warp_idx_stop_16k[4];
extern const LC3_INT  BW_warp_idx_start_24k[4];
extern const LC3_INT  BW_warp_idx_stop_24k[4];
extern const LC3_INT  BW_warp_idx_start_32k[4];
extern const LC3_INT  BW_warp_idx_stop_32k[4];
extern const LC3_INT  BW_warp_idx_start_48k[4];
extern const LC3_INT  BW_warp_idx_stop_48k[4];
extern const LC3_INT* BW_warp_idx_start_all[4];
extern const LC3_INT* BW_warp_idx_stop_all[4];

extern const LC3_INT  BW_warp_idx_start_16k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_16k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_24k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_24k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_32k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_32k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_start_48k_2_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_48k_2_5ms[4];
extern const LC3_INT* BW_warp_idx_start_all_2_5ms[4];
extern const LC3_INT* BW_warp_idx_stop_all_2_5ms[4];
extern const LC3_INT BW_cutoff_bin_all_2_5ms_HR[MAX_BW_BANDS_NUMBER];
extern const LC3_INT bands_number_2_5ms_HR[6];

extern const LC3_INT BW_cutoff_bin_all_2_5ms[MAX_BW_BANDS_NUMBER];
extern const LC3_INT bands_number_2_5ms[5];

extern const LC3_INT  BW_warp_idx_start_16k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_16k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_24k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_24k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_32k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_32k_5ms[4];
extern const LC3_INT  BW_warp_idx_start_48k_5ms[4];
extern const LC3_INT  BW_warp_idx_stop_48k_5ms[4];
extern const LC3_INT* BW_warp_idx_start_all_5ms[4];
extern const LC3_INT* BW_warp_idx_stop_all_5ms[4];
extern const LC3_INT  BW_cutoff_bin_all_5ms[MAX_BW_BANDS_NUMBER];
extern const LC3_INT bands_number_5ms[6];
extern const LC3_INT BW_cutoff_bin_all_HR[MAX_BW_BANDS_NUMBER];
extern const LC3_INT BW_cutoff_bin_all_5ms_HR[MAX_BW_BANDS_NUMBER];
extern const LC3_INT BW_cutoff_bin_all[MAX_BW_BANDS_NUMBER];
extern const LC3_INT BW_cutoff_bits_all[MAX_BW_BANDS_NUMBER];

extern const LC3_INT BW_cutoff_bin_all_7_5ms[MAX_BW_BANDS_NUMBER];
extern const LC3_INT bands_number_7_5ms[6];
extern const LC3_INT bands_number_7_5ms_HR[6];
extern const LC3_INT* BW_warp_idx_start_all_7_5ms[4];
extern const LC3_INT* BW_warp_idx_stop_all_7_5ms[4];
extern const LC3_INT brickwall_dist_7_5ms[4];
extern const LC3_INT* ACC_COEFF_PER_BAND_PLC_7_5ms[];

/* Arithmetic coding */
extern const LC3_INT tns_cf[8][18];
extern const LC3_INT tns_freq_cf[2][9];

/* MDCT Windows */
extern const LC3_FLOAT MDCT_WINDOW_80[160];
extern const LC3_FLOAT MDCT_WINDOW_160[320];
extern const LC3_FLOAT MDCT_WINDOW_240[480];
extern const LC3_FLOAT MDCT_WINDOW_320[640];
extern const LC3_FLOAT MDCT_WINDOW_480[960];
extern const LC3_FLOAT MDCT_WINDOW_960[1920];
extern const LC3_FLOAT* MDCT_WINS_10ms[2][6];
extern const LC3_INT    MDCT_la_zeroes[6];

extern const LC3_FLOAT MDCT_WINDOW_80_2_5ms[40];
extern const LC3_FLOAT MDCT_WINDOW_160_2_5ms[80];
extern const LC3_FLOAT MDCT_WINDOW_240_2_5ms[120];
extern const LC3_FLOAT MDCT_WINDOW_320_2_5ms[160];
extern const LC3_FLOAT MDCT_WINDOW_480_2_5ms[240];
extern const LC3_FLOAT* MDCT_WINS_2_5ms[2][6];
extern const LC3_INT    MDCT_la_zeroes_2_5ms[6];

extern const LC3_FLOAT MDCT_WINDOW_80_5ms[80];
extern const LC3_FLOAT MDCT_WINDOW_160_5ms[160];
extern const LC3_FLOAT MDCT_WINDOW_240_5ms[240];
extern const LC3_FLOAT MDCT_WINDOW_320_5ms[320];
extern const LC3_FLOAT MDCT_WINDOW_480_5ms[480];
extern const LC3_FLOAT* MDCT_WINS_5ms[2][6];
extern const LC3_INT    MDCT_la_zeroes_5ms[6];

extern const LC3_FLOAT* MDCT_WINS_7_5ms[2][6];
extern const LC3_INT32 MDCT_la_zeroes_7_5ms[6];

extern const LC3_INT MDCT_WINDOWS_LENGTHS_10ms[6];
extern const LC3_INT MDCT_WINDOWS_LENGTHS_7_5ms[6];
extern const LC3_INT MDCT_WINDOWS_LENGTHS_5ms[6];
extern const LC3_INT MDCT_WINDOWS_LENGTHS_2_5ms[6];

/* Per band energy */
extern const LC3_INT* ACC_COEFF_PER_BAND[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_HR[6];

extern const LC3_INT* ACC_COEFF_PER_BAND_2_5ms_HR[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_2_5ms[5];

extern const LC3_INT* ACC_COEFF_PER_BAND_7_5ms_HR[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_7_5ms[5];

extern const LC3_INT* ACC_COEFF_PER_BAND_5ms_HR[6];
extern const LC3_INT* ACC_COEFF_PER_BAND_5ms[5];

/* Near Nyquist detector */
extern const LC3_INT NN_thresh;
/* Tone detector */
extern const LC3_INT32 TD_HR_thresh_10ms;
extern const LC3_INT32 TD_HR_thresh_7_5ms;
extern const LC3_INT32 TD_HR_thresh_5ms;
extern const LC3_INT32 TD_HR_thresh_2_5ms; 

extern const LC3_INT32 xavg_N_grp[5];
extern const LC3_FLOAT *hannOla_wins[5];
extern const LC3_INT32 gwlpr[MAX_LGW+1];
extern const LC3_INT16 fade_scheme_tab[24 / 2][3];
extern const LC3_FLOAT scATHFx[MAX_LGW - 2];
extern const LC3_INT32 mdct_grp_bins[10];
extern const LC3_FLOAT* PhECU_whr16ms_wins[5];

extern const LC3_FLOAT plc_preemph_fac[];
extern const LC3_INT* ACC_COEFF_PER_BAND_PLC[];
extern const LC3_INT* ACC_COEFF_PER_BAND_PLC_2_5ms[];
extern const LC3_INT* ACC_COEFF_PER_BAND_PLC_5ms[];
extern const LC3_FLOAT *plc_tdc_lpc_all[6];
extern const LC3_FLOAT plc_tdc_lpc_8[17];
extern const LC3_FLOAT plc_tdc_lpc_16[17];
extern const LC3_FLOAT plc_tdc_lpc_24[17];
extern const LC3_FLOAT plc_tdc_lpc_32[17];
extern const LC3_FLOAT plc_tdc_lpc_48[17];
extern const LC3_FLOAT plc_tdc_lpc_96[17];
extern const LC3_FLOAT plc_tdc_lpc_8_25ms[9];

extern const LC3_INT16 plc_fadeout_param_maxlen[4];
extern const LC3_INT16 plc_fadeout_param_maxbytes[4];
extern const LC3_INT16 PLC_FADEOUT_TYPE_2_SELECTOR;
#endif /* CONSTANTS_H */
