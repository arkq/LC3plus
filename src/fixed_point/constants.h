/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "defines.h"

#include "basop_util.h"


#ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 BW_cutoff_bin_all_HR[];

extern RAM_ALIGN const Word32 *const LowDelayShapes_n960_HRA_2_5ms[2];
extern RAM_ALIGN const Word32 *const LowDelayShapes_n960_HRA_5ms[2];
extern RAM_ALIGN const Word32 *const LowDelayShapes_n960_HRA[2];
#endif

#  ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len[6];
#  else
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len[5];
#  endif
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes[NUM_SAMP_FREQ];
#    ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word32 *const LowDelayShapes_n960[6];
#    else
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960[6];
#    endif

#  ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_5ms[6];
#  else
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_5ms[5];
#  endif
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes_5ms[NUM_SAMP_FREQ];
#    ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word32 *const LowDelayShapes_n960_5ms[6];
#    else
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960_5ms[6];
#    endif
#ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_2_5ms[6];
#else
extern RAM_ALIGN const Word16 LowDelayShapes_n960_len_2_5ms[5];
#endif
extern RAM_ALIGN const Word16 LowDelayShapes_n960_la_zeroes_2_5ms[NUM_SAMP_FREQ];
#    ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word32 *const LowDelayShapes_n960_2_5ms[6];
#    else
extern RAM_ALIGN const Word16 *const LowDelayShapes_n960_2_5ms[6];
#    endif

extern RAM_ALIGN const Word16 NN_thresh;
extern RAM_ALIGN const Word16 NN_thresh_exp;
extern RAM_ALIGN const Word32 BW_thresh_quiet[4];
extern RAM_ALIGN const Word16 BW_thresh_quiet_exp;
extern RAM_ALIGN const Word16 BW_thresh_brickwall[4];

extern RAM_ALIGN const Word16 BW_cutoff_bin_all[];
extern RAM_ALIGN const Word16 BW_cutoff_bits_all[];

#ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 bands_number_2_5ms_HR[];
#endif

extern RAM_ALIGN const Word16 BW_brickwall_dist[4];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all[MAX_BW_BANDS_NUMBER - 1];

extern RAM_ALIGN const Word16 BW_brickwall_dist_5ms[4];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 BW_brickwall_dist_2_5ms[4];
extern RAM_ALIGN const Word16 *const BW_warp_idx_start_all_2_5ms[MAX_BW_BANDS_NUMBER - 1];
extern RAM_ALIGN const Word16 *const BW_warp_idx_stop_all_2_5ms[MAX_BW_BANDS_NUMBER - 1];

extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq[MAX_BW_BANDS_NUMBER];

#    ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_HR[2];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_HR[2];

extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_5ms_HR[2];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_5ms_HR[2];

extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_2_5ms_HR[2];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_2_5ms_HR[2];
#    endif

extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_5ms[MAX_BW_BANDS_NUMBER];


extern RAM_ALIGN const Word16 *const tns_subdiv_startfreq_2_5ms[MAX_BW_BANDS_NUMBER];
extern RAM_ALIGN const Word16 *const tns_subdiv_stopfreq_2_5ms[MAX_BW_BANDS_NUMBER];


extern RAM_ALIGN const Word16 Tab_esc_nb[4];

extern RAM_ALIGN const Word8 ari_spec_lookup[4096];
extern RAM_ALIGN const UWord16 ari_spec_cumfreq[64][17];
extern RAM_ALIGN const UWord16 ari_spec_freq[64][17];
extern RAM_ALIGN const UWord16 ari_spec_bits[64][17];

extern RAM_ALIGN const Word32 tnsAcfWindow[MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_bits[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_freq[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_order_cumfreq[2][MAXLAG];
extern RAM_ALIGN const Word16 ac_tns_coef_bits[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 ac_tns_coef_freq[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 ac_tns_coef_cumfreq[MAXLAG][TNS_COEF_RES];
extern RAM_ALIGN const Word16 tnsQuantPts[TNS_COEF_RES];
extern RAM_ALIGN const Word16 tnsQuantThr[TNS_COEF_RES - 1];

extern RAM_ALIGN const Word16 *const lpc_pre_emphasis[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_pre_emphasis_e[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_2_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_lin_pre_emphasis_e_2_5ms[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 *const lpc_warp_dee_emphasis[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const lpc_warp_dee_emphasis_e[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 bands_nrg_scale[32];

#  ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word16 *const bands_offset_2_5ms_HR[2];
extern RAM_ALIGN const Word16 *const bands_offset_5ms_HR[2];
extern RAM_ALIGN const Word16 *const bands_offset_HR[2];
#  endif

extern RAM_ALIGN const Word16 *const bands_offset[6];
extern RAM_ALIGN const Word16 bands_offset_with_one_max[NUM_OFFSETS];
extern RAM_ALIGN const Word16 bands_offset_with_two_max[NUM_OFFSETS];

extern RAM_ALIGN const Word16 bands_number_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const bands_offset_5ms[6];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_5ms[NUM_OFFSETS];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_5ms[NUM_OFFSETS];
extern RAM_ALIGN const Word16 bands_number_2_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const bands_offset_2_5ms[6];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_2_5ms[NUM_OFFSETS];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_2_5ms[NUM_OFFSETS];


extern RAM_ALIGN const Word16 pitch_max[5];
extern RAM_ALIGN const Word16 plc_preemph_fac[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 TDC_high_16[11];
extern RAM_ALIGN const Word16 TDC_high_32[11];
extern RAM_ALIGN const Word16 TDC_high_16_harm[11];
extern RAM_ALIGN const Word16 TDC_high_32_harm[11];

extern RAM_ALIGN const Word32 *const lag_win[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 *const bands_offset_lin[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const bands_offset_lin_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 *const bands_offset_lin_2_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_one_max_lin_2_5ms[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 bands_offset_with_two_max_lin_2_5ms[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word32 inv_odft_twiddle_80_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_80_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_60_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_60_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_40_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_40_im[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_20_re[M];
extern RAM_ALIGN const Word32 inv_odft_twiddle_20_im[M];

#ifdef SUBSET_WB
extern RAM_ALIGN const Word16 resamp_filt_16k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_16k[1];
#endif
#ifdef SUBSET_SSWB
extern RAM_ALIGN const Word16 resamp_filt_24k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_24k[1];
#endif
#ifdef SUBSET_SWB
extern RAM_ALIGN const Word16 resamp_filt_32k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_32k[1];
#endif
#ifdef SUBSET_FB
extern RAM_ALIGN const Word16 resamp_filt_48k[240];
#else
extern RAM_ALIGN const Word16 resamp_filt_48k[1];
#endif

extern RAM_ALIGN const Word16 resamp_params[NUM_SAMP_FREQ][4];
extern RAM_ALIGN const Word16 *const resamp_filts[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 highpass50_filt_num[3];
extern RAM_ALIGN const Word16 highpass50_filt_den[2];

extern RAM_ALIGN const Word16 olpa_ac_weighting[98];

extern RAM_ALIGN const Word16 ltpf_ac_interp_filt[7][9];
extern RAM_ALIGN const Word16 inter_filter[5][4][12];
extern RAM_ALIGN const Word16 inter_filter_shift[5];
extern RAM_ALIGN const Word16 inter_filter_len[5];
extern RAM_ALIGN const Word16 tilt_filter[5][4][11];
extern RAM_ALIGN const Word16 tilt_filter_len[5];
extern RAM_ALIGN const Word16 gain_scale_fac[4];
extern RAM_ALIGN const UWord16 pitch_scale[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 ltpf_overlap_len[NUM_SAMP_FREQ];

extern RAM_ALIGN const Word16 sns_vq_reg_adj_scf[2];
extern RAM_ALIGN const Word16 sns_vq_reg_lf_adj_scf[4];
extern RAM_ALIGN const Word16 sns_vq_near_adj_scf[4];
extern RAM_ALIGN const Word16 sns_vq_far_adj_scf[8];
extern RAM_ALIGN const Word16 *const sns_gaintabPtr[4];
extern RAM_ALIGN const Word16 sns_gainSz[4];
extern RAM_ALIGN const Word16 sns_gainMSBbits[4];
extern RAM_ALIGN const Word16 sns_gainLSBbits[4];
extern RAM_ALIGN const Word16 sns_Kval[4][2];
extern RAM_ALIGN const UWord32 sns_MPVQ_Sz[4][2];

extern RAM_ALIGN const Word16 st1SCF0_7_base5_32x8_Q11[256];
extern RAM_ALIGN const Word16 st1SCF8_15_base5_32x8_Q11[256];

#ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word32 st1SCF0_7_base5_32x8_Q27[256];
extern RAM_ALIGN const Word32 st1SCF8_15_base5_32x8_Q27[256];
#endif

/* PVQ deindexing tables */
extern RAM_ALIGN const UWord32 h_memN16K12[12 + 2];
extern RAM_ALIGN const UWord32 h_memN10K22[22 + 2];
extern RAM_ALIGN const UWord32 h_memN6K2[2 + 2];
extern RAM_ALIGN const Word16 tabledKMAX[16 + 1];
extern RAM_ALIGN const UWord32 *const MPVQ_offs_ptr[16 + 1];

extern RAM_ALIGN const Word16 isqrt_Q16tab[1 + SQRT_EN_MAX_FX];
#ifdef ENABLE_HR_MODE
extern RAM_ALIGN const Word32 isqrt_Q31tab[1 + SQRT_EN_MAX_FX];
#endif

extern RAM_ALIGN const Word16 adjust_global_gain_tables[5][NUM_SAMP_FREQ];


extern RAM_ALIGN const Word16 sqrt_table_phecu[];
extern RAM_ALIGN const Word16 POW_ATT_TABLE0[];
extern RAM_ALIGN const Word16 POW_ATT_TABLE1[];
#ifdef PLC2_FADEOUT_IN_MS
#if PLC2_FADEOUT_IN_MS == 0
extern RAM_ALIGN const Word16 *const POW_ATT_TABLES[3];
#else
extern RAM_ALIGN const Word16 *const POW_ATT_TABLES[11];
#endif
#else
extern RAM_ALIGN const Word16 *const POW_ATT_TABLES[3];
#endif

extern RAM_ALIGN const Word16 e_tot_headroom[];
extern RAM_ALIGN const Word16 xfp_wE_MDCT2FFTQ11[];

extern RAM_ALIGN const Word16 num_FsByResQ0[5];
extern RAM_ALIGN const Word16 *const LprotSzPtr;
extern RAM_ALIGN const Word16 InvLprot_Q22[5];
extern RAM_ALIGN const Word16 PhEcuFftScale[5];
extern RAM_ALIGN const    Word16  oneOverFrameQ15Tab[5];
extern RAM_ALIGN const Word16 PhEcu_Xsav_Flt2FxDnShift[];
extern RAM_ALIGN const Word16 PhEcu_Xsav_Flt2FxScaleQ15[];
extern RAM_ALIGN const Word16 PhEcu_frac_thr_rise_lin_Q15[];
extern RAM_ALIGN const Word16 PhEcu_frac_thr_decay_lin_Q15[];

extern RAM_ALIGN const Word16 mdct_grp_bins_fx[];
extern RAM_ALIGN const Word16 xavg_N_grp_fx[];
extern RAM_ALIGN const Word16 spec_shape_headroom[];
extern RAM_ALIGN const Word16 rectLengthTab[NUM_SAMP_FREQ];
extern RAM_ALIGN const Word16 hamm_len2Tab[];

extern RAM_ALIGN const Word16 gw_len_inv_shift_fx[];
extern RAM_ALIGN const Word16 gwlpr_fx[];

extern RAM_ALIGN const Word16 sin_quarterQ15_fx[];
extern RAM_ALIGN const Word16 sincos_lowres_tab_sinQ15_fx[];

extern RAM_ALIGN const Word16 *const PhECU_wins[5][3];

extern RAM_ALIGN const Word16 *const w_new[];
extern RAM_ALIGN const Word16 *const w_old[];

/* extern    RAM_ALIGN const  Word16 WORK_LEN[]; */
extern RAM_ALIGN const Word16 COPY_LEN[];
extern RAM_ALIGN const Word16 OLA_LEN[];


#endif
