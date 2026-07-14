/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_DEC_LC3_H
#define SETUP_DEC_LC3_H

#include "constants.h"
#include "defines.h"

#ifdef LL_INCL_HPVC
typedef struct
{
    Word16 Np;           /*one of 8,16,32,64,128 */
    Word16 Kp;           /*value := 0..36  */
    Word16 NsSafe;           /* {-1, 0,}  1...16  */
    Word16 NsHdrSafe;        /*  {0, 1, 3} */

    Word16 splitRule;
    Word16 Ns;           /* {-1, 0,}  1...16  */
    Word16 NsHdr;        /* 0, 1, 3  */

    /*************************************************************/
    Word16 LS;              /* leading sign 1= negative, 0 positive */

    Word16 hdrNdbg;           /*      1-16  or   3 */
    Word16 hdrKdbg;           /*      always Kp ? */
    Word32 hdrSz;
    Word32 hdrIdx;

    Word16 splitHdrNsDbg[3];
    Word16 splitHdrKsDbg[3];
    Word32 splitHdrLeafSz[3];
    Word32 splitHdrLeafIdx[3];

    Word32 flatLeafSz[16];   /*  N_MPVQ = Func(Nleaf,Kleaf)   */
    Word32 flatLeafIdx[16];  /*   only Ns in use */


    Word16 start_coeff_nb;   /*   start coeff for debugging */
    Word32 *Xqm_ptr;        /*  ptr to the start of the integers to be encoded using a HPVC tree, mainly used for verifcation */
    Word16 xDbg[LL_HPVC_NP_MAX]; /*input  to enum */
    Word16 xDeEnumDbg[LL_HPVC_NP_MAX]; /* output from deenum Stepd */
} HpvcTreeEnumCfg;

/* per channel setup  for HPVC lossless decoding */
typedef struct
{
    Word16 active_flag;
    Word16 mode;             /* -1==TCX_only,  0=early start of HPVC, 1 late start of HPVC */
    Word16 startCoefListNom[2];
    Word16 startCoefNom;         /* e.g. 14000/50  => 280,   18000/50 => 360  */
    Word16 startCoefList[2];  /* endpoint adjusted adjusted startcoefs */
    Word16 startCoef;         /*  selected lastnz and Nsignal adjusted startcoef  */

    Word16 nomTreeLim;             /* nominal number of HPVC trees allowed, WMOPS optimization */
#ifdef HPVC_MAXTREE_LIMIT
    Word16 maxTreeLim;             /* maximum number of HPVC trees allowed, hard limit */
#endif
    HpvcTreeEnumCfg* HpvcTreeEnumCfgPtr;  /* ptr to auxilliary singel tree   */
} HpvcDecCfg;
#endif /* LL_INCL_HPVC */

typedef struct
{
    Word16 *x_old_tot_fx;      /* MAX_LEN_PCM_PLC_TOT    */
    Word32 *PhECU_f0est;       /* MAX_PLOCS            interpolated plocs  */
    Word16 *PhECU_xfp_fx;      /* MAX_LPROT */
    Word16 *PhECU_X_sav_fx;    /* MAX_LPROT */
    Word16 *PhECU_plocs;       /* MAX_PLOCS */
    Word16 *PhECU_fg_wintaper; /* MDCT_MEM_LEN_MAX */
    Word16 *PhECU_win_pre_tda; /* MAX_WIN_PRE_TDA */
    Word32  tdc_gain_c;
    Word16  stab_fac;
    Word16  tdc_fract;
    Word16  tdc_seed;
    Word16  tdc_preemph_fac;
    Word16  tdc_lpc_order;
    Word16  cum_fflcAtten;
    Word16  harmonicBuf_fx[MAX_PITCH];
    Word16  harmonicBuf_Q;
    Word16  synthHist_fx[M];
    Word16  cum_fading_slow;
    Word16  cum_fading_fast;
    Word16  PhECU_LprotOrg_fx; /* needed to change the Prot size  adaptively  */
    Word16  PhECU_Lprot_fx;
    Word16  PhECU_fs_idx_fx;
    Word16  PhECU_frame_ms;  /* needed in PLC_Update and PLCMain functons*/
    Word16  PhECU_seed_fx;
    Word16  PhECU_xfp_exp_fx;
    Word16  PhECU_time_offs;
    Word16  PhECU_X_savQ_fx;
    Word16  PhECU_num_plocs;
    Word16  PhECU_f0hzLtpBinQ7; /*  ltp F0 in bins  if available  */
    Word16  PhECU_short_flag_prev;
    Word16  PhECU_whr_tot_taper;
    Word16  PhECU_whr_tot_flat;
    Word16  PhECU_LDWIN_OLAP;
    Word16  PhECU_LA;
    Word16  PhECU_t_adv;
    Word16  PhECU_beta_mute;
    Word16  norm_corrQ15_fx;
    Word16  q_fx_old_exp;
    Word16  max_len_pcm_plc;
    Word16  max_lprot;
    Word16  max_plocs;

    /* Word32 L_tot W_energy sum exponent */
    Word16  PhECU_oold_Ltot_exp_fx;
    Word16  PhECU_old_Ltot_exp_fx;
    Word32  PhECU_L_oold_xfp_w_E_fx;
    Word32  PhECU_L_old_xfp_w_E_fx;
    Word16  PhECU_oold_xfp_w_E_exp_fx;   /* input Word16 xfp exponnet  */
    Word16  PhECU_old_xfp_w_E_exp_fx;
    Word16  PhECU_oold_grp_shape_fx[MAX_LGW];
    Word16  PhECU_old_grp_shape_fx[MAX_LGW];
    Word16  PhECU_margin_xfp;
	Word16  PhECU_nonpure_tone_flag;         /*   non-pure single tone indicator state */
    Word16  PhECU_mag_chg_1st[MAX_LGW];
    Word16  PhECU_Xavg[MAX_LGW];
    Word16  old_scf_q[M];
    Word16  old_old_scf_q[M];
    Word16  tdc_A[M + 1];
    /* for now 20 ms saved Q14  or ptr to a combined ifft win and MDCT  preTDA synthesis window  16  ms */

    Word16 longterm_counter_plcTdc;
    Word16 longterm_counter_plcNsAdv;
    Word16 longterm_analysis_counter_max;  /* Maximum longterm frames number */
    Word16 longterm_analysis_counter_max_bytebuffer;  /* Same as above but reduced for circular bit-buffer */
    Word32 *plc_longterm_advc_tdc;
    Word32 *plc_longterm_advc_ns;
    UWord8 plc_fadeout_type;
    Word16 overall_counter;
    Word8  longterm_counter_byte_position;
    Word8  longterm_counter_bit_position;
#ifdef CR13_C_RESET_CLASSIFIER_AFTER_BAD_FRAMES
    Word16 numberOfGoodFrames;
#endif
} AplcSetup;

/* Channel state and bitrate-derived values go in this struct */
typedef struct
{
    Word16 *ltpf_mem_x;       /* LTPF_MEM_X_LEN */
    Word16 *ltpf_mem_y;       /* LTPF_MEM_Y_LEN */
#ifdef ENABLE_HR_MODE
    Word32 *stDec_ola_mem_fx; /* MDCT_MEM_LEN_MAX */
#else
    Word16 *stDec_ola_mem_fx; /* MDCT_MEM_LEN_MAX */
#endif
    AplcSetup *plcAd;
    Word16 *   q_old_d_fx; /* MAX_BW */
    Word16     q_old_fx_exp;
    Word16     ns_seed;
    Word16     ns_cum_alpha;
    Word16  pc_nbLostFramesInRow;
    Word16  pc_seed;
    Word16 *q_old_res_fx;
    Word16  q_old_res_fx_exp;
    Word16  prev_gg;
    Word16  prev_gg_e;
    Word16  prev_BW_cutoff_idx_nf;
    Word16 prev_fac_ns_fx;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 total_bits;
#else
    Word16 total_bits;
#endif
    Word16 enable_lpc_weighting;
    Word16 stDec_ola_mem_fx_exp;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 targetBytes;
#else
    Word16 targetBytes;
#endif
    Word16 ltpf_mem_e;
    Word16 ltpf_mem_pitch_int;
    Word16 ltpf_mem_pitch_fr;
    Word16 ltpf_mem_gain;
    Word16 ltpf_mem_active;
    Word16 ltpf_scale_fac_idx;
    Word16 ltpf_mem_scale_fac_idx;
    Word16 quantizedGainOff;
    Word16 prev_bfi;
    Word16 prev_prev_bfi;
    Word16 concealMethod;
    Word16 nbLostFramesInRow;
    Word16 plc_damping;
    Word16 last_size;
    Word32 rel_pitch_change;
#ifdef CR9_C_ADD_1p25MS
#ifdef FIX_TX_RX_STRUCT_STEREO
    Word16 ltpf_rx_status[2];
#endif
    Word16 ltpf_mem_continuation;
    Word16 ltpf_mem_active_prev;
    Word16 ltpf_mem_pitch_int_prev;
    Word16 ltpf_mem_pitch_fr_prev;
    Word16 ltpf_mem_beta_idx_prev;
    Word16 ltpf_mem_gain_prev;
    Word16 ltpf_pitch_stability_counter;
#ifdef NEW_SIGNALLING_SCHEME_1p25
    Word16 ltpfinfo_frame_cntr_fx;  /* individual cntr for each channel*/
#endif
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 quantizedGainOff_ll;
#ifdef LL_INCL_HPVC
    HpvcDecCfg hpvcDecCfg;
#endif
    Word16 ll_offQuant;
    Word16 ll_ari_bits;
    Word16 ll_adap_prev;
    Word16 ScaleSignal_Memory;
#endif
} DecSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3PLUS_Dec
{
    DecSetup *    channel_setup[MAX_CHANNELS];
#ifdef ENABLE_HR_MODE
    const Word32 *W_fx;
#else
    const Word16 *W_fx;
#endif
    const Word16 *bands_offset;
    Word32        fs;           /* sampling rate, 44.1 maps to 48 */
    Word32        fs_out;       /* output sampling rate */
    Word16        fs_idx;       /* sampling rate index */
    Word16        frame_length; /* sampling rate index */
    Word16        channels;     /* number of channels */
    Word16        plcMeth;      /* PLC method for all channels */
    LC3PLUS_FrameDuration        frame_dms;    /* frame length in dms (decimilliseconds, 10^-4)*/
    Word16        last_size;    /* size of last frame, without error protection */
    Word16        ep_enabled;   /* error protection enabled */
    Word16        error_report; /* corrected errors in last frame or -1 on error */

    Word16 n_pccw;
    Word16 be_bp_left;
    Word16 be_bp_right;
    Word16 n_pc;
    Word16 m_fec;
    Word16 epmr;
    Word16 combined_channel_coding;

    Word16 yLen;
    Word16 W_size;
    Word16 la_zeroes;
    Word16 stDec_ola_mem_fx_len;
    Word16 bands_number;
    Word16 ltpf_mem_x_len;
    Word16 ltpf_mem_y_len;
    Word16 BW_cutoff_bits;
    Word16 hrmode;
    Word16 alpha_type_2_table[160];/* PLC_FADEOUT_TYPE_1_IN_MS*100/125 */
#ifndef FIX_TX_RX_STRUCT_STEREO
#ifdef CR9_C_ADD_1p25MS
    Word16           ltpf_rx_status[2];
    Word16           ltpf_mem_continuation;
    Word16           ltpf_mem_active_prev;
    Word16           ltpf_mem_pitch_int_prev;
    Word16           ltpf_mem_pitch_fr_prev;
    Word16           ltpf_mem_beta_idx_prev;
    Word16           ltpf_mem_gain_prev;
    Word16           ltpf_pitch_stability_counter;
#endif
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 low_band_limit;
    Word16 lossless;
    Word16 ll_tns;
    Word16 ll_tns_remove;
    Word16 ll_adap_flag;
    Word16 wavFormat;
#endif
  
#ifdef DEBUG
    Word16 max_dec_scratch;
    Word16 max_dec_stack_index;
#endif
    Word16 lc3_scratch_initialized; /* Indicate if max. size has been calculated and scratch allocator has been initialized */
    UWord32 scratch_max_size;       /* Maximum scratch size used throughout decoder */
};

#endif /* SETUP_DEC_LC3_H */
