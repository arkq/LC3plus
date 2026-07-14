/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef SETUP_ENC_LC3_H
#define SETUP_ENC_LC3_H

#include "constants.h"
#include "defines.h"
#include "setup_dec_lc3plus.h"  /* HpvcTreeEnumCfg used below when LL_INCL_HPVC is set */

#ifdef LL_INCL_HPVC
/* per channel setup  for HPVC lossless encoding */
typedef struct
{
    Word16 active_flag;  /*disable for low frame sizes , and some Fs's */
    Word16 mode;        /* -1==TCX_only,  0=early start of HPVC, 1 late start of HPVC */
    Word16 startCoefListNom[2];
    Word16 startCoefNom;     /* nominal start coeff e.g. 14000/50  => 280,   18000/50 => 360  */
    Word16 startCoefList[2];  /* endpoint adjusted adjusted startcoefs */
    Word16 startCoef;        /* selected lastnz and N_signal adjusted start coef */

    Word16 nomTreeLim;               /* nominal number of PVQ trees allowed, WMOPS optimization */
#ifdef HPVC_MAXTREE_LIMIT
    Word16 maxTreeLim;               /* maximum number of PVQ trees allowed, hard limit */
#endif
    Word16* Tx_dec;               /*  vector values [2,8,16,32,64,128 ,] */
    Word16* Tx_splitRule;          /*  vector values [-1, 0,1,2,] */
    HpvcTreeEnumCfg* HpvcTreeEnumCfgPtr;        /* ptr to auxilliary coding_data_hpvc */
} HpvcEncCfg;
#endif /* LL_INCL_HPVC */

/* Channel state and bitrate-derived values go in this struct */
typedef struct
{
#ifdef ENABLE_HR_MODE
    Word32 *stEnc_mdct_mem; /* MDCT_MEM_LEN_MAX */
#else
    Word16 *stEnc_mdct_mem; /* MDCT_MEM_LEN_MAX */
#endif
    Word32 *mdct_mem32;     /* MDCT_MEM_LEN_MAX */
    Word32  targetBitsOff;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 targetBytes;
#else
    Word16 targetBytes;
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 total_bits;
    Word32 targetBitsInit;
    Word32 targetBitsAri;
#else
    Word16 total_bits;
    Word16 targetBitsInit;
    Word16 targetBitsAri;
#endif
    Word16  enable_lpc_weighting;
    Word16  ltpf_enable;
    Word16  quantizedGainOff;
    Word16  tns_bits;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 targetBitsQuant;
#else
    Word16 targetBitsQuant;
#endif
    Word16  olpa_mem_s6k4_exp;
    Word16  olpa_mem_pitch;
    Word16  pitch_flag;
    Word16  ltpf_mem_in_exp;
    Word16  ltpf_mem_normcorr[LEN_MEM_NORMCORR];
    Word16  ltpf_mem_mem_normcorr;
    Word16  ltpf_mem_ltpf_on;
    Word16  ltpf_mem_pitch;
#ifdef FIX_TX_RX_STRUCT_STEREO
    Word16  Tx_ltpf;
#endif
    Word16  mem_targetBits;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 mem_specBits;
#else
    Word16 mem_specBits;
#endif 
    Word16  x_exp;
    Word16  resamp_exp;
    Word16  attack_handling; /* flag to enable attack handling */
    Word16  attdec_filter_mem[2];
    Word16  attdec_detected;
    Word16  attdec_position;
    Word32  attdec_acc_energy;
    Word16  attdec_scaling;
#ifdef ENABLE_HR_MODE
    Word16 regBits;
    Word32 resamp_mem32[120];
#else
    Word32 resamp_mem32[60];
#endif
#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 r12k8_mem_in[240];
#else
    Word16 r12k8_mem_in[120];
#endif
#else
    Word16 r12k8_mem_in[60];
#endif
    Word32 r12k8_mem_50[2];
    Word16    r12k8_mem_out[44];
    Word16 olpa_mem_s12k8[3];
    Word16 olpa_mem_s6k4[LEN_6K4 + MAX_PITCH_6K4 + 16];
    Word16  ltpf_mem_in[LTPF_MEMIN_LEN + LEN_12K8 + 1];
    Word16 n_pccw;
    Word16 n_pc;
    Word16 lfe;
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 quantizedGainOff_ll;
    Word16 scaleSignal_Memory;
#ifdef LL_INCL_HPVC
    HpvcEncCfg hpvcEncCfg;
#endif
#endif
} EncSetup;

/* Constants and sampling rate derived values go in this struct */
struct LC3PLUS_Enc
{
    EncSetup *channel_setup[MAX_CHANNELS];
#ifdef ENABLE_HR_MODE
    const Word32 *W_fx;
#else
    const Word16 *W_fx;
#endif
    const Word16 *bands_offset;

    Word32 fs;           /* encoder sampling rate 44.1 -> 48 */
    Word32 fs_in;        /* input sampling rate */
    Word32 bitrate;      /* global bitrate */
    Word16 fs_idx;       /* sampling rate index */
    Word16 frame_length; /* audio samples / frame */
    Word16 channels;     /* number of channels */
    Word16 epmode;       /* error protection mode */
    LC3PLUS_FrameDuration frame_dms;    /* enum for frame length in steps of 12.5 dms   */
    Word8 lc3_br_set;    /* indicate if bitrate has been set */

    Word16 yLen;
    Word16 W_size;
    Word16 la_zeroes;
    Word16 stEnc_mdct_mem_len;
    Word16 bands_number;
    Word16 nSubdivisions;
    Word16 ltpf_mem_in_len;
    Word16 envelope_bits;
    Word16 global_gain_bits;
    Word16 noise_fac_bits;
    Word16 BW_cutoff_bits;
    Word16 r12k8_mem_in_len;
    Word16 r12k8_mem_out_len;
    Word16 near_nyquist_index;
    Word16 near_nyquist_flag;

    Word16 epmr;
    Word16 combined_channel_coding;
    Word32 bandwidth;
    Word32 bandwidth_preset;
    Word32 bw_ctrl_active;
    Word16 bw_ctrl_cutoff_bin;
    Word16 bw_index;
    Word16 attdec_nblocks;
    Word16 attdec_damping;
    Word16 attdec_hangover_thresh;
    Word16 hrmode;
    Word16 sns_damping;
#ifdef CR9_C_ADD_1p25MS
#ifndef FIX_TX_RX_STRUCT_STEREO
    Word16 Tx_ltpf;
#endif
    Word16 LT_normcorr;
#endif
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 low_band_limit;
    Word16 lossless;
    Word32 ll_ari_bits;
    Word32 ll_ari_bits_lb; // lb: low band
    Word16 ll_tns;
    Word16 ll_est_bit_usage;
    Word16 ll_bit_balance;
    Word16 ll_offQuant;
    Word16 ll_cbr;
    Word16 ll_adap_flag;
    Word16 b_relative;
    Word16 ll_tns_remove;
    Word16 ll_totalBytes;
    Word16 ll_tns_lsb_num_remove_limit;
    Word16 quantizedGainOff_ll;
    Word16 quantizedGainOff_ll_lb; 
    Word16 ll_carryOver;      // switch: carry over unused bytes to higher channel(s)
    Word16 ll_carryOverBytes; // num of bytes to carry over to higher channel(s)
    Word32 totalBytes;
    Word16 padding;
    Word16 wavFormat;
    Word16 scaleSignal;
    Word16 ll_shift;
#endif
  
#ifdef DEBUG
    Word16 max_enc_scratch;
    Word16 max_enc_stack_index;
#endif
    Word16 lc3_scratch_initialized; /* Indicate if max. size has been calculated and scratch allocator has been initialized */
    UWord32 scratch_max_size;       /* Maximum scratch size used throughout encoder */
};

#endif /* SETUP_ENC_LC3_H */
