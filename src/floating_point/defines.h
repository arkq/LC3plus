/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef DEFINES_H
#define DEFINES_H

#include "stdint.h"


typedef float    LC3_FLOAT;
typedef int32_t  LC3_INT;
typedef int16_t  LC3_INT16;
typedef uint16_t LC3_UINT16;
typedef short    LC3_SHORT;
typedef uint8_t  LC3_UINT8;
typedef int8_t   LC3_INT8;
typedef uint32_t LC3_UINT32;

/* Release defines */
#define ENABLE_ADVANCED_PLC_FL
#define ENABLE_ADVANCED_PLC_FL_DEFAULT
#define ENABLE_BW_CONTROLLER
#define ENABLE_HR_MODE_FL
#define ENABLE_PADDING
#define ENABLE_RFRAME_FL
#define ENABLE_PLC
/* flags */
#define ENABLE_PLC_MODE_FLAG
#define ENABLE_BANDWIDTH_FLAG
#define ENABLE_EP_MODE_FLAG
#define ENABLE_FRAME_MS_FLAG
#define ENABLE_HR_MODE_FL_FLAG

#ifndef NO_POST_REL_CHANGES
/* Post-release non-bitexact changes */

#define CR13_B_FIX_PC_BINS
#define CR13_C_RESET_CLASSIFIER_AFTER_BAD_FRAMES
#define CR12_D_FIX_BITRATE_LIMITS

#define CR9_C_ADD_1p25MS
#ifdef CR9_C_ADD_1p25MS
#  define ENABLE_12p5_DMS_MODE
#  define CR9_C_ADD_1p25MS_LRSNS
#  define FIX_LTPF_PITCH_1p25

#endif

#endif /* NO_POST_REL_CHANGES */

#ifdef CR9_C_ADD_1p25MS_LRSNS
#   define LRSNS_PC_SIGNAL_FIX              /*correct handling of incoming bfi==2 to DEC_ENTROPY for conformance */
#define LRSNS_10MS_BFISIGNAL_FIX          /* correct signaling of detected BER in SNS/LRSNS */
#define LRSNS_WMC_FIX
/*#define LRSNS_CBC_NO_LTPF_DEPENDENCY  */    /* turn off LRSNS CB_C dependency on LTPF activation flag */   
 
#define SNSLR_N_FIXENV   4              /* 4 fix envelopes multiplexed inside  the full codeword */
#define SNSLR_N_FIXENV_SHIFTS   4       /* 2 bits */
#define SNSLR_MAX_PVQ_CAND  6           /* splitLF(0), full(1), fixed_env 2+{0,1,2,3 },  */
#define SNSLR_MAX_PVQ_SEARCH_CAND (SNSLR_MAX_PVQ_CAND-SNSLR_N_FIXENV+1)     /* 3 = splitLF(0), full(1), fixed_envs(2), */



#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
#    define LC3_CONST_FLOATMAX FLT_MAX 
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
#  undef SNS_VQ_MAX_IDX
#  define SNS_VQ_MAX_IDX 8

#  define SNSLR_NST1         170    /*   3*170 + 2 = 512  */
#  define SNSLR_NPVQ_L5K6    1970   /*  mux into  11 bits 2048 pos , 78 slots remaining,  6 whole bits  */
#  define SNSLR_NPVQ_L5K8    5890   /*  additional split LF part muxed into 13 bits,  log2(5980)= 12.5241    */
#  define SNSLR_NPVQ_L5K5    1002   /*  mux into  10 bits 1024 pos  , 22 slots remainging 4 whole bits  */
#  define SNSLR_NPVQ_L8K2    128    /*  1+ 6 = 7  bits  */
#  define SNSLR_NPVQ_L15K5   207006 /*  1+ 16.6593 bits  */

#endif

#define MAX_UINT8 255 
#      define THRESH_100_DMS_TDC_CNT    9
#      define THRESH_100_DMS_NS_CNT     7
#      define THRESH_100_DMS_TDC_NS_CNT 73
#      define THRESH_075_DMS_TDC_CNT    7
#      define THRESH_075_DMS_NS_CNT     7
#      define THRESH_075_DMS_TDC_NS_CNT 87
#      define THRESH_050_DMS_TDC_CNT    22
#      define THRESH_050_DMS_NS_CNT     15
#      define THRESH_050_DMS_TDC_NS_CNT 141
#      define THRESH_025_DMS_TDC_CNT    20
#      define THRESH_025_DMS_NS_CNT     21
#      define THRESH_025_DMS_TDC_NS_CNT 278
#define REL_PITCH_THRESH 0.36
#define PLC_LONGTERM_ANALYSIS_MS 200                           /* Analysis window 2000 ms / 10 ms */
 
#define PLC_LONGTERM_ANALYSIS_STARTUP_FILL            0.5f     /* required buffer fill amount, set to 0.0 to not require any fill at all */
 

#define SNS_IDX_LF 0
#define SNS_IDX_HF 1
#define SNS_IDX_SHAPEJ 2
#define SNS_IDX_GAIN 3
#define SNS_IDX_LS_INDA 4
#define SNS_IDX_A 5
#define SNS_IDX_BORGAINLSB 6

#ifdef CR9_C_ADD_1p25MS

/* master integration fixes for 1p25  */
#  define FIX_FLOAT_ENC_QUANTIZE_1P25MS_512KBPS  /* add  two last MDCT coeffs into the last quadruple for global_gain _energy_ analysis   */
#  define FIX_FLOAT_LT_NORMCORR_INIT             /*align state to BASOP start value  of ~.5 as it has an effect on SNS_compute */

#  define FIX_BOTH_1p25_WB_GLOBGAINOFFSET_NONBE                /* 1p25 curve tilt calulation corrected and made into BASOP  */ 
#  define FIX_BOTH_1p25_WB_GLOBGAINOFFSET_LOWLIM_NONBE  -135   /* 1p25 curve tilt calulation limited to value   -135   kbps for NB&WB   */
#  define    FIX_BOTH_1p25_TEST_NEW_GG_EST2                    /*  GG_EST2_will use bands with two coeffs in each  instead of four, for 1p25 WB and 1p25 SSWB , note only active for  regular, hrmode==0  */

#  ifdef FIX_BOTH_1p25_TEST_NEW_GG_EST2 
#    define    FIX_BOTH_1p25_GG_EST_SWB_FB                /*  1.25ms  GG_EST  update for SWB/FB, better RD curve  float and BASOP , active for hrmode==0 */
#  endif

/* defines to activate 2 or 3  tuple 1.25ms loops for all or any of WB,SSWB,SWB,FB  */
#  define FIX_1p25_GG_EST_TUPLES   /*  1.25 ms  will use 2 or 3 tuples,  note only active for  regular, hrmode==0  */
#  ifdef FIX_1p25_GG_EST_TUPLES
#    define    GG_1p25_WB_TUPLES     2
#    define    GG_1p25_SSWB_TUPLES   2
#    define    GG_1p25_SWB_TUPLES    2
#    define    GG_1p25_FB_TUPLES     3
#   define    GG_1p25_MAX_TUPLES    MAX(MAX(GG_1p25_WB_TUPLES ,GG_1p25_SSWB_TUPLES  ), MAX(GG_1p25_SWB_TUPLES ,GG_1p25_FB_TUPLES))     /* used to control common energy loop */

#    ifndef     FIX_BOTH_1p25_ALLOC_SPECTRUM  
#      define    FIX_BOTH_1p25_ALLOC_SPECTRUM  
#     endif 

#    endif 


#  define LTPF_ADAPTIVE_GAIN_RATE 20 /* Number of frames it must take to reach maximum beta from the default value, provided the pitch remains constant */
#  define LTPF_ADAPTIVE_GAIN

#  ifdef LTPF_ADAPTIVE_GAIN
#    define LTPF_PITCH_STABILITY_THRESHOLD  5         /* Number of frames for which the pitch must be constant for adaptive gain and pitch correction to be applied */
#    define LTPF_ADAPTIVE_GAIN_RATE         20        /* Number of frames it must take to reach maximum beta from the default value, provided the pitch remains constant */

#    define LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
#  endif

#ifdef CR9_C_ADD_1p25MS
#define FIX_LTPF_MEM_CONTINUATION
#define FIX_LTPF_PITCH_MEM_LEN
#define FIX_PLC_CONFORM_ISSUES
#define FIX_TDC_BURST_ERROR
#define FIX_LTPF_DEC_FLFX_MISMATCH
#define FIX_TX_RX_STRUCT_STEREO
#define FIX_ADDITONAL_1p25_ISSUES
#define NEW_SIGNALLING_SCHEME_1p25
#define FIX_LTPF_1p25
#endif
#endif

/* Precision Defines */
#define LC3_FABS(x)    (fabsf(x))
#define LC3_POW(x, y)  (powf(x, y))
#define LC3_LOGTEN(x)  (log10f(x))
#define LC3_LOGTWO(x)  (log2f(x))
#define LC3_COS(x)     (cos(x))
#define LC3_SIN(x)     (sin(x))
#define LC3_SQRT(x)    (sqrtf(x))
#define LC3_EXP(x)     (expf(x))
#define LC3_FMIN(x, y) (fminf(x, y))
#define LC3_FMAX(x, y) (fmaxf(x, y))

#  define MAX_BR 320000              /*      400 * 800 */
#  define MIN_BR_100DMS   16000      /*       20 * 800 * 100/100  */
#  define MIN_BR_025DMS   64000      /*       20 * 800 * 100/ 25  */
#  define MIN_BR_050DMS   32000      /*       20 * 800 * 100/ 50  */
#  define MAX_BR_050DMS_NB   260800  /*      163 * 800 * 100/ 50  */
#  define MAX_BR_100DMS_NB   114400  /* for 100ms at  8kHz */
#  define MAX_BR_100DMS_WB   221600  /* for 100ms at 16kHz */
#  define MAX_BR_100DMS_SSWB 314400  /* for 100ms at 24kHz */

#  ifdef CR9_C_ADD_1p25MS
#    define MIN_BR_0125DMS   128000  /*       20 * 800 * 100/ 12.5  */
#    define MAX_BR_0125DMS   512000    
#  endif

#  define MIN_BR_075DMS_48KHZ_HR ((int)124800/ 800/2)* 800
#  define MIN_BR_075DMS_96KHZ_HR ((int)149600/ 800/2)* 800
#  define MIN_BR_075DMS   21334 /* ceil( 20 * 800 * 100/ 75) */
#  define MAX_BR_075DMS  426667 /* ceil(400 * 800 * 100/ 75) */
#  define MAX_BR_075DMS_NB   152534  /* ceil(143 * 800 * 100/ 75) */
#  define MAX_BR_075DMS_WB   295467  /* ceil(277 * 800 * 100/ 75) */
#  define MAX_BR_075DMS_SSWB 419200  /* ceil(393 * 800 * 100/ 75) */
typedef int32_t  LC3_INT32;

#  if defined(__xtensa__)
#    define ALIGNMENT_BALLOC 4
#    define ALIGNMENT_BALLOC_RED 3
#  else
#    define ALIGNMENT_BALLOC 8
#    define ALIGNMENT_BALLOC_RED 7
#  endif

/* PLC2/PhEcu fading settings */
/* PLC2/PHEcu muting Table setup settings */
#    define PLC2_FADEOUT_IN_MS_MIN        30          /* Table min */
#    define PLC2_FADEOUT_IN_MS_MAX        140         /* Table max  */
#    define PLC2_FADEOUT_RES   10                     /* 10 ms steps used in fadeout constant tables  */

/* current active settings */
#  define PLC2_FADEOUT_IN_MS        30         /* 30   P800 fadeout optimized  */
#  define PLC2_FADEOUT_LONG_IN_MS   120        /* 120  MUSHRA,  && stable tonal  fadeout optimized  */

#  define PHECU_FRES 62.5
#  define PHECU_C_JACOB 1.1429
#  define MAX_LGW 9 /*  LGW48K + 1 !! */
#  define QUOT_LPR_LTR 4
#  define MAX_PLC_LPROT ((512 * 48) / 32)
#  define MAX_PLC_NPLOCS ((MAX_PLC_LPROT / 4) + 1)
#  define MAX_PLC_LMSPEC ((MAX_PLC_LPROT / 2) + 1)
#  define MAX_PLC_LMEM (400) /*  "only"  up to 20kHz (400 MDCT bins at 10 ms) at 48 kHz supported by PhEcu    */

#  define POS_ONE_Q15 (32767.0 / 32768.0)
#  define PHECU_LTOT_MIN_MAN 1   /* lowest possible mantissa energy value */
#  define PHECU_LTOT_MIN_EXP -61 /* L_tot =  PHECU_LTOT_MIN_MAN*2^(PHECU_LTOT_MIN_EXP-31) */
#  define PHECU_LTOT_MIN
#  define PHECU_GRP_SHAPE_INIT 0 /* BASOP Q15 */
#  define PHECU_ENV_STAB_LOCAL POS_ONE_Q15
#  define PHECU_DELTA_CORR 5
#  define PHECU_PFIND_SENS 0.93
#  define PHECU_LA 0

#  define LC3_ROUND(x) (roundf(x))
#  define LC3_FLOOR(x) (floorf(x))

#  define LC3_CONST_POW_2_16 65536
#  define LC3_CONST_POW_2_M16 1.525878906250000e-05
#  define LC3_CONST_POW_2_100 1.267650600228229e+30

#  define MAX_LEN_PCM_PLC (MAX_PITCH + MAX_LEN)
#  define MAX_PITCH  CEILING((MAX_PITCH_12K8 * MAX_LEN * 100), 12800)
#  define TDC_L_FIR_HP 11
#  define PLC3_HPBLENDTHROTTLE 30                  /* higher numbers increase throttled blending from hp filtered to unfiltered uv excitation (0 is no throttle) */

#  define PLC_FADEOUT_TYPE_1_IN_MS 200
#  define PLC_FADEOUT_IN_MS 60                     /* fade-out to zero in ms for TD-PLC and NS, minimum value is 20 */
#  define PLC4_TRANSIT_START_IN_MS 20              /* begin of transition time for noise substitution for voiced signals */
#  define PLC4_TRANSIT_END_IN_MS PLC_FADEOUT_IN_MS /* end   of transition time for noise substitution */
#  define PLC34_ATTEN_FAC_100   0.5000           /* attenuation factor for NS and TDC @ 10  ms*/
#  define PLC34_ATTEN_FAC_075   0.5946           /* attenuation factor for NS and TDC @ 7.5 ms */
#  define PLC34_ATTEN_FAC_050   0.7071           /* attenuation factor for NS and TDC @ 5.0 ms*/
#  define PLC34_ATTEN_FAC_025   0.8409           /* attenuation factor for NS and TDC @ 2.5 ms*/

#ifdef CR9_C_ADD_1p25MS
#  define PLC34_ATTEN_FAC_0125  0.9995           /* attenuation factor for NS and TDC @ 1.25 ms*/
#endif

#  define FEC_SLOT_BYTES_MIN 40
#  define FEC_SLOT_BYTES_MAX 400
#  ifdef CR12_D_FIX_BITRATE_LIMITS
#    ifdef ENABLE_HR_MODE_FL
#      define FEC_SLOT_BYTES_MIN_025DMS_48KHZ_HR 54
#      define FEC_SLOT_BYTES_MIN_025DMS_96KHZ_HR 61
#      define FEC_SLOT_BYTES_MIN_050DMS_48KHZ_HR 87
#      define FEC_SLOT_BYTES_MIN_050DMS_96KHZ_HR 101
#      define FEC_SLOT_BYTES_MIN_075DMS_48KHZ_HR 110
#      define FEC_SLOT_BYTES_MIN_075DMS_96KHZ_HR 126
#      define FEC_SLOT_BYTES_MIN_100DMS_48KHZ_HR 140
#      define FEC_SLOT_BYTES_MIN_100DMS_96KHZ_HR 164
#    endif
#  endif

#  define LC3_CONST_POW_2_M15 3.051757812500000e-05
#  define LC3_CONST_POW_2_23 8388608
#  define LC3_CONST_POW_2_23_NEG -8388608
#  define LC3_CONST_POW_2_23_RED 8388607

#  define LC3_CONST_POW_2_100 1.267650600228229e+30

/* G192 bitstream writing/reading */
#define G192_REDUNDANCY_FRAME 0x6B22
#define G192_GOOD_FRAME 0x6B21
#define G192_BAD_FRAME 0x6B20
#define G192_ZERO 0x007F
#define G192_ONE 0x0081
#define READ_G192FER /* Allow C executable to also read G192 formatted FER files */

#ifdef  DEBUG 
#ifdef READ_G192FER
#  define  READ_G192_FER_BYTE  /* Allow C executable to also read G192 byte formatted FER files  0x20=BAD , 0x21=Good  */ 
#endif 
#endif 

#  define LC3_EPS (1.1e-7f)

#  define M_PI_LC3PLUS 3.14159265358979323846

/* FUNCTION MACROS */
#define CEILING(x, y) (((x) + (y)-1) / (y))
#define FRAME2FS_IDX_10MS(x)  (x<500 ? (x/100) : 5)        /*   80 -> 0, 160 -> 1, 240 -> 2, 320 -> 3, 480 -> 4 , 960  -> 5*/
#define FS2FS_IDX(x) ((x) == 96000 ? 5 : (x) / 10000) /* 8000 -> 0, 16000 -> 1, 24000 -> 2, 32000 -> 3, 48000 -> 4, 96000 -> 5 */

#define UNUSED(x) (void)(x) /* silence unused parameter warning */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define STATIC_ASSERTS(cond, s) typedef char assert_##s[(cond) ? 1 : -1]
#define STATIC_ASSERTI(cond, i) STATIC_ASSERTS(cond, i)
#define STATIC_ASSERT(cond) STATIC_ASSERTI(cond, __LINE__)

/* For dynamic memory calculations */
#define CODEC_FS(fs) ((fs) == 44100 ? 48000 : (fs))
#define DYN_MAX_LEN(fs) MAX(CODEC_FS(fs) / 100, 160)
#  define DYN_MAX_LEN_EXT(fs) MAX(CODEC_FS(fs) / 100, 160) /* extension to length 160 for NB(fs=8000)    */
#define DYN_MAX_MDCT_LEN(fs) (DYN_MAX_LEN(fs) - (180 * DYN_MAX_LEN(fs) / 480))

/* OPTIONS */

#define MAX_LEN_NR 480
#define MAX_SR 96000
#define EXT_RES_ITER_MAX 20
#define MAX_BW_BANDS_NUMBER 6
#define MAX_LEN MAX_SR/100 /* = 10ms at 96kHz */
#define MAX_RESBITS 5000
#define MAX_RESBITS_LEN ((MAX_RESBITS + 7)/8)

#define MAX_CHANNELS 2
#define MIN_NBYTES      20  /*  100dms:  16  kbps at !=44.1kHz,  14.7kbps at 44.1kHz
                                 50dms:  32  kbps at !=44.1kHz,  29.4kbps at 44.1kHz
                                 25dms:  64  kbps at !=44.1kHz,  58.8kbps at 44.1kHz */
#define MAX_NBYTES_025 100  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */
#define MAX_NBYTES_050 200  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */
#define MAX_NBYTES_100 400  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */

#ifdef ENABLE_HR_MODE_FL
#    define MIN_BR_25MS_48KHZ_HR ((int)172800/3200/2)*3200
#    define MIN_BR_25MS_96KHZ_HR ((int)198400/3200/2)*3200
#    define MIN_BR_50MS_48KHZ_HR ((int)148800/1600/2)*1600
#    define MIN_BR_50MS_96KHZ_HR ((int)174400/1600/2)*1600
#    define MIN_BR_100MS_48KHZ_HR ((int)124800/800/2)*800
#    define MIN_BR_100MS_96KHZ_HR ((int)149600/800/2)*800
#endif /* ENABLE_HR_MODE */

#define MAX_NBYTES2 625
#define BYTESBUFSIZE (MAX_NBYTES2 * MAX_CHANNELS)
#define MAX_BW_BIN 400
#if MAX_BW_BIN > MAX_LEN
#  define MAX_BW MAX_LEN
#else
#  define MAX_BW MAX_BW_BIN
#endif

#  ifdef ENABLE_HR_MODE_FL
#    define MAX_BW_HR 960
#  endif

/* SCF */
#define M 16
#define MAX_BANDS_NUMBER 64
#define MAX_BANDS_NUMBER_PLC 80
#define PVQ_MAX_VEC_SIZE M

/* PVQ VQ setup */
#define SCF_MAX_PARAM                                                  \
    7 /* (L+H) + submode_MSB +gain+(Ia_leads+Ia_mpvq)+(Ib_joint_mpvq), \
         submode-LSB */

#ifdef CR9_C_ADD_1p25MS
#  undef SNS_VQ_MAX_IDX 
#  define SNS_VQ_MAX_IDX 8

#  define SNSLR_NST1         170     /*,    3*170 + 2 = 512  */
#  define SNSLR_NPVQ_L5K6   1970     /*  mux into  11 bits 2048 pos , 78 slots remaining,  6 whole bits  */
#  define SNSLR_NPVQ_L5K5    1002    /*  mux into  10 bits 1024 pos  , 22 slots remainging 4 whole bits  */
#  define SNSLR_NPVQ_L8K2   128      /*  1+6 = 7  bits  */
#  define SNSLR_NPVQ_L15K5   207006  /*     1+ 16.6593 bits  ,  */

#  define SNSLR_ST1_SCALEC       (round(1.5*16384.0)/16384.0 )            /*  SNS_ST1_SCALEC 1.5  in  Q1.14  */
#  define SNSLR_ST1_INVSCALEC    (round((2.0/3.0)*(32768.0))/32768.0)     /*  SNSLR_ST1_INVSCALEC  0.66667  in  Q0.15  */

#endif /* CR9_C_ADD_1p25MS */

/* RESIDUAL CODING */
#define NPRM_RESQ 5 * MAX_LEN

/* MDCT */
#define MDCT_MEM_LEN_MAX (MAX_LEN - ((180 * MAX_LEN) / 480))

/* TNS */
#define TNS_NUMFILTERS_MAX 2
#define MAXLAG 8

/* OLPA/LTPF */
#define LEN_12K8 128
#define LEN_6K4 64
#define MIN_PITCH_6K4 17
#define MAX_PITCH_6K4 114
#define RANGE_PITCH_6K4 98
#define MIN_PITCH_12K8 32
#define MAX_PITCH_12K8 228
#define RES2_PITCH_12K8 157
#define RES4_PITCH_12K8 127
#define LTPF_MEMIN_LEN (MAX_PITCH_12K8 + 4)
#ifdef CR9_C_ADD_1p25MS
#define LEN_MEM_NORMCORR 5
#else
#define LEN_MEM_NORMCORR 2
#endif
/* Advanced PLC */



/* some configurations leave empty translation units. */
extern int fix_empty_translation_unit_warning;

#endif
