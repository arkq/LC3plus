/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef DEFINES_H
#define DEFINES_H

#ifndef DISABLE_HR_MODE
#  define ENABLE_HR_MODE
#endif

#ifdef ENABLE_HR_MODE
#  define MAX_BW_HR 960
#endif

#ifdef ENABLE_HR_MODE
#  define extractW16(a) extract_h(a)
#else
#  define extractW16(a) (a)
#endif

#ifdef ENABLE_HR_MODE
#  define MAX_BR 320000        /*      400 * 800 */
#endif

/* FUNCTION MACROS */
/* NB, divisions in some of these MACROs, use mainly for initial setup, do not use in loops */
#define CEILING(x, y) (((x) + (y)-1) / (y))

#define FRAME2FS_IDX(x) (x / 100) /*   80 -> 0, 160 -> 1, 240 -> 2, 320 -> 3, 480 -> 4 */
#ifdef ENABLE_HR_MODE
#  define FS2FS_IDX(x) ((x) == 96000 ? 5 : (x) / 10000) /* 8000 -> 0, 16000 -> 1, 24000 -> 2, 32000 -> 3, 48000 -> 4, 96000 -> 5 */
#else
#  define FS2FS_IDX(x) (x / 10000)  /*   8000 -> 0, 16000 -> 1, 24000 -> 2, 32000 -> 3, 48000 -> 4 */
#endif

#define UNUSED(x) (void)(x)       /* silence unused parameter warning */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define STATIC_ASSERTS(cond, s) typedef char assert_##s[(cond) ? 1 : -1]
#define STATIC_ASSERTI(cond, i) STATIC_ASSERTS(cond, i)
#define STATIC_ASSERT(cond) STATIC_ASSERTI(cond, __LINE__)
#define TRACE(x) PRINTF("%s:%i %s = %i\n", __FILE__, __LINE__, #x, (int)(x))

/* For dynamic memory calculations */
#define CODEC_FS(fs) ((fs) == 44100 ? 48000 : (fs))

#define DYN_MAX_LEN(fs) (CODEC_FS(fs) / 100)
#define DYN_MAX_LEN_EXT(fs) MAX(CODEC_FS(fs) / 100, 160) /* extension to length 160 for NB(fs=8000)    */
#define DYN_MAX_LPROT(fs) ((512 * (CODEC_FS(fs) / 100)) / 320)

#define DYN_MAX_PLOCS(fs) (DYN_MAX_LPROT(fs) / 4 + 1)
#define DYN_MAX_MDCT_LEN(fs) (DYN_MAX_LEN(fs) - (180 * DYN_MAX_LEN(fs) / 480))

#define MAX_PITCH_FS(fs) (CEILING((MAX_PITCH_12K8 * CODEC_FS(fs)), (12800)))

/*  get the maximum buffer size assuming 10 ms framing */
/*  MAX_PITCH was previously always  the highest fs in the current  SUBSET , i.e. typically 48 kHz  */
/*  now fs(from incoming wav file)  adaptive  MAX_PITCH_FS(fs) is used instead  */
#define DYN_MAX_LEN_PCM_PLC_CLASSIFIER(fs)                                                                             \
    (MAX_PITCH_FS(fs) + DYN_MAX_LEN(fs)) /* CLASSIFIER PCM memory requirement   */
#define DYN_MAX_LEN_PCM_PLC_TDCAPPLYFILTER(fs)                                                                         \
    ((M + 1) + MAX_PITCH_FS(fs) + (DYN_MAX_LEN(fs) / 2)) /* TDC filtering PCM memory requirement */

#ifdef ENABLE_HR_MODE
#  define DYN_MAX_LEN_PCM_PLC(fs) (MAX_PITCH_FS(fs) + DYN_MAX_LEN(fs))
#else
#  define DYN_MAX_LEN_PCM_PLC(fs) MAX(DYN_MAX_LEN_PCM_PLC_CLASSIFIER(fs), DYN_MAX_LEN_PCM_PLC_TDCAPPLYFILTER(fs))
#endif

#define FRAME_MS_BLOCK 25

/* OPTIONS */
#define ENABLE_2_5MS_MODE
#define ENABLE_5MS_MODE
#define ENABLE_10_MS_MODE
#define ENABLE_ADVANCED_PLC
#define ENABLE_ADVANCED_PLC_DEFAULT
#define ENABLE_BW_CONTROLLER
#define ENABLE_ERROR_PROTECTION
#define ENABLE_RFRAME
#define ENABLE_PC
#define ENABLE_PLC
#define ENABLE_PADDING
/* flags */
#define ENABLE_BANDWIDTH_FLAG
#define ENABLE_RBANDWIDTH_FLAG
#define ENABLE_EP_MODE_FLAG
#define ENABLE_FRAME_MS_FLAG

#ifndef NO_POST_REL_CHANGES
/* Post-release non-bitexact changes */

#endif /* NO_POST_REL_CHANGES Post-release changes */

/* G192 bitstream writing/reading */
#define G192_GOOD_FRAME 0x6B21
#define G192_BAD_FRAME 0x6B20
#define G192_REDUNDANCY_FRAME 0x6B22
#define G192_ZERO 0x007F
#define G192_ONE 0x0081
#define READ_G192FER      /* allow C executable to also read G192 formatted FER files */

#define DYNMEM_COUNT
#define STAMEM_COUNT


/*
do not change  __forceinline  for mex compilation using  gcc6.3.0 or larger
  gcc630 supported by MATLAB 2018b, via Mingw "app"
  */
#ifdef __GNUC__
#  define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#else
#  define GCC_VERSION 0
#endif

/* Define __forceinline as empty if ARM not activated to avoid any errors */
#undef __forceinline
/* The above undef is needed to compile using make_mex.m without the following warning in Matlab for windows */
#define __forceinline

/* SUBSETS */
#if !(defined(SUBSET_NB) || defined(SUBSET_WB) || defined(SUBSET_SSWB) || defined(SUBSET_SWB) || defined(SUBSET_FB) || defined(SUBSET_UB))
#  define SUBSET_NB
#  define SUBSET_WB
#  define SUBSET_SSWB
#  define SUBSET_SWB
#  define SUBSET_FB
#  ifdef ENABLE_HR_MODE
#    define SUBSET_UB
#  endif
#endif

#define PACK_RESBITS

/* FRAME/BUFFER */
#ifdef SUBSET_UB
#  define MAX_LEN 960 /* = 10ms at 96kHz */
#elif defined SUBSET_FB
#  define MAX_LEN 480 /* = 10ms at 48kHz */
#elif defined(SUBSET_SWB)
#  define MAX_LEN 320 /* = 10ms at 32kHz */
#elif defined(SUBSET_SSWB)
#  define MAX_LEN 240 /* = 10ms at 24kHz */
#elif defined(SUBSET_WB)
#  define MAX_LEN 160 /* = 10ms at 16kHz */
#elif defined(SUBSET_NB)
#  define MAX_LEN 80 /* = 10ms at 8kHz */
#endif

#ifdef ENABLE_HR_MODE
#  define ENABLE_FFT_RESCALE
#  define ENABLE_FFT_30X16
#  define ENABLE_DCTIV_RESCALE

#  define EXT_RES_ITER_MAX 20
#  define MAX_RESBITS 5000
#  define HR_MODE_SCRATCH_SIZE 60
#else /* ENABLE_HR_MODE */
#  define MAX_RESBITS MAX_LEN
#endif /* ENABLE_HR_MODE */
    
/* BW Cutoff-Detection */
#define MAX_BW_BANDS_NUMBER 5

#ifdef ENABLE_HR_MODE
#  define RESBITS_PACK_SHIFT 3
#  define RESBITS_PACK_N (1<<(RESBITS_PACK_SHIFT))
#  define RESBITS_PACK_MASK (RESBITS_PACK_N - 1)
#  define MAX_RESBITS_LEN ((MAX_RESBITS + RESBITS_PACK_MASK) >> RESBITS_PACK_SHIFT)
#  define MAX_RESBITS_LEN_32BIT_ALIGN (((MAX_RESBITS_LEN + 3)/4)*4)
#else /* ENABLE_HR_MODE */
#  define MAX_RESBITS_LEN (MAX_RESBITS)
#  define MAX_RESBITS_LEN_32BIT_ALIGN (((MAX_RESBITS_LEN + 3)/4)*4)
#endif /* ENABLE_HR_MODE */

#define MAX_CHANNELS 2
#define MIN_NBYTES      20  /*  100dms:  16  kbps at !=44.1kHz,  14.7kbps at 44.1kHz
                                 50dms:  32  kbps at !=44.1kHz,  29.4kbps at 44.1kHz
                                 25dms:  64  kbps at !=44.1kHz,  58.8kbps at 44.1kHz */
#define MAX_NBYTES_025 100  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */
#define MAX_NBYTES_050 200  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */
#define MAX_NBYTES_100 400  /* any dms: 320  kbps at !=44.1kHz, 294  kbps at 44.1kHz */
#ifdef ENABLE_HR_MODE
#    define MIN_BR_25MS_48KHZ_HR ((int)172800/3200/2)*3200
#    define MIN_BR_25MS_96KHZ_HR ((int)198400/3200/2)*3200
#    define MIN_BR_50MS_48KHZ_HR ((int)148800/1600/2)*1600
#    define MIN_BR_50MS_96KHZ_HR ((int)174400/1600/2)*1600
#    define MIN_BR_100MS_48KHZ_HR ((int)124800/800/2)*800
#    define MIN_BR_100MS_96KHZ_HR ((int)149600/800/2)*800
#  define MAX_NBYTES_RED_HR 625 /* From table 5.2 */
#  define BYTESBUFSIZE (MAX_NBYTES_RED_HR * LC3PLUS_MAX_CHANNELS)
#else
#  define BYTESBUFSIZE (MAX_NBYTES_100 * LC3PLUS_MAX_CHANNELS)
#endif
#define MAX_BW_BIN 400
#define FEC_SLOT_BYTES_MIN 40
#  define FEC_SLOT_BYTES_MAX 400
#if MAX_BW_BIN > MAX_LEN
#  define MAX_BW MAX_LEN
#else
#  define MAX_BW MAX_BW_BIN
#endif

#define NUM_OFFSETS 7
#define NUM_SAMP_FREQ 6

/* SCF */
#define M 16 /* LPC_ORDER */
#define MAX_BANDS_NUMBER 64
#define MAX_BANDS_NUMBER_PLC 80
#define PVQ_MAX_VEC_SIZE M
#define KMAX_FX 10 /* N=10  */
#define K_OUTL_FAR 6
#define MPVQ_SZ_OUTL_FAR (1549824U >> 1)
#define TOTBITS_SHAPE OUTL_FAR(1.0 + 19.5637)
#define SQRT_EN_MAX_FX 64 /*  table size, for fast  inv_sqrt  */
#define N_SETA 10
#define N_SETB (PVQ_MAX_VEC_SIZE - N_SETA)

#define SNS_DAMPING 27853        /* 0.85 in Q15 */
#ifdef ENABLE_HR_MODE
#  define SNS_DAMPING_HRMODE 19661 /* 0.6 in Q15 */
#endif

/* PVQ VQ setup */
#define VQMODES26                                                                                                      \
    ((0 << 4) + (1 << 3) + (1 << 2) + (1 << 1) +                                                                       \
     (1 << 0))                      /* select all shapes, including extended lf, and new far                           \
                                              shape for MODE26, add extreme far mode not enabled */
#define SCF_STAGE1_MAX_BITS (5 + 5) /* typically only (5+5) used */
#define SCF_STAGE1_NBCDKENTRIES 32
#define SCF_STAGE2_MAX_GAIN_BITS 2
#define SCF_STAGE2_NB_MODE_BITS 1
#define SCF_STAGE2_MAX_NB_SHAPES 3
#define SCF_MAX_PARAM 7 /* (L+H) + submode_MSB +gain+(Ia_leads+Ia_mpvq)+(Ib_joint_mpvq), submode-LSB */
#define SCF_MAX_PARAM_ST2 (SCF_MAX_PARAM - 2)
#define N_SCF_SHAPES_ST2 4
#define N_SCF_SHAPE_SETUPS 1

/* ARITHMETIC ENCODER */
#define NBITS_CONTEXT 8
#define NBITS_RATEQ 2
#define A_THRES_SHIFT 2
#define A_THRES (1 << A_THRES_SHIFT)
#define VAL_ESC 16
#define SYM_BITS_Q 11

/* RESIDUAL CODING */
#define NPRM_RESQ MAX_LEN

/* NOISE FILLING */
#define NOISEFILLWIDTH 3
#define NOISEFILLSTART 24
#define NOISEFILLWIDTH_5MS 1
#define NOISEFILLSTART_5MS 12
#define NOISEFILLWIDTH_2_5MS 1
#define NOISEFILLSTART_2_5MS 6

/* MDCT */
#ifdef ENABLE_HR_MODE
#  define TWIDDLE ((Word32)0x5a82799a)
#else
#  define TWIDDLE WORD322WORD16(0x5a82799a)
#endif
#define MDCT_MEM_LEN_MAX (MAX_LEN - ((180 * MAX_LEN) / 480))

/* TNS */
#define TNS_NUMFILTERS_MAX 2
#define TNS_COEF_RES 17
#define INDEX_SHIFT 8
#define NSUBDIV 3
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
#define M_LTPF 24
#define LTPF_MEM_X_LEN (MAX_LEN + MAX_LEN / 40 - 2)
#define LTPF_MEM_Y_LEN (MAX_LEN + CEILING(MAX_PITCH_12K8 * MAX_LEN, 128) + (MAX_LEN / 80))

/* PLC */
#define PLC_DEBUG

#define PLC_FADEOUT_IN_MS         60    /* fade-out to zero in ms for TD-PLC and NS, minimum value is 20 */
#define PLC2_FADEOUT_IN_MS        30    /*  0 uses original fixed values for PLC2
                                           -1 uses PLC_FADEOUT_IN_MS as basis for a PLC2 macro-re-calculation
                                           30..100 uses separate setting for PLC2 */
#define PLC4_TRANSIT_END_IN_MS PLC_FADEOUT_IN_MS    /* end   of transition time for noise substitution */
#define PLC4_TRANSIT_START_IN_MS  20    /* begin of transition time for noise substitution for voiced signals */
#define PLC34_ATTEN_FAC_100_FX   0x4000 /* attenuation factor for NS and TDC @ 10  ms (0.5000)*/
#define PLC34_ATTEN_FAC_050_FX   0x5A83 /* attenuation factor for NS and TDC @ 5.0 ms (0.7071)*/
#define PLC34_ATTEN_FAC_025_FX   0x6BA3 /* attenuation factor for NS and TDC @ 2.5 ms (0.8409)*/
#define PLC3_HPBLENDTHROTTLE      30    /* higher numbers increase throttled blending from hp filtered to unfiltered uv excitation (0 is no throttle) */

#define MAX_PITCH_8K (CEILING((MAX_PITCH_12K8 * 8000), (12800))) /*NB  was a risky MACRO   at 0.5 border !!,    */
#define MAX_PITCH_16K ((MAX_PITCH_12K8 * 16000) / (12800))       /* exact integer */
#define MAX_PITCH_24K                                                                                                  \
    (CEILING((MAX_PITCH_12K8 * 24000), (12800)))           /* was a risky MACRO truncation at 0.5 border !! ,   */
#define MAX_PITCH_32K ((MAX_PITCH_12K8 * 32000) / (12800)) /* exact integer */
#define MAX_PITCH_48K ((MAX_PITCH_12K8 * 48000) / (12800)) /* exact integer */

#ifdef ENABLE_HR_MODE
#  define MAX_PITCH_96K ((MAX_PITCH_12K8 * 96000) / (12800)) /* exact integer */
#endif

#ifdef ENABLE_HR_MODE
#  define MAX_PITCH MAX_PITCH_96K
#elif defined(SUBSET_FB)
#  define MAX_PITCH MAX_PITCH_48K
#elif defined(SUBSET_SWB)
#  define MAX_PITCH MAX_PITCH_32K
#elif defined(SUBSET_SSWB)
#  define MAX_PITCH MAX_PITCH_24K
#elif defined(SUBSET_WB)
#  define MAX_PITCH MAX_PITCH_16K
#elif defined(SUBSET_NB)
#  define MAX_PITCH MAX_PITCH_8K
#endif


/*  macri MAX_LEN_PCM_PLC  is used for WC memory  size checking
   DYN_MAX_LEN_PCM_PLC(fs) is used by static memory buffers
 */
#define MAX_LEN_PCM_PLC (MAX_PITCH + MAX_LEN)

#define TDC_L_FIR_HP 11
#define TDC_L_FIR_LP 11

#define BASE_LPROT 512 /* BASE Lprot set to 512 for 32 kHz sampling */

#define LPROT48K (BASE_LPROT * 48 / 32) /* Fs based Lprot */
#define LPROT40K (BASE_LPROT * 40 / 32) /* For Fs == 48000 only use 0 - 20000 Hz in 8 bands */
#define LPROT32K BASE_LPROT
#define LPROT24K (BASE_LPROT * 24 / 32)
#define LPROT16K (BASE_LPROT * 16 / 32)
#define LPROT8K (BASE_LPROT * 8 / 32)

/* #define PHECU_XFP_LA  1 */ /* Original soluiton xfp generated with 1ms look ahead into MDCT memory
                               ctrl.PhECU.LA = 0.001 * ctrl.fs; */
/* #define PHECU_XFP_LA  4 */ /* Option with some lookahead 1/4 ms look ahead into MDCT memory
                               ctrl.PhECU.LA = 0.001 * ctrl.fs / 4; */
#define PHECU_XFP_LA                                                                                                   \
    0 /* Option without lookahead 0 ms look ahead into MDCT memory                                                     \
       ctrl.PhECU.LA = 0; */

#if (PHECU_XFP_LA == 0)

#  define PHECU_LA_48K 0 /* 0 ms */
#  define PHECU_LA_32K 0 /* 0 ms */
#  define PHECU_LA_24K 0 /* 0 ms */
#  define PHECU_LA_16K 0 /* 0 ms */
#  define PHECU_LA_8K 0  /* 0 ms */

#  define MAX_PHECU_LA (0)

#else
#  if (PHECU_XFP_LA == 4)

#    define PHECU_LA_48K (48 / 4) /* 0.25 ms */
#    define PHECU_LA_32K (32 / 4) /* 0.25 ms */
#    define PHECU_LA_24K (24 / 4) /* 0.25 ms */
#    define PHECU_LA_16K (16 / 4) /* 0.25 ms */
#    define PHECU_LA_8K (8 / 4)   /* 0.25 ms */

#    define MAX_PHECU_LA (MAX_LEN / 10 / 4)

#  else

#    define PHECU_LA_48K 48 /* 1 ms */
#    define PHECU_LA_32K 32 /* 1 ms */
#    define PHECU_LA_24K 24 /* 1 ms */
#    define PHECU_LA_16K 16 /* 1 ms */
#    define PHECU_LA_8K 8   /* 1 ms */

#    define MAX_PHECU_LA (MAX_LEN / 10)

#  endif
#endif


/* (PHECU_LA_48K == 0) */

#define COPY_LEN_8K 16
#define COPY_LEN_16K 32
#define COPY_LEN_24K 48
#define COPY_LEN_32K 64
#define COPY_LEN_48K 96

#define OLA_LEN_8K 14
#define OLA_LEN_16K 28
#define OLA_LEN_24K 42
#define OLA_LEN_32K 56
#define OLA_LEN_48K 84

#define LPROT48K_RED LPROT40K /* limit peak searched part of spectrum for FB */
#define LPROT32K_RED LPROT32K /* limit peak searched part of spectrum for SWB */
#define LPROT24K_RED LPROT24K /* limit peak searched part of spectrum for HQ */
#define LPROT16K_RED LPROT16K /* limit peak searched part of spectrum for SQ */
#define LPROT8K_RED LPROT8K   /* limit peak searched part of spectrum for NB */

#define INV_LPROT48K_Q22 5461  /* round(2^22/LPROT48K) , 1/6 */
#define INV_LPROT32K_Q22 8192  /* round(2^22/LPROT32K) , 1/4 */
#define INV_LPROT24K_Q22 10923 /* round(2^22/LPROT24K) , 1/3 */
#define INV_LPROT16K_Q22 16384 /* round(2^22/LPROT16K) , 1/2 */
#define INV_LPROT8K_Q22 32767  /* round(2^22/LPROT8K)  , 1/1 */

#define LGW48K 8 /* 8 bands defined , but may the same frequency groups as for SWB(32k)   */
#define LGW32K 7
#define LGW24K 6
#define LGW16K 5
#define LGW8K 4

#define MAX_LGW 9 /*  LGW48K + 1 !! */

#define UNINIT_OR_UNSAFE_OOLD_SENTINEL -32768
#define LTOT_INIT_FLAG -32768
#define LTOT_MIN_MAN 1   /* lowest possible energy value */
#define LTOT_MIN_EXP -61 /* L_tot=  LTOT_MIN_MAN*2^(LTOT_MIN_EXP-31) */
#define GRP_SHAPE_INIT 0 /* Q15 */

#define TRANA_TIME 4 /* Transient analysis length in ms  */
#define LTRANA48K (48000 * TRANA_TIME / 1000)
#define LTRANA32K (32000 * TRANA_TIME / 1000)
#define LTRANA24K (24000 * TRANA_TIME / 1000)
#define LTRANA16K (16000 * TRANA_TIME / 1000)
#define LTRANA8K (8000 * TRANA_TIME / 1000)

#define MAX_PLOCS ((MAX_LPROT / 4) + 1) /* maximum number of spectral peaks to be searched */
#define QUOT_LPR_LTR 4

#define BETA_MUTE_FAC_INI 16384 /* Q15, initial noise attenuation factor */

#define OFF_FRAMES_LIMIT 30 /* Hard limit  to a maximum of 300ms of  burst decay frames   */

#define Lprot_hamm_len2_48k (3 * 48)
#define Lprot_hamm_len2_32k (3 * 32)
#define Lprot_hamm_len2_24k (3 * 24)
#define Lprot_hamm_len2_16k (3 * 16)
#define Lprot_hamm_len2_8k (3 * 8)

#define FRAME_TIME 10 /* Frame length in ms */

#define L_FRAME48K (48000 * FRAME_TIME / 1000)
#define L_FRAME32K (32000 * FRAME_TIME / 1000)
#define L_FRAME24K (24000 * FRAME_TIME / 1000)
#define L_FRAME16K (16000 * FRAME_TIME / 1000)
#define L_FRAME8K (8000 * FRAME_TIME / 1000)

#ifdef DISABLE_PLC
#  define MAX_LTRANA 0
#  define MAX_LPROT LPROT48K
#  define MAX_L_FRAME L_FRAME48K
#  define MAX_LPROT_RED LPROT48K_RED
#else

#  ifdef SUBSET_UB
#    define MAX_LPROT LPROT48K         /* Max length of protype frame for buffer allocation */
#    define MAX_LPROT_RED LPROT48K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA48K
#    define MAX_L_FRAME L_FRAME48K
#  elif defined(SUBSET_FB)
#    define MAX_LPROT LPROT48K         /* Max length of protype frame for buffer allocation */
#    define MAX_LPROT_RED LPROT48K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA48K
#    define MAX_L_FRAME L_FRAME48K
#  elif defined(SUBSET_SWB)
#    define MAX_LPROT LPROT32K
#    define MAX_LPROT_RED LPROT32K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA32K
#    define MAX_L_FRAME L_FRAME32K
#  elif defined(SUBSET_SSWB)
#    define MAX_LPROT LPROT24K
#    define MAX_LPROT_RED LPROT24K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA24K
#    define MAX_L_FRAME L_FRAME24K
#  elif defined(SUBSET_WB)
#    define MAX_LPROT LPROT16K
#    define MAX_LPROT_RED LPROT16K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA16K
#    define MAX_L_FRAME L_FRAME16K
#  elif defined(SUBSET_NB)
#    define MAX_LPROT LPROT8K
#    define MAX_LPROT_RED LPROT8K_RED /* stack alloc peak searched part */
#    define MAX_LTRANA LTRANA8K
#    define MAX_L_FRAME L_FRAME8K
#  endif

#endif /* #ifdef DISABLE_PLC */

#define INV_L_FRAME48K_Q15 32768 / L_FRAME48K
#define INV_L_FRAME32K_Q15 32768 / L_FRAME32K
#define INV_L_FRAME24K_Q15 32768 / L_FRAME24K
#define INV_L_FRAME16K_Q15 32768 / L_FRAME16K
#define INV_L_FRAME8K_Q15 32768 / L_FRAME8K

#define MDCT_MEM_LEN_48 (L_FRAME48K - ((180 * L_FRAME48K) / 480))

#define R1_48 690
#define R2_48 420
#define R1_16 230
#define R2_16 140
#define R1_25 368
#define R2_25 224

#define MAX_LEN_PCM_PLC_TOT (MAX_LEN_PCM_PLC) /* TDC_pcm_buff   */
#define MAX_WIN_PRE_TDA (0) /* not in RAM any longer , now stored in ROM only as PhECU_wins[fs_idx][1] */

#define MAX_PLCMETH 1

/* Scratch buffer defines */
#define scratchBuffer_ACTIVE
#ifdef ENABLE_HR_MODE
#  define SCRATCH_BUF_LEN_ENC (4 * MAX_LEN + 32 + 32 + 4 * MAX_LEN + 3 * MAX_LEN + MAX_RESBITS_LEN_32BIT_ALIGN)
#  define SCRATCH_BUF_LEN_ENC_CURRENT_SCRATCH (4 * MAX_LEN + 8 * 60 + MAX_LEN)
#else
#  define SCRATCH_BUF_LEN_ENC (4 * MAX_LEN + 32 + 32 + 2 * MAX_LEN + 3 * MAX_LEN + MAX_RESBITS_LEN_32BIT_ALIGN)
#  define SCRATCH_BUF_LEN_ENC_CURRENT_SCRATCH (4 * MAX_LEN)
#endif

#define SCRATCH_BUF_LEN_ENC_TOT (SCRATCH_BUF_LEN_ENC + SCRATCH_BUF_LEN_ENC_CURRENT_SCRATCH)

#ifdef ENABLE_HR_MODE
#  define SCRATCH_BUF_LEN_DEC (4 * MAX_LEN + MAX_RESBITS_LEN_32BIT_ALIGN + 32 + 32 + 4 * MAX_LEN + 32 + 128 + 128)
#else
#  define SCRATCH_BUF_LEN_DEC (4 * MAX_LEN + 2 * MAX_LEN + 32 + 32 + 2 * MAX_LEN + 32 + 128 + 128)
#endif
#define SCRATCH_BUF_LEN_DEC_CURRENT_SCRATCH (2 * MAX_LGW + 8 * MAX_LPROT + 12 * MAX_L_FRAME)
#ifdef ENABLE_HR_MODE
#  define SCRATCH_BUF_LEN_DEC_TOT (SCRATCH_BUF_LEN_DEC + (4 * HR_MODE_SCRATCH_SIZE) + SCRATCH_BUF_LEN_DEC_CURRENT_SCRATCH)
#else
#  define SCRATCH_BUF_LEN_DEC_TOT (SCRATCH_BUF_LEN_DEC + SCRATCH_BUF_LEN_DEC_CURRENT_SCRATCH)
#endif

#define ADVACED_PLC_SIZE                                                                                               \
    (sizeof(AplcSetup) + (4 + 2) * MAX_PLOCS + 2 * MAX_LEN_PCM_PLC_TOT + 0 * MAX_LPROT + 2 * MDCT_MEM_LEN_MAX +        \
     2 * MAX_WIN_PRE_TDA)
/* ERI : Correct this sum after ROM/RAM , joint RAM optimizations */


#define ENC_MAX_SIZE (sizeof(LC3PLUS_Enc) + MAX_CHANNELS * (sizeof(EncSetup) + 6 * MDCT_MEM_LEN_MAX))
#define DEC_MAX_SIZE                                                                                                   \
    (sizeof(LC3PLUS_Dec) + MAX_CHANNELS * (sizeof(DecSetup) + ADVACED_PLC_SIZE + 2 * LTPF_MEM_X_LEN + 2 * LTPF_MEM_Y_LEN + \
                                       2 * MDCT_MEM_LEN_MAX + 2 * MAX_BW))

#if (defined(_M_ARM) || defined(__CC_ARM) || defined(__TMS470__) || defined(__arm__) || defined(__aarch64__)) &&       \
    (defined(__TARGET_FEATURE_NEON) || defined(__ARM_NEON) || defined(__ARM_NEON__) || defined(__ARM_AARCH64_NEON__))
#  define ALIGNMENT_DEFAULT 16
#elif defined(__bfin__)
#  define ALIGNMENT_DEFAULT 4
#else
#  define ALIGNMENT_DEFAULT 8
#endif

/* RAM_ALIGN keyword causes memory alignment of global variables. */
#if defined(__GNUC__)
#  define RAM_ALIGN __attribute__((aligned(ALIGNMENT_DEFAULT)))
#elif defined(__CC_ARM)
#  define RAM_ALIGN __align(ALIGNMENT_DEFAULT)
#elif defined(__ANALOG_EXTENSIONS__)
#  define RAM_ALIGN
#  pragma pack(ALIGNMENT_DEFAULT)
#elif defined(_MSC_VER)
#  define RAM_ALIGN __declspec(align(ALIGNMENT_DEFAULT))
#else
#  define RAM_ALIGN
#endif

#define scratchAlign(ptr, offset) (void *)(((uintptr_t)(ptr) + (offset) + 0x3) & ~0x3)
#define ALIGN_BUFFER_STRUCT


/* some configurations leave empty translation units. */
extern int fix_empty_translation_unit_warning;

#endif
