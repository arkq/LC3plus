/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef DEFINES_H
#define DEFINES_H

#include "stdint.h"

/* Precision Defines */
/* LC3_DOUBLE_PRECISION can be used for doubles */

#ifndef LC3_DOUBLE_PRECISION
#define LC3_SINGLE_PRECISION
#endif



#ifdef LC3_SINGLE_PRECISION
#define LC3_FABS(x) (fabsf(x))
#define LC3_POW(x, y) (powf(x, y))
#define LC3_LOG10(x) (log10f(x))
#define LC3_LOG2(x) (log2f(x))
#define LC3_COS(x) (cosf(x))
#define LC3_SIN(x) (sinf(x))
#define LC3_SQRT(x) (sqrtf(x))
#define LC3_EXP(x) (expf(x))
typedef float LC3_FLOAT;
#endif

#ifdef LC3_DOUBLE_PRECISION
typedef double LC3_FLOAT;
#define LC3_FABS(x) (fabs(x))
#define LC3_POW(x, y) (pow(x, y))
#define LC3_LOG10(x) (log10(x))
#define LC3_LOG2(x) (log2(x))
#define LC3_COS(x) (cos(x))
#define LC3_SIN(x) (sin(x))
#define LC3_SQRT(x) (sqrt(x))
#define LC3_EXP(x) (exp(x))
#define kiss_fft_scalar double
#endif

typedef int32_t LC3_INT;
typedef int16_t LC3_INT16;
typedef uint16_t LC3_UINT16;
typedef short LC3_SHORT;
typedef uint8_t LC3_UINT8;
typedef int8_t LC3_INT8;
typedef uint32_t LC3_UINT32;


/* Release defines */
#define ENABLE_2_5MS_MODE
#define ENABLE_5MS_MODE
#define ENABLE_10_MS_MODE
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

#endif /* NO_POST_REL_CHANGES */

/* G192 bitstream writing/reading */
#define G192_REDUNDANCY_FRAME 0x6B22
#define G192_GOOD_FRAME 0x6B21
#define G192_BAD_FRAME 0x6B20
#define G192_ZERO 0x007F
#define G192_ONE 0x0081
#define READ_G192FER /* Allow C executable to also read G192 formatted FER files */



#define M_PI 3.14159265358979323846

/* FUNCTION MACROS */
#define CEILING(x, y) (((x) + (y)-1) / (y))
#define FRAME2FS_IDX(x) (x / 100) /*   80 -> 0, 160 -> 1, 240 -> 2, 320 -> 3, 480 -> 4*/
#define FS2FS_IDX(x)                                                                       \
    (x / 10000)             /*   8000 -> 0, 16000 -> 1, 24000 -> 2, 32000 -> 3, 48000 -> 4 \
                             */
#define UNUSED(x) (void)(x) /* silence unused parameter warning */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define STATIC_ASSERTS(cond, s) typedef char assert_##s[(cond) ? 1 : -1]
#define STATIC_ASSERTI(cond, i) STATIC_ASSERTS(cond, i)
#define STATIC_ASSERT(cond) STATIC_ASSERTI(cond, __LINE__)

/* For dynamic memory calculations */
#define CODEC_FS(fs) ((fs) == 44100 ? 48000 : (fs))
#define DYN_MAX_LEN(fs) MAX(CODEC_FS(fs) / 100, 160)
#define DYN_MAX_MDCT_LEN(fs) (DYN_MAX_LEN(fs) - (180 * DYN_MAX_LEN(fs) / 480))

/* OPTIONS */

#define MAX_SR 96000
#define EXT_RES_ITER_MAX 20
#define MAX_BW_BANDS_NUMBER 6
#define MAX_LEN MAX_SR/100 /* = 10ms at 96kHz */
#define MAX_RESBITS 5000
#define MAX_RESBITS_LEN ((MAX_RESBITS + 7)/8)

#define MAX_CHANNELS 2
#define MIN_NBYTES 20       /* 16kbps at 8/16/24/32/48kHz */
#define MAX_NBYTES 400      /* 320kbps at 48kHz */
#define MAX_NBYTES2 625
#define BYTESBUFSIZE (MAX_NBYTES2 * MAX_CHANNELS)
#define MAX_BW_BIN 400
#if MAX_BW_BIN > MAX_LEN
#define MAX_BW MAX_LEN
#else
#define MAX_BW MAX_BW_BIN
#endif

/* SCF */
#define M 16
#define MAX_BANDS_NUMBER 64
#define MAX_BANDS_NUMBER_PLC 80
#define PVQ_MAX_VEC_SIZE M

/* PVQ VQ setup */
#define SCF_MAX_PARAM                                                  \
    7 /* (L+H) + submode_MSB +gain+(Ia_leads+Ia_mpvq)+(Ib_joint_mpvq), \
         submode-LSB */

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

/* Advanced PLC */



/* some configurations leave empty translation units. */
extern int fix_empty_translation_unit_warning;

#endif
