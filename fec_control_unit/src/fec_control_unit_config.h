/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __fec_control_unit_config_h__
#define __fec_control_unit_config_h__


/* configure fer buffer (controls bufConfig defined in fec_control_unit.c) */

#define BFI_HISTORY_LENGTH_S  15
#define BFI_HISTORY_LENGTH_M  15*15    // 225
#define BFI_HISTORY_LENGTH_L  15*15*15 //3375
#define BFI_HISTORY_RATE_S     1
#define BFI_HISTORY_RATE_M    15   
#define BFI_HISTORY_RATE_L    15*15


/* configure fer reduction buffer */

#define FER_REDUTCION_LENGTH  100
#define FER_REDUCTION_RATE     10
#define FER_REDUCTION_DELTA     1.5f


/* configure initial fers for each ep mode */

#define FER_INIT_EPM4  0.00
#define FER_INIT_EPM3  0.00
#define FER_INIT_EPM2  0.00
#define FER_INIT_EPM1  0.00


/* configure thresholds */

#define STABILITY_THRESHOLD  0.2f
#define FER_SM  30
#define FER_SL  30
#define FER_ML   4
#define EPMR_CONFIDENCE_LEN  10 // threshold for low confidence component in epmode decision
#define FER_BUF_NUM  3  // three histories: short, mid-zise, long

/* definition of content and setup. Select one of the following */

/* the following setups are currently sopported:
 * 0 : NB with 40 byte slots
 * 1 : WB with 40 byte slots
 * 2 : SSWB with 40 byte slots
 * 3 : SWB with 80 byte slots
 */

#ifndef SETUP
#  define SETUP 1
#endif

#if SETUP == 0
#  define NB_40BYTE_SLOT
#  define FER_MAX 30
#elif SETUP == 1
#  define WB_40BYTE_SLOT
#  define FER_MAX 51
#elif SETUP == 2
#  define SSWB_40BYTE_SLOT
#  define FER_MAX 30
#elif SETUP == 3
#  define SWB_80BYTE_SLOT
#  define FER_MAX 30
#endif

/* declarations of constants */

extern const int    bufConfig[3][2];
extern const float  mosLqo[4][FER_MAX+1];
extern const float  default_reduction_factors[3];


/* some helper macros */
#define EP2EPMR(x) (LC3PLUS_EpModeRequest) ((x)-1)
#define EPMR2EP(x) (LC3PLUS_EpMode)(((x) & 3) + 1)
#ifdef MIN
#undef MIN
#endif
#define MIN(x, y) (x) < (y) ? (x) : (y)

#endif // __fec_control_unit_config_h__
