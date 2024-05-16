/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __BASOP_UTIL_ROM_H__
#define __BASOP_UTIL_ROM_H__

#include "basop_util.h"
#include "functions.h"

#define LD_INT_TAB_LEN 120
#define INV_TABLE_SIZE 256
#define SQRT_TABLE_SIZE 256

#ifndef CHEAP_NORM_SIZE
#define CHEAP_NORM_SIZE 161
#endif

#define MINSFTAB 7
#define MAXSFTAB 25

#ifdef ENABLE_HR_MODE
#define SHC(x) ((Word32)x)
#else
#define SHC(x) WORD322WORD16((Word32)x)
#endif

/**
 * \brief  Lookup-Table for binary logarithm
 */
extern const Word16 ldCoeff[7];

/**
  \brief     Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2_tab_long[32];

/**
  \brief     Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2w_tab_long[32];

/**
  \brief     Lookup-Table for binary power algorithm
*/
extern const UWord32 exp2x_tab_long[32];

/**
 * \brief 1/x, x=[0,1,2,3...]  table
 */
#ifdef ENABLE_HR_MODE
extern const Word16 InvIntTable[74];
#else
extern const Word16 InvIntTable[32];
#endif

/**
 * \ brief Lookup for sine tables and windows.
 */
void BASOP_getTables(
#ifdef ENABLE_HR_MODE
                     const PWord32 **ptwiddle,
                     const PWord32 **sin_twiddle,
#else
                     const PWord16 **ptwiddle,
                     const PWord16 **sin_twiddle,
#endif                     
                     Word16        *sin_step,
                     Word16        length);

extern const Word32 RealFFT20_twid[6];
extern const Word32 RealFFT32_twid[10];
extern const Word32 RealFFT40_twid[12];
extern const Word32 RealFFT60_twid[17];
extern const Word32 RealFFT64_twid[18];
extern const Word32 RealFFT80_twid[22];
extern const Word32 RealFFT96_twid[26];
extern const Word32 RealFFT128_twid[34];
extern const Word32 RealFFT192_twid[50];
extern const Word32 RealFFT256_twid[66];
extern const Word32 RealFFT384_twid[98];
extern const Word32 RealFFT512_twid[130];
extern const Word32 RealFFT768_twid[194];

#ifdef ENABLE_HR_MODE
extern const PWord32 SineTable480[241];
extern const PWord32 SineTable320[161];
extern const PWord32 SineTable960[481];
#else
extern const PWord16 SineTable480[241];
extern const PWord16 SineTable320[161];
extern const PWord16 SineTable960[481];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_15_6[2 * (90 - 15)];
#else
extern const Word16 RotVector_15_6[2 * (90 - 15)];
#endif

extern const Word32 RotVector_32_32[2 * 20];
extern const Word32 RotVector_40_32[2 * 28];

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_320[2 * (320 - 20)];
#else
extern const Word16 RotVector_320[2 * (320 - 20)];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_360[2 * (360 - 30)];
#else
extern const Word16 RotVector_360[2 * (360 - 30)];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_480[2 * (480 - 30)];
#else
extern const Word16 RotVector_480[2 * (480 - 30)];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_720[2 * (720 - 30)];
extern const Word32 RotVector_960[2 * (480 - 60)];
#else
extern const Word16 RotVector_960[2 * (480 - 60)];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_30_16[2 * (480 - 30)];
#endif

#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_32_8[2 * (256 - 32)];
#else
extern const Word16 RotVector_32_8[2 * (256 - 32)];
#endif

#if defined(SUBSET_SSWB) || defined(SUBSET_SWB) || defined(SUBSET_FB)
#ifdef ENABLE_HR_MODE
extern const Word32 RotVector_32_12[2 * (384 - 32)];
#else
extern const Word16 RotVector_32_12[2 * (384 - 32)];
#endif
#else
#define RotVector_32_12 NULL
#endif

extern const Word32 isqrt_table[128 + 2];

extern const Word32 Log2_16_table1[16];
extern const Word16 Log2_16_table2[16];

extern const Word32  InvLog2_16_table1[64];
extern const Word16  InvLog2_16_table2[64];

extern const UWord8  gf16_mult_table[256];
extern const UWord8  rs16_elp_deg2_table[256];
extern const UWord16 rs16_elp_deg3_table[256];

#ifdef ENABLE_HR_MODE
extern const Word32 invSqrtTab[(128 + 2)];
#endif

#endif
