/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "rom_basop_util.h"


void dct_IV(Word32 *pDat,       /* i/o: pointer to data buffer */
            Word16 *pDat_e,     /* i/o: pointer to data exponent */
            Word16  L,          /* i  : length of block */
#  ifdef ENABLE_HR_MODE
            Word16  hrmode,     /* indicate high precision usage */
#  endif
            Word32 *workBuffer) /* : size of L */
            
{
    Word16 sin_step;
    Word16 idx;
    Word16 M_var;
    Word16 M2;

    Word32 *pDat_0;
    Word32 *pDat_1;

    Word32 accu1;
    Word32 accu2;
    Word32 accu3;
    Word32 accu4;

    Counter i;

#ifdef ENABLE_HR_MODE
    const PWord32 *twiddle;
    const PWord32 *sin_twiddle;
#else
    const PWord16 *twiddle;
    const PWord16 *sin_twiddle;
#endif

#ifdef ENABLE_DCTIV_RESCALE
    Word16 scale;
#endif

#ifdef DYNMEM_COUNT
#ifdef ENABLE_HR_MODE
    Dyn_Mem_In("dct_IV", sizeof(struct {
                   Word16  sin_step;
                   Word16  idx;
                   Counter i;
                   Word16  M_var;
                   Word16  M2;

                   Word32 *pDat_0;
                   Word32 *pDat_1;

                   Word32 accu1;
                   Word32 accu2;
                   Word32 accu3;
                   Word32 accu4;

                   const PWord32 *twiddle;
                   const PWord32 *sin_twiddle;
               }));

#else
    Dyn_Mem_In("dct_IV", sizeof(struct {
                   Word16  sin_step;
                   Word16  idx;
                   Counter i;
                   Word16  M_var;
                   Word16  M2;

                   Word32 *pDat_0;
                   Word32 *pDat_1;

                   Word32 accu1;
                   Word32 accu2;
                   Word32 accu3;
                   Word32 accu4;

                   const PWord16 *twiddle;
                   const PWord16 *sin_twiddle;
               }));
#endif /* ENABLE_HR_MODE */
#endif /* DYNMEM_COUNT */

    M_var = shr_pos_pos(L, 1);
    M2    = sub(M_var, 1);

    BASOP_getTables(&twiddle, &sin_twiddle, &sin_step, L);

    pDat_0 = &pDat[0];
    pDat_1 = &pDat[L - 2];

    FOR (i = 0; i < M2; i += 2)
    {
#ifdef ENABLE_HR_MODE
      if (hrmode) {
        cplxMpy32_32_32_2(accu1, accu2, pDat_1[1], pDat_0[0], twiddle[i].v.re, twiddle[i].v.im);
        cplxMpy32_32_32_2(accu3, accu4, pDat_1[0], pDat_0[1], twiddle[i + 1].v.re, twiddle[i + 1].v.im);
      } else {
        cplxMpy32_32_16_2(accu1, accu2, pDat_1[1], pDat_0[0], round_fx_sat(twiddle[i].v.re), round_fx_sat(twiddle[i].v.im));
        cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_0[1], round_fx_sat(twiddle[i + 1].v.re), round_fx_sat(twiddle[i + 1].v.im));
      }
#else
        cplxMpy32_32_16_2(accu1, accu2, pDat_1[1], pDat_0[0], twiddle[i].v.re, twiddle[i].v.im);
        cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_0[1], twiddle[i + 1].v.re, twiddle[i + 1].v.im);
#endif

        pDat_0[0] = accu2;           move32();
        pDat_0[1] = accu1;           move32();
        pDat_1[0] = accu4;           move32();
        pDat_1[1] = L_negate(accu3); move32();

        pDat_0 = pDat_0 + 2;
        pDat_1 = pDat_1 - 2;
    }

#ifdef ENABLE_DCTIV_RESCALE
    if (hrmode) 
    {
        
        scale = s_max(getScaleFactor32(pDat, L), 0); move16();
    
        FOR (i = 0; i < L; i++)
        {
            pDat[i] = L_shl_pos(pDat[i], scale); move32();
        }

        *pDat_e = sub(*pDat_e, scale); move16();
    }
#endif

    BASOP_cfft(&pDat[0], &pDat[1], M_var, 2, pDat_e, workBuffer);

    pDat_0 = &pDat[0];
    pDat_1 = &pDat[L - 2];

    idx = sin_step;
    M2  = sub(shr_pos_pos(add(M_var, 1), 1), 1);

    /* Sin and Cos values are 0.0f and 1.0f */
#ifdef ENABLE_HR_MODE
    cplxMpy32_32_32_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#else
    cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#endif

    pDat_1[1] = L_negate(L_shr_pos(pDat_0[1], 1)); move32();
    pDat_0[0] = L_shr_pos(pDat_0[0], 1);           move32();

    FOR (i = 1; i < M2; i++)
    {
        pDat_0[1] = accu3; move32();
        pDat_1[0] = accu4; move32();

        pDat_0 = pDat_0 + 2;
        pDat_1 = pDat_1 - 2;

#ifdef ENABLE_HR_MODE
        cplxMpy32_32_32_2(accu1, accu2, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#else
        cplxMpy32_32_16_2(accu1, accu2, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#endif

        idx += sin_step;

#ifdef ENABLE_HR_MODE
        cplxMpy32_32_32_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#else
        cplxMpy32_32_16_2(accu3, accu4, pDat_1[0], pDat_1[1], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#endif

        pDat_1[1] = L_negate(accu1); move32();
        pDat_0[0] = accu2;           move32();
    }

    pDat_0[1] = accu3; move32();
    pDat_1[0] = accu4; move32();

    pDat_0 = pDat_0 + 2;
    pDat_1 = pDat_1 - 2;

#ifdef ENABLE_HR_MODE
    cplxMpy32_32_32_2(accu3, accu4, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#else
    cplxMpy32_32_16_2(accu3, accu4, pDat_0[1], pDat_0[0], sin_twiddle[idx].v.re, sin_twiddle[idx].v.im);
#endif

/* Last Sin and Cos value pair are the same */
#ifdef ENABLE_HR_MODE
    accu1  = L_shr_pos(Mpy_32_32(pDat_1[0], TWIDDLE), 1);
    accu2  = L_shr_pos(Mpy_32_32(pDat_1[1], TWIDDLE), 1);
#else
    accu1  = L_shr_pos(Mpy_32_16(pDat_1[0], TWIDDLE), 1);
    accu2  = L_shr_pos(Mpy_32_16(pDat_1[1], TWIDDLE), 1);
#endif

    pDat_1[0] = L_add(accu1, accu2); move32();
    pDat_0[1] = L_sub(accu1, accu2); move32();

    pDat_1[1] = L_negate(accu3); move32();
    pDat_0[0] = accu4;           move32();

    /* twiddeling scale is 2 */
    *pDat_e = add(*pDat_e, 2); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

