/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"



void processAdjustGlobalGain_fx(Word16 *gg_idx, Word16 gg_idx_min, Word16 gg_idx_off,
#ifdef ENABLE_HR_MODE
                                Word32 *gain,
#else
                                Word16 *gain,
#endif
                                Word16 *gain_e,
                                Word16 target, Word16 nBits, Word16 *gainChange, Word16 fs_idx
#ifdef ENABLE_HR_MODE
                                , Word16 hrmode, Word16 frame_dms
#endif
                                )
{

    Word32 L_tmp;
    Word16 delta, delta2;
#ifdef ENABLE_HR_MODE
    Word16 gg_idx_inc;
    Word16 gg_idx_inc_max;
    Word16 gg_idx_inc_s;
    Word32 factor;
#endif

#ifdef DYNMEM_COUNT
#if defined(ENABLE_HR_MODE)
    Dyn_Mem_In("processAdjustGlobalGain_fx", sizeof(struct {
                   Word32 L_tmp;
                   Word16 delta, delta2;
                   Word16 gg_idx_inc;
                   Word16 gg_idx_inc_max;
                   Word16 gg_idx_inc_s;
                   Word16 factor; 
               }));
#else
    Dyn_Mem_In("processAdjustGlobalGain_fx", sizeof(struct {
                   Word32 L_tmp;
                   Word16 delta, delta2;
               }));
#endif /* ENABLE_HR_MODE */
#endif /* DYNMEM_COUNT */

#ifdef ENABLE_HR_MODE
    IF (sub(frame_dms, 25) == 0)
    {
        IF (sub(target, 520) < 0)
        {
            factor = 3; move16();
            gg_idx_inc_max = 30; move16();
        } ELSE {
            factor = 4; move16();
            gg_idx_inc_max = 40; move16();
        }
    }
    ELSE IF (sub(frame_dms, 50) == 0)
    {
        factor = 2; move16();
        gg_idx_inc_max = 20; move16();
    }
    ELSE IF (sub(frame_dms, 75) == 0)
    {
        factor = 40265318; move16(); // factor = 1.2 * 2^25
        gg_idx_inc_max = 12 ; move16();
    }
    ELSE
    {
        factor = 1; move16();
        gg_idx_inc_max = 10; move16();
    }
#endif

    IF (sub(nBits, adjust_global_gain_tables[0][fs_idx]) < 0)
    {
        delta = mult_r(add(nBits, 48), 2048);
    }
    ELSE IF (sub(nBits, adjust_global_gain_tables[1][fs_idx]) < 0)
    {
        delta = mult_r(add(nBits, adjust_global_gain_tables[4][fs_idx]), adjust_global_gain_tables[3][fs_idx]);
    }
    ELSE IF (sub(nBits, adjust_global_gain_tables[2][fs_idx]) < 0)
    {
        delta = mult_r(nBits, 683);
    }
    ELSE
    {
        delta = mult_r(adjust_global_gain_tables[2][fs_idx], 683);
    }
    delta2 = add(delta, 2);

    *gainChange = 0; move16();

    test();
    IF (sub(*gg_idx, 255) == 0 && sub(nBits, target) > 0)
    {
        *gainChange = 1; move16();
    }

    test(); test(); test();
    IF ((sub(*gg_idx, 255) < 0 && sub(nBits, target) > 0) || (*gg_idx > 0 && sub(nBits, sub(target, delta2)) < 0))
    {
#ifdef ENABLE_HR_MODE
        IF (hrmode)
        {
            IF (sub(nBits, target) > 0)
            {
                gg_idx_inc = sub(nBits, target);
                IF (sub(frame_dms, 75) == 0)
                {
                    gg_idx_inc = extract_l(L_shr_pos(Mpy_32_16(factor, gg_idx_inc), 10)); // Mpy_32_16(1.2*2^25, gg_idx_inc), 25 - 15)
                    gg_idx_inc = BASOP_Util_Divide1616_Scale(gg_idx_inc, delta, &gg_idx_inc_s);
                    gg_idx_inc = shr_sat(gg_idx_inc, sub(15, gg_idx_inc_s));
                    gg_idx_inc = add(gg_idx_inc, 1); // adding 1 instead of 1.2
                }
                ELSE
                {
                    gg_idx_inc = extract_l(L_mult0(gg_idx_inc, factor));
                    gg_idx_inc = BASOP_Util_Divide1616_Scale(gg_idx_inc, delta, &gg_idx_inc_s);
                    gg_idx_inc = shr_sat(gg_idx_inc, sub(15, gg_idx_inc_s));
                    gg_idx_inc = add(gg_idx_inc, factor);
                }
                gg_idx_inc = s_min(gg_idx_inc, gg_idx_inc_max);
                
                *gg_idx = add(*gg_idx, gg_idx_inc); move16();
            }

            *gg_idx = s_min(*gg_idx, 255); move16();
        }
        ELSE
#endif
        {
            test();
            IF (sub(nBits, sub(target, delta2)) < 0)
            {
                *gg_idx = sub(*gg_idx, 1); move16();
            }
            ELSE IF (sub(*gg_idx, 254) == 0 || sub(nBits, add(target, delta)) < 0)
            {
                *gg_idx = add(*gg_idx, 1); move16();
            }
            ELSE
            {
                *gg_idx = add(*gg_idx, 2); move16();
            }
        }

        *gg_idx = s_max(*gg_idx, sub(gg_idx_min, gg_idx_off)); move16();

#ifdef ENABLE_HR_MODE
        L_tmp = Mpy_32_16(0x3CBE6B83, L_shl_pos(add(*gg_idx, gg_idx_off), 7));
#else
        L_tmp       = L_shl_pos(L_mult0(add(*gg_idx, gg_idx_off), 0x797D), 7); /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
#endif
        *gain_e     = add(extract_l(L_shr_pos(L_tmp, 25)), 1);                 /* get exponent */
#ifdef ENABLE_HR_MODE
        *gain       = BASOP_Util_InvLog2(L_or(L_tmp, (Word32)0xFE000000));
#else
        *gain       = round_fx(BASOP_Util_InvLog2(L_or(L_tmp, (Word32)0xFE000000)));
#endif
        *gainChange = 1; move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

