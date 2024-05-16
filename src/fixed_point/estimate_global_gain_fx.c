/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processEstimateGlobalGain_fx(Word32 x[], Word16 x_e, Word16 lg, Word16 nbitsSQ,
#ifdef ENABLE_HR_MODE
                                  Word32 *gain,
#else
                                  Word16 *gain,
#endif
                                  Word16 *gain_e,
                                  Word16 *quantizedGain, Word16 *quantizedGainMin, Word16 quantizedGainOff,
                                  Word32 *targetBitsOff, Word16 *old_targetBits, Word16 old_specBits,
                                  Word8 *scratchBuffer
#ifdef ENABLE_HR_MODE
                                  , Word16 hrmode, Word16 regBits, Word16 frame_dms
#endif
)
{

    Word16  lg_4, tmp16, iszero, s;
    Word32  ener, tmp32, x_max;
    Word32  target, fac, offset;
    Word32 *en;
    Counter iter, i;
    Word32 diff, diff2;
#ifdef ENABLE_HR_MODE
    Word16 *en_exp = NULL;
    Word32 M0, M1;
#endif

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processEstimateGlobalGain_fx", sizeof(struct {
                   Word16  lg_4, s, tmp16, iszero;
                   Word32  ener, tmp32, x_max;
                   Word32  target, fac, offset;
                   Word32 *en;
                   Counter i, iter;
                   Word32  diff, diff2;
               }));
#endif

    en = (Word32 *)scratchAlign(scratchBuffer,
                                0); /* Size = MAX_LEN bytes */

#ifdef ENABLE_HR_MODE
    if (hrmode)
    {
        M0 = 1;
        M1 = 1; /* Regularization factor; needs´to be 1e-5, but 1e-5 is 0 in Q15 */
        en_exp = (Word16 *) scratchAlign(en, sizeof(*en) * MAX_LEN);
    }
#endif
    IF (*old_targetBits < 0)
    {
        *targetBitsOff = 0;
        move16();
    }
    ELSE
    {
        tmp32          = L_add(*targetBitsOff, L_deposit_h(sub(*old_targetBits, old_specBits)));
        tmp32          = L_min((40 << 16), L_max(-(40 << 16), tmp32));
        *targetBitsOff = L_add(Mpy_32_16(*targetBitsOff, 26214), Mpy_32_16(tmp32, 6554));
        move16();
    }

    *old_targetBits = nbitsSQ;
    move16();
    nbitsSQ = add(nbitsSQ, round_fx(*targetBitsOff));

    lg_4  = shr_pos(lg, 2);
    x_max = 0;
    move32();

/* energy of quadruples with 9dB offset */
#ifdef ENABLE_HR_MODE
    IF (hrmode)
    {
        FOR (i = 0; i < lg_4; i++)
        {
            Word32 absval;
            Word16 idx;
            /* normalization */
            s = 31;
            move16();

            /* M1 requires a 32x16 mult with Q0 i, resulting in Q15. Keeping both M0 and M1 in same Q */
            /* Use Q15 for M0 and M1 calculation */
            idx = shl(i, 2);

            tmp32  = L_abs(x[0]);
            absval = L_shr(tmp32, 16);
            M0     = L_add(M0, absval);              /* M0 += fabs(x[idx])*/
            M1     = L_add(M1, L_mult(absval, idx)); /* M1 += i*fabs(x[idx])*/
            idx    = add(idx, 1);

            absval = L_abs(x[1]);
            tmp32  = L_max(tmp32, absval);
            absval = L_shr(tmp32, 16);
            M0     = L_add(M0, absval);              /* M0 += fabs(x[idx])*/
            M1     = L_add(M1, L_mult(absval, idx)); /* M1 += idx*fabs(x[idx])*/
            idx    = add(idx, 1);

            absval = L_abs(x[2]);
            tmp32  = L_max(tmp32, absval);
            absval = L_shr(tmp32, 16);
            M0     = L_add(M0, absval);              /* M0 += fabs(x[idx])*/
            M1     = L_add(M1, L_mult(absval, idx)); /* M1 += idx*fabs(x[idx])*/
            idx    = add(idx, 1);

            absval = L_abs(x[3]);
            tmp32  = L_max(tmp32, absval);
            absval = L_shr(tmp32, 16);
            M0     = L_add(M0, absval);              /* M0 += fabs(x[idx])*/
            M1     = L_add(M1, L_mult(absval, idx)); /* M1 += idx*fabs(x[idx])*/

            x_max = L_max(x_max, tmp32);

            if (tmp32 != 0)
                s = norm_l(tmp32);

            s = sub(s, 2); /* 2 bits headroom */

            /* calc quadruple energy */
            ener = L_deposit_l(1);

            tmp16 = round_fx(L_shl(x[0], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[1], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[2], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[3], s));
            ener  = L_mac(ener, tmp16, tmp16);

            s = shl_pos(sub(x_e, s), 1);
            if (ener == 1 && s < 0)
                s = 0;
            IF (regBits > 0)
            {
                en[i]     = ener;
                en_exp[i] = s;
                move32();
            }
            ELSE
            {
                /* log */
                tmp32 = L_add(BASOP_Util_Log2(ener), L_shl_pos(L_deposit_l(s), 25)); /* log2, 6Q25 */
                tmp32 = L_add(L_shr_pos(Mpy_32_16(tmp32, 0x436E), 7),
                              0x4E666); /* -> (28/20)*(7+10*tmp32/log2(10)), 16Q15 */
                en[i] = tmp32;
                move32();
            }

            x += 4;
        }
    }
    ELSE
#endif
    {
        FOR (i = 0; i < lg_4; i++)
        {
            /* normalization */
            s = 31;
            move16();

            tmp32 = L_abs(x[0]);
            tmp32 = L_max(tmp32, L_abs(x[1]));
            tmp32 = L_max(tmp32, L_abs(x[2]));
            tmp32 = L_max(tmp32, L_abs(x[3]));
            x_max = L_max(x_max, tmp32);

            if (tmp32 != 0)
                s = norm_l(tmp32);

            s = sub(s, 2); /* 2 bits headroom */

            /* calc quadruple energy */
            ener = L_deposit_l(1);

            tmp16 = round_fx(L_shl(x[0], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[1], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[2], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[3], s));
            ener  = L_mac(ener, tmp16, tmp16);

            s = shl_pos(sub(x_e, s), 1);
            if (ener == 1 && s < 0)
                s = 0;

            /* log */
            tmp32 = L_add_sat(BASOP_Util_Log2(ener), L_shl_sat(L_deposit_l(s), 25)); /* log2, 6Q25 */
            tmp32 =
                L_add(L_shr_pos(Mpy_32_16(tmp32, 0x436E), 6), 0x9CCCD); /* -> (28/20)*(7+10*tmp32/log2(10)), 15Q16 */
            en[i] = tmp32;
            move32();
            x += 4;
        }
    }
    
    IF (x_max == 0)
    {
        *quantizedGainMin = quantizedGainOff;
        move16();
        *quantizedGain = 0;
        move16();
        *old_targetBits = -1;
        move16();
    }
    ELSE
    {
        Word32 sub_val = 0xFCDD38F;
        /*28 * log10(32767 - 0.375) * (1 - 1e-7) in Q21 */

#ifdef ENABLE_HR_MODE
        if (hrmode)
        {
            /*
            Original float code :
                rB_offset = 8 * (1 - MIN(M1/M0, 2*frame_ms)/(2*frame_ms)
            */
            Word16 ratio;
            Word16 ratio_exp = 0;
            Word32 regterm   = MAX_32; /* 1 in Q31 */
            Word32 rB_offset, reg_val;
            Word32 ratio_prod;
            Word16 n_reg_val;

            if (M0 <= 0x7fff)
            {
                Word16 inv_M0 = Inv16(M0, &ratio_exp);                        /* Inverse in Q(15 - ratio_exp) */
                ratio         = L_shr(Mpy_32_16(M1, inv_M0), 15 - ratio_exp); /* Integer ratio */
            }
            else
            {
                Word16 M0_h = L_shr(M0, 15);
                Word32 M1_h = L_shr(M1, 15);

                Word16 inv_M0 = Inv16(M0_h, &ratio_exp);                        /* Inverse in Q(15 - ratio_exp) */
                ratio         = extract_l(L_shr(Mpy_32_16(M1_h, inv_M0), 15 - ratio_exp)); /* Integer ratio */
            }

            /*
            regterm = MIN (M1 / M0, 2 *frame_ms) / (2 * frame_ms)
            regterm = MIN (M1 * (10 / 2) / M0, frame_dms) / frame_dms
            regterm = MIN ( ratio * 5, frame_dms) / frame_dms
            */

            ratio_prod = L_mult(ratio, 5);

            if (ratio_prod < frame_dms)
            {
                Word16 mult_factor = 0;

                SWITCH (frame_dms) /* 1 / frame_dms in Q15 */
                {
                case 25: mult_factor = 1311; break;
                case 50: mult_factor = 655; break;
                case  75: mult_factor =  437; break;
                case 100: mult_factor = 328; break;
                }

                /* ratio_prod < frame_dms. Hence Word16 can be used */

                regterm = L_shl(L_mult(extract_l(ratio_prod), mult_factor), 15); /* result in Q31 */
            }

            rB_offset = L_sub(MAX_32, regterm);
            /* Calculation in Q28 to prevent overflow. Subtraction result in Q31, downshift by 3 results in Q28.
               Multiplication by 8 is implemented as upshift by 3. 
            */            

            /*
            FLOAT code : reg_val = x_max * LC3_POW(2,-regBits - rB_offset);
            */
            Word16 reg_val_e = x_e;

            IF(rB_offset > 0)
            {
                Word32 reg_exp = L_negate(L_add(L_shl(regBits, 25), L_shr(rB_offset, 3)));
                reg_val = Mpy_32_32(x_max, BASOP_Util_InvLog2(reg_exp)); /* Product is in Q31 */
                /* reg_val is in Q(31-x_e) */
            }
            ELSE
            {
                reg_val = x_max;
                reg_val_e = sub(x_e, regBits);
            }
            
            sub_val = 0x183BA045;
            move16();
            /*28 * log10(32768*256 - 2) in Q21 */

            /*
            Adding LC3_POW(2, -31) to reg_val.2^-31 * 2^(31-x_e) = 2^-x_e.
            If x_e is positive, this is below precision requirements to be used.
            */

            if (reg_val_e < 0)
            {
                reg_val = L_add_sat(reg_val, L_shl_sat(1, negate(reg_val_e)));
            }
            n_reg_val = norm_l(reg_val);

            FOR (i = 0; i < lg_4; i++)
            {
                ener = en[i];
                move16();
                s = en_exp[i];
                move16();

                Word16 shift_val = sub(reg_val_e, s);

                IF (sub(n_reg_val, shift_val) > 0)
                {
                    IF(shift_val > -32)
                    {
                        ener = L_add(ener, L_shl(reg_val, shift_val)); /* Match q formats */
                    }                    
                }
                ELSE
                {
                    IF (sub(shift_val, 32) >= 0 )
                    {
                        ener = reg_val;
                    }
                    ELSE
                    {
                        ener = L_add_sat(reg_val, L_shr(ener, shift_val));
                    }
                    s    = reg_val_e;
                }

                tmp32 = L_add(BASOP_Util_Log2(ener), L_shl_pos(L_deposit_l(s), 25)); /* log2, 6Q25 */
                tmp32 = L_add(L_shr_pos(Mpy_32_32(tmp32, 0x436E439A), 7), 0x4E666); /* -> (28/20)*(7+10*tmp32/log2(10)), 15Q16 */
                en[i] = tmp32;
                move32();
            }
        }
#endif
        x_max = BASOP_Util_Log2(x_max);
        /* Minimum gain */
        x_max = L_add(x_max, L_shl_pos(L_deposit_l(x_e), 25)); /* log2(x_max) in 6Q25 */
        x_max = L_sub(
            Mpy_32_32(x_max, 0x436E439A),
            sub_val); /* 28*log10(x_max/(32768-0.375)) = log2(x_max)*(28/log2(10))-28*log10(32768-0.375) in 10Q21 */
        /* 28/log1(10) is in Q27
        Mpy_32_32 : Q25(x_max) + Q27 + Q1(Mpy_32_32_ss) - Q32 = Q21  */
        *quantizedGainMin = extract_l(L_shr_pos(L_add(x_max, (1 << 21) + (1 << 11)), 21));
        move16();
        ASSERT(*quantizedGainMin <= 255 + quantizedGainOff);
        *quantizedGainMin = s_max(quantizedGainOff, s_min(add(255, quantizedGainOff), *quantizedGainMin));
        
        offset = L_deposit_h(add(255, quantizedGainOff)); /* -> 127 */
        
#ifdef ENABLE_HR_MODE
        IF(hrmode)
        {
            offset = L_shr_pos(offset, 2);
            /* SQ scale: 4 bits / 6 dB per quadruple */
            target = L_mult(0x3EB8, nbitsSQ); /* -> (28/20) * (1.4) * nbitsSQ */

            fac = L_add(0x400000, 0); /* -> 256 */
            /* find offset (0 to 127) */
            FOR (iter = 0; iter < 8; iter++)
            {
                fac    = L_shr_pos(fac, 1);
                offset = L_sub(offset, fac);

                ener   = L_deposit_l(0);
                iszero = 1;
                move16();

                FOR (i = lg_4 - 1; i >= 0; i--)
                {
                    tmp32 = L_sub(L_shr_pos(en[i], 1), offset);
                    diff  = L_sub(tmp32, 0x27333); /* 0x4E666 -> (28/20)*(7) in Q15 */
                    if (diff < 0)
                    {
                        if (iszero == 0)
                        {
                            ener = L_add(ener, 0xF1EC); /* 0x1E3D7 -> (28/20)*(2.7) in Q15 */
                        }
                    }
                    else
                    {
                        ener   = L_add(ener, tmp32);
                        iszero = 0;
                        move16();

                        diff2 = L_sub(tmp32, 0x118000); /* 0x230000 -> (28/20)*(50) */
                        if (diff2 >= 0)
                        {
                            ener = L_add(ener, diff2);
                        }
                    }
                }

                /* if ener is above target -> increase offset */
                test();
                if (L_sub(ener, target) > 0 && iszero == 0)
                {
                    offset = L_add(offset, fac);
                }
            }
            tmp16 = extract_h(L_shl(offset, 2));
        }
        ELSE
#endif
        {
            /* SQ scale: 4 bits / 6 dB per quadruple */
            target = L_shl_pos(L_mult(0x7D71, nbitsSQ), 1); /* -> (28/20) * (1.4) * nbitsSQ */
            fac = L_add(0x1000000, 0); /* -> 256 */

            /* find offset (0 to 127) */
            FOR (iter = 0; iter < 8; iter++)
            {
                fac    = L_shr_pos(fac, 1);
                offset = L_sub(offset, fac);

                ener   = L_deposit_l(0);
                iszero = 1;
                move16();
                FOR (i = lg_4 - 1; i >= 0; i--)
                {
                    tmp32 = L_sub(en[i], offset);
                    diff  = L_sub(tmp32, 0x9CCCD); /* 0x9CCCD -> (28/20)*(7) in Q16 */
                    if (diff < 0)
                    {
                        if (iszero == 0)
                        {
                            ener = L_add_sat(ener, 0x3C7AE); /* 0x3C7AE -> (28/20)*(2.7) in Q16 */
                        }
                    }
                    else 
                    {
                        ener   = L_add_sat(ener, tmp32);
                        iszero = 0;
                        move16();

                        diff2 = L_sub(tmp32, 0x460000); /* 0x460000 -> (28/20)*(50) in Q16 */
                        if (diff2 >= 0)
                        {
                            ener = L_add_sat(ener, diff2);
                        }
                    }
                
                }

                /* if ener is above target -> increase offset */
                test();
                if (L_sub(ener, target) > 0 && iszero == 0)
                {
                    offset = L_add(offset, fac);
                }
            }
            tmp16 = extract_h(offset);

        }

        if (sub(tmp16, *quantizedGainMin) < 0)
        {
            *old_targetBits = -1;
            move16();
        }
        *quantizedGain = sub(s_max(*quantizedGainMin, tmp16), quantizedGainOff);
        move16();
    }

#ifdef ENABLE_HR_MODE
    tmp32 =
        Mpy_32_16(0x797CD707, L_shl_pos(add(*quantizedGain, quantizedGainOff), 6));
#else
    tmp32 =
        L_shl_pos(L_mult0(add(*quantizedGain, quantizedGainOff), 0x797D), 7); /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
#endif
    *gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);                        /* get exponent */
#ifdef ENABLE_HR_MODE
    *gain = BASOP_Util_InvLog2(L_or(tmp32, (Word32)0xFE000000));
#else
    *gain   = round_fx(BASOP_Util_InvLog2(L_or(tmp32,(Word32) 0xFE000000)));
#endif

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

