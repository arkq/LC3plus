/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
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
                                  , Word16 hrmode, Word16 regBits, LC3PLUS_FrameDuration frame_dms
#else
#if defined(FIX_BOTH_1p25_TEST_NEW_GG_EST2) ||defined (FIX_1p25_GG_EST_TUPLES) 
                                  , LC3PLUS_FrameDuration frame_dms
#endif 
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
#ifdef   FIX_1p25_GG_EST_TUPLES 
    /* NB (-1)  WB(20), SSWB(30), SWB(40), FB(50) */
    Word16  nCoeffTab1p25[5] = { -1/*NB*/, GG_1p25_WB_TUPLES, GG_1p25_SSWB_TUPLES, GG_1p25_SWB_TUPLES, GG_1p25_FB_TUPLES };
    Word16  bwIdx;
    Word16 divTabQ15[5] = { 0, -32768, 16384, 10923, 8192 };
#endif 
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
    Word16 nCoeffPerBand = 4;    /* 4 coeff per band  */
    Word16 twoCoeffFlag = 0;   /*  boolean */
#endif 
#if defined( FIX_1p25_FLEX_ITER_TUPLE_LOOP ) || defined( FIX_1p25_FLEX_BISECT_LOOP ) 
    Word16 twoCoeffOutShiftA;
    Word16 twoCoeffOutShiftB;
    Word32 L_twoCoeffSuppressMaskA;
    Word32 L_twoCoeffSuppressMaskB;
    Word32 L_dB_scale_offset, L_dB_scale_offset2, L_dB_scale_offset3;

    Word16 twoCoeffOutShift;
    Word32 L_twoCoeffSuppressMask;

#elif defined(     FIX_BASOP_1p25_NEW_GG_EST3) 
    Word16 twoCoeffOutShift;
    Word32 L_twoCoeffSuppressMask;
    Word32 L_dB_scale_offset, L_dB_scale_offset2, L_dB_scale_offset3;
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
    UNUSED(iszero);

#ifdef FIX_1p25_FLEX_ITER_TUPLE_LOOP 
    UNUSED(twoCoeffOutShift);
    UNUSED(L_twoCoeffSuppressMask);
#endif 
#ifdef FIX_1p25_FLEX_BISECT_LOOP
    UNUSED(twoCoeffOutShiftA);
    UNUSED(L_twoCoeffSuppressMaskA);
    UNUSED(twoCoeffOutShiftB);
    UNUSED(L_twoCoeffSuppressMaskB);
#endif

#ifdef  FIX_1p25_FLEX_ITER_TUPLE_LOOP 
    UNUSED(twoCoeffOutShiftA);
    UNUSED(L_twoCoeffSuppressMaskA);
    UNUSED(twoCoeffOutShiftB);
    UNUSED(L_twoCoeffSuppressMaskB);
#endif

#ifndef   FIX_1p25_GG_EST_TUPLES 
#  ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && (sub(lg, 30) <= 0))  /* WB and SSWB */
    {
        nCoeffPerBand = 2; move16();  /* 2 coeff per band  */
        twoCoeffFlag = 1;  move16();
    }
#  endif 
#endif 


#ifdef   FIX_1p25_GG_EST_TUPLES 
    twoCoeffFlag = 0;     move16();
    nCoeffPerBand = 4;    move16();    /* default assume 4 coeff per band(4-tuples block)  */
    lg_4 = shr_pos_pos(lg, 2);

    test();
    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {   
        bwIdx = mult(lg, 3276);
        assert(bwIdx == ((lg - 1) / 10));
        nCoeffPerBand = nCoeffTab1p25[bwIdx]; move16();

        if (sub(nCoeffPerBand, 2) == 0)
        {
            twoCoeffFlag = 1; move16();
        }
        s = mult(add(lg, sub(nCoeffPerBand, 1)), divTabQ15[nCoeffPerBand]);

        tmp16 = i_mult(s, nCoeffPerBand);

#ifdef FIX_1p25_32kHz_CLANG_WARNING_EST_GAIN  
       /* common 2,3,4 tuple energy loop may access tail values in  blocks of GG_1p25_MAX_TUPLES  */
       /* even if they are actually unused as values and shifted out, they should not be accessed uninitialized */
        tmp16 = add(tmp16, sub(GG_1p25_MAX_TUPLES,nCoeffPerBand));
#endif

        FOR(i = lg; i < tmp16; i++)
        {
            x[i] = L_deposit_l(0);  /* zero the  top virtual MDCT lines in the spectrum signal x */
        }
        assert(MAX_LEN >= s);
        lg_4 = s;    /* lg_4 also used for 3_tuples, 2 tuples*/
    }
  
#  elif  defined( FIX_BASOP_ENC_QUANTIZE_1P25MS_512KBPS)
    /* handle non-eaxactly quadruple lengths, in case  frame_length is (e.g 1p25ms 24kHz --> 30 = 7*4 + 2 ) */
    /*                                                                  e.g 1p25ms 48kHz --> 50 = 12*4 + 2 )         */
    /* quadruple extension is needed provide a correct x_max analysis as req. for scaling analysis   */
    lg_4 = shr_pos(lg, 2); /*1/4*/

    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) 
    {
        tmp16 = sub(lg, shl_pos(lg_4, 2));
        test();
        if (tmp16 != 0)
        {
            lg_4 = add(lg_4, 1);
        }
        FOR(i = lg; i < lg_4 * 4; i++)
        {
            x[i] = L_deposit_l(0);  /* zero these top virtual MDCT lines */
        }
    }
#  endif
  

    en = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = MAX_LEN bytes */


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

    *old_targetBits = nbitsSQ;     move16();
    nbitsSQ = add(nbitsSQ, round_fx(*targetBitsOff));

#ifndef  FIX_BASOP_ENC_QUANTIZE_1P25MS_512KBPS
    lg_4  = shr_pos(lg, 2);
#endif

    x_max = 0;    move32();

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
    {  /* regular  hrmode==0  */
#if defined    ( FIX_1p25_FLEX_ITER_TUPLE_LOOP )  
        if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
        {
            /* 4 or 2 tuples */
            L_dB_scale_offset = L_shr_pos(0x9CCCD, twoCoeffFlag);    /* S15Q16   7.5(4 tuples) or 3.75(2_tuples)   offset */
            L_dB_scale_offset2 = L_shr_pos(0x460000, twoCoeffFlag);
            L_dB_scale_offset3 = L_shr_pos(0x3C7AE, twoCoeffFlag);

            /* set up constants for blocks of 4   (or block of 3 or 2)*/

            /*init for 2-tuples energy loop */
            twoCoeffOutShiftA = 31;       move16();   /*  shift out */
            twoCoeffOutShiftB = 31;       move16();   /*  shift out */
            L_twoCoeffSuppressMaskA = 0;  move32(); /* suppress value*/
            L_twoCoeffSuppressMaskB = 0;  move32(); /* suppress value*/


#if   GG_1p25_MAX_TUPLES>=3 
             IF(sub(nCoeffPerBand, 3) == 0)
             {
                 L_dB_scale_offset = L_add(Mpy_32_16(0x9CCCD, 24576), 1L); /* S15Q16/  round(7.5 * 3/4),  481689.75 -> 481690 */
                 L_dB_scale_offset2 = Mpy_32_16(0x460000, 24576); /*3440640.0 */
                 L_dB_scale_offset3 = Mpy_32_16(0x3C7AE, 24576); /* 185794.5*/

                 twoCoeffOutShiftA = 0;          move16();              /* do not shift out */
                 L_twoCoeffSuppressMaskA = -1L;  move32();              /* all ones,  do not suppress*/
                 twoCoeffOutShiftB = 31;          move16();    /* shift out        x[3]  */
                 L_twoCoeffSuppressMaskB = 0;    move32();     /* suppress value*/
             }
#endif 

#if   GG_1p25_MAX_TUPLES>=4 
             IF(sub(nCoeffPerBand, 4) == 0)
             {
                 twoCoeffOutShiftA = 0;          move16();              /* do not shift out */
                 twoCoeffOutShiftB = 0;          move16();              /* do not shift out */
                 L_twoCoeffSuppressMaskA = -1L;  move32();              /* all ones,  do not suppress*/
                 L_twoCoeffSuppressMaskB = -1L;  move32();              /* all ones,  do not suppress*/
             }
#endif 
#endif 

#ifdef     FIX_1p25_FLEX_ITER_TUPLE_LOOP 

             FOR(i = 0; i < lg_4; i++)   /*  1.25 ms energy loop  used for all    N/4 , (N)/3 and N/2 */
             {
                 /* normalization */
                 s = 31;     move16();

                 tmp32 = L_abs(x[0]);
                 tmp32 = L_max(tmp32, L_abs(x[1]));

                 /* 0 or 31 shift, 0 for blocks of 4,  31 to accumulate a 0 value (for block #2,#3)  for blocks of 2 or blocks of 3  */
#if   GG_1p25_MAX_TUPLES>=3 
                 tmp32 = L_max(tmp32, L_shr_pos_pos(L_abs(x[2]), twoCoeffOutShiftA));
#endif 
#if   GG_1p25_MAX_TUPLES>=4  
                 tmp32 = L_max(tmp32, L_shr_pos_pos(L_abs(x[3]), twoCoeffOutShiftB));
#endif 
                 x_max = L_max(x_max, tmp32);  /* global x_max */


                 if (tmp32 != 0)
                 {
                     s = norm_l(tmp32);
                 }

                 s = sub(s, 2); /* 2 bits headroom . for up to four values accumulation */
                                  /* NB only 1 bit headroom needed with L_mac0            */
                                  /* NB only 2tuples will only need 1 bit even with L_mac */
                 /* calc quadruple or , triple or tuple energy */
                 ener = L_deposit_l(1);


                 tmp16 = round_fx(L_shl(x[0], s)); /* s is in range [29...-2], i.e. can become a downshift of two  */
                 ener = L_mac(ener, tmp16, tmp16);

                 tmp16 = round_fx(L_shl(x[1], s));
                 ener = L_mac(ener, tmp16, tmp16);

#if   GG_1p25_MAX_TUPLES >= 3 
                 tmp32 = L_and(x[2], L_twoCoeffSuppressMaskA);  /* mask is  0 for 2blocks, all ones  for 3 or 4 blocks  */
                 tmp16 = round_fx(L_shl(tmp32, s));
                 ener = L_mac(ener, tmp16, tmp16);
#endif 

#if   GG_1p25_MAX_TUPLES >= 4 
                 tmp32 = L_and(x[3], L_twoCoeffSuppressMaskB);   /* 0 for 2 or 3  blocks, all ones for 4 blocks  */
                 tmp16 = round_fx(L_shl(tmp32, s));
                 ener = L_mac(ener, tmp16, tmp16);
#endif 
                 s = shl_pos(sub(x_e, s), 1);

#  ifdef NONBE_GAIN_EST_FIX
                 test(); test();
                 if ((L_sub(ener, 1L) == 0) && (s < 0))
                 {
                     s = 0; move16();
                 }
#  endif

#ifndef  LOG2_LC_APPROX  
                 /* log2 */
                 tmp32 = L_add_sat(BASOP_Util_Log2(ener), L_shl_sat(L_deposit_l(s), 25)); /* log2, 6Q25 */ 
#endif 
#ifdef LOG2_LC_APPROX    
                 tmp32 = L_add_sat(BASOP_Util_Log2_LC(ener), L_shl_sat(L_deposit_l(s), 25)); /* log2LC , 6Q25 */
#endif 
                 tmp32 = L_add(L_shr_pos(Mpy_32_16(tmp32, 0x436E), 6), L_dB_scale_offset);            /* -> (28/20)*( (7or5.25or3.5 ) + 10*tmp32/log2(10)), 15Q16 */

                 en[i] = tmp32;   move32();   

                 x += nCoeffPerBand;
                 /* N/4  or  N/3  or N/2 */
             }   /* for loop*/
         }  /* end of IF 1p25 */
         else
         { /*  4 or 2 merged  loop  for non-1p25ms */
#endif

#ifndef   FIX_1p25_GG_EST_TUPLES 
 #  ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
        test();
        if (twoCoeffFlag != 0)
        {
            lg_4 = add(lg_4, lg_4);  /*lg_4(extended) for N/4  doubled and used as "lg_2" for N/2  */
        }
#  endif 
#endif 

 #ifdef     FIX_BASOP_1p25_NEW_GG_EST3 
        /*set up constants for blocks of two  (or block of 4)*/
        L_dB_scale_offset = L_shr_pos(0x9CCCD, twoCoeffFlag); /* S15Q16/  7.5 or 3.5  offset */
#ifdef     FIX_1p25_FLEX_ITER_TUPLE_LOOP 
        Word16 twoCoeffOutShift;
        Word32   L_twoCoeffSuppressMask;
#endif 
        twoCoeffOutShift = 0;           move16();              /* do not shift out */
        L_twoCoeffSuppressMask = -1L;   move32();              /* all ones-mask,  do not suppress*/

        IF(twoCoeffFlag != 0)
        {
            twoCoeffOutShift = 31; move16();
            L_twoCoeffSuppressMask = 0;  move16();/* supress value*/
            x[lg_4 * 2] = 0; move32();     /*make sure extra two coeffs at tail are  zeroed */
            x[lg_4 * 2 + 1] = 0; move32();
        };
#endif

        FOR (i = 0; i < lg_4; i++)
        {
            /* normalization */
            s = 31;        move16();

            tmp32 = L_abs(x[0]);
            tmp32 = L_max(tmp32, L_abs(x[1]));
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3 
            /* 0 or 31 shift, 0 for blocks of 4,  31 to accumulate a 0 value (for block #2,#3)  for blocks of 2  */
            tmp32 = L_max(tmp32, L_shr_pos(L_abs(x[2]), twoCoeffOutShift));
            tmp32 = L_max(tmp32, L_shr_pos(L_abs(x[3]), twoCoeffOutShift));
#  endif
#else 
            tmp32 = L_max(tmp32, L_abs(x[2]));
            tmp32 = L_max(tmp32, L_abs(x[3]));
#endif 
            x_max = L_max(x_max, tmp32);

            if (tmp32 != 0) {
                s = norm_l(tmp32);
            }
            s = sub(s, 2); /* 2 bits headroom */ /* note: slightly suboptimal 2 bit headroom  kept for pairs  */

            /* calc quadruple energy */
            ener = L_deposit_l(1);

            tmp16 = round_fx(L_shl(x[0], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[1], s));
            ener  = L_mac(ener, tmp16, tmp16);

#ifdef      FIX_BOTH_1p25_TEST_NEW_GG_EST2 
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3 
            tmp32 = L_and(x[2], L_twoCoeffSuppressMask);  /* 0 for 2blocks, all ones for 4 blocks  */
            tmp16 = round_fx(L_shl(tmp32, s));
            ener = L_mac(ener, tmp16, tmp16);

            tmp32 = L_and(x[3], L_twoCoeffSuppressMask);   /* 0 for 2blocks, all ones for 4 blocks  */
            tmp16 = round_fx(L_shl(tmp32, s));
            ener = L_mac(ener, tmp16, tmp16);
#  endif
#else 
            tmp16 = round_fx(L_shl(x[2], s));
            ener  = L_mac(ener, tmp16, tmp16);

            tmp16 = round_fx(L_shl(x[3], s));
            ener  = L_mac(ener, tmp16, tmp16);
#endif 
            s = shl_pos(sub(x_e, s), 1);
            if (ener == 1 && s < 0) {
                s = 0;
            }

            /* log */
            tmp32 = L_add_sat(BASOP_Util_Log2(ener), L_shl_sat(L_deposit_l(s), 25)); /* log2, 6Q25 */
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3  
            tmp32 = L_add(L_shr_pos(Mpy_32_16(tmp32, 0x436E), 6), L_dB_scale_offset);            /* -> (28/20)*( (7or3.5) + 10*tmp32/log2(10)), 15Q16 */
#  endif 
#else 
            tmp32 = L_add(L_shr_pos(Mpy_32_16(tmp32, 0x436E), 6), 0x9CCCD); /* -> (28/20)*(7+10*tmp32/log2(10)), 15Q16 */
#endif 
            en[i] = tmp32;         move32();

#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
            x += nCoeffPerBand;
#else 
            x += 4;
#endif
        } /* N/4 or  N/2 energy  loop */
#ifdef  FIX_1p25_FLEX_ITER_TUPLE_LOOP 
        }  /* end of non-1.25ms  4 or 2 tuple  energy loop */
#endif 
    }/* hrmode==0 */
    
    IF (x_max == 0)
    {
        *quantizedGainMin = quantizedGainOff;         move16();
        *quantizedGain = 0;         move16();
        *old_targetBits = -1;         move16();
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

            int frame_dms_val = 0;
            SWITCH (frame_dms) /* 1 / frame_dms in Q15 */
            {
#ifdef CR9_C_ADD_1p25MS
                case LC3PLUS_FRAME_DURATION_1p25MS: assert(0); frame_dms_val = 13; BREAK;
#endif
                case LC3PLUS_FRAME_DURATION_2p5MS: frame_dms_val = 25; BREAK;
                case LC3PLUS_FRAME_DURATION_5MS: frame_dms_val = 50; BREAK;
                case  LC3PLUS_FRAME_DURATION_7p5MS: frame_dms_val =  75; BREAK;
                case LC3PLUS_FRAME_DURATION_10MS: frame_dms_val = 100; BREAK;
                case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
            }
            
            if (ratio_prod < frame_dms_val)
            {
                Word16 mult_factor = 0;

                SWITCH (frame_dms) /* 1 / frame_dms in Q15 */
                {
#ifdef CR9_C_ADD_1p25MS
                    case LC3PLUS_FRAME_DURATION_1p25MS: mult_factor = 15721; BREAK;
#endif
                case LC3PLUS_FRAME_DURATION_2p5MS: mult_factor = 1311; BREAK;
                    case LC3PLUS_FRAME_DURATION_5MS: mult_factor = 655; BREAK;
                    case  LC3PLUS_FRAME_DURATION_7p5MS: mult_factor =  437; BREAK;
                    case LC3PLUS_FRAME_DURATION_10MS: mult_factor = 328; BREAK;
                    case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
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
            
            sub_val = 0x183BA045;             move16();
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
                ener = en[i];                 move16();
                s = en_exp[i];                 move16();

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
                en[i] = tmp32;                 move32();
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
                        iszero = 0;                         move16();

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
        /* regular mode */ 
        /* SQ scale: 4 bits / 6 dB per quadruple */
        {

 #ifndef FIX_1p25_FLEX_BISECT_LOOP 
#  ifdef  FIX_BASOP_1p25_NEW_GG_EST3 
            /*these limits were constant for blocks of 4 , now we adaptively adjust them for blocks of two */
            L_dB_scale_offset = L_shr_pos(0x9CCCD, twoCoeffFlag);
            L_dB_scale_offset2 = L_shr_pos(0x460000, twoCoeffFlag);
            L_dB_scale_offset3 = L_shr_pos(0x3C7AE, twoCoeffFlag);
#  endif 
#endif
            target = L_shl_pos(L_mult(0x7D71, nbitsSQ), 1); /* -> (28/20) * (1.4) * nbitsSQ */
            fac = L_add(0x1000000, 0); /* -> 256 */

#ifdef FIX_1p25_FLEX_BISECT_LOOP 

          /*  1p25 iteration code, for 4-tuples  3-tuples and 2 tuples */
            test();
            IF((frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) && (sub(nCoeffPerBand, 4) <= 0))
            {
                Word16 idx;    /* 0, 1, 2 , section index */
                Word32 L_tmp32;
                Word32 L_offset2;

                Word16 shiftl_fac[3] = { -31, 0, 1 };        /* ROM float mult_fact[3] = { 0.0, 1.0, 2.0 }; */  /* -31 shifts out all energy in an signed Word32 */

                /* 4 or 2 tuples */
                L_dB_scale_offset = L_shr_pos(0x9CCCD, twoCoeffFlag);    /* S15Q16   7.5(4 tuples) or 3.75(2_tuples)   offset */
                L_dB_scale_offset2 = L_shr_pos(0x460000, twoCoeffFlag);
                L_dB_scale_offset3 = L_shr_pos(0x3C7AE, twoCoeffFlag);

                /* set up constants for blocks of 4   (or block of 3 or 2)*/
#if   GG_1p25_MAX_TUPLES>=3 
                IF(sub(nCoeffPerBand, 3) == 0)
                {
                    L_dB_scale_offset = L_add(Mpy_32_16(0x9CCCD, 24576), 1L); /* S15Q16/  round(7.5 * 3/4),  481689.75 -> 481690 */
                    L_dB_scale_offset2 = Mpy_32_16(0x460000, 24576); /*3440640.0 */
                    L_dB_scale_offset3 = Mpy_32_16(0x3C7AE, 24576); /* 185794.5*/
                }
#endif 
                /* tuple-size dynamic constant tables, treated as registers */
                Word32 L_en_lim[2];
                Word32 L_add_fac[3];

                L_en_lim[0] = L_add(L_dB_scale_offset, 0);  /*28/20*7.0*/;
                L_en_lim[1] = L_add(L_dB_scale_offset2, 0); /*28/20*50.0*/

                L_add_fac[0] = L_add(L_dB_scale_offset3, 0);
                L_add_fac[1] = L_add(0, 0);
                L_add_fac[2] = L_negate(L_dB_scale_offset2);


                /* find offset (0 to 255) */
                FOR(iter = 0; iter < 8; iter++)
                {
                    fac = L_shr_pos_pos(fac, 1);
                    offset = L_sub(offset, fac);
                    iszero = 1;     move16();
                    ener = L_deposit_l(0);
                    /* skip iszero section,  go to first band with enough energy    */
                    /* find first index where we start  summing up adding bits. */

                    L_offset2 = L_add(offset, L_dB_scale_offset);
                    FOR(i = sub(lg_4, 1); i >= 0; i--)
                    {
                        L_tmp32 = L_sub(en[i], L_offset2);
                        IF(L_tmp32 > 0)
                        {
                            iszero = 0;  move16(); /* now also allow upward update */
                            BREAK;
                        }
                    }
                    /*  exiting counter value "i" is remembered , for the next FOR loop initialization */


                    /* Equations   realized   acc = {    0*d_en + 2.7   ,  1*dEn + 0.0,    2*dEn + (-50.0)  }; */
                    FOR(; i >= 0; i--)
                    {
                        L_tmp32 = L_sub(en[i], offset);

                        /* pre_quantize energy "L_tmp32" into a  idx 0, 1, 2;  to address tabled acc. values     */
                        idx = 0;  move16();       /* "ptr"  init to low section */
                        if (L_sub(L_tmp32, L_en_lim[0]) > 0)
                        {
                            idx = add(idx, 1);     /* single op */  /* middle section */
                        }

                        if (L_sub(L_tmp32, L_en_lim[1]) > 0)
                        {
                            idx = add(idx, 1);     /* single op */  /* high  section*/
                        }

                        ener = L_add_sat(ener, L_shl_sat(L_tmp32, shiftl_fac[idx]));     /* *0.0, *1.0, *2.0 */  /* NB! this shifts L_tmp32, (the ener contribution)  both up and down */
                        ener = L_add_sat(ener, L_add_fac[idx]);                          /*  2.7,  0,  -50  */
                    }
                  

                    /* if ener is above target -> increase offset */
                    /* if only low level signal , iszero was never set to  0  --> offset not moved */
                    assert(target > 0);
                    test(); test();
                    if (L_sub(ener, target) > 0 && iszero == 0)
                    {
                        offset = L_add(offset, fac);
                    }
                }  /* iter 1--8 ,1p25  loop */
            }
            /* if 1p25 */
            ELSE
            {
 #endif  /* FLEX BISECT*/ 
                /* find offset (0 to 127) */
#ifdef     FIX_BASOP_1p25_NEW_GG_EST3  
                 /* 4 ,or  2 tuples */
                 L_dB_scale_offset = L_shr_pos(0x9CCCD, twoCoeffFlag);
            /* S15Q16   7.5(4 tuples) or 3.75(2_tuples)   offset */
                  L_dB_scale_offset2 = L_shr_pos(0x460000, twoCoeffFlag);
                  L_dB_scale_offset3 = L_shr_pos(0x3C7AE, twoCoeffFlag);
#endif 
                  FOR(iter = 0; iter < 8; iter++)
                     {
                         fac = L_shr_pos(fac, 1);
                         offset = L_sub(offset, fac);

                         ener = L_deposit_l(0);
                         iszero = 1;  move16();


                         FOR(i = lg_4 - 1; i >= 0; i--)
                         {
                             tmp32 = L_sub(en[i], offset);
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3              
                             diff = L_sub(tmp32, L_dB_scale_offset); /* 0x9CCCD -> (28/20)*( 7or3.5) in Q16 */
#  endif 
#else 
                             diff = L_sub(tmp32, 0x9CCCD); /* 0x9CCCD -> (28/20)*(7) in Q16 */
#endif 
                             if (diff < 0)
                             {
                                 if (iszero == 0)
                                 {
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3   
                                   ener = L_add_sat(ener, L_dB_scale_offset3);   /* 0x3C7AE -> (28/20)*(2.7or 3.5) in Q16  */
#  endif 
#else 
                                   ener = L_add_sat(ener, 0x3C7AE); /* 0x3C7AE -> (28/20)*(2.7) in Q16 */
#endif 
                                 }
                             }
                             else
                             {
                                 ener = L_add_sat(ener, tmp32);
                                 iszero = 0;               move16();
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2
#  ifdef     FIX_BASOP_1p25_NEW_GG_EST3
                                 diff2 = L_sub(tmp32, L_dB_scale_offset2); /* 0x460000 -> (28/20)*(50or25) in Q16 */
#  endif 
#else 
                                 diff2 = L_sub(tmp32, 0x460000); /* 0x460000 -> (28/20)*(50) in Q16 */
#endif 
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
                     } /* legacy 1-8 iter loop , FOR */
#ifdef FIX_1p25_FLEX_BISECT_LOOP 
            } /* 1p25, non 1p25*/
#endif 
            tmp16 = extract_h(offset);
        }

            if (sub(tmp16, *quantizedGainMin) < 0)
            {
                *old_targetBits = -1;            move16();
            }
        *quantizedGain = sub(s_max(*quantizedGainMin, tmp16), quantizedGainOff);         move16();
    }

#ifdef ENABLE_HR_MODE
    tmp32 = Mpy_32_16(0x797CD707, L_shl_pos(add(*quantizedGain, quantizedGainOff), 6));
#else
    tmp32 =  L_shl_pos(L_mult0(add(*quantizedGain, quantizedGainOff), 0x797D), 7); /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
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

