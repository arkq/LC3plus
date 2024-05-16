/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#  ifdef ENABLE_HR_MODE
void processQuantizeSpec_fx(Word32 x[], Word16 x_e, Word32 gain, Word16 gain_e, Word32 xq[], Word16 nt, Word16 target,
                            Word16 totalBits, Word16 *nBits, Word16 *nBits2, Word16 fs_idx, Word16 *lastnzout,
                            Word16 *codingdata, Word16 *lsbMode, Word16 mode, Word16 hrmode)
{

    Word32  a1, b1, a1_i, b1_i;
    Word16  t, lev1;
    Word16  lastnz, lastnz2;
    Word16  rateFlag;
    Word32  nbits32, nbits232, target32;
    Word16  nt_half;
    Word32  c, ab_max, msb, a1_msb, b1_msb;
    Word16  levmax;
    Word16  s;
    Word16  totBits, nbits_lsb;
    Counter k, lev;
    Word16  maxlevs;
#    ifndef FUNCTION_quantizeSpec_func1
    Word16  tmp16;
    Word32  offs32;
    Counter i;
#    else
    Word32 ARM_params[3];
#    endif

#    ifdef DYNMEM_COUNT
    Dyn_Mem_In("processQuantizeSpec_fx", sizeof(struct {
                   Word32  a1, b1, a1_i, b1_i, ab_max, c;
                   Word16  t, lev1;
                   Word16  lastnz, lastnz2;
                   Word16  rateFlag;
                   Word32  nbits32, nbits232;
                   Word16  nt_half;
                   Word16  msb, a1_msb, b1_msb, levmax;
                   Counter k, lev, i;
                   Word16  s;
                   Word16  tmp16;
                   Word32  offs32, target32;
                   Word16  totBits, nbits_lsb;
                   Word16  maxlevs;
               }));
#    endif

    assert(target >= 0);
    
    /* Quantization */
    gain = invFixp(gain, &gain_e);

    maxlevs = 21;
    IF (hrmode)
    {
        s = sub(add(x_e, gain_e), 23);
        s = s_min(s, 31);
#    ifdef FUNCTION_quantizeSpec_func1
            ARM_params[0] = gain;
            ARM_params[1] = s;

            quantizeSpec_func1_hr_ip(xq, x, nt, ARM_params);
#    else
            IF(s > 4)
            {
                s = sub(s, 4); /* Use extra bits of precision for fine calculation */
                FOR(i = 0; i < nt; i++)
                {
                    offs32 = Mpy_32_32(L_shl_pos(x[i], s), gain); /* multiply */
                    xq[i] =     L_shr_r_pos(offs32, 4);    /* Convert to Q0 with rounding */
                    /* rounding is the equivalent of adding 0.5, which is the offset in hrmode */

                    move32();
                }
            }
            ELSE
            {
                FOR(i = 0; i < nt; i++)
                {
                    offs32 = Mpy_32_32(x[i], gain); /* multiply */
                    offs32 = L_shl(offs32, s);      /* convert to 23Q8 */
                    xq[i] = L_shr_r_pos(offs32, 8);    /* Convert to Q0 with rounding */
                    /* rounding is the equivalent of adding 0.5, which is the offset in hrmode */

                    move32();
                }
            }
#    endif /* FUNCTION_quantizeSpec_func1 */
    }
    ELSE
    {
        s = sub(add(x_e, gain_e), 15);
        s = s_max(s_min(s, 15), -15);
#    ifdef FUNCTION_quantizeSpec_func1
            ARM_params[0] = gain;
            ARM_params[1] = s;
            ARM_params[2] = -4096;
            quantizeSpec_func1_ip(xq, x, nt, ARM_params);
#    else
        FOR (i = 0; i < nt; i++)
        {
                offs32 = Mpy_32_32(L_abs(x[i]), gain); /* multiply */
                offs32 = L_shl(offs32, s);             /* convert to 15Q16 */
                tmp16  = mac_r(offs32, -4096, 1);      /* add offset and truncate */

                if (x[i] < 0)
                    tmp16 = negate(tmp16); /* restore sign */

                /* Normal quantization: xq[i] =  x[i] / gg + sign(x[i]) * 0.375
                   quant_offset is -0.125 in Q15 and round adds 0.5 in Q16. Hence
                   mac_r results in abs(x[i])/gain - 0.125 + 0.5 = abs(x[i])/gain + 0.375.
                   Due to the abs and negate combination this achieves the same result
                   as spec.
                */

                xq[i] = tmp16;
                move16();
            }
#    endif /* FUNCTION_quantizeSpec_func1 */
    }
    /* Rate flag */
    rateFlag = 0;
    move16();
    if (fs_idx != 5)
    {
        if (sub(totalBits, add(160, i_mult(fs_idx, 160))) > 0)
        {
            rateFlag = 2 << NBITS_CONTEXT;
            move16();
        }
    }

    /* Init */
    nt_half = shr_pos(nt, 1);
    c       = 0;
    move16();
    t = 0;
    move16();
    a1_i = 0;
    move16();
    b1_i = 1;
    move16();
    target32 = L_shl_pos(L_deposit_l(target), SYM_BITS_Q);
    nbits32  = L_negate(target32);
    nbits232 = 0;
    move32();
    nbits_lsb = 0;
    move16();

    if (fs_idx != 5)
    {
        IF (mode == 0 && sub(totalBits, add(480, i_mult(fs_idx, 160))) >= 0)
        {
            mode = 1;
            move16();
        }
    }

    /* Find last non-zero tuple */
    lastnz = find_last_nz_pair(xq, nt);
    IF (mode >= 0)
    {
        lastnz2 = 2;
    }
    ELSE
    {
        lastnz2 = lastnz;
    }

    IF (mode < 0)
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }
            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = L_abs(xq[a1_i]);
            b1     = L_abs(xq[b1_i]);
            ab_max = L_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (L_sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = L_add(a1, L_shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = L_add(shl_pos(s_and(c, 0xf), 4), L_add(L_add(a1, b1), 1));
            }
            ELSE IF (L_sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                nbits32       = L_add(nbits32, 2 << SYM_BITS_Q);
                a1_msb        = L_shr_pos_pos(a1, 1);
                b1_msb        = L_shr_pos_pos(b1, 1);
                msb           = L_add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = L_add(shl_pos(s_and(c, 0xf), 4), L_add(L_shl_pos(L_add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(maxlevs, sub(norm_l(ab_max), 8));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(levmax), SYM_BITS_Q + 1));
                a1_msb        = L_shr(a1, levmax);
                b1_msb        = L_shr(b1, levmax);
                msb           = L_add(a1_msb, L_shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }
    ELSE IF (mode == 0)
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }

            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = L_abs(xq[a1_i]);
            b1     = L_abs(xq[b1_i]);
            ab_max = L_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (L_sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = L_add(a1, L_shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = L_add(shl_pos(s_and(c, 0xf), 4), L_add(L_add(a1, b1), 1));
            }
            ELSE IF (L_sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                nbits32       = L_add(nbits32, 2 << SYM_BITS_Q);
                a1_msb        = L_shr_pos_pos(a1, 1);
                b1_msb        = L_shr_pos_pos(b1, 1);
                msb           = L_add(a1_msb, L_shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = add(shl_pos(s_and(c, 0xf), 4), L_add(L_shl_pos(L_add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(maxlevs, sub(norm_l(ab_max), 8));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(levmax), SYM_BITS_Q + 1));
                a1_msb        = L_shr(a1, levmax);
                b1_msb        = L_shr(b1, levmax);
                msb           = L_add(a1_msb, L_shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }
    ELSE
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);

            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }

            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = L_abs(xq[a1_i]);
            b1     = L_abs(xq[b1_i]);
            ab_max = L_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (L_sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = L_add(a1, L_shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a1, b1), 1));
            }
            ELSE IF (L_sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                a1_msb        = L_shr_pos_pos(a1, 1);
                b1_msb        = L_shr_pos_pos(b1, 1);
                msb           = L_add(a1_msb, L_shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                nbits_lsb = add(nbits_lsb, 2);
                if (L_sub(a1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (L_sub(b1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(maxlevs, sub(norm_l(ab_max), 8));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(sub(levmax, 1)), SYM_BITS_Q + 1));
                a1_msb        = L_shr(a1, levmax);
                b1_msb        = L_shr(b1, levmax);
                msb           = L_add(a1_msb, L_shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                a1_msb  = L_shr_pos(a1, 1);
                b1_msb  = L_shr_pos(b1, 1);
                if (a1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                nbits_lsb = add(nbits_lsb, 2);
                if (L_sub(a1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (L_sub(b1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }

                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }

                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }

    /* Number of consumed bits */
    nbits32 = L_add(nbits32, target32);
    totBits = add(extract_l(L_shr_pos_pos(L_sub(nbits32, 1), SYM_BITS_Q)), 1);
    IF (mode > 0)
    {
        totBits = add(totBits, nbits_lsb);
    }
    IF (nBits != NULL)
    {
        *nBits = totBits;
    }
    IF (mode >= 0)
    {
        nbits232 = L_add(nbits232, target32);
        *nBits2  = add(extract_l(L_shr_pos(L_sub(nbits232, 1), SYM_BITS_Q)), 1);
    }
    ELSE
    {
        *nBits2 = *nBits;
        move16();
    }
    IF (mode > 0)
    {
        *nBits2 = add(*nBits2, nbits_lsb);
    }
    *lastnzout = lastnz2;

    /* Truncation of high frequency coefficients */
    IF (lastnz > lastnz2)
    {
        basop_memset(&xq[lastnz2], 0, (lastnz - lastnz2) * sizeof(*xq));
    }

    /* Truncation of LSBs */
    test();
    IF (mode > 0 && sub(totBits, target) > 0)
    {
        *lsbMode = 1;
        move16();
    }
    ELSE
    {
        *lsbMode = 0;
        move16();
    }

#    ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#    endif
}

#  else /* ENABLE_HR_MODE */

void processQuantizeSpec_fx(Word32 x[], Word16 x_e, Word16 gain, Word16 gain_e, Word16 xq[], Word16 nt, Word16 target,
                            Word16 totalBits, Word16 *nBits, Word16 *nBits2, Word16 fs_idx, Word16 *lastnzout,
                            Word16 *codingdata, Word16 *lsbMode, Word16 mode
)
{

    Word16  a1, b1, a1_i, b1_i;
    Word16  t, lev1;
    Word16  lastnz, lastnz2;
    Word16  rateFlag;
    Word32  nbits32, nbits232, target32;
    Word16  nt_half;
    Word16  c, ab_max, msb, a1_msb, b1_msb, levmax;
    Word16  s;
    Word16  totBits, nbits_lsb;
    Counter k, lev;
    Word16  tmp16;
    Word32  offs32;
    Counter i;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processQuantizeSpec_fx", sizeof(struct {
                   Word16  a1, b1, a1_i, b1_i;
                   Word16  t, lev1;
                   Word16  lastnz, lastnz2;
                   Word16  rateFlag;
                   Word32  nbits32, nbits232;
                   Word16  nt_half;
                   Word16  c, ab_max, msb, a1_msb, b1_msb, levmax;
                   Counter k, lev, i;
                   Word16  s;
                   Word16  tmp16;
                   Word32  offs32, target32;
                   Word16  totBits, nbits_lsb;
               }));
#endif

    /* Quantization */
    gain = Inv16(gain, &gain_e);
    s    = sub(add(x_e, gain_e), 15);
    s = s_max(s_min(s, 15), -15);
     
    FOR (i = 0; i < nt; i++)
    {
        offs32 = Mpy_32_16(L_abs(x[i]), gain);          /* multiply */
        offs32 = L_shl(offs32, s);                      /* convert to 15Q16 */
        tmp16 = extract_h(L_add(offs32, 0x00006000));   /* add offset and truncate */
        Word16 x_sign = (Word16) L_shr(x[i], 31);
        xq[i] = sub(s_xor(tmp16, x_sign), x_sign);
        move16();
        /*
           Normal quantization: xq[i] =  x[i] / gg + sign(x[i]) * 0.375
           -> 0.375 = 0x00006000 in 15Q16
        */
    }

    /* Rate flag */
    rateFlag = 0;
    move16();
    if (sub(totalBits, add(160, i_mult(fs_idx, 160))) > 0)
    {
        rateFlag = 2 << NBITS_CONTEXT;
        move16();
    }

    /* Init */
    nt_half = shr_pos(nt, 1);
    c       = 0;
    move16();
    t = 0;
    move16();
    a1_i = 0;
    move16();
    b1_i = 1;
    move16();
    target32 = L_shl_pos(L_deposit_l(target), SYM_BITS_Q);
    nbits32  = L_negate(target32);
    nbits232 = 0;
    move32();
    nbits_lsb = 0;
    move16();
    IF (mode == 0 && sub(totalBits, add(480, i_mult(fs_idx, 160))) >= 0)
    {
        mode = 1;
        move16();
    }

    /* Find last non-zero tuple */
    lastnz = find_last_nz_pair(xq, nt);
    IF (mode >= 0)
    {
        lastnz2 = 2;
    }
    ELSE
    {
        lastnz2 = lastnz;
    }

    IF (mode < 0)
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }
            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = abs_s(xq[a1_i]);
            b1     = abs_s(xq[b1_i]);
            ab_max = s_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = add(a1, shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a1, b1), 1));
            }
            ELSE IF (sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                nbits32       = L_add(nbits32, 2 << SYM_BITS_Q);
                a1_msb        = shr_pos_pos(a1, 1);
                b1_msb        = shr_pos_pos(b1, 1);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(13, norm_s(ab_max));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(levmax), SYM_BITS_Q + 1));
                a1_msb        = shr(a1, levmax);
                b1_msb        = shr(b1, levmax);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }
    ELSE IF (mode == 0)
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }
            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = abs_s(xq[a1_i]);
            b1     = abs_s(xq[b1_i]);
            ab_max = s_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = add(a1, shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a1, b1), 1));
            }
            ELSE IF (sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                nbits32       = L_add(nbits32, 2 << SYM_BITS_Q);
                a1_msb        = shr_pos_pos(a1, 1);
                b1_msb        = shr_pos_pos(b1, 1);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(13, norm_s(ab_max));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(levmax), SYM_BITS_Q + 1));
                a1_msb        = shr(a1, levmax);
                b1_msb        = shr(b1, levmax);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }
    ELSE
    {
        /* Main Loop through the 2-tuples */
        FOR (k = 0; k < lastnz; k += 2)
        {

            /* Get context */
            t = add(c, rateFlag);
            if (sub(k, nt_half) > 0)
            {
                t = add(t, 1 << NBITS_CONTEXT);
            }
            codingdata[0] = t;
            move16();

            /* Init current 2-tuple encoding */
            a1     = abs_s(xq[a1_i]);
            b1     = abs_s(xq[b1_i]);
            ab_max = s_max(a1, b1);

            IF (ab_max == 0)
            {
                codingdata[1] = -1;
                move16();
                codingdata[2] = 0;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][0]);
                c       = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(ab_max, A_THRES) < 0)
            {
                codingdata[1] = 0;
                move16();
                msb           = add(a1, shl_pos(b1, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][msb]);
                if (a1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1 != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a1, b1), 1));
            }
            ELSE IF (sub(ab_max, 2 * A_THRES) < 0)
            {
                codingdata[1] = 1;
                move16();
                nbits32       = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t]][VAL_ESC]);
                a1_msb        = shr_pos_pos(a1, 1);
                b1_msb        = shr_pos_pos(b1, 1);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[1]]][msb]);
                if (a1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                nbits_lsb = add(nbits_lsb, 2);
                if (sub(a1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (sub(b1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1_msb, b1_msb), 1), 1));
            }
            ELSE
            {
                levmax        = sub(13, norm_s(ab_max));
                codingdata[1] = levmax;
                move16();
                FOR (lev = 0; lev < levmax; lev++)
                {
                    lev1    = s_min(lev, 3);
                    nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][VAL_ESC]);
                }
                nbits32       = L_add(nbits32, L_shl_pos(L_deposit_l(sub(levmax, 1)), SYM_BITS_Q + 1));
                a1_msb        = shr(a1, levmax);
                b1_msb        = shr(b1, levmax);
                msb           = add(a1_msb, shl_pos(b1_msb, A_THRES_SHIFT));
                codingdata[2] = msb;
                move16();
                lev1    = s_min(levmax, 3);
                nbits32 = L_add(nbits32, ari_spec_bits[ari_spec_lookup[t + Tab_esc_nb[lev1]]][msb]);
                a1_msb  = shr_pos(a1, 1);
                b1_msb  = shr_pos(b1, 1);
                if (a1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                if (b1_msb != 0)
                {
                    nbits32 = L_add(nbits32, 1 << SYM_BITS_Q);
                }
                nbits_lsb = add(nbits_lsb, 2);
                if (sub(a1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (sub(b1, 1) == 0)
                {
                    nbits_lsb = add(nbits_lsb, 1);
                }
                if (nbits32 <= 0)
                {
                    lastnz2 = add(k, 2);
                }
                if (nbits32 <= 0)
                {
                    nbits232 = nbits32;
                    move32();
                }
                c = add(shl_pos(s_and(c, 0xf), 4), add(12, s_min(levmax, 3)));
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /* end of the 2-tuples loop */
    }

    /* Number of consumed bits */
    nbits32 = L_add(nbits32, target32);
    totBits = add(extract_l(L_shr_pos_pos(L_sub(nbits32, 1), SYM_BITS_Q)), 1);
    IF (mode > 0)
    {
        totBits = add(totBits, nbits_lsb);
    }
    IF (nBits != NULL)
    {
        *nBits = totBits;
    }
    IF (mode >= 0)
    {
        nbits232 = L_add(nbits232, target32);
        *nBits2  = add(extract_l(L_shr_pos(L_sub(nbits232, 1), SYM_BITS_Q)), 1);
    }
    ELSE
    {
        *nBits2 = *nBits;
        move16();
    }
    IF (mode > 0)
    {
        *nBits2 = add(*nBits2, nbits_lsb);
    }
    *lastnzout = lastnz2;

    /* Truncation of high frequency coefficients */
    IF (lastnz > lastnz2)
    {
        basop_memset(&xq[lastnz2], 0, (lastnz - lastnz2) * sizeof(*xq));
    }

    /* Truncation of LSBs */
    test();
    IF (mode > 0 && sub(totBits, target) > 0)
    {
        *lsbMode = 1;
        move16();
    }
    ELSE
    {
        *lsbMode = 0;
        move16();
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#  endif /* ENABLE_HR_MODE */
