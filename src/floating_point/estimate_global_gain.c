/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"


void processEstimateGlobalGain_fl(LC3_FLOAT x[], LC3_INT lg, LC3_INT nbitsSQ, LC3_FLOAT* gain, LC3_INT* quantizedGain,
                                  LC3_INT* quantizedGainMin, LC3_INT quantizedGainOff, LC3_FLOAT* targetBitsOff,
                                  LC3_INT* old_targetBits, LC3_INT old_specBits
                                  , LC3_INT hrmode , LC3_INT regBits, LC3PLUS_FrameDuration frame_ms
)
{

    LC3_INT   i, N, offset, j, iszero, fac;
    LC3_FLOAT g_min, x_max, tmp, ind, ind_min, target, ener;
    LC3_FLOAT en[MAX_LEN / 4];
    LC3_FLOAT reg_val = 4.656612873077393e-10;
#ifdef    FIX_1p25_GG_EST_TUPLES
    LC3_INT tuples[1 + 4] = { -1, GG_1p25_WB_TUPLES, GG_1p25_SSWB_TUPLES, GG_1p25_SWB_TUPLES, GG_1p25_FB_TUPLES };
    LC3_INT bw_idx;
#endif 

#ifdef   FIX_1p25_GG_EST_TUPLES  
    LC3_INT32 lg_extra;
    if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS) 
    {
        assert((MAX_LEN) >= ((lg / GG_1p25_MAX_TUPLES) + 1)*GG_1p25_MAX_TUPLES);   /* en size check, for 1p25ms   max tuple size*/
    }
#else 
#  ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
    assert( (MAX_LEN / 4) >= ((8*4 + 2) / 2));   /* en size check, for 1p25ms  SSWB */
#  endif 
#endif 

#ifdef    FIX_1p25_GG_EST_TUPLES
    lg_extra = 0;    
    bw_idx = MIN((lg / 10) - 1, 4);
 
    if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS)
    {  
        lg_extra = tuples[bw_idx] * (( lg + (tuples[bw_idx] - 1)) / tuples[bw_idx]) - lg; 
        /*x buffer is assert checked  in enc_lc3_fl.c() for 1.25ms  energy calculation */
        for (i = lg; i < (lg + lg_extra); i++)
        {
            x[i] = 0.0; /* zero extended  tail if a truncated tuple-block exists */
        }
    }
#else 
#  ifdef  FIX_FLOAT_ENC_QUANTIZE_1P25MS_512KBPS
    LC3_INT32 lg_extra;
    
    lg_extra = lg - 4 * (lg / 4); 
    for ( i=lg; i< (lg+lg_extra); i++)
    {
         x[i] = 0.0; /* zero tail if a truncated quadruple exists */
    }
#  endif
#endif 
    if (*old_targetBits < 0) {
        *targetBitsOff = 0;
    } else {
        tmp            = MIN(40, MAX(-40, *targetBitsOff + *old_targetBits - old_specBits));
        *targetBitsOff = 0.8 * *targetBitsOff + 0.2 * tmp;
    }

    *old_targetBits = nbitsSQ;
    nbitsSQ         = nbitsSQ + round(*targetBitsOff);

    x_max = array_max_abs(x, lg);

    if (hrmode && regBits > 0)
    {
        LC3_FLOAT M0 = 1e-5, M1 = 1e-5, rB_offset;
        LC3_FLOAT thresh = 2*frame_ms*1.25;
        for (i = 0; i < lg; i++)
        {
                M0 += fabs(x[i]);
                M1 += i*fabs(x[i]);
        }

        rB_offset = 8 * (1 - MIN(M1/M0, thresh) / thresh);
        reg_val = x_max * LC3_POW(2,-regBits - rB_offset);
    }

    if (x_max < LC3_EPS)
    {
        ind_min         = quantizedGainOff;
        ind             = 0;
        *old_targetBits = -1;
    }
    else {
        if (hrmode == 1) {
            g_min = x_max / (32768 * 256 - 2);
        }
        else {
            g_min = x_max / (32767 - 0.375);
        }
        /* Prevent positive rounding errors from LC3_LOG10 function */
        ind_min = 28.0 * LC3_LOGTEN(g_min);

        ind_min = ceil(ind_min + LC3_FABS(ind_min) * LC3_EPS);

        assert(LC3_POW(10, ind_min / 28.0) >= g_min);
        assert(ind_min <= (255 + quantizedGainOff));

#ifdef  FIX_FLOAT_ENC_QUANTIZE_1P25MS_512KBPS
        N = lg + lg_extra;
#else 
        N = lg;
#endif


#ifndef FIX_1p25_GG_EST_TUPLES
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2
        /*  increase to at least 10 analysis bands for  WB, SSWB 1p25ms */
        if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS && lg <= 30)
        {
            j = 0;
            for (i = 0; i < N; i = i + 2) /* steps of 2  for 1.25ms frame lengths */
            {
                tmp = x[i] * x[i];
                tmp += x[i + 1] * x[i + 1];
                en[j] = (28.0 / 20.0) * (7 * 0.5 + 10.0 * LC3_LOGTEN(tmp + reg_val)); /*  offset of 7 per 4 coeff band   now changed to 3.5 per 2 coeff band */
                j++;
            }

            target = (28.0 / 20.0) * (1.4) * nbitsSQ;  /* global index domain  sum over all 10 bands  */

            fac = 256;
            offset = 255 + quantizedGainOff;

            for (i = 0; i < 8; i++)
            {
                fac = fac >> 1;
                offset = offset - fac;
                ener = 0;
                iszero = 1;

                /* we sum up energy from the top,  to not add up noisefilled coeffs */
                for (j = N / 2 - 1; j >= 0; j--)
                {
                    tmp = en[j] - offset;

                    if (tmp < ((7.0) * (28.0 / 20.0) * 0.5))
                    {
                        if (iszero == 0)
                        {
                            ener = ener + ((2.7) * (28.0 / 20.0)* 0.5);  /* low cost in coded band zero */
                        }
                    }
                    else
                    {
                        if (tmp > ((50.0) * (28.0 / 20.0)*0.5))
                        {
                            /* high value with many escapes */

                            ener = ener + 2.0 * tmp - (50.0) * (28.0 / 20.0)*0.5;

                        }
                        else
                        {
                            ener = ener + (tmp*1.0);
                        }
                        iszero = 0;
                    }
                } /* loop over over band  N/2-1 ...  0   */

                if (ener > target && iszero == 0)
                {
                    offset = offset + fac;
                }

            } /* over 8 splits/iterations to handle all 256 possible shifts  */

            if (offset < ind_min)
            {
                *old_targetBits = -1;
            }

            ind = MAX(ind_min, offset) - quantizedGainOff;
        }
        else
#endif
#endif 



#ifdef FIX_1p25_GG_EST_TUPLES

#if  GG_1p25_MAX_TUPLES == 2   
            /*  the tuple/2-block  with halved limits, is here separated from the  4-block loop ,
                 optionally they can be parametrized into one function */
            if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS)
#else 
            if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS && tuples[bw_idx] == 2)
#endif 
            {
                assert(((lg / 2) * 2) == lg);

                j = 0;
                for (i = 0; i < N; i = i + 2)
                {
                    tmp = x[i] * x[i];
                    tmp += x[i + 1] * x[i + 1];
                    en[j] = (28.0 / 20.0) * (7 * 0.5 + 10.0 * LC3_LOGTEN(tmp + reg_val)); /*  offset of 7 per 4 coeff band   now changed to 3.5 per 2 coeff band */
                    j++;
                }

                target = (28.0 / 20.0) * (1.4) * nbitsSQ;  /* global index domain  sum over all bands  */

                fac = 256;
                offset = 255 + quantizedGainOff;

                for (i = 0; i < 8; i++)
                {
                    fac = fac >> 1;
                    offset = offset - fac;
                    ener = 0;
                    iszero = 1;

                    /* we sum up energy from the top,  to not add up noisefilled coeffs */
                    for (j = N / 2 - 1; j >= 0; j--)
                    {
                        tmp = en[j] - offset;

                        if (tmp < ((7.0) * (28.0 / 20.0) * 0.5))
                        {
                            if (iszero == 0)
                            {
                                ener = ener + ((2.7) * (28.0 / 20.0)* 0.5);  /* low cost in coded band zero */
                            }
                        }
                        else
                        {
                            if (tmp > ((50.0) * (28.0 / 20.0)*0.5))
                            {
                                /* high value with many escapes */

                                ener = ener + 2.0 * tmp - (50.0) * (28.0 / 20.0)*0.5;

                            }
                            else
                            {
                                ener = ener + (tmp*1.0);
                            }
                            iszero = 0;
                        }
                    } /* loop over over band  N/2-1 ...  0   */

                    if (ener > target && iszero == 0)
                    {
                        offset = offset + fac;
                    }

                } /* over 8 splits/iterations to handle all 256 possible shifts  */

                if (offset < ind_min)
                {
                    *old_targetBits = -1;
                }

                ind = MAX(ind_min, offset) - quantizedGainOff;
            }
            else

#if  GG_1p25_MAX_TUPLES == 3   
                /*  the 3tuple -block  with modified  limits, is here separated from the quadruple/4-block loop ,
                     optionally they can be parametrized into one function */
                if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS)
#else 

                if (frame_ms == LC3PLUS_FRAME_DURATION_1p25MS && tuples[bw_idx] == 3)
#endif 
                {   /* 3 tuple loop */
                    N = 3 * (int)ceil(lg / 3.0);

                    assert(lg == 50 || lg == 40 || lg == 30 || lg == 20);
                    assert(N == (17 * 3) || N == (14 * 3) || N == (10 * 3) || N == (7 * 3));

                    j = 0;
                    for (i = 0; i < N; i = i + 3)
                    {
                        tmp = x[i] * x[i];
                        tmp += x[i + 1] * x[i + 1];
                        tmp += x[i + 2] * x[i + 2];
                        en[j] = (28.0 / 20.0) * (7 * 0.75 + 10.0 * LC3_LOGTEN(tmp + reg_val)); /*  offset of 7*0.75  per  3 coeff band   n */
                        j++;
                    }

                    target = (28.0 / 20.0) * (1.4) * nbitsSQ;  /* global index domain  sum over all bands  */

                    fac = 256;
                    offset = 255 + quantizedGainOff;

                    for (i = 0; i < 8; i++)
                    {
                        fac = fac >> 1;
                        offset = offset - fac;
                        ener = 0;
                        iszero = 1;

                        /* we sum up energy from the top,  to not add up noisefilled coeffs */
                        for (j = N / 3 - 1; j >= 0; j--)
                        {
                            tmp = en[j] - offset;

                            if (tmp < ((7.0) * (28.0 / 20.0) * 0.75))
                            {
                                if (iszero == 0)
                                {
                                    ener = ener + ((2.7) * (28.0 / 20.0)* 0.75);  /* low cost in coded band zero */
                                }
                            }
                            else
                            {
                                if (tmp > ((50.0) * (28.0 / 20.0)*0.75))
                                {
                                    /* high value with many escapes */

                                    ener = ener + 2.0 * tmp - (50.0) * (28.0 / 20.0)*0.75;

                                }
                                else
                                {
                                    ener = ener + (tmp*1.0);
                                }
                                iszero = 0;
                            }
                        } /* loop over over band  N/3-1 ...  0   */

                        if (ener > target && iszero == 0)
                        {
                            offset = offset + fac;
                        }

                    } /* over 8 splits/iterations to handle all 256 possible shifts  */

                    if (offset < ind_min)
                    {
                        *old_targetBits = -1;
                    }

                    ind = MAX(ind_min, offset) - quantizedGainOff;
                }
                else

#endif  /* end  2, 3  tuple loops */

#ifdef    FIX_1p25_GG_EST_TUPLES
                {  /* brackets for 4 tuple loop */
                    assert((((lg + lg_extra) / 4) * 4) == (lg + lg_extra));
#endif
#ifndef    FIX_1p25_GG_EST_TUPLES
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
                    {
#endif
#endif 

                        j = 0;
                        for (i = 0; i < N; i = i + 4) {
                            tmp = x[i] * x[i];
                            tmp += x[i + 1] * x[i + 1];
                            tmp += x[i + 2] * x[i + 2];
                            tmp += x[i + 3] * x[i + 3];
                            en[j] = (28.0 / 20.0) * (7 + 10.0 * LC3_LOGTEN(tmp + reg_val));
                            j++;
                        }

                        target = (28.0 / 20.0) * (1.4) * nbitsSQ;
                        fac = 256;
                        offset = 255 + quantizedGainOff;

                        for (i = 0; i < 8; i++)
                        {
                            fac = fac >> 1;
                            offset = offset - fac;
                            ener = 0;
                            iszero = 1;

                            for (j = N / 4 - 1; j >= 0; j--) {
                                tmp = en[j] - offset;

                                if (tmp < (7.0) * (28.0 / 20.0)) {
                                    if (iszero == 0) {
                                        ener = ener + (2.7) * (28.0 / 20.0);
                                    }
                                }
                                else {
                                    if (tmp > (50.0) * (28.0 / 20.0)) {
                                        ener = ener + 2.0 * tmp - (50.0) * (28.0 / 20.0);
                                    }
                                    else {
                                        ener = ener + tmp;
                                    }

                                    iszero = 0;
                                }
                            }

                            if (ener > target && iszero == 0) {
                                offset = offset + fac;
                            }
                        }

                        if (offset < ind_min) {
                            *old_targetBits = -1;
                        }

                        ind = MAX(ind_min, offset) - quantizedGainOff;
                }
#ifdef    FIX_1p25_GG_EST_TUPLES
    }
  /* end bracket for 4 tuple loop */
#else 
#ifdef     FIX_BOTH_1p25_TEST_NEW_GG_EST2 
    }
#endif
#endif

    *quantizedGainMin = ind_min;
    *quantizedGain = ind;

    *gain = LC3_POW(10.0, ((ind + quantizedGainOff) / 28.0));
}
