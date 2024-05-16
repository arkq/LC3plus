/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void scale_signal24_fx(Word32 x[], /* i:   time input signal */
#ifdef ENABLE_HR_MODE
                       Word32 x_scaled[],
#else
                       Word16 x_scaled[],
#endif
                       Word16 *x_exp,
#ifdef ENABLE_HR_MODE
                       Word32 mdct_mem[],
#else
                       Word16 mdct_mem[],
#endif
                       Word16 mdct_mem_len,
                       Word16 resample_mem_in[], Word16 resample_mem_in_len, Word32 resample_mem_in50[],
                       Word16 resample_mem_out[], Word16 resample_mem_out_len, Word32 mdct_mem32[], Word16 N,
                       Word32 resamp_mem32[], Word16 mem_s12k8[], Word16 *resamp_scale)
{
    Word16 i;
    Word16 s;
    Word16 scales[6];

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("scale_signal24_fx", sizeof(struct {
                   Word16 i;
                   Word16 s;
                   Word16 scales[6];
               }));
#endif

    /* Scale input for 24 bit case */
    /* assure 24 bit input */
    FOR (i = 0; i < N; i++)
    {
        IF (x[i] >= 0)
        {
            x[i] = L_and(x[i], 0x007fffff);
        }
        ELSE
        {
            x[i] = (Word32)L_or( (UWord32)x[i], 0xff800000);
        }
    }

    /* Find maximum exponent */
    scales[0] = sub(15 + 8, getScaleFactor32_0(x, N));
    scales[1] = sub(15 + 8, getScaleFactor32_0(mdct_mem32, mdct_mem_len));
    scales[2] = sub(15 + 8, getScaleFactor32_0(resamp_mem32, resample_mem_in_len));
    scales[3] = sub(sub(*resamp_scale, 2), getScaleFactor32_0(resample_mem_in50, 2));
    scales[4] = sub(sub(*resamp_scale, 2), getScaleFactor16_0(resample_mem_out, resample_mem_out_len));
    scales[5] = sub(sub(*resamp_scale, 2), getScaleFactor16_0(mem_s12k8, 3));
    *x_exp    = 7; move16();
    FOR (i = 0; i < 6; i++)
    {
        *x_exp = s_max(*x_exp, scales[i]); move16();
    }

    /* Shift input buffers */
    s = sub(15 + 8, *x_exp);
    FOR (i = 0; i < N; i++)
    {
#ifdef ENABLE_HR_MODE
        x_scaled[i] = L_shl(x[i], s);
#else
        x_scaled[i] = round_fx_sat(L_shl(x[i], s));
#endif
    }

    FOR (i = 0; i < mdct_mem_len; i++)
    {
#ifdef ENABLE_HR_MODE
        mdct_mem[i] = L_shl(mdct_mem32[i], s);
#else
        mdct_mem[i] = round_fx_sat(L_shl(mdct_mem32[i], s));
#endif
    }

    FOR (i = 0; i < resample_mem_in_len; i++)
    {
        resample_mem_in[i] = round_fx_sat(L_shl(resamp_mem32[i], s));
    }

    /* Adjust resampler filter and output buffers */
    s             = sub(sub(*resamp_scale, 2), *x_exp);
    *resamp_scale = add(*x_exp, 2);

    IF (s)
    {
        FOR (i = 0; i < 2; i++)
        {
            resample_mem_in50[i] = L_shl(resample_mem_in50[i], s);
        }
        FOR (i = 0; i < resample_mem_out_len; i++)
        {
            resample_mem_out[i] = shl(resample_mem_out[i], s);
        }

        FOR (i = 0; i < 3; i++)
        {
            mem_s12k8[i] = shl(mem_s12k8[i], s);
        }
    }
    /* Store part of current frame as mdct memory buffer and resampler input buffer for next frame */
    basop_memcpy(mdct_mem32, &x[N - mdct_mem_len], mdct_mem_len * sizeof(Word32));
    basop_memmove(resamp_mem32, &x[N - resample_mem_in_len], resample_mem_in_len * sizeof(Word32));

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

