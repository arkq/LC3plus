/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef CR9_C_ADD_1p25MS
static float limitShaping (LC3_FLOAT* xl4 ) {
    LC3_FLOAT fac;
    LC3_FLOAT score;
    LC3_FLOAT min_fac;
    LC3_FLOAT max_fac;
    LC3_FLOAT start;
    LC3_FLOAT stop;
    
    min_fac = 1.f;
    max_fac = 0.3f;
    start   = 5.f;
    stop    = 8.f;

    score = ((2.f*(xl4[0]-xl4[1])) + (xl4[0]-xl4[2])) / 2.f;
    score = fmin(fmax(score, start), stop);

    fac = (stop-score)/(stop-start);
    return ((min_fac - max_fac) * fac + max_fac);
}
#endif

void processSnsComputeScf_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* gains, LC3_INT smooth, LC3_FLOAT sns_damping, LC3_FLOAT attdec_damping_factor, LC3_INT fs_idx
#ifdef CR9_C_ADD_1p25MS
                            , LC3PLUS_FrameDuration frame_dms, LC3_FLOAT *LT_normcorr, LC3_FLOAT normcorr
#endif
                            )
{
    LC3_INT   bands_number, d, i, j, n, n2, n4, mapping[64];
    LC3_FLOAT x_tmp1[MAX_LEN], sum = 0, mean, nf, gains_smooth[M], ratio;
    LC3_FLOAT sum_gains_smooth = 0;
#ifdef CR9_C_ADD_1p25MS
    LC3_FLOAT fac;
    LC3_FLOAT start;
    LC3_FLOAT limiterGain;
#endif
    const LC3_FLOAT *sns_preemph_adapt, *sns_preemph;
    bands_number = xLen;
    
    sns_preemph = sns_preemph_all[fs_idx];
    
#ifdef CR9_C_ADD_1p25MS
    limiterGain = 1.f;
    sns_preemph_adapt = sns_preemph_adaptMaxTilt_all[fs_idx];
#else
    (void) sns_preemph_adapt;
#endif

    /* 5 ms */
    if (bands_number < 64) {
        d = 64 - bands_number;

        if (d < xLen)
        {
            j = 0;
            for (i = 0; i < 2 * d; i = i + 2) {
                x_tmp1[i]     = x[j];
                x_tmp1[i + 1] = x[j];
                j++;
            }

            move_float(&x_tmp1[2 * d], &x[d], 64 - 2 * d);
        }
        else if (bands_number < 32)
        {
            ratio = LC3_FABS((LC3_FLOAT) (1.0 - 32.0 / (LC3_FLOAT) xLen));
            n4 = round(ratio * xLen);
            n2 = xLen - n4;
            
            j = 0;
            for(i = 1; i <= n4; i++)
            {
                mapping[j] = i;
                mapping[j + 1] = i;
                mapping[j + 2] = i;
                mapping[j + 3] = i;
                j += 4;
            }
            
            for (i = n4 + 1; i <= n4 + n2; i++)
            {
                mapping[j] = i;
                mapping[j + 1] = i;
                j += 2;
            }
            
            
            for (i = 0; i < 64; i++)
            {
                x_tmp1[i] = x[mapping[i] - 1];
            }
        } else {
            assert(0 && "Unsupported number of bands!");
        }

        move_float(x, x_tmp1, 64);

        bands_number = 64;
        xLen         = bands_number;
    }


    /* Smoothing */

    x_tmp1[0] = x[0];
    move_float(&x_tmp1[1], &x[0], 63);

    for (i = 0; i < 63; i++) {
        x[i] = 0.5 * x[i] + 0.25 * (x_tmp1[i] + x[i + 1]);
    }
    
    x[63] = 0.5 * x[63] + 0.25 * (x_tmp1[63] + x[63]);

    /* Pre-emphasis */
#ifdef CR9_C_ADD_1p25MS
    if (sns_preemph_adapt != NULL && frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        *LT_normcorr = normcorr * 0.125f + *LT_normcorr * 0.875f;

        start = 0.8f;   /* adaptive preemphasis active from start to 1.0 */
        fac = (fmax(*LT_normcorr - start, 0) * (1. / (1.-start)));

        for (i = 0; i < 64; i++) {
            x[i] = x[i] * (sns_preemph[i] + fac*sns_preemph_adapt[i]);
        }
    } 
    else
    {
        for (i = 0; i < 64; i++) {
            x[i] = x[i] * sns_preemph[i];
        }
    }
#else  
    for (i = 0; i < 64; i++) {
        x[i] = x[i] * sns_preemph[i];
    }
#endif

    /* Noise floor at -40dB */
    for (i = 0; i < 64; i++) {
        sum += x[i];
    }

    mean = sum * 0.015625; /* 1/64 */

    nf = mean * 1.00e-04;
    nf = MAX(nf, 2.328306436538696e-10);

    for (i = 0; i < 64; i++) {
        if (x[i] < nf) {
            x[i] = nf;
        }
    }

    /* Log-domain */
    for (i = 0; i < 64; i++) {
        x[i] = LC3_LOGTWO(x[i]) * 0.5;
    }

    /* Downsampling */
    for (n = 0; n < 16; n++) {
        if (n == 0) {
            x_tmp1[0] = x[0];

            move_float(&x_tmp1[1], &x[0], 5);

        } else if (n == 15) {
            move_float(x_tmp1, &x[59], 5);

            x_tmp1[5] = x[63];

        } else {
            move_float(x_tmp1, &x[n * 4 - 1], ((n * 4 + 5 - 1) - (n * 4 - 1) + 1));
        }

        sum = 0;
        for (i = 0; i < 6; i++) {
            sum += x_tmp1[i] * sns_W[i];
        }

        gains_smooth[n] = sum;
        sum_gains_smooth += sum;
    }
#ifdef CR9_C_ADD_1p25MS
    /* limit shaping */
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
        limiterGain *= limitShaping(gains_smooth);
    }
#endif
    /* Remove mean and scaling */
    mean = sum_gains_smooth / 16.0;

    for (i = 0; i < 16; i++) {
#ifdef CR9_C_ADD_1p25MS
        gains[i] = limiterGain * sns_damping * (gains_smooth[i] - mean);
#else
        gains[i] = sns_damping * (gains_smooth[i] - mean);
#endif
    }

    /* Smoothing */
#ifdef CR9_C_ADD_1p25MS_LRSNS 
     if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {   /* smoothing loop for 1.25 ms  */

        const LC3_FLOAT A0 = 3.0 / 16.0;  /* 2/16= 0.125,   3/16 = 0.1875  */
        const LC3_FLOAT A1 = (1.0 - 2 * A0);
        const LC3_FLOAT A2 = A0;

        gains_smooth[0] = A0 * (gains[0] + gains[1]) * 0.5 + A1 * gains[0] + A2 * gains[1];
        /* BASOP-loop:: preload gains[-1] with  0.5*(gains[0]+gains[1]) */
        for (i = 1; i < (M - 1); i++)
        {
            gains_smooth[i] = A0 * gains[i - 1] + A1 * gains[i] + A2 * gains[i + 1];
        }
        gains_smooth[M - 1] = A0 * gains[M - 2] + A1 * gains[M - 1] + A2 * 0.5 * (gains[M - 2] + gains[M - 1]);
        /* BASOP-loop :: preload gains[M] with  0.5*(gains[M-2]+gains[M-1]) */
        sum = 0;
        for (i = 0; i < M; i++)
        {
            sum += gains_smooth[i];
        }

        mean = sum / (LC3_FLOAT)M;
        for (i = 0; i < M; i++)
        {
            gains[i] = attdec_damping_factor * (gains_smooth[i] - mean);
        }
    }
    else if (smooth == 1) /* original attack smoothing */
#    else
    if (smooth)
#    endif
    {
        gains_smooth[0] = (gains[0] + gains[1] + gains[2]) / 3.0;
        gains_smooth[1] = (gains[0] + gains[1] + gains[2] + gains[3]) / 4.0;

        for (i = 2; i < 14; i++) {
            gains_smooth[i] = (gains[i - 2] + gains[i - 1] + gains[i] + gains[i + 1] + gains[i + 2]) / 5.0;
        }

        gains_smooth[M - 2] = (gains[M - 4] + gains[M - 3] + gains[M - 2] + gains[M - 1]) / 4.0;
        gains_smooth[M - 1] = (gains[M - 3] + gains[M - 2] + gains[M - 1]) / 3.0;

        sum = 0;
        for (i = 0; i < M; i++) {
            sum += gains_smooth[i];
        }

        mean = sum / (LC3_FLOAT)M;

        for (i = 0; i < M; i++) {
            gains[i] = attdec_damping_factor * (gains_smooth[i] - mean);
        }
    }
}
