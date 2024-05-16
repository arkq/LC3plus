/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processSnsComputeScf_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* gains, LC3_INT smooth, LC3_FLOAT sns_damping, LC3_FLOAT attdec_damping_factor, LC3_INT fs_idx)
{
    LC3_INT   bands_number, d, i, j, n, n2, n4, mapping[64];
    LC3_FLOAT x_tmp1[MAX_LEN], sum, mean, nf, gains_smooth[M], ratio;
    LC3_FLOAT sum_gains_smooth;
    const LC3_FLOAT* sns_preemph;
    
    sum_gains_smooth = 0; sum = 0;
    sns_preemph = sns_preemph_all[fs_idx];
    
    bands_number = xLen;
    assert(bands_number <= 64);

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
    for (i = 0; i < 64; i++) {
        x[i] = x[i] * sns_preemph[i];
    }

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


    /* Remove mean and scaling */
    mean = sum_gains_smooth / 16.0;

    for (i = 0; i < 16; i++) {
        gains[i] = sns_damping * (gains_smooth[i] - mean);
    }

    /* Smoothing */
    if (smooth) {
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
