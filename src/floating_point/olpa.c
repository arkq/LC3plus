/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_INT32 len_input);
static LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT32 len);

void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_INT32 len_input)
{
    /* a = 1, so denominator == 1, nothing to do here */
    LC3_INT32 i, j;


    j = 0;
    for (i = 4; i < len_input; i += 2) {
        out[j++] = (buf[0] * in[i]) + (buf[1] * in[i - 1]) + (buf[2] * in[i - 2]) + (buf[3] * in[i - 3]) + (buf[4] * in[i - 4]);
    }
}

LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT len)
{
    LC3_INT   max_i = 0, i;
    LC3_FLOAT max = in[0];

    if (len <= 0) {
        return -128;
    }

    for (i = 0; i < len; i++) {
        if (in[i] > max) {
            max   = in[i];
            max_i = i;
        }
    }

    return max_i;
}

void processOlpa_fl(LC3_FLOAT* s_12k8, LC3_FLOAT* mem_s12k8, LC3_FLOAT* mem_s6k4, LC3_INT* mem_old_T0, 
                    LC3_INT* pitch_flag, 
                    LC3_INT* T0_out, LC3_FLOAT* normcorr_out, LC3_INT len, LC3_INT frame_dms)
{
    LC3_FLOAT norm_corr = 0, sum = 0, sum0 = 0, sum1 = 0, sum2 = 0, norm_corr2 = 0, *s6k4;
    LC3_FLOAT buf[LEN_6K4 + MAX_PITCH_6K4 + MAX_LEN], R0[RANGE_PITCH_6K4]; /* constant length */
    LC3_INT   i = 0, len2 = 0, T0 = 0, T02 = 0, min_pitch = 0, max_pitch = 0, L = 0, mem_in_len = 0, acflen = 0, delta = 0;

    len2       = len / 2;
    switch(frame_dms)
    {
        case 50:
            delta = len / 2;
            acflen = len2 * 2;
            break;

        case 25:
            delta = 3*(len /2);
            acflen = len2*4;
            break;

        default:
    delta      = 0;
    acflen     = len2;
    }

    mem_in_len = MAX_PITCH_6K4 + delta;

    /* Downsampling */
    move_float(buf, mem_s12k8, 3);
    move_float(&buf[3], s_12k8, len);
    move_float(mem_s12k8, &buf[len], 3);
    filter_olpa(buf, R0, olpa_down2, len + 3);

    /* Compute autocorrelation */
    s6k4 = &buf[mem_in_len - delta];
    move_float(&buf[mem_in_len], R0, len2);
    move_float(buf, mem_s6k4, mem_in_len);
    move_float(mem_s6k4, &buf[len2], mem_in_len);
    for (i = MIN_PITCH_6K4; i <= MAX_PITCH_6K4; i++) {
        sum = mac_loop(s6k4, &s6k4[-i], acflen);
        R0[i - MIN_PITCH_6K4] = sum;
    }

    /* Weight autocorrelation and find maximum */
    
    /* Second try in the neighborhood of the previous pitch */
    min_pitch = MAX(MIN_PITCH_6K4, *mem_old_T0 - 4);
    max_pitch = MIN(MAX_PITCH_6K4, *mem_old_T0 + 4);

    L = searchMaxIndice(&R0[min_pitch - MIN_PITCH_6K4], max_pitch - min_pitch + 1 );
    T02 = L + min_pitch;
    
    for (i = 0; i < RANGE_PITCH_6K4; i++) {
        R0[i] = R0[i] * olpa_acw[i];
    }
    L  = searchMaxIndice(R0, RANGE_PITCH_6K4);
    T0 = L + MIN_PITCH_6K4;

    /* Compute normalized correlation */
    sum0 = sum1 = sum2 = 0;

    for (i = 0; i < acflen; i++) {
        sum0 += s6k4[i] * s6k4[i - T0];
        sum1 += s6k4[i - T0] * s6k4[i - T0];
        sum2 += s6k4[i] * s6k4[i];
    }
    sum1 = sum1 * sum2;
    sum1 = LC3_SQRT(sum1) + 1.00e-05;
    norm_corr = sum0 / sum1;
    norm_corr = MAX(0, norm_corr);

    if (T02 != T0) {
        sum0 = sum1 = sum2 = 0;
        for (i = 0; i < acflen; i++) {
            sum0 += s6k4[i] * s6k4[i - T02];
            sum1 += s6k4[i - T02] * s6k4[i - T02];
            sum2 += s6k4[i] * s6k4[i];
        }
        sum1 = sum1 * sum2;
        sum1 = LC3_SQRT(sum1) + 1.00e-05;
        norm_corr2 = sum0 / sum1;
        norm_corr2 = MAX(0, norm_corr2);

        if (norm_corr2 > (norm_corr * 0.85)) {
            T0        = T02;
            norm_corr = norm_corr2;
        }
    }

    switch(frame_dms)
    {
        case 50:
            if (*pitch_flag == 1)
            {
                *mem_old_T0   = T0;
                *pitch_flag = 0;
            }
            else
            {
                *pitch_flag += 1;
            }
            break;

        case 25:
            if (*pitch_flag == 3)
            {
                *mem_old_T0   = T0;
                *pitch_flag = 0;
            }
            else
            {
                *pitch_flag += 1;
            }
            break;

        default:
    *mem_old_T0   = T0;
    }

    *T0_out       = T0 * 2.0;
    *normcorr_out = norm_corr;

}
