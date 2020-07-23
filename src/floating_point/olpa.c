/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

static void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_FLOAT len_buf, LC3_INT len_input);
static LC3_INT  searchMaxIndice(LC3_FLOAT* in, LC3_INT len);

void filter_olpa(LC3_FLOAT* in, LC3_FLOAT* out, const LC3_FLOAT* buf, LC3_FLOAT len_buf, LC3_INT len_input)
{
    LC3_INT   i = 0, j = 0;
    LC3_FLOAT sum = 0;
    /* a = 1, so denominator == 1, nothing to do here */

    for (i = 0; i < len_input; i++) {
        j   = 0;
        sum = 0;
        for (j = 0; (j < len_buf) && (j <= i); j++) {
            sum += buf[j] * in[i - j];
        }

        out[i] = sum;
    }
}

LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT len)
{
    LC3_INT   max_i = 0, i = 0;
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

void processOlpa_fl(LC3_FLOAT* s_12k8, LC3_FLOAT* mem_s12k8, LC3_FLOAT* mem_s6k4, LC3_INT* mem_old_T0, LC3_INT* T0_out,
                    LC3_FLOAT* normcorr_out, LC3_INT len, LC3_INT frame_dms)
{
    LC3_FLOAT norm_corr = 0, sum = 0, sum0 = 0, sum1 = 0, sum2 = 0, norm_corr2 = 0, *s6k4;
    LC3_FLOAT buf[LEN_6K4 + MAX_PITCH_6K4] = {0}, filt_out[LEN_12K8 + 3] = {0}, d_wsp[LEN_6K4] = {0}, R0[RANGE_PITCH_6K4] = {0}, R[RANGE_PITCH_6K4] = {0}; /* constant length */
    LC3_INT   i = 0, j = 0, len2 = 0, T0 = 0, T02 = 0, min_pitch = 0, max_pitch = 0, L = 0, mem_in_len = 0, acflen = 0;
    

    mem_in_len = MAX_PITCH_6K4;
    len2       = len / 2;
    acflen     = len2;
    if (frame_dms == 25)
    {
        mem_in_len += 16;
        acflen     += 16;
    }

    /* Downsampling */
    move_float(buf, mem_s12k8, 3);
    move_float(&buf[3], s_12k8, len);
    move_float(mem_s12k8, &buf[len], 3);
    filter_olpa(buf, filt_out, olpa_down2, 5, len + 3);
    for (i = 4, j = 0; i < len + 3; i = i + 2) {
        d_wsp[j] = filt_out[i];
        j++;
    }

    /* Compute autocorrelation */
    s6k4 = &buf[mem_in_len];
    move_float(buf, mem_s6k4, mem_in_len);
    move_float(s6k4, d_wsp, len2);
    move_float(mem_s6k4, &buf[len2], mem_in_len);
    if (frame_dms == 25)
    {
        s6k4 = s6k4 - 16;
    }
    for (i = MIN_PITCH_6K4; i <= MAX_PITCH_6K4; i++) {
        sum = 0;
        for (j = 0; j < acflen; j++) {
            sum += s6k4[j] * s6k4[j - i];
        }
        R0[i - MIN_PITCH_6K4] = sum;
    }

    /* Weight autocorrelation and find maximum */
    move_float(R, R0, RANGE_PITCH_6K4);
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
    sum1 = LC3_SQRT(sum1) + LC3_POW(10.0, -5.0);
    norm_corr = sum0 / sum1;
    norm_corr = MAX(0, norm_corr);

    /* Second try in the neighborhood of the previous pitch */
    min_pitch = MAX(MIN_PITCH_6K4, *mem_old_T0 - 4);
    max_pitch = MIN(MAX_PITCH_6K4, *mem_old_T0 + 4);
    L = searchMaxIndice(&R[min_pitch - MIN_PITCH_6K4], max_pitch - min_pitch + 1 );
    T02 = L + min_pitch;

    if (T02 != T0) {
        sum0 = sum1 = sum2 = 0;
        for (i = 0; i < acflen; i++) {
            sum0 += s6k4[i] * s6k4[i - T02];
            sum1 += s6k4[i - T02] * s6k4[i - T02];
            sum2 += s6k4[i] * s6k4[i];
        }
        sum1 = sum1 * sum2;
        sum1 = LC3_SQRT(sum1) + LC3_POW(10.0, -5.0);
        norm_corr2 = sum0 / sum1;
        norm_corr2 = MAX(0, norm_corr2);

        if (norm_corr2 > (norm_corr * 0.85)) {
            T0        = T02;
            norm_corr = norm_corr2;
        }
    }

    *mem_old_T0   = T0;
    *T0_out       = T0 * 2.0;
    *normcorr_out = norm_corr;
}
