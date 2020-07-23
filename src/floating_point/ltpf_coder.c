/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

static LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT len);

LC3_INT searchMaxIndice(LC3_FLOAT* in, LC3_INT len)
{
    LC3_INT   max_i = 0, i = 0;
    LC3_FLOAT max = 0;

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

void process_ltpf_coder_fl(LC3_FLOAT* xin, LC3_INT xLen, LC3_INT ltpf_enable, LC3_INT pitch_ol, LC3_FLOAT pitch_ol_norm_corr, LC3_INT frame_dms,
                           LC3_FLOAT* mem_old_x, LC3_INT memLen, LC3_FLOAT* mem_norm_corr_past, LC3_INT* mem_on, LC3_FLOAT* mem_pitch,
                           LC3_INT* param, LC3_FLOAT* mem_norm_corr_past_past, LC3_INT* bits)
{
    LC3_FLOAT buffer[LTPF_MEMIN_LEN + LEN_12K8 + 1 + (LEN_12K8 >> 2)] = {0}, sum = 0, buf_tmp[MAX_LEN] = {0}, cor_up[MAX_LEN] = {0}, *x;
    LC3_INT   i = 0, j = 0, k = 0, n = 0, step = 0, N = 0, ltpf_active = 0, pitch_search_delta = 0,
        pitch_search_upsamp = 0, pitch_search_L_interpol1 = 0,
        t0_min = 0, t0_max = 0, t_min = 0, t_max = 0, temp2 = 0, t1 = 0, pitch_int = 0, pitch_fr = 0, midpoint = 0,
        delta_up = 0, delta_down = 0, pitch_index = 0, gain = 0, acflen = 0;
    LC3_FLOAT norm_corr = 0, cor[MAX_LEN] = {0}, cor_int[MAX_LEN] = {0}, currFrame[MAX_LEN] = {0}, predFrame[MAX_LEN] = {0}, sum1 = 0, sum2 = 0, sum3 = 0;
    LC3_FLOAT pitch = 0;

    /* Signal Buffer */
    N = xLen - 1;
    x = &buffer[memLen];

    move_float(buffer, mem_old_x, memLen);
    move_float(x, xin, xLen);
    move_float(mem_old_x, &buffer[N], xLen + memLen - N);

    ltpf_active = 0;
    norm_corr   = 0;

    pitch_search_delta       = 4;
    pitch_search_upsamp      = 4;
    pitch_search_L_interpol1 = 4;

    if (pitch_ol_norm_corr > 0.6) {
        /* Search Bounds */
        t0_min = pitch_ol - pitch_search_delta;
        t0_max = pitch_ol + pitch_search_delta;
        t0_min = MAX(t0_min, MIN_PITCH_12K8);
        t0_max = MIN(t0_max, MAX_PITCH_12K8);
        acflen = N;
        
        if (frame_dms == 25)
        {
            acflen = 2 * N;
            x = x - N;
        }

        /* Cross-Correlation Bounds */
        t_min = t0_min - pitch_search_L_interpol1;
        t_max = t0_max + pitch_search_L_interpol1;

        /* Compute norm */
        sum1 = sum2 = 0;
        for (j = 0; j < acflen; j++) {
            sum1 += x[j] * x[j];
            sum2 += x[j - t_min] * x[j - t_min];
        }

        /* Compute Cross-Correlation */
        for (i = t_min; i <= t_max; i++) {
            sum = 0;
            for (j = 0; j < acflen; j++) {
                sum += x[j] * x[j - i];
            }

            if (i > t_min) {
                sum2 = sum2 + x[-i]*x[-i]
                            - x[acflen - 1 - ( i - 1 )]*x[acflen - 1 - ( i - 1 )];
            }
            sum3      = LC3_SQRT(sum1 * sum2) + LC3_POW(10, -5);
            norm_corr = sum / sum3;
            
            norm_corr = MAX(0, norm_corr);
            cor[i - t_min] = norm_corr;

        }

        /* Find Integer Pitch-Lag */
        j = 0;
        for (i = pitch_search_L_interpol1; i <= t_max - t_min - pitch_search_L_interpol1; i++) {
            buf_tmp[j] = cor[i];
            j++;
        }

        temp2 = searchMaxIndice(buf_tmp, j);

        t1 = temp2 + t0_min;
        assert(t1 >= t0_min && t1 <= t0_max);

        /* Find Fractional Pitch-Lag */
        if (t1 >= RES2_PITCH_12K8) {
            pitch_int = t1;
            pitch_fr  = 0;
        } else {
            j = 0;
            
            for (i = 0; i < pitch_search_upsamp * (t_max - t_min) + 1; i = i + pitch_search_upsamp) {
                cor_up[i] = cor[j];
                j++;
            }

            for (i = 0; i < pitch_search_upsamp * (t0_max - t0_min + 1); i++) {
                sum = 0;

                k = 0;
                for (j = i; j < i + 32; j++) {
                    sum += cor_up[j] * inter4_1[k];
                    k++;
                }

                cor_int[i] = sum;
            }

            if (t1 >= RES4_PITCH_12K8) {
                step = 2;
            } else {
                step = 1;
            }

            midpoint = pitch_search_upsamp * (t1 - t0_min) + 1;
            delta_up = pitch_search_upsamp - step;

            if (t1 == t0_min) {
                delta_down = 0;
            } else {
                delta_down = pitch_search_upsamp - step;
            }

            j = 0;
            for (i = midpoint - delta_down - 1; i <= midpoint + delta_up; i = i + step) {
                buf_tmp[j] = cor_int[i];
                j++;
            }

            temp2    = searchMaxIndice(buf_tmp, ((midpoint + delta_up) - (midpoint - delta_down)) / step + 1);
            pitch_fr = temp2 * step - delta_down;

            if (pitch_fr >= 0) {
                pitch_int = t1;
            } else {
                pitch_int = t1 - 1;
                pitch_fr  = pitch_search_upsamp + pitch_fr;
            }
        }

        assert((pitch_int <= MAX_PITCH_12K8 && pitch_int >= RES2_PITCH_12K8 && pitch_fr == 0) ||
               (pitch_int < RES2_PITCH_12K8 && pitch_int >= RES4_PITCH_12K8 && (pitch_fr == 0 || pitch_fr == 2)) ||
               (pitch_int < RES4_PITCH_12K8 && pitch_int >= MIN_PITCH_12K8 &&
                (pitch_fr == 0 || pitch_fr == 1 || pitch_fr == 2 || pitch_fr == 3)));
        
        if (pitch_int < RES4_PITCH_12K8) {
            pitch_index = pitch_int * 4 + pitch_fr - (MIN_PITCH_12K8 * 4);
        } else if (pitch_int < RES2_PITCH_12K8) {
            pitch_index = pitch_int * 2 + floor(pitch_fr / 2) - (RES4_PITCH_12K8 * 2) + ((RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4);
        } else {
            pitch_index = pitch_int - RES2_PITCH_12K8 + ((RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4) + ((RES2_PITCH_12K8 - RES4_PITCH_12K8) * 2);
        }

        assert(pitch_index >= 0 && pitch_index < 512);
        pitch = (LC3_FLOAT) pitch_int + (LC3_FLOAT) pitch_fr / 4.0;
        

        for (n = 0; n < acflen; n++)
        {
            currFrame[n] = x[n + 1] * enc_inter_filter[0][0] +
                           x[n]     * enc_inter_filter[0][1] +
                           x[n - 1] * enc_inter_filter[0][2];

            predFrame[n] = x[n - pitch_int + 1] * enc_inter_filter[pitch_fr][0] +
                           x[n - pitch_int]     * enc_inter_filter[pitch_fr][1] +
                           x[n - pitch_int - 1] * enc_inter_filter[pitch_fr][2] +
                           x[n - pitch_int - 2] * enc_inter_filter[pitch_fr][3];
        }

        /* Normalized Correlation */
        sum1 = sum2 = sum3 = 0;

        for (i = 0; i < acflen; i++) {
            sum1 += currFrame[i] * predFrame[i];
        }

        for (i = 0; i < acflen; i++) {
            sum2 += currFrame[i] * currFrame[i];
        }

        for (i = 0; i < acflen; i++) {
            sum3 += predFrame[i] * predFrame[i];
        }

        sum2      = LC3_SQRT(sum2 * sum3) + LC3_POW(10, -5);
        norm_corr = sum1 / sum2;

        assert(norm_corr >= -1.00001 && norm_corr <= 1.00001);
        norm_corr = MIN(1, MAX(-1, norm_corr));
        if (norm_corr < 0) {
            norm_corr = 0;
        }

        if (ltpf_enable == 1) {
            /* Decision if ltpf active */
            if ((*mem_on == 0 && (frame_dms == 100 || *mem_norm_corr_past_past > 0.94) && *mem_norm_corr_past > 0.94 &&
                 norm_corr > 0.94) ||
                (*mem_on == 1 && norm_corr > 0.9) ||
                (*mem_on == 1 && LC3_FABS(pitch - *mem_pitch) < 2 && (norm_corr - *mem_norm_corr_past) > -0.1 &&
                 norm_corr > 0.84)) {
                ltpf_active = 1;
            }
        }

        gain = 4;

    } else {
        gain      = 0;
        norm_corr = pitch_ol_norm_corr;
        pitch     = 0;
    }

    if (gain > 0) {
        param[0] = 1;
        param[1] = ltpf_active;
        param[2] = pitch_index;
        *bits    = 11;
    } else {
        zero_int(param, 3);

        *bits = 1;
    }

    if (frame_dms < 100) {
        *mem_norm_corr_past_past = *mem_norm_corr_past;
    }

    *mem_norm_corr_past = norm_corr;
    *mem_on             = ltpf_active;
    *mem_pitch          = pitch;
}
