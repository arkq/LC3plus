/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void process_ltpf_decoder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* y, LC3_INT fs, LC3_FLOAT* mem_old_x, LC3_FLOAT* mem_old_y,
                             LC3_INT* mem_pitch_int, LC3_INT* mem_pitch_fr, LC3_FLOAT* mem_gain, LC3_INT* mem_beta_idx, LC3_INT bfi,
                             LC3_INT* param, LC3_INT* mem_param, LC3_INT conf_beta_idx, LC3_FLOAT conf_beta, LC3_INT concealMethod,
                             LC3_FLOAT damping
                             , LC3_INT *mem_ltpf_active                         
                             , LC3_FLOAT *rel_pitch_change, LC3_INT hrmode, LC3_INT frame_dms
)
{
    LC3_INT i, j, n, N, L_past_x, N4, N34,
        pitch_int, pitch_fr, p1, p2, L_past_y, inter_len, tilt_len = 0,
        tilt_len_r, inter_len_r, old_x_len, old_y_len;

    LC3_FLOAT conf_alpha, gain, a1[12], a2[12], b1[11], b2[11],
          buf_x[4 * MAX_LEN], buf_y[4 * MAX_LEN], buf_z[4 * MAX_LEN], pitch, sum1, sum2;
	LC3_FLOAT *p_x, *p_y, *p_y2, *p_x_init, *p_y_init, *p_a1, *p_b1, *p_a2, *p_b2, fade_fac, current_fade_fac_up, current_fade_fac_down;
    LC3_FLOAT pitch_fl_c_old, pitch_delta;
    const LC3_FLOAT *inter_filter[4], *tilt_filter[4];

#ifdef WMOPS
    push_wmops("process_ltpf_decoder_fl");
#endif
    pitch_fl_c_old = (LC3_FLOAT) *mem_pitch_int + (LC3_FLOAT)*mem_pitch_fr / 4.0;
    conf_alpha = 0.85;

    if (bfi != 1) {
        /* Decode pitch */
        if (param[0] == 1) {
            if (param[2] < (RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4) {
                pitch_int = MIN_PITCH_12K8 + floor(param[2] / 4);
                pitch_fr  = param[2] - ((pitch_int - MIN_PITCH_12K8) * 4);
            } else if (param[2] < ((RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4) + ((RES2_PITCH_12K8 - RES4_PITCH_12K8) * 2)) {
                param[2]  = param[2] - ((RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4);
                pitch_int = RES4_PITCH_12K8 + floor(param[2] / 2);
                pitch_fr  = param[2] - ((pitch_int - RES4_PITCH_12K8) * 2);
                pitch_fr  = pitch_fr * 2;
            } else {
                pitch_int =
                    param[2] + (RES2_PITCH_12K8 - ((RES4_PITCH_12K8 - MIN_PITCH_12K8) * 4) - ((RES2_PITCH_12K8 - RES4_PITCH_12K8) * 2));
                pitch_fr = 0;
            }

            pitch     = ((LC3_FLOAT)pitch_int + (LC3_FLOAT)pitch_fr / 4.0) * (LC3_FLOAT)fs / 12800.0;
            pitch     = round(pitch * 4.0) / 4.0;
            pitch_int = floor(pitch);
            pitch_fr  = (LC3_INT)((pitch - (LC3_FLOAT)pitch_int) * 4.0);
        } else {
            pitch_int = 0;
            pitch_fr  = 0;
        }

        /* Decode gain */
        if (conf_beta_idx < 0) {
            param[1] = 0;
        }

        if (param[1] == 1) {
            gain = conf_beta;
        } else {
            gain = 0;
        }
    }
    else if (concealMethod > 0) {
        if (conf_beta_idx < 0) {
            if (mem_param[1] && *mem_beta_idx >= 0)
            {
                conf_beta_idx = *mem_beta_idx;
            }
        }

        memmove(param, mem_param, sizeof(LC3_INT32) * 3);
        if (concealMethod == 2)
        {
            /* cause the ltpf to "fade_out" and only filter during initial 2.5 ms and then its buffer during 7.5 ms */
            assert(bfi == 1);
            param[1] = 0; /* ltpf_active = 0 */
        }

        pitch_int = *mem_pitch_int;
        pitch_fr  = *mem_pitch_fr;
        gain      = (LC3_FLOAT) *mem_gain * damping;
    }

        if ((conf_beta <= 0) && (*mem_ltpf_active == 0))
        {
            if (fs == 8000 || fs == 16000) {
                tilt_len = 4 - 2;
            }
            else if (fs == 24000) {
                tilt_len = 6 - 2;
            }
            else if (fs == 32000) {
                tilt_len = 8 - 2;
            }
            else if (fs == 44100 || fs == 48000) {
                tilt_len = 12 - 2;
            }

            N = xLen;
            old_x_len = tilt_len;
            inter_len = MAX(fs, 16000) / 8000;
            old_y_len = ceilf((LC3_FLOAT)228.0 * fs / 12800.0) + inter_len;     /* 228.0 needed to make use of ceil */

            move_float(mem_old_y, &mem_old_y[N], (old_y_len - N));
            move_float(&mem_old_y[old_y_len - N], x, N);
            move_float(mem_old_x, &x[N - old_x_len], old_x_len);

            *mem_ltpf_active = 0;
        }
        else
        {
            inter_len_r = 0; tilt_len_r = 0;
        if (fs == 8000 || fs == 16000) {
            inter_filter[0] = conf_inter_filter_16[0];
            inter_filter[1] = conf_inter_filter_16[1];
            inter_filter[2] = conf_inter_filter_16[2];
            inter_filter[3] = conf_inter_filter_16[3];
            inter_len_r     = 4;

            tilt_filter[0] = conf_tilt_filter_16[0];
            tilt_filter[1] = conf_tilt_filter_16[1];
            tilt_filter[2] = conf_tilt_filter_16[2];
            tilt_filter[3] = conf_tilt_filter_16[3];
            tilt_len       = 4 - 2;
            tilt_len_r     = 3;
        } else if (fs == 24000) {
            inter_filter[0] = conf_inter_filter_24[0];
            inter_filter[1] = conf_inter_filter_24[1];
            inter_filter[2] = conf_inter_filter_24[2];
            inter_filter[3] = conf_inter_filter_24[3];
            inter_len_r     = 6;

            tilt_filter[0] = conf_tilt_filter_24[0];
            tilt_filter[1] = conf_tilt_filter_24[1];
            tilt_filter[2] = conf_tilt_filter_24[2];
            tilt_filter[3] = conf_tilt_filter_24[3];
            tilt_len       = 6 - 2;
            tilt_len_r     = 5;
        } else if (fs == 32000) {
            inter_filter[0] = conf_inter_filter_32[0];
            inter_filter[1] = conf_inter_filter_32[1];
            inter_filter[2] = conf_inter_filter_32[2];
            inter_filter[3] = conf_inter_filter_32[3];
            inter_len_r     = 8;

            tilt_filter[0] = conf_tilt_filter_32[0];
            tilt_filter[1] = conf_tilt_filter_32[1];
            tilt_filter[2] = conf_tilt_filter_32[2];
            tilt_filter[3] = conf_tilt_filter_32[3];
            tilt_len       = 8 - 2;
            tilt_len_r     = 7;
        } else if (fs == 44100 || fs == 48000) {
            inter_filter[0] = conf_inter_filter_48[0];
            inter_filter[1] = conf_inter_filter_48[1];
            inter_filter[2] = conf_inter_filter_48[2];
            inter_filter[3] = conf_inter_filter_48[3];
            inter_len_r     = 12;

            tilt_filter[0] = conf_tilt_filter_48[0];
            tilt_filter[1] = conf_tilt_filter_48[1];
            tilt_filter[2] = conf_tilt_filter_48[2];
            tilt_filter[3] = conf_tilt_filter_48[3];
            tilt_len       = 12 - 2;
            tilt_len_r     = 11;
        }

        inter_len = MAX(fs, 16000) / 8000;

        /* Init buffers */
        N         = xLen;
        old_x_len = tilt_len;
        old_y_len = ceilf(228.0 * fs / 12800.0) + inter_len;
        L_past_x = old_x_len;
        move_float(buf_x, mem_old_x, old_x_len);
        move_float(&buf_x[old_x_len], x, xLen);
        L_past_y = old_y_len;
        move_float(buf_y, mem_old_y, old_y_len);
        zero_float(&buf_y[old_y_len], xLen);

        N4  = fs * 0.0025;
        N34 = N - N4;

        /* Init filter parameters */
        if (mem_param[1] == 1) {
            for (i = 0; i < inter_len_r; i++) {
                a1[i] = *mem_gain * inter_filter[*mem_pitch_fr][i];
            }

            for (i = 0; i < tilt_len_r; i++) {
                b1[i] = conf_alpha * (*mem_gain) * tilt_filter[*mem_beta_idx][i];
            }

            p1 = *mem_pitch_int;
        }

        if (param[1] == 1) {
            for (i = 0; i < tilt_len_r; i++) {
                b2[i] = conf_alpha * gain * tilt_filter[conf_beta_idx][i];
            }

            for (i = 0; i < inter_len_r; i++) {
                a2[i] = gain * inter_filter[pitch_fr][i];
            }

            p2 = pitch_int;
        }

        /* First quarter of the current frame: cross-fading */
        fade_fac = 1. / (LC3_FLOAT) N4;
        current_fade_fac_up = 0.f;
        current_fade_fac_down = 1.f;
        (void) p_x; (void) p_y; (void) p_a1; (void) p_b1;
        
        if (mem_param[1] == 0 && param[1] == 0) {
            memmove(&buf_y[L_past_y], &buf_x[L_past_x], sizeof(LC3_FLOAT) * N4);

        } else if (mem_param[1] == 1 && param[1] == 0) {
            for (n = 0; n < N4; n++) {
                sum1 = 0;
                sum2 = 0;
                j    = 0;
                for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                    sum1 += b1[j] * buf_x[i];
                    j++;
                }

                j = 0;
                for (i = L_past_y + n - p1 + inter_len - 1; i >= L_past_y + n - p1 - inter_len; i--) {
                    sum2 += a1[j] * buf_y[i];
                    j++;
                }

                buf_y[L_past_y + n] = buf_x[L_past_x + n] - current_fade_fac_down * sum1 +
                                      current_fade_fac_down * sum2;
                current_fade_fac_down -= fade_fac;
            }

        } else if (mem_param[1] == 0 && param[1] == 1) {
            for (n = 0; n < N4; n++) {
                sum1 = 0;
                sum2 = 0;
                j    = 0;
                for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                    sum1 += b2[j] * buf_x[i];
                    j++;
                }

                j = 0;
                for (i = L_past_y + n - p2 + inter_len - 1; i >= L_past_y + n - p2 - inter_len; i--) {
                    sum2 += a2[j] * buf_y[i];
                    j++;
                }

                buf_y[L_past_y + n] = buf_x[L_past_x + n] - current_fade_fac_up * sum1 + current_fade_fac_up * sum2;
                current_fade_fac_up += fade_fac;
            }
        } else if (*mem_pitch_int == pitch_int && *mem_pitch_fr == pitch_fr) {
            for (n = 0; n < N4; n++) {
                sum1 = 0;
                sum2 = 0;
                j    = 0;
                for (i = L_past_x + n; i >= L_past_x + n - tilt_len; i--) {
                    sum1 += b2[j] * buf_x[i];
                    j++;
                }

                j = 0;
                for (i = L_past_y + n - p2 + inter_len - 1; i >= L_past_y + n - p2 - inter_len; i--) {
                    sum2 += a2[j] * buf_y[i];
                    j++;
                }

                buf_y[L_past_y + n] = buf_x[L_past_x + n] - sum1 + sum2;
            }
        } else {
            p_x_init = &buf_x[L_past_x];
            p_y_init = &buf_y[L_past_y - p1 + inter_len - 1];
            p_y2 = &buf_y[L_past_y];
            for (n = 0; n < N4; n++) {
                sum1 = 0;
                sum2 = 0;
                p_b1 = b1;
                p_x  = p_x_init;
                for (i = tilt_len; i >= 0; i--) {
                    sum1 += *p_b1 * *p_x;
                    p_b1++;
                    p_x--;
                }

                p_y  = p_y_init;
                p_a1 = a1;
                for (i = 2*inter_len - 1; i >= 0; i--) {
                    sum2 += *p_a1 * *p_y;
                    p_a1++;
                    p_y--;
                }

                *p_y2 = *p_x_init - current_fade_fac_down * sum1 +
                                      current_fade_fac_down * sum2;
                current_fade_fac_down -= fade_fac;
                p_x_init++;
                p_y_init++;
                p_y2++;
            }

            move_float(buf_z, buf_y, (old_y_len + xLen));
            p_x_init = &buf_z[L_past_y];  /* buf z in this case */
            p_y_init = &buf_y[L_past_y - p2 + inter_len - 1];
            p_y2 = &buf_y[L_past_y];

            for (n = 0; n < N4; n++) {
                sum1 = 0;
                sum2 = 0;
                j    = 0;
                p_x  = p_x_init;
                p_b2 = b2;
                for (i = tilt_len; i >= 0; i--) {
                    sum1 += *p_b2 * *p_x;
                    p_b2++;
                    p_x--;
                }

                p_y = p_y_init;
                p_a2 = a2;
                for (i = 2*inter_len - 1; i >= 0; i--) {
                    sum2 += *p_a2 * *p_y;
                    p_a2++;
                    p_y--;
                }

                *p_y2 = *p_x_init - current_fade_fac_up * sum1 + current_fade_fac_up * sum2;
                current_fade_fac_up += fade_fac;
                p_x_init++;
                p_y_init++;
                p_y2++;
            }
        }

        /* Second quarter of the current frame */
        if (param[1] == 0) {
            move_float(&buf_y[L_past_y + N4], &buf_x[L_past_x + N4],
                    ((L_past_x + N4 + N34) - (L_past_x + N4)));
        } else {
            p_x_init = &buf_x[L_past_x + N4];
            p_y_init = &buf_y[L_past_y + N4 - p2 + inter_len - 1];
            p_y2 = &buf_y[L_past_y + N4];
            for (n = 0; n < N34; n++) {
                sum1 = 0;
                sum2 = 0;
                p_b2  = b2;
                p_x   = p_x_init;

                for (i = 0; i <= tilt_len; i++) {
                    sum1 += *p_b2 * *p_x;
                    p_b2++;
                    p_x--;
                }
                
                p_a2 = a2;
                p_y  = p_y_init;
                 
                for (i = 2*inter_len - 1; i >= 0; i--) {
                    sum2 += *p_a2 * *p_y;
                    p_a2++;
                    p_y--;
                }
                p_y_init++;
                *p_y2 = *p_x_init - sum1 + sum2;
                p_x_init++;
                p_y2++;
            }
        }
        /* Output */
        move_float(y, &buf_y[L_past_y], N);

        /* Update memory */
        move_float(mem_old_x, &buf_x[N], old_x_len);
        move_float(mem_old_y, &buf_y[N], old_y_len);

        *mem_ltpf_active = (conf_beta > 0);
    }

    /* Update ltpf param memory */
    move_int(mem_param, param, 3);
    *mem_pitch_int = pitch_int;
    *mem_pitch_fr  = pitch_fr;
    *mem_gain      = gain;
    *mem_beta_idx  = conf_beta_idx;                         
    if (bfi == 0 && hrmode == 1 && (frame_dms == 50 || frame_dms == 25)){
        pitch_delta = LC3_FABS(pitch_fl_c_old - (LC3_FLOAT)pitch_int - (LC3_FLOAT)(pitch_fr / 4.0));
        *rel_pitch_change = pitch_delta / MAX(pitch_fl_c_old, 1);
    }

#ifdef WMOPS
    pop_wmops();
#endif
}
