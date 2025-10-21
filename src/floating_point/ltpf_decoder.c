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
static LC3_INT16 get_continuation (LC3_INT32 fading_case, LC3PLUS_FrameDuration frame_dms, LC3_INT32 pos, LC3_INT32 total) 
{
    if ( frame_dms != LC3PLUS_FRAME_DURATION_1p25MS )
    {
        return 0;
    }
    else
    {
        if ( pos == total )
        {
            return 0;
        }
        else
        {
            return fading_case;
        }
    }
}
#endif

#    ifdef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
#      ifdef LTPF_PRINT_PARAMS
static bool compare_normalized_corrs(LC3_FLOAT const *sig, LC3_INT32 len, LC3_INT32 pitch_int, LC3_INT32 mem_pitch_int, LC3_FLOAT *corr, LC3_FLOAT *corr_prev)
#      else
static bool compare_normalized_corrs(LC3_FLOAT const *sig, LC3_INT32 len, LC3_INT32 pitch_int, LC3_INT32 mem_pitch_int)
#      endif
{
    LC3_FLOAT norm_0, norm_t, norm_t_prev, xcorr, xcorr_prev;
    norm_0 = norm_t = norm_t_prev = xcorr = xcorr_prev = 0;

    for ( int i=0; i < len; i++ )
    {
        norm_0 += sig[i] * sig[i];
        norm_t += sig[i - pitch_int] * sig[i - pitch_int];
        xcorr  += sig[i] * sig[i - pitch_int];

        norm_t_prev += sig[i - mem_pitch_int] * sig[i - mem_pitch_int];
        xcorr_prev  += sig[i] * sig[i - mem_pitch_int];
    }

    xcorr = MIN( 1.f, MAX( 0.f, xcorr / ( LC3_SQRT( norm_0 * norm_t ) + 1.00e-05f ) ) );
    xcorr_prev = MIN( 1.f, MAX( 0.f, xcorr_prev / ( LC3_SQRT( norm_0 * norm_t_prev ) + 1.00e-05f ) ) );

#      ifdef LTPF_PRINT_PARAMS
    *corr = xcorr;
    *corr_prev = xcorr_prev;
#      endif

    if ( xcorr_prev - xcorr > 5e-4 )
    {
        return true;
    }

    return false;
}
#    endif

void process_ltpf_decoder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* y, LC3_INT fs, LC3_FLOAT* mem_old_x, LC3_FLOAT* mem_old_y,
                             LC3_INT* mem_pitch_int, LC3_INT* mem_pitch_fr, LC3_FLOAT* mem_gain, LC3_INT* mem_beta_idx, LC3_INT bfi,
                             LC3_INT* param, LC3_INT* mem_param, LC3_INT conf_beta_idx, LC3_FLOAT *conf_beta, LC3_INT concealMethod, LC3_FLOAT damping
                             , LC3_INT *mem_ltpf_active
                             , LC3_FLOAT *rel_pitch_change, LC3_INT hrmode, LC3PLUS_FrameDuration frame_dms
#ifdef CR9_C_ADD_1p25MS 
                             , LC3_INT16* mem_continuation, LC3_INT32* mem_param_prev, LC3_INT16* mem_pitch_int_prev, 
                             LC3_INT16* mem_pitch_fr_prev, LC3_INT32* mem_beta_idx_prev, LC3_FLOAT* mem_gain_prev,
                             LC3_INT16* pitch_stability_counter, 
                              LC3_FLOAT* gain_step,
                              LC3_FLOAT conf_beta_max
#endif 
 )
{
    LC3_INT i, j, n, N, L_past_x, N4, N34,
        pitch_int, pitch_fr, p1, p2, L_past_y, inter_len, tilt_len = 0,
        tilt_len_r, inter_len_r, old_x_len, old_y_len, fading_case, N4_D, N4_S;

    LC3_FLOAT conf_alpha, gain, a1[12], a2[12], b1[11], b2[11],
          buf_x[4 * MAX_LEN], buf_y[4 * MAX_LEN], buf_z[4 * MAX_LEN], pitch, sum1, sum2;
	LC3_FLOAT *p_x, *p_y, *p_y2, *p_x_init, *p_y_init, *p_a1, *p_b1, *p_a2, *p_b2, fade_fac, current_fade_fac_up, current_fade_fac_down;
    LC3_FLOAT pitch_delta;
    const LC3_FLOAT *inter_filter[4], *tilt_filter[4];
#ifdef LTPF_ADAPTIVE_GAIN
    bool ltpf_active = false;
    bool ltpf_active_prev = false;
    bool pitch_changed = false;
    bool pitch_was_stable = false;
#endif

#  ifdef LTPF_ADAPTIVE_GAIN
    LC3_INT32 original_param[3];
    LC3_INT32 original_pitch_int, original_pitch_fr;
    LC3_FLOAT original_gain;

#    ifdef LTPF_PRINT_PARAMS
    static int local_frame = 0;
#      ifdef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
    LC3_FLOAT xcorr, xcorr_prev;
#      endif
#    endif
#  endif

#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
    UNUSED(mem_ltpf_active);
    fading_case = 0;
#endif
    tilt_len = 0;
    conf_alpha = 0.85;
    p1 = 0;
    p2 = 0;
    
#ifdef WMOPS
    push_wmops("process_ltpf_decoder_fl");
#endif
 
#ifdef FIX_LTPF_1p25
    if (param[2] == -1) {
        param[1] = 0;
        param[0] = 0;
    }
#endif
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
            gain = *conf_beta;
        } else {
            gain = 0;
        }
    }
    else if ((concealMethod > 0) 
#ifdef CR9_C_ADD_1p25MS    
    && (*mem_continuation == 0)
#endif  
    ) {
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
    
#ifdef CR9_C_ADD_1p25MS
    if (*mem_continuation) 
    {
#ifdef LTPF_ADAPTIVE_GAIN
        /* Save original LTPF parameters */
        move_int(original_param, param, 3);
        original_pitch_int   = pitch_int;
        original_pitch_fr    = pitch_fr;
        original_gain        = gain;
#endif

        fading_case = *mem_continuation;
        move_int(param, mem_param, 3);
        pitch_int   = *mem_pitch_int;
        pitch_fr    = *mem_pitch_fr;
        gain        = *mem_gain;
#ifdef FIX_LTPF_1p25
        conf_beta_idx = *mem_beta_idx;
#endif

        move_int(mem_param, mem_param_prev, 3);
        *mem_pitch_int   = *mem_pitch_int_prev;
        *mem_pitch_fr    = *mem_pitch_fr_prev;
        *mem_gain        = *mem_gain_prev;
#ifdef FIX_LTPF_1p25
        *mem_beta_idx = *mem_beta_idx_prev;
#endif

    }
#endif

    if ( fs <= 48000 )
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

        inter_len = MAX(fs, 16000) / 8000;

        /* Init buffers */
        N         = xLen;
        old_x_len = tilt_len;
        old_y_len = ceil(228.0 * fs / 12800.0) + inter_len;
        L_past_x = old_x_len;

        move_float(buf_x, mem_old_x, old_x_len);
        move_float(&buf_x[old_x_len], x, xLen);
        L_past_y = old_y_len;
        move_float(buf_y, mem_old_y, old_y_len);
        move_float(&buf_y[old_y_len], x, xLen);
    }
#ifdef LTPF_ADAPTIVE_GAIN
    if ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
    {
        conf_alpha = 0.98;

        /* Control variables */
        ltpf_active = param[1];
        ltpf_active_prev = mem_param[1];
        pitch_changed = !( ( pitch_int == *mem_pitch_int ) && ( pitch_fr == *mem_pitch_fr ) );
        pitch_was_stable = ( ( *pitch_stability_counter >= LTPF_PITCH_STABILITY_THRESHOLD ) );
        
#  ifdef CR9_C_ADD_1p25MS
    if ( *mem_continuation == 0 )
    {
#  endif
#  ifdef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
        if ( !pitch_was_stable && pitch_changed && pitch_int != 0 && *mem_pitch_int != 0 )
        {
#    ifdef LTPF_PRINT_PARAMS
            pitch_was_stable = compare_normalized_corrs(buf_y + L_past_y, xLen, pitch_int, *mem_pitch_int, &xcorr, &xcorr_prev);
#    else
            pitch_was_stable = compare_normalized_corrs(buf_y + L_past_y, xLen, pitch_int, *mem_pitch_int);
#    endif
        }
#  endif
        
        if ( ltpf_active && !pitch_changed )
        {
            /* increment gain and increment pitch stability counter */
            gain = pitch_was_stable ? MIN( conf_beta_max, MAX( gain, *mem_gain ) + *gain_step ) : MAX( gain, *mem_gain );
            ( *pitch_stability_counter )++;
        }
        else if ( ltpf_active && pitch_changed && !pitch_was_stable )
        {
            /* decrement gain and reset pitch stability counter */
            gain = ( *mem_gain > gain ) ? MAX( gain, *mem_gain - ( conf_beta_max / LTPF_ADAPTIVE_GAIN_RATE ) ) : gain;
            *pitch_stability_counter = 0;
        }
        else if ( !ltpf_active && !pitch_was_stable && ltpf_active_prev && pitch_changed )
        {
            /* decrement gain, use previous pitch and reset pitch stability counter */
            gain = *mem_gain - *gain_step;

#ifdef FIX_LTPF_1p25
            if (*conf_beta > 0 && (gain - *conf_beta) > -(20.f/(1<<15)))
#else
            if ( (gain - *conf_beta) > -(20.f/(1<<15)))
#endif
            {
                move_int( param, mem_param, 3 );
                pitch_int = *mem_pitch_int;
                pitch_fr = *mem_pitch_fr;
            }
            else
            {
                gain = 0.f;
            }
            *pitch_stability_counter = 0;
        }
        else if ( ( ltpf_active && pitch_changed && pitch_was_stable ) || ( !ltpf_active && pitch_was_stable ) || ( !ltpf_active && !pitch_was_stable && ltpf_active_prev && !pitch_changed ) )
        {
            /* use previous pitch and gain and reset pitch stability counter */
            move_int( param, mem_param, 3 );
            pitch_int = *mem_pitch_int;
            pitch_fr = *mem_pitch_fr;
            gain        = *mem_gain;
            *pitch_stability_counter = 0;
        }
#  ifdef CR9_C_ADD_1p25MS       /* This is added here for the future when adaptive gain will be enabled for other frame sizes. */
    }
    else if ( *mem_continuation != 0 && original_param[1] == 1 
#ifdef FIX_PLC_CONFORM_ISSUES
    && bfi == 0 
#endif
    )
    {
        /* Code enters this block if LTPF is reenabled when adaptive gain is being applied.    */
        /* In this case, use new LTPF parameters but with a smaller gain than in the prev frame.*/
        fading_case = 0;
        *mem_continuation = 0;

        move_int( mem_param_prev, mem_param, 3 );
        *mem_pitch_int_prev   = *mem_pitch_int;
        *mem_pitch_fr_prev    = *mem_pitch_fr;
        *mem_gain_prev        = *mem_gain;

        move_int( mem_param, param, 3 );
        *mem_pitch_int   = pitch_int;
        *mem_pitch_fr    = pitch_fr;
        *mem_gain        = gain;

        move_int( param, original_param, 3 );
        pitch_int   = original_pitch_int;
        pitch_fr    = original_pitch_fr;
        gain        = original_gain;

        /* if prev gain > curr gain, then decrease gain slowly. */
        gain = ( *mem_gain > gain) ? MAX( gain, *mem_gain - ( conf_beta_max / LTPF_ADAPTIVE_GAIN_RATE ) ) : gain;

        *pitch_stability_counter = 0;
    }
#  endif
    }
#endif /* LTPF_ADAPTIVE_GAIN */
    if ( mem_param[1] && *mem_beta_idx < 0 )
    {
        mem_param[1] = 0;
    }

    if ( param[1] && conf_beta_idx < 0 )
    {
        param[1] = 0;
    }

    if ( fs <= 48000 )
    {
#ifdef FIX_LTPF_MEM_CONTINUATION
#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        if ( (( param[1] == 0 ) && ( mem_param[1] == 0 )) || fading_case == 1 )
#else
        if ( (( *conf_beta <= 0 ) && ( *mem_ltpf_active == 0 )) || fading_case == 1 )
#endif
#else
        if ( ( *conf_beta <= 0 ) && ( *mem_ltpf_active == 0 ) )
#endif
        {
#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
            if ( fs == 8000 || fs == 16000 )
            {
                tilt_len = 4 - 2;
            }
            else if ( fs == 24000 )
            {
                tilt_len = 6 - 2;
            }
            else if ( fs == 32000 )
            {
                tilt_len = 8 - 2;
            }
            else if ( fs == 44100 || fs == 48000 )
            {
                tilt_len = 12 - 2;
            }
            N = xLen;
            old_x_len = tilt_len;
            inter_len = MAX(fs, 16000) / 8000;
            old_y_len = ceilf((LC3_FLOAT)228.0 * fs / 12800.0) + inter_len;     /* 228.0 needed to make use of ceil */
#    endif

#    ifdef FIX_LTPF_MEM_CONTINUATION
#      ifdef CR9_C_ADD_1p25MS
            fading_case = 1;
            N4 = fs * 0.0025;
            N4_S = 0;
            N4_D = N4;
            if ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
            {
                N4_D = 2 * N4;
                if ( *mem_continuation )
                {
                    N4_S = N4;
                }
            }
            *mem_continuation = get_continuation( fading_case, frame_dms, ( N4 + N4_S ), N4_D );
#      endif
#    endif

            move_float(mem_old_y, &mem_old_y[N], (old_y_len - N));
            move_float(&mem_old_y[old_y_len - N], x, N);
            move_float(mem_old_x, &x[N - old_x_len], old_x_len);

#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
            mem_param[1] = 0;
#else
            *mem_ltpf_active = 0;
#endif
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
#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
                tilt_len = 4 - 2;
#    endif
                tilt_len_r = 3;
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
#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
                tilt_len = 6 - 2;
#    endif
                tilt_len_r = 5;
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
#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
                tilt_len = 8 - 2;
#    endif
                tilt_len_r = 7;
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
#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
                tilt_len = 12 - 2;
#    endif
                tilt_len_r = 11;
        }

#    ifndef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
            inter_len = MAX( fs, 16000 ) / 8000;

        /* Init buffers */
        N         = xLen;
        old_x_len = tilt_len;
        old_y_len = ceilf(228.0 * fs / 12800.0) + inter_len;
        L_past_x = old_x_len;
        move_float(buf_x, mem_old_x, old_x_len);
        move_float(&buf_x[old_x_len], x, xLen);
        L_past_y = old_y_len;
        move_float(buf_y, mem_old_y, old_y_len);
            zero_float( &buf_y[old_y_len], xLen );
#    endif

        N4  = fs * 0.0025;
#ifdef CR9_C_ADD_1p25MS
        N4 = MIN(N4, xLen);
#endif
        N34 = N - N4;

        /* Init filter parameters */
        if (mem_param[1] == 1) {
            for (i = 0; i < inter_len_r; i++) {
                a1[i] = *mem_gain * inter_filter[*mem_pitch_fr][i];
            }

            assert( *mem_beta_idx >= 0 );
            for (i = 0; i < tilt_len_r; i++) {
                b1[i] = conf_alpha * (*mem_gain) * tilt_filter[*mem_beta_idx][i];
            }

            p1 = *mem_pitch_int;
        }

        if (param[1] == 1) {
			assert( conf_beta_idx >= 0 );
            for (i = 0; i < tilt_len_r; i++) {
                b2[i] = conf_alpha * gain * tilt_filter[conf_beta_idx][i];
            }

            for (i = 0; i < inter_len_r; i++) {
                a2[i] = gain * inter_filter[pitch_fr][i];
            }

            p2 = pitch_int;
        }
        
#ifdef CR9_C_ADD_1p25MS
            if (*mem_continuation == 0) {
#endif
                /* check fading case */
                if (mem_param[1] == 0 && param[1] == 0) {
                    fading_case = 1;
                } else if (mem_param[1] == 1 && param[1] == 0) {
                    fading_case = 3;
                } else if (mem_param[1] == 0 && param[1] == 1) {
                    fading_case = 2;
                } else if (*mem_pitch_int == pitch_int && *mem_pitch_fr == pitch_fr) {
                    fading_case = 4;
                } else {
                    fading_case = 5;
                }
#ifdef CR9_C_ADD_1p25MS
            }
#endif

            N4_D = N4;
            UNUSED(N4_D); UNUSED(N4_S);
            N4_S = 0;
#ifdef CR9_C_ADD_1p25MS
            if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
            {
                N4_D = 2*N4;
                if (*mem_continuation)
                {
                    N4_S = N4;
                }
            }
#endif

        /* First quarter of the current frame: cross-fading */
        fade_fac = 1. / (LC3_FLOAT) N4_D;
        current_fade_fac_up = N4_S*fade_fac;
        current_fade_fac_down = 1.f - current_fade_fac_up;
        (void) p_x; (void) p_y; (void) p_a1; (void) p_b1;
        
        if (fading_case == 1) {
            
            memmove(&buf_y[L_past_y], &buf_x[L_past_x], sizeof(LC3_FLOAT) * N4);
#    if defined( CR9_C_ADD_1p25MS ) && defined( FIX_LTPF_MEM_CONTINUATION )
            *mem_continuation = get_continuation(fading_case, frame_dms, (N4+N4_S), N4_D);
#endif
        } else if (fading_case == 3) {
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
#ifdef CR9_C_ADD_1p25MS
            *mem_continuation = get_continuation(fading_case, frame_dms, (n+N4_S), N4_D);
#endif
        } else if (fading_case == 2) {
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
#ifdef CR9_C_ADD_1p25MS                
            *mem_continuation = get_continuation(fading_case, frame_dms, (n+N4_S), N4_D);
#endif
        } else if (fading_case == 4) {
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
#ifdef CR9_C_ADD_1p25MS                
            *mem_continuation = get_continuation(fading_case, frame_dms, (n+N4_S), N4_D);
#endif
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
#ifdef CR9_C_ADD_1p25MS
            *mem_continuation = get_continuation(fading_case, frame_dms, (n+N4_S), N4_D);
#endif
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

#ifndef FIX_LTPF_DEC_FLFX_MISMATCH
        *mem_ltpf_active = ( *conf_beta > 0 );
#endif
        }
    }


    if ( bfi == 0 && hrmode == 1 && ( frame_dms == LC3PLUS_FRAME_DURATION_5MS || frame_dms == LC3PLUS_FRAME_DURATION_2p5MS ) )
    {
        pitch_delta = LC3_FABS( (LC3_FLOAT) *mem_pitch_int + (LC3_FLOAT) ( *mem_pitch_fr / 4.0 ) - (LC3_FLOAT) pitch_int - (LC3_FLOAT) ( pitch_fr / 4.0 ) );
        *rel_pitch_change = pitch_delta / MAX( (LC3_FLOAT) *mem_pitch_int + (LC3_FLOAT) ( *mem_pitch_fr / 4.0 ), 1 );
    }

#ifdef CR9_C_ADD_1p25MS
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) 
    {
        move_int(mem_param_prev, mem_param, 3);
        *mem_pitch_int_prev = *mem_pitch_int;
        *mem_pitch_fr_prev  = *mem_pitch_fr;
        *mem_gain_prev     = *mem_gain;
        *mem_beta_idx_prev = *mem_beta_idx;
    }
#endif

    /* Update ltpf param memory */
    move_int(mem_param, param, 3);
    *mem_pitch_int = pitch_int;
    *mem_pitch_fr  = pitch_fr;
    *mem_gain      = gain;
    *mem_beta_idx  = conf_beta_idx;                         


#ifdef WMOPS
    pop_wmops();
#endif
}
