/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void ltpf_synth_filter(Word16 *synth_ltp, Word16 *synth, Word16 length, Word16 pitch_int, Word16 pitch_fr,
                              Word16 gain, Word16 scale_fac_idx, Word16 fs_idx,
                              Word16 fade /* 0=normal, +1=fade-in, -1=fade-out */
                              , Word16 continuation, Word16 splitFading
#ifdef LTPF_ADAPTIVE_GAIN
                              , LC3PLUS_FrameDuration frame_dms
#endif
                              );


/*************************************************************************/

#ifdef CR9_C_ADD_1p25MS
static void get_continuation (Word16 fading_case, LC3PLUS_FrameDuration frame_dms, Word16 *continuation) 
{
    if (frame_dms != LC3PLUS_FRAME_DURATION_1p25MS) {
        *continuation = 0;
    }
    else
    {
        if ( *continuation > 0 )
        {
            *continuation = 0;
        }
        else
        {
            *continuation = fading_case;
        }
    }
}
#endif

#    ifdef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
static Word16 compare_normalized_corrs(Word16 *sig, Word16 len, Word16 pitch_int, Word16 old_pitch_int)
{
    Word32 sum0, sum1, sum2, prod, inv;
    Word32 sum0_prev, sum2_prev, prod_prev;
    Word16 scale0, scale1, scale2, shift, prod_exp, prod_exp_prev;
    Word32 norm_corr, norm_corr_prev; 

    sum0 = L_mult0( sig[0], sig[-pitch_int] );
    sum1 = L_mac0( 1, sig[0], sig[0] );
    sum2 = L_mac0( 1, sig[-pitch_int], sig[-pitch_int] );

    sum0_prev = L_mult0( sig[0], sig[-old_pitch_int] );
    sum2_prev = L_mac0( 1, sig[-old_pitch_int], sig[-old_pitch_int] );

    FOR (int i=1; i < len; i++)
    {
        sum0 = L_mac0( sum0, sig[i], sig[i - pitch_int] );
        sum1 = L_mac0( sum1, sig[i], sig[i] );
        sum2 = L_mac0( sum2, sig[i - pitch_int], sig[i - pitch_int] );

        sum0_prev = L_mac0( sum0_prev, sig[i], sig[i - old_pitch_int] );
        sum2_prev = L_mac0( sum2_prev, sig[i - old_pitch_int], sig[i - old_pitch_int] );
    }

    scale1 = norm_l( sum1 );
    scale2 = norm_l( sum2 );
    sum1 = L_shl_pos( sum1, scale1 );
    sum2 = L_shl_pos( sum2, scale2 );
    prod = Mpy_32_32( sum1, sum2 );
    shift = norm_l( prod );
    prod = L_shl_pos( prod, shift );
    prod_exp = sub( 62, add( add( scale1, scale2 ), shift ) );
    inv = Isqrt( prod, &prod_exp );
    scale0 = norm_l( sum0 );
    sum0 = L_shl_pos( sum0, scale0 );
    prod = Mpy_32_32( sum0, inv );
    prod_exp = add( sub( 31, scale0 ), prod_exp );
    
#ifdef FIX_LTPFDEC_BASOP 
    norm_corr = INT_MAX;         move32();
    test();
    if(prod == 0 || sub(norm_l(prod), prod_exp) >= 0)
    {
        norm_corr = L_shl_sat(prod, prod_exp);
    }
    norm_corr = L_max(0, norm_corr);
#else 
    test();
    IF( prod == 0 || sub( norm_l( prod ), prod_exp ) >= 0 )
    {
        norm_corr = L_max( 0, L_shl_sat( prod, prod_exp ) );       
    }
    ELSE
    {
        norm_corr = 2147483647;         move32();
    }
 
    if ( norm_corr < 0 )
    {
        norm_corr = 0;  move32();
    }
#endif 
 

    scale2 = norm_l( sum2_prev );
    sum2_prev = L_shl_pos( sum2_prev, scale2 );
    prod_prev = Mpy_32_32( sum1, sum2_prev );
    shift = norm_l( prod_prev );
    prod_prev = L_shl_pos( prod_prev, shift );
    prod_exp_prev = sub( 62, add( add( scale1, scale2 ), shift ) );
    inv = Isqrt( prod_prev, &prod_exp_prev );
    scale0 = norm_l( sum0_prev );
    sum0_prev = L_shl_pos( sum0_prev, scale0 );
    prod_prev = Mpy_32_32( sum0_prev, inv );
    prod_exp_prev = add( sub( 31, scale0 ), prod_exp_prev );

#ifdef FIX_LTPFDEC_BASOP 
    norm_corr_prev = INT_MAX;         move32(); 
    test();
    if( prod_prev == 0 || sub( norm_l( prod_prev ), prod_exp_prev ) >= 0 )
    {   
        norm_corr_prev = L_shl_sat( prod_prev, prod_exp_prev );
    }
    norm_corr_prev = L_max( 0L, norm_corr_prev);
#else 
    test();
    IF( prod_prev == 0 || sub(norm_l(prod_prev), prod_exp_prev) >= 0)
    {
        norm_corr_prev = L_max(0, L_shl_sat(prod_prev, prod_exp_prev)); move32();
    }
    ELSE
    {
        norm_corr_prev = 2147483647;        move32();
    }
 
    IF(norm_corr_prev < 0)
    {
        norm_corr_prev = 0;
    }
#endif 

#ifdef FIX_LTPFDEC_BASOP  
    IF(L_sub_sat(L_sub(norm_corr_prev, norm_corr), 1073742L) > 0)
#else 
    IF( L_sub(L_sub(norm_corr_prev, norm_corr), 1073742L )   > 0)
#endif 
    {
        return 1;
    }

    return 0;
}
#  endif

void process_ltpf_decoder_fx(Word16 *x_e, Word16 L_frame, Word16 old_x_len, Word16 fs_idx, Word16 old_y_len,
                             Word16 *old_e, Word16 *x_in, Word16 *old_x, Word16 *y_out, Word16 *old_y, Word16 ltpf,
                             Word16 ltpf_active, Word16 pitch_index, Word16 *old_pitch_int, Word16 *old_pitch_fr,
                             Word16 *old_gain, Word16 *mem_ltpf_active, Word16 scale_fac_idx, Word16 bfi,
                             Word16 concealMethod,
                             Word16 damping, Word16 *old_scale_fac_idx,                      
                             Word32 *rel_pitch_change, Word16 hrmode, LC3PLUS_FrameDuration frame_dms,
                             Word8 *scratchBuffer
#ifdef CR9_C_ADD_1p25MS
                             , Word16* mem_continuation, Word16* mem_pitch_int_prev, 
                             Word16* mem_pitch_fr_prev, Word16* mem_beta_idx_prev, Word16* mem_gain_prev,  Word16 *ltpf_mem_active_prev, Word16* pitch_stability_counter
#endif
                             )
{
    Counter i;
    Word16 gain, s, s0, s1, pitch, pitch_int, pitch_fr, N4, N34, fading_case, split_fading;
    Word16 *x, *y;
    Word16* z;
    Word32 tmp32, pitch_delta;

#  ifdef LTPF_ADAPTIVE_GAIN
    Word16 bkp_ltpf_active = 0;
    Word16 bkp_pitch_int = 0, bkp_pitch_fr = 0;
    Word16 bkp_gain = 0;
#  endif

#  ifdef DYNMEM_COUNT
#    ifndef CR9_N_SHORT_FADE_FOR_UNSTABLE_PITCH
    Dyn_Mem_In( "process_ltpf_decoder_fx", sizeof( struct {
                    Counter i;
                    Word16 gain, s, s0, s1, pitch, pitch_int, pitch_fr, N4, N34;
                    Word16 *x, *y;
                    Word16* z;
                } ) );
#    else 
    Dyn_Mem_In("process_ltpf_decoder_fx", sizeof(struct {
        Counter i;
        Word16 gain, s, s0, s1, pitch, pitch_int, pitch_fr, N4, N34;
        Word16 *x, *y;
        Word16* z;
        Word32 tmp32, pitch_delta;
    }));
#    endif 

#  endif

    z = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = MAX_LEN / 4 + 10 */


#  ifdef CR9_C_ADD_1p25MS
    UNUSED( frame_dms );
    UNUSED( mem_continuation );
    UNUSED( mem_pitch_int_prev );
    UNUSED( mem_pitch_fr_prev );
    UNUSED( mem_beta_idx_prev );
    UNUSED( mem_gain_prev );
#  else
    Word16 mem_cont = 0;
    Word16* mem_continuation = &mem_cont;
#    ifndef FIX_LTPF_MEM_CONTINUATION
    fading_case = 0;
#    endif
    UNUSED( mem_continuation );
#  endif

#  ifdef FIX_LTPF_MEM_CONTINUATION
    fading_case = 0;
#  endif

    split_fading = 0;
#  ifdef CR9_C_ADD_1p25MS
    if ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
    {
        split_fading = 1;
    }
#  endif
    test();
    IF ((sub(bfi, 1) == 0) && (sub(concealMethod, LC3_CON_TEC_NS_STD) == 0))
    {
        ltpf        = 0; move16();
        ltpf_active = 0; move16();
        pitch_int = 0; move16();
        pitch_fr  = 0; move16();
        gain      = 0; move16();
    }

#  ifdef FIX_LTPF_1p25
    IF (pitch_index == -1) {
        ltpf = 0;
        ltpf_active = 0;
    }
#  endif

    /* Filter parameters */
    IF (sub(bfi, 1) != 0)
    {
        IF (ltpf == 0)
        {
            pitch_int = 0; move16();
            pitch_fr  = 0; move16();
        }
        ELSE
        {
            /* Decode pitch */
            IF (sub(pitch_index, 380) < 0)
            {
                pitch_int = shr_pos(add(pitch_index, 64), 2);
                pitch_fr  = add(sub(pitch_index, shl_pos(pitch_int, 2)), 128);
            }
            ELSE IF (sub(pitch_index, 440) < 0)
            {
                pitch_int = shr_pos(sub(pitch_index, 126), 1);
                pitch_fr  = sub(sub(shl_pos(pitch_index, 1), shl_pos(pitch_int, 2)), 252);
            }
            ELSE
            {
                pitch_int = sub(pitch_index, 283);
                pitch_fr  = 0; move16();
            }
            pitch     = add(shl_pos(pitch_int, 2), pitch_fr);
#ifdef ENABLE_HR_MODE
            IF (sub(fs_idx, 5) == 0)
            {
                pitch = round_fx(L_shl_pos(L_mult(shl_pos(pitch, 2), pitch_scale[4]), 1));
            }
            ELSE
#endif
            {
                pitch = mult_r(shl_pos(pitch, 2), pitch_scale[fs_idx]);
            }
            pitch_int = shr_pos(pitch, 2);
            pitch_fr  = sub(pitch, shl_pos(pitch_int, 2));
        }

        /* Decode gain */
        if (scale_fac_idx < 0)
        {
            ltpf_active = 0;
            ASSERT(!(*old_scale_fac_idx < 0 && *mem_ltpf_active == 1));
        }
        IF (ltpf_active == 0)
        {
            gain = 0; move16();
        }
        ELSE
        {
            ASSERT(scale_fac_idx >= 0);
            gain = gain_scale_fac[scale_fac_idx]; move16();
        }
    }
#ifdef CR9_C_ADD_1p25MS 
    ELSE IF ( (sub(concealMethod, LC3_CON_TEC_NS_STD) != 0)   
            && (*mem_continuation == 0) )
#else
    ELSE IF (sub(concealMethod, LC3_CON_TEC_NS_STD) != 0)
#endif
    {
        /* fix to avoid not initialized filtering for concelament 
           might be necessary in case of bit errors or rate switching */
        if (scale_fac_idx < 0) {
            if (*mem_ltpf_active && *old_scale_fac_idx>=0)
            {
                scale_fac_idx = *old_scale_fac_idx;
            }
        }

        ltpf_active = *mem_ltpf_active; move16();

        if ((sub(concealMethod, LC3_CON_TEC_PHASE_ECU) == 0))
        { /* always start fade off to save filtering WMOPS for the remaining 7.5 ms  */
            assert(bfi == 1);
            ltpf_active = 0; move16(); /*always start fade off , still maintain  *mem_ltpf_active */
        }

        pitch_int = *old_pitch_int;
        pitch_fr  = *old_pitch_fr;
        gain      = mult_r(*old_gain, damping);
    }
    
#ifdef CR9_C_ADD_1p25MS
    IF (*mem_continuation > 0) 
    {
#ifdef LTPF_ADAPTIVE_GAIN
        /* Backup LTPF parameters */
        bkp_ltpf_active = ltpf_active;
        bkp_pitch_int   = pitch_int;
        bkp_pitch_fr    = pitch_fr;
        bkp_gain        = gain;
#endif

        fading_case = *mem_continuation;
        ltpf_active = *mem_ltpf_active;
        pitch_int   = *old_pitch_int;
        pitch_fr    = *old_pitch_fr;
        gain        = *old_gain;
#ifdef FIX_LTPF_1p25
        scale_fac_idx = *old_scale_fac_idx;
#endif
        *mem_ltpf_active = *ltpf_mem_active_prev;
        *old_pitch_int   = *mem_pitch_int_prev;
        *old_pitch_fr    = *mem_pitch_fr_prev;
        *old_gain        = *mem_gain_prev;
#ifdef FIX_LTPF_1p25
        *old_scale_fac_idx = *mem_beta_idx_prev;
#endif
    }
#endif

#  ifdef LTPF_ADAPTIVE_GAIN
    IF (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        /* Control variables */
        Word8 ltpf_active_ctrl = ltpf_active;
        Word8 ltpf_active_prev_ctrl = *mem_ltpf_active;
#endif
        Word8 pitch_changed = !( ( pitch_int == *old_pitch_int ) && ( pitch_fr == *old_pitch_fr ) );
        Word8 pitch_was_stable = ( ( *pitch_stability_counter >= LTPF_PITCH_STABILITY_THRESHOLD ) );
        
#    ifdef CR9_C_ADD_1p25MS          
    IF ( *mem_continuation == 0 )
    {
#    endif

#  ifdef LTPF_ADAPTIVE_GAIN_WITH_NORM_CORR
        Word16 tmp_y[LTPF_MEM_Y_LEN] = {0};
        basop_memmove( tmp_y, old_y, ( old_y_len ) * sizeof( Word16 ) );
        basop_memmove( tmp_y + old_y_len, x_in, ( L_frame ) * sizeof( Word16 ) );
       
        Word16 scale1 = sub( getScaleFactor16_0( tmp_y, old_y_len + L_frame ), 3 );
        Scale_sig( tmp_y, old_y_len + L_frame, scale1 );

        IF ( !pitch_was_stable && pitch_changed && pitch_int != 0 && *old_pitch_int != 0 )
        {
            pitch_was_stable = compare_normalized_corrs(tmp_y + old_y_len, L_frame, pitch_int, *old_pitch_int);
        }
#  endif

#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        IF ( ltpf_active_ctrl && !pitch_changed )
#else
        IF ( ltpf_active && !pitch_changed )
#endif
        {
            /* increment gain and increment pitch stability counter */
            gain = pitch_was_stable ? MIN( max_adaptive_gain[scale_fac_idx], MAX( gain, *old_gain ) + adaptive_gain_step ) : MAX( gain, *old_gain );
            (*pitch_stability_counter)++;
        }
#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        ELSE IF ( ltpf_active_ctrl && pitch_changed && !pitch_was_stable )
#else
        ELSE IF ( ltpf_active && pitch_changed && !pitch_was_stable )
#endif
        {
            /* decrement gain and reset pitch stability counter */
            gain = ( *old_gain > gain) ? MAX( gain, *old_gain - max_adaptive_gain_step[scale_fac_idx] ) : gain;
            *pitch_stability_counter = 0;
        }
#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        ELSE IF ( !ltpf_active_ctrl && !pitch_was_stable && ltpf_active_prev_ctrl && pitch_changed )
#else
        ELSE IF ( !ltpf_active && !pitch_was_stable && *mem_ltpf_active && pitch_changed )
#endif
        {
            /* decrement gain, use previous pitch and reset pitch stability counter */
            gain = *old_gain - adaptive_gain_step;

            if (scale_fac_idx>=0 && (gain-gain_scale_fac[scale_fac_idx]) > -20)
            {
                ltpf_active = *mem_ltpf_active;
                pitch_int = *old_pitch_int;
                pitch_fr = *old_pitch_fr;
            }
            else
            {
                gain = 0;
            }
            *pitch_stability_counter = 0;
        }
#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
        ELSE IF ( ( ltpf_active_ctrl && pitch_changed && pitch_was_stable ) || ( !ltpf_active_ctrl && pitch_was_stable ) || ( !ltpf_active_ctrl && !pitch_was_stable && ltpf_active_prev_ctrl && !pitch_changed ) )
#else
        ELSE IF ( ( ltpf_active && pitch_changed && pitch_was_stable ) || ( !ltpf_active && pitch_was_stable ) || ( !ltpf_active && !pitch_was_stable && *mem_ltpf_active && !pitch_changed ) )
#endif
        {
            /* use previous pitch and gain and reset pitch stability counter */
            ltpf_active = *mem_ltpf_active;
            pitch_int = *old_pitch_int;
            pitch_fr = *old_pitch_fr;
            gain        = *old_gain;
            *pitch_stability_counter = 0;
        }
#    ifdef CR9_C_ADD_1p25MS       /* This preprocessor is added here for the future when adaptive gain will be enabled for other frame sizes. */
    }
#ifdef FIX_PLC_CONFORM_ISSUES
    ELSE IF ( *mem_continuation != 0 && bkp_ltpf_active == 1 && bfi == 0 )
#else
    ELSE IF ( *mem_continuation != 0 && bkp_ltpf_active == 1 )
#endif
    {
        /* Code enters this block if LTPF is reenabled when adaptive gain is being applied.     */
        /* In this case, use new LTPF parameters but with a smaller gain than in the prev frame. */
        fading_case = 0;
        *mem_continuation = 0;

        *ltpf_mem_active_prev = *mem_ltpf_active;
        *mem_pitch_int_prev   = *old_pitch_int;
        *mem_pitch_fr_prev    = *old_pitch_fr;
        *mem_gain_prev        = *old_gain;

        *mem_ltpf_active = ltpf_active;
        *old_pitch_int   = pitch_int;
        *old_pitch_fr    = pitch_fr;
        *old_gain        = gain;

        ltpf_active = bkp_ltpf_active;
        pitch_int   = bkp_pitch_int;
        pitch_fr    = bkp_pitch_fr;
        gain        = bkp_gain;

        /* if prev gain > curr gain, then decrease gain slowly.*/
        gain = ( *old_gain > gain) ? MAX( gain, *old_gain - max_adaptive_gain_step[scale_fac_idx] ) : gain;

        *pitch_stability_counter = 0;
    }
#    endif
    }
#  endif /* LTPF_ADAPTIVE_GAIN */

#ifdef FIX_LTPF_1p25
    if ( *mem_ltpf_active && *old_scale_fac_idx < 0 )
    {
        *mem_ltpf_active = 0;
    }

    if ( ltpf_active && scale_fac_idx < 0 )
    {
        ltpf_active = 0;
    }
#endif

    IF( sub( fs_idx, 5 ) < 0 )
    {
        test();
        test();

#    ifdef FIX_LTPF_MEM_CONTINUATION
        IF( ( ltpf_active == 0 && *mem_ltpf_active == 0 ) || fading_case == 1 )
#    else
        IF( ltpf_active == 0 && *mem_ltpf_active == 0 && *mem_continuation == 0 )
#    endif
        {
            /* LTPF inactive */
            basop_memmove( y_out, x_in, L_frame * sizeof( Word16 ) );
#    ifdef FIX_LTPF_MEM_CONTINUATION
#      ifdef CR9_C_ADD_1p25MS
            fading_case = 1;
            get_continuation( fading_case, frame_dms, mem_continuation );
#      endif
#    endif


        /* Update */
        s = sub(*old_e, *x_e);
        IF (s > 0)
        {
            basop_memmove(old_y, &old_y[L_frame], (old_y_len - L_frame) * sizeof(Word16));

            IF (sub(s, 15) > 0)
            {
                basop_memset(&old_y[old_y_len - L_frame], 0, (L_frame) * sizeof(Word16));

                basop_memset(old_x, 0, (old_x_len) * sizeof(Word16));
            }
            ELSE
            {
                FOR (i = 0; i < L_frame; i++)
                {
                    old_y[i + old_y_len - L_frame] = shr(x_in[i], s); move16();
                }
                FOR (i = 0; i < old_x_len; i++)
                {
                    old_x[i] = shr(x_in[i + L_frame - old_x_len], s); move16();
                }
            }
        }
        ELSE
        {
            IF (sub(s, -15) < 0)
            {
                basop_memset(old_y, 0, (old_y_len - L_frame) * sizeof(Word16));
            }
            ELSE
            {
                FOR (i = 0; i < old_y_len - L_frame; i++)
                {
                    old_y[i] = shl(old_y[i + L_frame], s); move16();
                }
            }

            basop_memmove(&old_y[old_y_len - L_frame], x_in, (L_frame) * sizeof(Word16));

            basop_memmove(old_x, &x_in[L_frame - old_x_len], (old_x_len) * sizeof(Word16));

                *old_e = *x_e;
                move16();
            }
            *old_gain = 0;
            move16();
            *mem_ltpf_active = 0;
            move16();
        }
        ELSE
        {
            /* Input/Output buffers */
            x = old_x + old_x_len;
            y = old_y + old_y_len;

#    ifdef ENABLE_HR_MODE
            assert( fs_idx < 5 && "Ltpf not supported for 96kHz!\n" );
#    endif

            N4 = ltpf_overlap_len[fs_idx];
            move16();
#    ifdef CR9_C_ADD_1p25MS
            N4 = N4 >> split_fading;
#    endif
            N34 = sub( L_frame, N4 );
            move16();

        /* Input */
        basop_memmove(x, x_in, (L_frame) * sizeof(Word16));

        /* Scaling */
        s0     = sub(s_min(getScaleFactor16_0(old_x, old_x_len), getScaleFactor16_0(old_y, old_y_len)), 1);
        *old_e = sub(*old_e, s0); move16();
        s1     = sub(getScaleFactor16(x, L_frame), 1);
        *x_e   = sub(*x_e, s1); move16();
        s      = sub(*old_e, *x_e);
        IF (s > 0)
        {
            Scale_sig(x, L_frame, sub(s1, s));
            Scale_sig(old_x, old_x_len, s0);
            Scale_sig(old_y, old_y_len, s0);
            *x_e = *old_e; move16();
        }
        ELSE
        {
            Scale_sig(x, L_frame, s1);
            Scale_sig(old_x, old_x_len, add(s0, s));
            Scale_sig(old_y, old_y_len, add(s0, s));
            *old_e = *x_e; move16();
        }

            test();
#ifndef FIX_LTPF_1p25
            IF (sub(*mem_ltpf_active, 1) == 0 && *old_scale_fac_idx < 0)
            {
                *mem_ltpf_active = 0;
            }
#endif
#ifdef CR9_C_ADD_1p25MS
            IF (*mem_continuation == 0)
            {
#endif
                /* fading case */
                IF (ltpf_active == 0)
                {
                    fading_case = 3;
                }
                ELSE IF (*mem_ltpf_active == 0)
                {
                    fading_case = 2;
                }
                ELSE IF (sub(pitch_int, *old_pitch_int) == 0 && sub(*old_pitch_fr, pitch_fr) == 0)
                {
                    fading_case = 4;
                }
                ELSE
                {
                    fading_case = 5;
                }
#ifdef CR9_C_ADD_1p25MS
            }
#endif

            /* Filtering */
            IF( sub( fading_case, 3 ) == 0 )
            {
                ltpf_synth_filter( y, x, N4, *old_pitch_int, *old_pitch_fr, *old_gain, *old_scale_fac_idx, fs_idx, -1, *mem_continuation, split_fading 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );
            }
            ELSE IF( sub( fading_case, 2 ) == 0 )
            {
                ltpf_synth_filter( y, x, N4, pitch_int, pitch_fr, gain, scale_fac_idx, fs_idx, 1, *mem_continuation, split_fading 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );
            }
            ELSE IF( sub( fading_case, 4 ) == 0 )
            {
                ltpf_synth_filter( y, x, N4, pitch_int, pitch_fr, gain, scale_fac_idx, fs_idx, 0, *mem_continuation, split_fading 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );
            }
            ELSE
            {
                ltpf_synth_filter( y, x, N4, *old_pitch_int, *old_pitch_fr, *old_gain, *old_scale_fac_idx, fs_idx,
                                   -1, *mem_continuation, split_fading 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );

                basop_memmove( z, y - tilt_filter_len[fs_idx], ( N4 + tilt_filter_len[fs_idx] ) * sizeof( Word16 ) );

                ltpf_synth_filter( y, z + tilt_filter_len[fs_idx], N4, pitch_int, pitch_fr, gain, scale_fac_idx,
                                   fs_idx, 1, *mem_continuation, split_fading 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );
            }
#    ifdef CR9_C_ADD_1p25MS
            IF( sub( fading_case, 1 ) > 0 )
            {
                get_continuation( fading_case, frame_dms, mem_continuation );
            }
#    endif

#    ifdef CR9_C_ADD_1p25MS
            IF( ltpf_active > 0 && frame_dms > LC3PLUS_FRAME_DURATION_1p25MS )
#    else
            IF( ltpf_active > 0 )
#    endif
            {
                ltpf_synth_filter( y + N4, x + N4, N34, pitch_int, pitch_fr, gain,
                                   scale_fac_idx, fs_idx, 0, 0, 0 
#    ifdef LTPF_ADAPTIVE_GAIN
                                   , 
                                   frame_dms
#    endif
                                 );
            }
            ELSE
            {
                basop_memmove( &y[N4], &x[N4], N34 * sizeof( Word16 ) );
            }

            /* Output */
            basop_memmove( y_out, y, ( L_frame ) * sizeof( Word16 ) );

            /* Update */
            basop_memmove( old_x, &old_x[L_frame], ( old_x_len ) * sizeof( Word16 ) );
            basop_memmove( old_y, &old_y[L_frame], ( old_y_len ) * sizeof( Word16 ) );

#ifndef FIX_LTPF_DEC_FLFX_MISMATCH
#    ifdef CR9_C_ADD_1p25MS
            if ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
            {
                *mem_pitch_int_prev = *old_pitch_int;
                *mem_pitch_fr_prev = *old_pitch_fr;
                *mem_gain_prev = *old_gain;
                *ltpf_mem_active_prev = *mem_ltpf_active;
            }
#    endif

            *old_gain = gain;
            move16();
            *mem_ltpf_active = ltpf_active;
            move16();
#endif
        }
    }




    IF( bfi == 0 && sub( hrmode, 1 ) == 0 && ( sub( frame_dms, LC3PLUS_FRAME_DURATION_5MS ) == 0 || sub( frame_dms, LC3PLUS_FRAME_DURATION_2p5MS ) == 0 ) )
    {
        pitch_delta = abs_s( add( sub( *old_pitch_int, pitch_int ), shr_pos( sub( *old_pitch_fr, pitch_fr ), 2 ) ) );
        tmp32 = BASOP_Util_Divide3216_Scale( pitch_delta, MAX( add( *old_pitch_int, shr_pos( *old_pitch_fr, 2 ) ), 1 ), &s0 );
        IF( s0 + 16 < 0 )
        {
            *rel_pitch_change = L_shr_pos( tmp32, -( s0 + 16 ) );
        }
        ELSE
        {
            *rel_pitch_change = L_shl_pos( tmp32, s0 + 16 );
        }
    }

#ifdef FIX_LTPF_DEC_FLFX_MISMATCH
#    ifdef CR9_C_ADD_1p25MS
        if ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
        {
            *mem_pitch_int_prev = *old_pitch_int;
            *mem_pitch_fr_prev = *old_pitch_fr;
            *mem_gain_prev = *old_gain;
            *ltpf_mem_active_prev = *mem_ltpf_active;
            *mem_beta_idx_prev = *old_scale_fac_idx;
        }
#    endif

        *old_gain = gain;
        move16();
        *mem_ltpf_active = ltpf_active;
        move16();
#endif

    *old_pitch_int = pitch_int;
    move16();
    *old_pitch_fr = pitch_fr;
    move16();
    *old_scale_fac_idx = scale_fac_idx;
    move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


static void ltpf_synth_filter(Word16 *synth_ltp, Word16 *synth, Word16 length, Word16 pitch_int, Word16 pitch_fr,
                              Word16 gain, Word16 scale_fac_idx, Word16 fs_idx,
                              Word16 fade /* 0=normal, +1=fade-in, -1=fade-out */ 
                              ,Word16 continuation, Word16 splitFading
#ifdef LTPF_ADAPTIVE_GAIN
                              , LC3PLUS_FrameDuration frame_dms
#endif
                              )
{
    Word16 *x0;
    Word16 *y0;
    Word32  s;
    Word16  alpha, step;
    Word16  i, k;
    Counter j, l;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("ltpf_synth_filter", sizeof(struct {
                   Word16 *x0;
                   Word16 *y0;
                   Word32  s;
                   Word16  alpha, step;
                   Word16  i, k;
                   Counter j, l;
               }));
#endif

    ASSERT(scale_fac_idx >= 0);

    step  = 0; /* initialize just to avoid compiler warning */
    alpha = 0; /* initialize just to avoid compiler warning */
    x0    = &synth_ltp[-pitch_int + inter_filter_shift[fs_idx]];
    y0    = synth;

    alpha = 0; move16();
    
    IF (fade != 0)
    {
        if (fade < 0)
        {
            alpha = 0x7FFF; move16();
        }
        if (continuation > 0) {
            alpha = shr(0x7FFF,1);
        }

/* step = 1.f/(float)(length); */
        if (sub(length, 20) == 0)
        {
            step = 1638 /*1.f/20.f Q15*/; move16();
        }
        if (sub(length, 30) == 0)
        {
            step = 1092 /*1.f/30.f Q15*/; move16();
        }
        if (sub(length, 40) == 0)
        {
            step = 819 /*1.f/40.f Q15*/; move16();
        }
        if (sub(length, 60) == 0)
        {
            step = 546 /*1.f/60.f Q15*/; move16();
        }
        if (sub(length, 80) == 0)
        {
            step = 409 /*1.f/80.f Q15*/; move16();
        }
        if (sub(length, 120) == 0)
        {
            step = 273 /*1.f/120.f Q15*/; move16();
        }

        if (fade < 0)
            step = negate(step);

        if (splitFading != 0)
            step = shr(step,1);
    }

    FOR (j = 0; j < length; j++)
    {
        const Word16* tilt_filter_ptr = tilt_filter[fs_idx][scale_fac_idx];
#      ifdef LTPF_ADAPTIVE_GAIN
        IF ( frame_dms == LC3PLUS_FRAME_DURATION_1p25MS )
        {
            tilt_filter_ptr = tilt_filter_1p25ms[fs_idx][scale_fac_idx];
        }
#      endif

        s = L_mult( x0[0], inter_filter[fs_idx][pitch_fr][0] );
        FOR( l = 1; l < inter_filter_len[fs_idx]; l++ )
        {
            s = L_mac(s, x0[-l], inter_filter[fs_idx][pitch_fr][l]);
        }
        FOR (l = 0; l < tilt_filter_len[fs_idx]; l++)
        {
            s = L_msu( s, y0[-l], tilt_filter_ptr[l] );
        }

        i = msu_r( s, y0[-l], tilt_filter_ptr[l] );

        k = mult_r(gain, i);

        if (fade != 0)
            k = mult_r(k, alpha);

        synth_ltp[j] = add(synth[j], k); move16();

        if (fade != 0)
            alpha = add(alpha, step);

        x0++;
        y0++;
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

