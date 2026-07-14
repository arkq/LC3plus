/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __LC3_UTIL_H__
#define __LC3_UTIL_H__

#ifdef CR14_A_ADD_LOSSLESS_MODE

static inline void scale_signal24_fx( Word32 x[], /* i:   time input signal */
#    ifdef ENABLE_HR_MODE
                                      Word32 x_scaled[],
#    else
                                      Word16 x_scaled[],
#    endif
                                      Word16* x_exp,
#    ifdef ENABLE_HR_MODE
                                      Word32 mdct_mem[],
#    else
                                      Word16 mdct_mem[],
#    endif
                                      Word16 mdct_mem_len,
                                      Word16 resample_mem_in[],
                                      Word16 resample_mem_in_len,
                                      Word32 resample_mem_in50[],
                                      Word16 resample_mem_out[],
                                      Word16 resample_mem_out_len,
                                      Word32 mdct_mem32[],
                                      Word16 N,
                                      Word32 resamp_mem32[],
                                      Word16 mem_s12k8[],
                                      Word16* resamp_scale )
{
    Word16 i;
    Word16 s;
    Word16 scales[6];

    Dyn_Mem_In( "scale_signal24_fx", sizeof( struct {
                    Word16 i;
                    Word16 s;
                    Word16 scales[6];
                } ) );

    /* Scale input for 24 bit case */

    /* assure 24 bit input */
    FOR( i = 0; i < N; i++ )
    {
        IF( x[i] >= 0 )
        {
            x[i] = L_and( x[i], 0x007fffff );
        }
        ELSE
        {
            x[i] = (Word32) L_or( (UWord32) x[i], 0xff800000 );
        }
    }

    /* Find maximum exponent */
    scales[0] = sub( 15 + 8, getScaleFactor32_0( x, N ) );
    scales[1] = sub( 15 + 8, getScaleFactor32_0( mdct_mem32, mdct_mem_len ) );
    scales[2] = sub( 15 + 8, getScaleFactor32_0( resamp_mem32, resample_mem_in_len ) );
    scales[3] = sub( sub( *resamp_scale, 2 ), getScaleFactor32_0( resample_mem_in50, 2 ) );
    scales[4] = sub( sub( *resamp_scale, 2 ), getScaleFactor16_0( resample_mem_out, resample_mem_out_len ) );
    scales[5] = sub( sub( *resamp_scale, 2 ), getScaleFactor16_0( mem_s12k8, 3 ) );
    *x_exp = 7;
    move16();
    FOR( i = 0; i < 6; i++ )
    {
        *x_exp = s_max( *x_exp, scales[i] );
    }

    /* Shift input buffers */
    s = sub( 15 + 8, *x_exp );
    FOR( i = 0; i < N; i++ )
    {
#    ifdef ENABLE_HR_MODE
        x_scaled[i] = L_shl( x[i], s );
#    else
        x_scaled[i] = round_fx_sat( L_shl( x[i], s ) );
#    endif
    }

    FOR( i = 0; i < mdct_mem_len; i++ )
    {
#    ifdef ENABLE_HR_MODE
        mdct_mem[i] = L_shl( mdct_mem32[i], s );
#    else
        mdct_mem[i] = round_fx_sat( L_shl( mdct_mem32[i], s ) );
#    endif
    }

    FOR( i = 0; i < resample_mem_in_len; i++ )
    {
        resample_mem_in[i] = round_fx_sat( L_shl( resamp_mem32[i], s ) );
    }

    /* Adjust resampler filter and output buffers */
    s = sub( sub( *resamp_scale, 2 ), *x_exp );
    *resamp_scale = add( *x_exp, 2 );

    IF( s )
    {
        FOR( i = 0; i < 2; i++ )
        {
            resample_mem_in50[i] = L_shl( resample_mem_in50[i], s );
        }
        FOR( i = 0; i < resample_mem_out_len; i++ )
        {
            resample_mem_out[i] = shl( resample_mem_out[i], s );
        }

        FOR( i = 0; i < 3; i++ )
        {
            mem_s12k8[i] = shl( mem_s12k8[i], s );
        }
    }
    /* Store part of current frame as mdct memory buffer and resampler input buffer for next frame */
    basop_memcpy( mdct_mem32, &x[N - mdct_mem_len], mdct_mem_len * sizeof( Word32 ) );
    basop_memmove( resamp_mem32, &x[N - resample_mem_in_len], resample_mem_in_len * sizeof( Word32 ) );

    Dyn_Mem_Out();
}

static inline void format_in_pcm_fx( Word16 wavFormat,
                                     void* s_in,
#  ifdef ENABLE_HR_MODE
                                     Word32* s_in_scaled,
#  else
                                     Word16* s_in_scaled,
#  endif
                                     Word16* x_exp,
#    ifdef ENABLE_HR_MODE
                                     Word32* mdct_mem,
#    else
                                     Word16* mdct_mem,
#    endif
                                     Word16 mdct_mem_len,
                                     Word16* resample_mem_in,
                                     Word16 resample_mem_in_len,
                                     Word32* resample_mem_50,
                                     Word16* resample_mem_out,
                                     Word16 resample_mem_out_len,
                                     Word32* mdct_mem32,
                                     Word32* resamp_mem32,
                                     Word16* mem_s12k8,
                                     Word16* resamp_exp,
                                     Word16 len_input,
                                     Word16 lossless
                                     )
{
    UNUSED(mdct_mem);
    UNUSED(mdct_mem_len);
    UNUSED(resample_mem_in);
    UNUSED(resample_mem_in_len);
    UNUSED(resample_mem_50);
    UNUSED(resample_mem_out);
    UNUSED(resample_mem_out_len);
    UNUSED(mdct_mem32);
    UNUSED(resamp_mem32);
    UNUSED(mem_s12k8);
    UNUSED(resamp_exp);
    
    SWITCH( wavFormat )
    {
    case 16:
    {
#  ifdef ENABLE_HR_MODE
        Word16* ip_buf = (Word16*) s_in;
        Word32 i;
        FOR( i = 0; i < len_input; i++ )
        {
            s_in_scaled[i] = L_deposit_h( ip_buf[i] );
        }
        *x_exp = 15;  move16();
        
        if (lossless) { 
            const int headroom = 6;             
            FOR( i = 0; i < len_input; i++ )
            {
                s_in_scaled[i] = L_shr( s_in_scaled[i], headroom);
            }
            *x_exp = 15+headroom;  move16();
        }
#  else
        memcpy( s_in_scaled, s_in, len_input * sizeof( *s_in_scaled ) );
#  endif
    }
    break;
    case 24:
    {
        Word32* ip_buf = (Word32*) s_in;
        FOR( Word32 i = 0; i < len_input; i++ )
        {
            s_in_scaled[i] = L_shl(ip_buf[i],3);
        }
        *x_exp = 31-8-3;
    }
    break;
    default:
        assert( 0 );
        break;
    }
}

static inline void format_in_pcm( Word16 wavFormat,
                                  void* s_in,
#  ifdef ENABLE_HR_MODE
                                  Word32* s_in_scaled,
#  else
                                  Word16* s_in_scaled,
#  endif
                                  EncSetup* h_EncSetup,
                                  LC3PLUS_Enc* encoder )
{

    Word16 hrmode = encoder->hrmode;
      
    if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        hrmode = 0;
    }
  
    format_in_pcm_fx( wavFormat,
                      s_in,
                      s_in_scaled,
                      &h_EncSetup->x_exp,
                      h_EncSetup->stEnc_mdct_mem,
                      encoder->stEnc_mdct_mem_len,
                      h_EncSetup->r12k8_mem_in,
                      encoder->r12k8_mem_in_len,
                      h_EncSetup->r12k8_mem_50,
                      h_EncSetup->r12k8_mem_out,
                      encoder->r12k8_mem_out_len,
                      h_EncSetup->mdct_mem32,
                      h_EncSetup->resamp_mem32,
                      h_EncSetup->olpa_mem_s12k8,
                      &h_EncSetup->resamp_exp,
                      encoder->frame_length,
                      hrmode                  
                       );
}

static inline void format_out_pcm( Word16 wavFormat,
#  ifdef ENABLE_HR_MODE
                                   Word32* x_fx_ip,
#  else
                                   Word16* x_fx,
#  endif
                                   Word16 q_fx_exp,
                                   void* s_out,
                                   int frame_length,
                                   Word16 lossless,
                                   Word16 ll_adap_flag
                                   )
{
    Word32 offset;
    Word16 scale;
    Counter i;

    IF( wavFormat == 16 )
    {
        scale = sub( 15, q_fx_exp );

        IF (lossless && ll_adap_flag)
        {
            scale = abs_s(q_fx_exp);

            FOR( i = 0; i < frame_length; i++ )
            {
                ( (Word16*) s_out )[i] = (Word16) x_fx_ip[i];
            }
            
            goto end;
        }

        FOR( i = 0; i < frame_length; i++ )
        {
#    ifdef ENABLE_HR_MODE
            ( (Word16*) s_out )[i] = round_fx_sat( L_shr_sat( x_fx_ip[i], scale ) );
#    else
            ( (Word16*) s_out )[i] = round_fx_sat( L_shr_sat( L_deposit_h( x_fx[i] ), scale ) );
#    endif
        }
    }
    ELSE IF( wavFormat == 24 )
    {
        scale = sub( sub( 31 + 16, 24 ), q_fx_exp );
        offset = L_shr_sat( 32768, sub( 16, scale ) );

        IF (lossless && ll_adap_flag)
        {
            FOR( i = 0; i < frame_length; i++ )
            {
                ( (Word32*) s_out )[i] = x_fx_ip[i];
            }
            
            goto end;
        }

        FOR( i = 0; i < frame_length; i++ )
        {
#      ifdef ENABLE_HR_MODE
            ( (Word32*) s_out )[i] = L_shr_sat( L_add_sat( x_fx_ip[i], offset ), scale );
#      else
            ( (Word32*) s_out )[i] = L_shr_sat( L_add_sat( L_deposit_h( x_fx[i] ), offset ), scale );
#      endif
        }
    }
end: ;
}

#endif /* CR14_A_ADD_LOSSLESS_MODE */

#endif /* __LC3_UTIL_H__ */
