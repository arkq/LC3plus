/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"
#include "util.h"

static void Enc_LC3PLUS_Channel(LC3PLUS_Enc *encoder, int channel, int bits_per_sample, Word32 *s_in, UWord8 *bytes, lc3_scratch_t scratch, int bfi_ext)
{
#ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(Word16 d_fx_exp;
                      Word16 gain_e, gain, quantizedGain, quantizedGainMin;
                      Word16 ener_fx_exp;
                      Word16 pitch, normcorr;
                      Word16 ltpf_bits;
                      Word16 tns_numfilters;
                      Word16 lsbMode, lastnz, BW_cutoff_idx;
                      Word16 gainChange, fac_ns_idx;
                      Word16 bp_side, mask_side;
                      Word16 s_12k8_len;
                      Word16 b_left;
                      Word32 * L_scf_idx;
                      Word32 * d_fx, *ener_fx;
                      Word16 * s_12k8, *int_scf_fx_exp, tns_order[TNS_NUMFILTERS_MAX], *indexes;
                      Word32 * q_d_fx24;
                      Word16 * scf;
                      Word16 * codingdata;
                      Word32 * scf_q, *int_scf_fx;
                      Word32 * s_in_scaled;
                      Word16 * s_in_scaled_lp;
                      UWord8 * resBits;
                      Word16 ltpf_idx[3];
                      EncSetup * h_EncSetup;
                      Word16 hrmode;
                     );
#else
    Dyn_Mem_Deluxe_In(Word16 d_fx_exp; Word16 gain_e, gain, quantizedGain, quantizedGainMin; Word16 ener_fx_exp;
                      Word16 pitch, normcorr; Word16 ltpf_bits; Word16 tns_numfilters;
                      Word16 lsbMode, lastnz, BW_cutoff_idx; Word16 gainChange, fac_ns_idx;
                      Word16 bp_side, mask_side; Word16 s_12k8_len; Word16 b_left;

                      Word32 * L_scf_idx; Word32 * d_fx, *ener_fx;
                      Word16 * s_12k8, *int_scf_fx_exp, *q_d_fx16, *int_scf_fx, tns_order[TNS_NUMFILTERS_MAX], *indexes;
                      Word16 * scf, *scf_q; Word16 * codingdata; Word16 * s_in_scaled; UWord8 * resBits;
                      Word16 ltpf_idx[3]; EncSetup * h_EncSetup;);
#endif /* ENABLE_HR_MODE */
  
    d_fx_exp = 0;
    gain_e = 0; gain = 0; quantizedGain = 0; quantizedGainMin = 0;
    ener_fx_exp = 0;
    pitch = 0; normcorr = 0;
    ltpf_bits = 0;
    tns_numfilters = 0;
    lsbMode = 0; lastnz = 0; BW_cutoff_idx = 0;
    gainChange = 0; fac_ns_idx = 0;
    bp_side = 0; mask_side = 0;
    s_12k8_len = 0;
    b_left = 0;
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
  Word32 nBits;
  Word32 numResBits;
#else
  Word16 nBits;
  Word16 numResBits;
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
    Word16 envelope_bits_fx;
    Word16 pitch_rx_fx;
    Word16 ltpf_rx_fx;
    Word16 lrsns_st1C_in_use ;
    Word16 ltptx_lowest_bit_lim;
#endif
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    fac_ns_idx = 0;
    gain = 0;
    ltpf_bits = 0;
    normcorr = 0;

    Word16 ll_side[2];
    Word16 off_idx = 0;
    Word32 gg_idx_off_ll_adap;
    UNUSED(gg_idx_off_ll_adap);
    Word16 tns_lsb_num_remove = 0;
    UWord8* deterministic_curve;

    UWord8* tns_lsb_remove;
    Word32* d_fx_orig;
    Word16 ll_adap_flag = 0;
    Word16 ll_adap_flag_2nd = 0;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 scaleSignal = 0;
#endif
    Word32 active_bits =0;
    Word32 entropy_bits =0;
    Word16 fallback_bit_planes = 0;
    Word16 fallback_num_bytes = 0;
    
    Word32* d_tda_fx;
    Word32* d_fx_tmp;

    Word16 d_fx_exp_tmp = 0;
    Word32 max_resBits_len = 0;
    Word16 bitsPerSample = 0;

    Word32 ll_ari_bits_lb = 0;
    Word32 active_bits_lb = 0;
    UWord8 ll_deltaCodesBits[HIGH_BANDS_NUMBER];
    basop_memset(ll_deltaCodesBits, 0, HIGH_BANDS_NUMBER);

#ifdef LL_INCL_HPVC
    /* HPVC analysis needs truncated_data and L_bits_tcxQ9 in scope across
     * Quant. 2 and the analysis block. Allocate at function scope.
     * hpvc_trees must live through processAriEncoder_fx (which is called
     * later in this function) since hpvcEncCfg.HpvcTreeEnumCfgPtr aliases it. */
    Word32 hpvc_truncated_data[MAX_LEN];
    Word32 hpvc_L_bits_tcxQ9[MAX_LEN / 8];
    HpvcTreeEnumCfg hpvc_trees[4 * HPVC_NOMTREE_COUNT_FB];
    Word16 Tx_dec[MAX_LEN / LL_HPVC_N_SIGNAL];
    Word16 Tx_splitRule[MAX_LEN / LL_HPVC_N_SIGNAL];
#endif

    IF (encoder->lossless)
    {
        ll_adap_flag = 1;
        ll_adap_flag_2nd = 1;
    }
    ELSE
    {
        ll_adap_flag = 0;
        ll_adap_flag_2nd = 0;
    }
    
    UNUSED(ll_adap_flag_2nd);

    if (encoder->ll_cbr){
        ll_adap_flag = 0;
    }
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
    envelope_bits_fx = -1;  /*later move up to Dynmem() struct area */
#endif
    h_EncSetup = encoder->channel_setup[channel];

#ifdef ENABLE_HR_MODE
    hrmode = encoder->hrmode;
#endif
 
#ifndef FIX_1p25_GG_EST_TUPLES 
   
#if defined (FIX_BASOP_ENC_QUANTIZE_1P25MS_512KBPS)   
        /* make sure the initial allocated d_fx buffer size is a multiple of 4 */
        /* ms_mode   Fs    frame_length  frame_ylen     misc             */
        /* 1.25ms  16 kHz  20/2=10         same          2-tuples      */
        /* 1.25ms  24kHz   30/2=15        same           2-tuples action in estimate_global_gain() */
        /* 1.25ms  32kHz   40/2=10        same           2-tuples     */
        /* 1.25ms  48kHz   60/3=20        50/3=16.66      3-tuples, 60/3 =15 allocated , separate action  in estimate_global_gain(),  50==60*40/48 */
        /* 2.5 ms  8 kHz   20/4=5         same            4-tuples    */

        Word16 lg_4  = shr(encoder->frame_length, 2);
        Word16 rem_4 = sub(encoder->frame_length, shl(lg_4, 2));
        if (rem_4 != 0 && sub(encoder->frame_length, encoder->yLen) == 0)
        {
            lg_4 = add(lg_4, 1); /* add one quadruple */
            ASSERT(encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS);
            ASSERT(MAX_LEN >= lg_4 * 4 ); /* "en" size in estimate_global_gain_fx() */
        }
#  ifdef         FIX_BASOP_1p25_NEW_GG_EST3 
        if(encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && (sub(encoder->frame_length, 30) <= 0))   /* WB and SSWB */
        {
            ASSERT(MAX_LEN >= 2 + lg_4 * 4); /*make sure that there are 2 extra tail coeffs  in d_fx*/
        }
#  endif 
#endif
#endif

    BASOP_sub_start("Encoder");
  
#ifdef FIX_1p25_GG_EST_TUPLES
#ifndef CR14_A_ADD_1p25MS_HR
    if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        ASSERT(MAX(80, encoder->frame_length) >= ((GG_1p25_MAX_TUPLES - 1) + encoder->yLen));
        /*make sure that there are  extra tail coeffs in d_fx for Global gain estimation routine */
    }
#endif
#endif

    d_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *d_fx ) * encoder->frame_length );
    L_scf_idx = (Word32*) lc3_scratch_push( scratch, sizeof( *L_scf_idx ) * SCF_MAX_PARAM );
    indexes = (Word16*) lc3_scratch_push( scratch, sizeof( *indexes ) * TNS_NUMFILTERS_MAX * MAXLAG );
#ifdef CR14_A_ADD_LOSSLESS_MODE
    d_tda_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *d_tda_fx) * encoder->frame_length );
#endif 
#  ifdef ENABLE_HR_MODE
    q_d_fx24 = (Word32*) lc3_scratch_push( scratch, sizeof( *q_d_fx24 ) * encoder->frame_length );
#  else
    q_d_fx16 = (Word16*) lc3_scratch_push( scratch, sizeof( *q_d_fx16 ) * encoder->frame_length );
#  endif
    codingdata = (Word16*) lc3_scratch_push( scratch, sizeof( *codingdata ) * ( 3 * encoder->frame_length / 2 ) );
#ifdef CR14_A_ADD_LOSSLESS_MODE
    deterministic_curve = (UWord8*) lc3_scratch_push( scratch, sizeof( *deterministic_curve ) * encoder->frame_length );
    d_fx_orig           = (Word32*) lc3_scratch_push( scratch, sizeof( *d_fx_orig ) * encoder->frame_length );
    tns_lsb_remove      = (UWord8*) lc3_scratch_push( scratch, sizeof( *tns_lsb_remove ) * encoder->frame_length );
#endif
#  ifdef ENABLE_HR_MODE
    s_in_scaled = (Word32*) lc3_scratch_push( scratch, sizeof( *s_in_scaled ) * encoder->frame_length );
#  else
    s_in_scaled = (Word16*) lc3_scratch_push( scratch, sizeof( *s_in_scaled ) * encoder->frame_length );
#  endif

#ifdef ENABLE_HR_MODE
    s_in_scaled_lp = (Word16 *)s_in_scaled;
#endif

    resBits = (UWord8*)lc3_scratch_push( scratch, sizeof( *resBits ) * 48001 );

    s_12k8 = (Word16*) lc3_scratch_push( scratch, sizeof( *s_12k8 ) * ( LEN_12K8 + 1 ) );
  
#  ifdef ENABLE_HR_MODE
    scf_q = (Word32*) lc3_scratch_push( scratch, sizeof( *scf_q ) * M );
#  else
    scf_q = (Word16*) lc3_scratch_push( scratch, sizeof( *scf_q ) * M );
#  endif
  
    scf = (Word16*) lc3_scratch_push( scratch, sizeof( *scf ) * M );
  
#  ifdef ENABLE_HR_MODE
    int_scf_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *int_scf_fx ) * MAX_BANDS_NUMBER );
#  else
    int_scf_fx = (Word16*) lc3_scratch_push( scratch, sizeof( *int_scf_fx ) * MAX_BANDS_NUMBER );
#  endif
    int_scf_fx_exp = (Word16*) lc3_scratch_push( scratch, sizeof( *int_scf_fx_exp ) * MAX_BANDS_NUMBER );

#ifndef CR14_A_ADD_LOSSLESS_MODE
    /* Scale 24-bit input data */
    IF (sub(bits_per_sample, 24) == 0)
    {
        BASOP_sub_start("Scale_signal24");
        scale_signal24_fx(s_in, s_in_scaled, &h_EncSetup->x_exp, h_EncSetup->stEnc_mdct_mem,
                          encoder->stEnc_mdct_mem_len, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len,
                          h_EncSetup->mdct_mem32, encoder->frame_length, h_EncSetup->resamp_mem32,
                          h_EncSetup->olpa_mem_s12k8, &h_EncSetup->resamp_exp);
        BASOP_sub_end();
    }
    ELSE
    {
#ifdef ENABLE_HR_MODE
        Word16 *ip_buf = (Word16*)s_in;
        Word32 i;
        FOR(i = 0; i < encoder->frame_length; i++)
        {
            s_in_scaled[i] = L_deposit_h(ip_buf[i]);
        }
        h_EncSetup->x_exp = 15; move16();
#else
        memcpy(s_in_scaled, s_in, encoder->frame_length * sizeof(*s_in_scaled));
#endif
    }
#endif
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    BASOP_sub_start("Scale_signal24");
  
    bitsPerSample = (bits_per_sample == 16) ? 16 : 24;
    IF (encoder->lossless)
    {        
        /* No input data scaling in VBR lossless mode */
        IF (bits_per_sample == 16)
        {
            FOR (int k = 0; k < encoder->frame_length; k++)
            {
                s_in_scaled[k] = ((Word16*) s_in)[k];
                h_EncSetup->x_exp = 31;
            }
        } ELSE IF (bits_per_sample == 24)
        {
            FOR (int k = 0; k < encoder->frame_length; k++)
            {
                s_in_scaled[k] = ((Word32*) s_in)[k];
                h_EncSetup->x_exp = 23;
            }
        } ELSE {
            assert(0);
        }
    } ELSE IF(hrmode)
    {
        /* Scale 24-bit input data */
    IF (sub(bits_per_sample, 24) == 0)
    {
        BASOP_sub_start("Scale_signal24");
        scale_signal24_fx(s_in, s_in_scaled, &h_EncSetup->x_exp, h_EncSetup->stEnc_mdct_mem,
                          encoder->stEnc_mdct_mem_len, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len,
                          h_EncSetup->mdct_mem32, encoder->frame_length, h_EncSetup->resamp_mem32,
                          h_EncSetup->olpa_mem_s12k8, &h_EncSetup->resamp_exp);
        BASOP_sub_end();
    }
    ELSE
    {
#ifdef ENABLE_HR_MODE
        Word16 *ip_buf = (Word16*)s_in;
        Word32 i;
        FOR(i = 0; i < encoder->frame_length; i++)
        {
            s_in_scaled[i] = L_deposit_h(ip_buf[i]);
        }
        h_EncSetup->x_exp = 15; move16();
#else
        memcpy(s_in_scaled, s_in, encoder->frame_length * sizeof(*s_in_scaled));
#endif
    }
    } ELSE
    {
        format_in_pcm( bits_per_sample, s_in, s_in_scaled, h_EncSetup, encoder );
    }
  
    BASOP_sub_end();

    BASOP_sub_start("Mdct");

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (encoder->lossless && encoder->ll_shift > 0)
    {
        scaleSignal = encoder->ll_shift;
        FOR (Word16 k = 0; k < encoder->frame_length; k++)
        {
            s_in_scaled[k] = L_shr(s_in_scaled[k], scaleSignal);
        }
    }
#endif

#ifdef CR15_A_LOSSLESS_1p25MS
    if (encoder->lossless) {
#else
    if (encoder->lossless && encoder->frame_dms != LC3PLUS_FRAME_DURATION_1p25MS) {
#endif
        Word16 applyMdctRounding = 0;
        if (encoder->lossless) {
            applyMdctRounding = 1;
        }
#ifdef CR14_A_ADD_LOSSLESS_MODE
        IF( (scaleSignal - h_EncSetup->scaleSignal_Memory) > 0)
        {
            FOR(int n = 0; n < encoder->stEnc_mdct_mem_len; n++)
            {
                h_EncSetup->stEnc_mdct_mem[n] = L_shr(h_EncSetup->stEnc_mdct_mem[n], scaleSignal - h_EncSetup->scaleSignal_Memory);
            }
        }
        IF( (h_EncSetup->scaleSignal_Memory - scaleSignal) > 0)
        {
            FOR(int n = 0; n < encoder->stEnc_mdct_mem_len; n++)
            {
                h_EncSetup->stEnc_mdct_mem[n] = L_shl(h_EncSetup->stEnc_mdct_mem[n], h_EncSetup->scaleSignal_Memory - scaleSignal);
            }
        }
        h_EncSetup->scaleSignal_Memory = scaleSignal;
#endif

        const int32_t *lift[4];
        getLiftingCoeffs2_fx(encoder->frame_dms, encoder->fs, lift);
        intMdct2_fx( s_in_scaled, d_fx, d_tda_fx, encoder->frame_length, lift, encoder->la_zeroes, h_EncSetup->stEnc_mdct_mem, &h_EncSetup->stEnc_mdct_mem[encoder->frame_length>>1], applyMdctRounding, scratch );
    }
    else
#endif
    
    {
        processMdct_fx(s_in_scaled, h_EncSetup->x_exp, encoder->frame_length,
#ifdef ENABLE_HR_MODE
                       hrmode,
#endif
                       encoder->W_fx, encoder->W_size,
                       h_EncSetup->stEnc_mdct_mem, encoder->stEnc_mdct_mem_len, d_fx, &d_fx_exp, scratch);
    }
  
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }

    BASOP_sub_end();


#ifdef CR14_A_ADD_LOSSLESS_MODE
    BASOP_sub_start( "PerBandEnergy_NearNyquist" );
    if(encoder->lossless)
    {
        ener_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *ener_fx ) * MAX_BANDS_NUMBER );
        d_fx_tmp = (Word32*) lc3_scratch_push( scratch, sizeof( *d_fx_tmp ) * encoder->frame_length );
        d_fx_exp = 31;

        memcpy( d_fx_tmp, d_fx, sizeof(Word32)*encoder->frame_length);

        Word16 s = getScaleFactor32( d_fx_tmp, encoder->frame_length );

        s = s - 4;
        FOR(int i=0;i<encoder->frame_length;i++)
        {
            d_fx_tmp[i] = L_shl(d_fx_tmp[i], s);
        }
        d_fx_exp_tmp = d_fx_exp - s;

        processPerBandEnergy_fx( ener_fx, &ener_fx_exp, d_fx_tmp, d_fx_exp_tmp - (31-h_EncSetup->x_exp), encoder->bands_offset,
                encoder->fs_idx, encoder->bands_number, 0, encoder->frame_dms, scratch
#ifdef ENABLE_HR_MODE
                , encoder->hrmode
#endif
        );

        encoder->near_nyquist_index = 0;
        processNearNyquistdetector_fx( &encoder->near_nyquist_flag, encoder->fs_idx, encoder->near_nyquist_index,
                                   encoder->bands_number, ener_fx, ener_fx_exp
#ifdef ENABLE_HR_MODE
                                   , encoder->frame_dms, encoder->hrmode );
#else
        );
#endif
        d_fx_tmp = (Word32*) lc3_scratch_pop( scratch, d_fx_tmp );
        ener_fx = (Word32*) lc3_scratch_pop( scratch, ener_fx );
    }

    BASOP_sub_end(); // "PerBandEnergy_NearNyquist"
#endif
    tns_order[0] = 0;
    tns_order[1] = 0;

#ifdef CR14_A_ADD_LOSSLESS_MODE
    tns_numfilters = 0;
    if (encoder->lossless)
    {
    tns_order[0] = 0;
    tns_order[1] = 0;
        if (encoder->ll_cbr && encoder->frame_dms != LC3PLUS_FRAME_DURATION_2p5MS)
        {
            BW_cutoff_idx = encoder->fs_idx;
            memcpy( d_fx_orig, d_fx, sizeof(Word32)*encoder->frame_length);
            BASOP_sub_start( "Tns_enc" );
            processAdaptiveTns_fx( &( h_EncSetup->tns_bits ), indexes, d_fx, BW_cutoff_idx, tns_order, &tns_numfilters, h_EncSetup->enable_lpc_weighting, encoder->nSubdivisions, encoder->frame_dms,
                                encoder->frame_length, scratch, encoder->hrmode, encoder->near_nyquist_flag, &ll_adap_flag, encoder->fs_idx, h_EncSetup->total_bits, encoder->ll_tns_lsb_num_remove_limit,
                                &tns_lsb_num_remove);
            BASOP_sub_end(); // "Tns_enc"
        }
        else
        {
        BASOP_sub_start( "Tns_enc" );

        BW_cutoff_idx = encoder->fs_idx;
        Word32 predictionGain[2] = {0};
        Word16 RC_1[MAXLAG];
        Word16 RC_2[MAXLAG];
        Word16 *RC[2];
        RC[0] = RC_1;
        RC[1] = RC_2;
        TnsStartStopFreqs startstopfreqs = {{0,0},0}; 
        processTnsCoder_fx( &( h_EncSetup->tns_bits ), indexes, d_fx, BW_cutoff_idx, tns_order, &tns_numfilters,
                            h_EncSetup->enable_lpc_weighting, encoder->nSubdivisions, encoder->frame_dms,
                            encoder->frame_length, scratch, encoder->hrmode,
                            encoder->near_nyquist_flag, predictionGain, 0, RC, &startstopfreqs, encoder->lossless);
        
        BASOP_sub_end(); // "Tns_enc"
    }
#  endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    }
  
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }
  
    if(encoder->ll_est_bit_usage && ll_adap_flag == 0 && encoder->ll_cbr)
    {
        BASOP_sub_start( "BitUsage_enc" );
        Word32 bit_usage_estimate = 0;

        #ifdef LOSSLESS_192kHz
        IF(encoder->fs_idx == 6)
        {
            active_bits  = calculate_active_bits(d_fx, encoder->frame_length / 2);
            entropy_bits = calculate_active_bits(&d_fx[encoder->frame_length / 2], encoder->frame_length / 2);
        }
        ELSE
        {
            active_bits = calculate_active_bits(d_fx, encoder->frame_length);
        }
        #else
             active_bits = calculate_active_bits(d_fx, encoder->frame_length);
        #endif 
        #ifdef LOSSLESS_192kHz
        bit_usage_estimate = estimate_bit_usage( encoder->fs_idx, active_bits, entropy_bits);
        #else
        bit_usage_estimate = estimate_bit_usage( encoder->fs_idx, active_bits);
        #endif 
        ll_adap_flag = get_ll_adap_flag( bit_usage_estimate, h_EncSetup->total_bits, encoder->ll_bit_balance );

        BASOP_sub_end();// BASOP_sub_start( "BitUsage_enc" );
    }
  
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }

    // Count WMOPS for lossless vs. lossy code path
    if ( 1 == ll_adap_flag )
    {
        BASOP_sub_sub_start( "Enc(lossless)" );
    }
    else
    {
        BASOP_sub_sub_start( "Enc(lossy)" );
    }
    
    IF(ll_adap_flag)
    {
        d_fx_exp = 31;
    }
    ELSE
    {
        if ( encoder->lossless )
        {
            /* scaling for lossy mode; do not apply for lossless*/
            Word16 s = getScaleFactor32( d_fx, encoder->frame_length );
            s = s - 4;
            d_fx_exp = h_EncSetup->x_exp - s;
            FOR( int i = 0; i < encoder->frame_length; i++ )
            {
                d_fx[i] = L_shl( d_fx[i], s );        move32();
            }
        }
    }
#endif
  
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }

    /* begin s_12k8 */
    BASOP_sub_start("Resamp12k8");
#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 LTPF_shift = 16;
  
    if (encoder->lossless) {
        LTPF_shift = s_min(16, 31 - h_EncSetup->x_exp);
    }
  
    downshift_w32_arr(s_in_scaled, s_in_scaled_lp, LTPF_shift, encoder->frame_length);
#else
    downshift_w32_arr(s_in_scaled, s_in_scaled_lp, 16, encoder->frame_length);
#endif
    /* s_in_scaled is no longer required */

    process_resamp12k8_fx(s_in_scaled_lp, encoder->frame_length, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len, s_12k8,
                          &s_12k8_len, encoder->fs_idx, encoder->frame_dms, scratch
                          , bits_per_sample
                         );
#else
    process_resamp12k8_fx(s_in_scaled, encoder->frame_length, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len, s_12k8,
                          &s_12k8_len, encoder->fs_idx, encoder->frame_dms, scratch
                          , bits_per_sample
                          );
#endif /* ENABLE_HR_MODE */
    BASOP_sub_end();
  
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }

    BASOP_sub_start("Olpa");
    process_olpa_fx(&h_EncSetup->olpa_mem_s6k4_exp, h_EncSetup->olpa_mem_s12k8, h_EncSetup->olpa_mem_s6k4, &pitch,
                    s_12k8, s_12k8_len, &normcorr, &h_EncSetup->olpa_mem_pitch,
                    &h_EncSetup->pitch_flag,
                    h_EncSetup->resamp_exp, encoder->frame_dms, scratch);
    BASOP_sub_end();
  
//    printf("pitch = %d\n", pitch);
//    printf("h_EncSetup->ltpf_enable = %d\n", h_EncSetup->ltpf_enable);
//    printf("h_EncSetup->ltpf_mem_normcorr = %d\n", *h_EncSetup->ltpf_mem_normcorr);
//    printf("normcorr = %d\n", normcorr);
//    printf("h_EncSetup->resamp_exp = %d\n", h_EncSetup->resamp_exp);
//  
//    printf("encoder->ltpf_mem_in_len = %d\n", encoder->ltpf_mem_in_len);
//    printf("s_12k8_len = %d\n", s_12k8_len);
//    printf("h_EncSetup->ltpf_mem_pitch = %d\n", h_EncSetup->ltpf_mem_pitch);

    BASOP_sub_start("LtpfEnc");
    process_ltpf_coder_fx(&ltpf_bits, pitch, h_EncSetup->ltpf_enable, &h_EncSetup->ltpf_mem_in_exp,
                          h_EncSetup->ltpf_mem_in, encoder->ltpf_mem_in_len, ltpf_idx, s_12k8, s_12k8_len,
                          h_EncSetup->ltpf_mem_normcorr, &h_EncSetup->ltpf_mem_mem_normcorr, normcorr,
                          &h_EncSetup->ltpf_mem_ltpf_on, &h_EncSetup->ltpf_mem_pitch, h_EncSetup->resamp_exp,
                          encoder->frame_dms, scratch
                          , encoder->hrmode
 
#ifdef CR9_C_ADD_1p25MS
#ifdef FIX_TX_RX_STRUCT_STEREO
#    ifdef NEW_SIGNALLING_SCHEME_1p25
                          ,&h_EncSetup->Tx_ltpf
#    else
                          ,h_EncSetup->Tx_ltpf
#    endif
#else
                          , encoder->Tx_ltpf
#endif
#endif
);
    BASOP_sub_end();
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag)
    {
        //set the ltpf active flag to zero in case of lossless frames
        ltpf_idx[1] = 0;
    }
#endif

//    printf("ltpf_idx[0] = %d\n", ltpf_idx[0]);
//    printf("ltpf_idx[1] = %d\n", ltpf_idx[1]);
//    printf("ltpf_idx[2] = %d\n", ltpf_idx[2]);

    /* end s_12k8 */
    BASOP_sub_start("AttackDetector");
#ifdef ENABLE_HR_MODE
    attack_detector_fx(encoder, h_EncSetup, s_in_scaled_lp, sub(h_EncSetup->x_exp, 15), scratch);
#else
    attack_detector_fx(encoder, h_EncSetup, s_in_scaled, sub(h_EncSetup->x_exp, 15), scratch);
#endif
    BASOP_sub_end();
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (encoder->lossless)
    {
        h_EncSetup->attdec_detected = 0;
    }

    IF (ll_adap_flag)
    {
        active_bits = calculate_active_bits( d_fx, encoder->frame_length );
        entropy_bits = calculate_active_bits(&d_fx[encoder->frame_length / 2], encoder->frame_length / 2);
  
        IF(encoder->fs_idx == 4 && bits_per_sample == 16)
        {
            encoder->ll_ari_bits = L_add(Mpy_32_32(-241591910,active_bits),2100);
            encoder->ll_ari_bits =  MAX(0, MIN(encoder->ll_ari_bits,h_EncSetup->total_bits));
            encoder->ll_ari_bits =  Mpy_32_32_0(L_mult0(encoder->ll_ari_bits, encoder->frame_length), 8947849);
        }
        IF(encoder->fs_idx == 5 && bits_per_sample == 24)
        {
            encoder->ll_ari_bits = L_add(Mpy_32_32(-858993459,active_bits),7800);
            encoder->ll_ari_bits =  MAX(0, MIN(encoder->ll_ari_bits,6000));
            encoder->ll_ari_bits =  Mpy_32_32_0(L_mult0(encoder->ll_ari_bits, encoder->frame_length), 4473924);
        }
  
        IF(encoder->fs_idx == 6)
        {
            IF( bits_per_sample == 24)
            {
                encoder->ll_ari_bits = -0.1906 * active_bits + 7978;
                encoder->ll_ari_bits =  MAX( 0, MIN(encoder->ll_ari_bits, 6000) );
            }
            ELSE
            {
                encoder->ll_ari_bits = -0.175 * active_bits + 7038;
                encoder->ll_ari_bits =  MAX( 0, MIN(encoder->ll_ari_bits, 6000) );
            }
            encoder->ll_ari_bits =  encoder->ll_ari_bits *  (1000.0f * encoder->frame_length / 192000.0f )  / 10;
        }
  
        encoder->ll_ari_bits = MIN(encoder->ll_ari_bits,h_EncSetup->total_bits);

        IF( encoder->ll_offQuant )
        {
            IF( encoder->fs_idx >= 5 && bits_per_sample == 24 )
            {
                IF( encoder->fs_idx >= 5 && bits_per_sample == 24 )
                {
                    Word32 fact[7] = {214748365, 107374183, 71582789, 53687092, 42949673, 35791395, 30678338}; 
                    gg_idx_off_ll_adap = L_negate( L_add(s_min( 115, Mpy_32_32(encoder->ll_ari_bits, fact[encoder->fs_idx])), L_add(105, L_mult0( 5, ( encoder->fs_idx + 1 ))) ));

                    off_idx = extract_l(L_negate(Mpy_32_32_0( L_add(135,gg_idx_off_ll_adap), 300647711)));
                    encoder->quantizedGainOff_ll = extract_l(L_sub(L_negate((Mpy_32_32_0(L_mult0(off_idx,100), 613566757))),135));
                }
                IF( bits_per_sample == 24 )
                {
                    encoder->quantizedGainOff_ll = encoder->quantizedGainOff_ll + 48;
                }
            }
        }
    }

    IF( bits_per_sample == 24 && encoder->fs_idx == 6 && ll_adap_flag )
    {
        active_bits_lb = calculate_active_bits(d_fx, encoder->bands_offset[encoder->low_band_limit] - 1);
        ll_ari_bits_lb = L_add(Mpy_32_32(-858993459,active_bits_lb),7800);
        ll_ari_bits_lb =  MAX(0, MIN(ll_ari_bits_lb,6000)); 

        Word16 tmp2 = extract_l(L_shr(L_mult(encoder->bands_offset[encoder->low_band_limit], 34),1));
        ll_ari_bits_lb = L_shr(L_mult(ll_ari_bits_lb,tmp2),16); 
                        
        Word32 gg_idx_off_ll_adap = L_negate( L_add(s_min( 115, Mpy_32_32(ll_ari_bits_lb, 35791395 )), L_add(105, L_mult0( 5, ( encoder->fs_idx ))) ));
        off_idx = extract_l(L_negate(Mpy_32_32_0( L_add(135,gg_idx_off_ll_adap), 300647711)));            
        encoder->quantizedGainOff_ll = extract_l(L_sub(L_negate((Mpy_32_32_0(L_mult0(off_idx,100), 613566757))),135)) + 48;
    }

#endif

    /* begin ener_fx */
    #ifdef CR14_A_ADD_LOSSLESS_MODE
        d_fx_tmp = (Word32*) lc3_scratch_push( scratch, sizeof( *d_fx_tmp ) * encoder->frame_length );
    #endif 
    ener_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *ener_fx ) * MAX_BANDS_NUMBER );
    BASOP_sub_start("PerBandEnergy");
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag)
    {
        memcpy( d_fx_tmp, d_fx, sizeof(Word32)*encoder->frame_length);

        Word16 s = getScaleFactor32( d_fx_tmp, encoder->frame_length );

        s = s - 4;
        FOR(int i=0;i<encoder->frame_length;i++)
        {
            d_fx_tmp[i] = L_shl(d_fx_tmp[i], s);
        }
        d_fx_exp_tmp = d_fx_exp - s;

        IF( encoder->fs_idx == 6 && bits_per_sample == 24)
        {
            processPerBandEnergy_fx( ener_fx, &ener_fx_exp, d_fx_tmp, d_fx_exp_tmp - (31-h_EncSetup->x_exp), encoder->bands_offset, encoder->fs_idx,
                                    encoder->low_band_limit, 0, encoder->frame_dms, scratch, encoder->hrmode ); 
        }
        ELSE
        {  
            processPerBandEnergy_fx(ener_fx, &ener_fx_exp, d_fx_tmp, d_fx_exp_tmp - (31-h_EncSetup->x_exp), encoder->bands_offset, encoder->fs_idx,
                                    encoder->bands_number, 0, encoder->frame_dms, scratch
    #ifdef ENABLE_HR_MODE
                                    , encoder->hrmode
    #endif
            );
        }
    }
    ELSE
    {
        processPerBandEnergy_fx(ener_fx, &ener_fx_exp, d_fx, d_fx_exp, encoder->bands_offset, encoder->fs_idx,
                                encoder->bands_number, 0, encoder->frame_dms, scratch
#ifdef ENABLE_HR_MODE
                                , encoder->hrmode
#endif
        );
    }

#else 
  
    processPerBandEnergy_fx(ener_fx, &ener_fx_exp, d_fx, d_fx_exp, encoder->bands_offset, encoder->fs_idx,
                            encoder->bands_number, 0, encoder->frame_dms, scratch
#ifdef ENABLE_HR_MODE
                            , encoder->hrmode
#endif
    );
#endif
    BASOP_sub_end();

#ifdef CR14_A_ADD_LOSSLESS_MODE
IF( !(encoder->fs_idx == 6 && bits_per_sample == 24) )
{
#endif 

    BASOP_sub_start("Near Nyquist Detector");
        /* Near Nyquist Detector */
        processNearNyquistdetector_fx(&encoder->near_nyquist_flag, encoder->fs_idx, encoder->near_nyquist_index,
                                      encoder->bands_number, ener_fx, ener_fx_exp
#ifdef ENABLE_HR_MODE
                                  , encoder->frame_dms, encoder->hrmode );
#else
                                  );
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
}
#endif 
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    if (encoder->lossless)
    {
        encoder->near_nyquist_flag = 0;
    }
#endif
  
        /* Disable LTPF if nyquist detector triggers */
        IF (encoder->near_nyquist_flag != 0 || sub(h_EncSetup->lfe, 1) == 0)
        {
            h_EncSetup->ltpf_mem_ltpf_on = 0;  move16();
            ltpf_idx[1] = 0;  move16();
        }
    BASOP_sub_end();
    BASOP_sub_start("BW Cutoff-Detection");
    IF (h_EncSetup->lfe == 0)
    {
#ifdef ENABLE_HR_MODE
    /* No BW Cutoff for 8 kHz and 96 kHz */
    IF (encoder->fs_idx > 0 && encoder->hrmode == 0 && encoder->bw_ctrl_active == 0)
    {
#else  /* ENABLE_HR_MODE */
    IF (encoder->fs_idx > 0 && encoder->bw_ctrl_active == 0)
    {
#endif /* ENABLE_HR_MODE */
        processDetectCutoffWarped_fx(&BW_cutoff_idx, ener_fx, ener_fx_exp, encoder->fs_idx, encoder->frame_dms);
    }
    ELSE
    {
        BW_cutoff_idx = encoder->fs_idx;
        move16();
    }
    }
    ELSE
    {
        BW_cutoff_idx = 0;
    }
    BASOP_sub_end();
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    if (encoder->lossless)
    {
        BW_cutoff_idx = encoder->fs_idx;
    }
#endif

    BASOP_sub_start("SnsCompScf");
      
#if defined( CR14_A_ADD_LOSSLESS_MODE)
    IF( ll_adap_flag && bits_per_sample == 24 )
    {
        ener_fx_exp = encoder->fs_idx == 6 ? ener_fx_exp  : ener_fx_exp + 8;
        if( ener_fx_exp > 43 )
            ener_fx_exp = 43;
    }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE  
    IF( encoder->fs_idx == 6 && ll_adap_flag && bits_per_sample == 24 )
    {
        processSnsComputeScf_fx( ener_fx, ener_fx_exp, encoder->fs_idx, encoder->low_band_limit, scf,
                                    h_EncSetup->attdec_detected, encoder->attdec_damping,
                                    scratch, encoder->sns_damping
#ifdef CR9_C_ADD_1p25MS
                            , encoder->frame_dms, normcorr, &encoder->LT_normcorr
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                             , ll_adap_flag, hrmode
#endif
        );
    }
    ELSE
    {
#endif

        processSnsComputeScf_fx(ener_fx, ener_fx_exp, encoder->fs_idx, encoder->bands_number, scf,
                            h_EncSetup->attdec_detected, encoder->attdec_damping, scratch, encoder->sns_damping
#ifdef CR9_C_ADD_1p25MS
                            , encoder->frame_dms, normcorr, &encoder->LT_normcorr
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                             , ll_adap_flag, hrmode
#endif
        );

#ifdef CR14_A_ADD_LOSSLESS_MODE
    }
#endif


    BASOP_sub_end();

    ener_fx = (Word32*) lc3_scratch_pop( scratch, ener_fx );

#ifdef CR9_C_ADD_1p25MS_LRSNS
    IF(sub(encoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)
    {
        pitch_rx_fx = ltpf_idx[0];  /* pitch_rx status flag   */
        ltpf_rx_fx = ltpf_idx[1];    /* ltpf_activation  used in  snslr_st1C mode   */
      
        BASOP_sub_start("EncLC3_SnsQuantScfEncLR_fx");
        envelope_bits_fx = snsQuantScfEncLR_fx(
            scf,             /*  input scf is always Word16 Q11 (both for ENABLE_HR and DISABLE_HR) */
            L_scf_idx,       /* output:  to send to enc_entropy_fx() */
            scf_q,           /* output:  Word32(for ENABLE_HR)   or Word16 (for DISABLE_HR ) */
            pitch_rx_fx, ltpf_rx_fx, /* input: pitch information to st1C */
            scratch
        );
        BASOP_sub_end();
    }
    ELSE
    {
        BASOP_sub_start("SnsQuantScfEnc38bit");
        processSnsQuantizeScfEncoder_fx(scf, L_scf_idx, scf_q, scratch);
        BASOP_sub_end();
    }
#else
    BASOP_sub_start("SnsQuantScfEnc");
    processSnsQuantizeScfEncoder_fx(scf, L_scf_idx, scf_q, scratch);
    BASOP_sub_end();
#endif /*  CR9_C_ADD_1p25MS_LRSNS */

    BASOP_sub_start("SnsInterpScfEnc");
    #ifdef CR14_A_ADD_LOSSLESS_MODE
    IF(encoder->fs_idx == 6 && bits_per_sample == 24 && ll_adap_flag)
    {
        processSnsInterpolateScf_fx(scf_q, int_scf_fx, int_scf_fx_exp, 1, encoder->low_band_limit, scratch);
    }
    ELSE
    {
    #endif 
        processSnsInterpolateScf_fx(scf_q, int_scf_fx, int_scf_fx_exp, 1, encoder->bands_number, scratch);
    #ifdef CR14_A_ADD_LOSSLESS_MODE
    }
    #endif 
    BASOP_sub_end();

    BASOP_sub_start("Mdct shaping_enc");

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF(ll_adap_flag)
    {
        Word16 max_scf_exp = 0;
        Word16 bands_limit = ( encoder->fs_idx == 6 && bits_per_sample == 24 ) ? encoder->low_band_limit : encoder->bands_number;
        FOR( int i = 0; i < bands_limit; i++ )
            max_scf_exp = s_max( max_scf_exp, int_scf_fx_exp[i] );

        Word16 extra_margin = s_max( 0, max_scf_exp - 5 );
        IF( extra_margin > 0 )
        {
            FOR( int i = 0; i < encoder->frame_length; i++ )
            {
                d_fx_tmp[i] = L_shr( d_fx_tmp[i], extra_margin );
                move32();
            }
            d_fx_exp_tmp += extra_margin;
        }

        IF( encoder->fs_idx == 6 && bits_per_sample == 24 )
        {
            processMdctShaping_fx( d_fx_tmp, int_scf_fx, int_scf_fx_exp, encoder->bands_offset, encoder->low_band_limit );
        }
        ELSE
        {
            processMdctShaping_fx( d_fx_tmp, int_scf_fx, int_scf_fx_exp, encoder->bands_offset, encoder->bands_number );
        }
    }
    ELSE
    {
        processMdctShaping_fx( d_fx, int_scf_fx, int_scf_fx_exp, encoder->bands_offset, encoder->bands_number );
    }
#else
    processMdctShaping_fx(d_fx, int_scf_fx, int_scf_fx_exp, encoder->bands_offset, encoder->bands_number);
#endif

    BASOP_sub_end();
      
    /* end int_scf_fx_exp */
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    //only for lossless path 
    Word16 tmp_int_scf_fx_exp[MAX_BANDS_NUMBER];
    Word16 tmp_int_scf_fx_exp_requantized[MAX_BANDS_NUMBER];
    Word32 tmp_int_scf_fx[MAX_BANDS_NUMBER];
        
    memcpy(tmp_int_scf_fx_exp, int_scf_fx_exp, sizeof(*int_scf_fx_exp) * encoder->bands_number);
    memcpy(tmp_int_scf_fx, int_scf_fx, sizeof(*int_scf_fx) * encoder->bands_number);
#endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (encoder->lossless == 0)
    {
#endif
    BASOP_sub_start("BandwidthControl_enc");
    if (encoder->bandwidth < L_shr_pos(encoder->fs, 1))
    {
        process_cutoff_bandwidth(d_fx, encoder->yLen, encoder->bw_ctrl_cutoff_bin);
        BW_cutoff_idx = s_min(BW_cutoff_idx, encoder->bw_index);
    }
    BASOP_sub_end();
    BASOP_sub_start("Tns_enc");

#ifdef CR9_C_ADD_1p25MS
    IF (h_EncSetup->lfe == 0 && encoder->frame_dms > LC3PLUS_FRAME_DURATION_1p25MS)
#else
    IF (h_EncSetup->lfe == 0)
#endif
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 predictionGain[2] = {0};
            
    Word16 RC_1[MAXLAG];
    Word16 RC_2[MAXLAG];
    Word16 *RC[2];
    RC[0] = RC_1;
    RC[1] = RC_2;
    TnsStartStopFreqs startstopfreqs = {{0,0},0};
#endif
      
    processTnsCoder_fx(&(h_EncSetup->tns_bits), indexes, d_fx, BW_cutoff_idx, tns_order, &tns_numfilters,
                       h_EncSetup->enable_lpc_weighting, encoder->nSubdivisions, encoder->frame_dms,
                       encoder->frame_length, scratch
#ifdef ENABLE_HR_MODE
                       , encoder->hrmode
#endif
                       , encoder->near_nyquist_flag
#ifdef CR14_A_ADD_LOSSLESS_MODE
                       , predictionGain, 0, RC, &startstopfreqs, encoder->lossless
#endif
    );
        }
    ELSE
    {
        tns_numfilters = 1;
        move16();
#ifdef CR9_C_ADD_1p25MS
        IF (encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
            tns_numfilters = 0; move16();
        }
#endif
        tns_order[0] = 0;
        move16();
        h_EncSetup->tns_bits = tns_numfilters;
        move16();
    }
    BASOP_sub_end();

#ifdef CR14_A_ADD_LOSSLESS_MODE
    }
#endif
      
    BASOP_sub_start("Est. Global Gain");
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsInit, (Word32)(add( h_EncSetup->tns_bits, ltpf_bits )) );
#else
    h_EncSetup->targetBitsQuant = sub(h_EncSetup->targetBitsInit, add(h_EncSetup->tns_bits, ltpf_bits));
#endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (encoder->lossless)
    {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->total_bits, (Word32)(add( encoder->envelope_bits, encoder->global_gain_bits) + getLastNzBits_fx( encoder->frame_length ) + 2 + 1) ) - ltpf_bits ;

        IF (encoder->ll_tns)
        {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, (Word32)h_EncSetup->tns_bits);
        }

        IF (encoder->ll_tns_remove && (tns_order[0] + tns_order[1] > 0) && encoder->ll_cbr)
        {
            IF (tns_lsb_num_remove == 0)
            {
                h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 1;
            }
            ELSE
            {
                IF (bits_per_sample == 24)
                {
                    h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 4;
                }
                IF (bits_per_sample == 16)
                {
                    h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 3;
                }
            }
        }

        h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 1; /* ll_adap_flag */
        IF (ll_adap_flag)
        {
            h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 1; /* b_relative */
        }
        IF (encoder->ll_offQuant && ll_adap_flag)
        {
            h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 3; /* quantized gain offset */
        }

        IF (ll_adap_flag == 0)
        {
            h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - 3; /* NOISE_FAC_BITS */

            #ifdef LOSSLESS_192kHz
#ifdef CR15_C_VARIOUS_LOSSLESS_FIXES
                h_EncSetup->targetBitsQuant = L_sub(h_EncSetup->targetBitsQuant, L_max(L_shr(h_EncSetup->total_bits, 11), 2)); /*  AC Finalization overhead margin   */
#else
                h_EncSetup->targetBitsQuant = L_sub(h_EncSetup->targetBitsQuant, L_max(L_shr(h_EncSetup->total_bits, 11), 1)); /*  AC Finalization overhead margin   */
#endif
            #else
                h_EncSetup->targetBitsQuant = sub(h_EncSetup->targetBitsQuant, shr_pos(h_EncSetup->total_bits, 11)); /* targetBitsQuant -= total_bits / 2048 */
            #endif
        }
        
        IF (encoder->ll_tns == 0)
        {
            h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsQuant - encoder->BW_cutoff_bits; /* BW_CUTOFF_BITS */
        }
        
        IF ( h_EncSetup->total_bits > 1280 )
        {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, 1 );
        }
        IF ( h_EncSetup->total_bits > 2560 )
        {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, 1 );
        }
        IF ( hrmode )
        {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, 1 );
        }
#ifdef CR14_A_ADD_LOSSLESS_MODE
        h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, 4 );
#endif
#ifdef LL_INCL_HPVC
        IF (h_EncSetup->hpvcEncCfg.active_flag != 0 && encoder->lossless != 0 && ll_adap_flag != 0)
        {
            h_EncSetup->targetBitsQuant = L_sub( h_EncSetup->targetBitsQuant, 1 );
        }
#endif
    }

    IF (ll_adap_flag)
    {
        IF ( bits_per_sample == 24 && encoder->fs_idx == 6)
        {
#ifdef CR15_C_VARIOUS_LOSSLESS_FIXES
            h_EncSetup->targetBitsQuant = L_sub(h_EncSetup->targetBitsQuant,16);
#else
            h_EncSetup->targetBitsQuant = L_sub(h_EncSetup->targetBitsQuant,12);
#endif
        }
        encoder->ll_ari_bits = MIN(encoder->ll_ari_bits, h_EncSetup->targetBitsQuant);
    }
#endif /* CR14_A_ADD_LOSSLESS_MODE */

#  ifdef CR9_C_ADD_1p25MS_LRSNS
    IF( sub(encoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0 )
    {    /* adjust(==reduce) target bits based on the selected LRSNS-VQ  bitrate  */

        h_EncSetup->targetBitsQuant = add(h_EncSetup->targetBitsQuant, 38);    /*  legacy 38bits was already pre-subtracted  in h_EncSetup->targetBitsInit setup  */
        ASSERT(envelope_bits_fx >= 9 && envelope_bits_fx <= 30);
        h_EncSetup->targetBitsQuant = sub(h_EncSetup->targetBitsQuant, envelope_bits_fx); /*  9,10, 29,30, subtract  actual LRSNS bitrate  */
    }
#  endif

    test();
#ifdef  CR9_C_ADD_1p25MS_LRSNS
    /* incoming (state based)  ltpf_bits for tranmission set in function ltpf_coder_fx()  */
    test(); test();
    lrsns_st1C_in_use = 0;      move16();
    /* LRSNS stage1C use 10 bits  and ltp/ltpf active flags ,  i.e. one can thus not always disable/steal the ltp/ltpf flags  */
    if ( (sub(encoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0) && (sub(envelope_bits_fx, 10) == 0) && (L_sub(L_scf_idx[0], 170L) > 0)  )
    {
        ASSERT(L_scf_idx[0] < (2 * 170) && (L_scf_idx[1] != 0) );
        lrsns_st1C_in_use = 1;             move16();
    }

    ltptx_lowest_bit_lim = 1;   move16();
    if ( sub(encoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0 )
    {
        ltptx_lowest_bit_lim = lrsns_ltp_bits_fx[0];     move16();
    }
    test(); test(); test(); test();

    /* allow cut away of LTP-active, ltpf-active, lag-index (including phase A, B) info when possible and not actually in use by LRSNS-VQ */

    IF((h_EncSetup->targetBitsQuant < 0) && (sub(ltpf_bits, ltptx_lowest_bit_lim) > 0) && (lrsns_st1C_in_use == 0))
#else
    IF (h_EncSetup->targetBitsQuant < 0 && sub(ltpf_bits, 1) > 0)
#endif
    {
        /* Disable LTPF */
        h_EncSetup->ltpf_mem_ltpf_on = 0;  move16();
        ltpf_idx[1]                  = 0;  move16();
#ifdef CR9_C_ADD_1p25MS_LRSNS
        ASSERT( (ltpf_bits-ltptx_lowest_bit_lim) > 0 );
        h_EncSetup->targetBitsQuant = add(h_EncSetup->targetBitsQuant, sub(ltpf_bits, ltptx_lowest_bit_lim)); move16(); /* subtract saving */
        ltpf_bits = ltptx_lowest_bit_lim;  move16();
#else
        ltpf_bits                    = 1;  move16();
        h_EncSetup->targetBitsQuant  = sub(h_EncSetup->targetBitsInit, add(h_EncSetup->tns_bits, ltpf_bits));
#endif
    }

#ifdef ENABLE_HR_MODE
    Word32 gain32;
#endif
      
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx[%d] = %d\n", k, d_fx[k]);
//        }
      
//        for (int k = 0; k < encoder->frame_length; k++)
//        {
//            printf("d_fx_tmp[%d] = %d\n", k, d_fx_tmp[k]);
//        }
//      
//        printf("d_fx_exp_tmp = %d\n", d_fx_exp_tmp);
//        printf("encoder->yLen = %d\n", encoder->yLen);
//        printf("ll_adap_flag = %d\n", ll_adap_flag);
//        printf("h_EncSetup->targetBitsQuant = %d\n", h_EncSetup->targetBitsQuant);
//        printf("h_EncSetup->quantizedGainOff = %d\n", h_EncSetup->quantizedGainOff);
//        printf("h_EncSetup->mem_specBits = %d\n", h_EncSetup->mem_specBits);
//        printf("h_EncSetup->mem_targetBits = %d\n", h_EncSetup->mem_targetBits);
//        printf("encoder->hrmode = %d\n", encoder->hrmode);
//        printf("h_EncSetup->regBits = %d\n", h_EncSetup->regBits);
//        printf("encoder->frame_dms = %d\n", encoder->frame_dms);
//        printf("encoder->ll_ari_bits = %d\n", encoder->ll_ari_bits);
//        printf("encoder->quantizedGainOff_ll = %d\n", encoder->quantizedGainOff_ll);
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (!ll_adap_flag)
    {
            processEstimateGlobalGain_fx( d_fx, d_fx_exp, encoder->yLen, h_EncSetup->targetBitsQuant,
                                  &gain32,
                                  &gain_e,
                                  &quantizedGain, &quantizedGainMin, h_EncSetup->quantizedGainOff,
                                  &h_EncSetup->targetBitsOff, &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits,
                                  scratch, encoder->hrmode, h_EncSetup->regBits, encoder->frame_dms);
    } ELSE {

        IF( encoder->fs_idx == 6 && bits_per_sample == 24 )
        {
            processEstimateGlobalGain_fx( d_fx_tmp, d_fx_exp_tmp, encoder->bands_offset[encoder->low_band_limit], ll_ari_bits_lb,
                                   &gain32,
                                   &gain_e,
                                   &quantizedGain, &quantizedGainMin, encoder->quantizedGainOff_ll,
                                   &h_EncSetup->targetBitsOff, &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits,
                                   scratch, encoder->hrmode, h_EncSetup->regBits, encoder->frame_dms);
        }
        ELSE
        {
            processEstimateGlobalGain_fx( d_fx_tmp, d_fx_exp_tmp, encoder->yLen, encoder->ll_ari_bits,
                                  &gain32,
                                  &gain_e,
                                  &quantizedGain, &quantizedGainMin, encoder->quantizedGainOff_ll,
                                  &h_EncSetup->targetBitsOff, &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits,
                                  scratch, encoder->hrmode, h_EncSetup->regBits, encoder->frame_dms);
        }
    }
    d_fx_tmp = (Word32*) lc3_scratch_pop( scratch, d_fx_tmp );
#else /* CR14_A_ADD_LOSSLESS_MODE */
    processEstimateGlobalGain_fx(d_fx, d_fx_exp, encoder->yLen, h_EncSetup->targetBitsQuant,
#ifdef ENABLE_HR_MODE
                                 &gain32,
#else
                                 &gain,
#endif
                                 &gain_e,
                                 &quantizedGain, &quantizedGainMin, h_EncSetup->quantizedGainOff,
                                 &h_EncSetup->targetBitsOff, &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits,
                                 scratch
#ifdef ENABLE_HR_MODE
                                 , encoder->hrmode, h_EncSetup->regBits, encoder->frame_dms
#else
#if defined(FIX_BOTH_1p25_TEST_NEW_GG_EST2) || defined (FIX_1p25_GG_EST_TUPLES) 
                                 ,  encoder->frame_dms
#endif 
#endif
    );
#endif /* CR14_A_ADD_LOSSLESS_MODE */
      
//    printf("gain = %d\n", gain);
//    printf("quantizedGain = %d\n", quantizedGain);
      
    BASOP_sub_end();
    /* begin q_d_fx16 */
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag)
    {
#ifdef LL_INCL_HPVC
        /* Reuse outer-scope buffer so contents persist for HPVC analysis. */
        Word32 *truncated_data = hpvc_truncated_data;
#else
        Word32* truncated_data = (Word32*) lc3_scratch_push( scratch, sizeof( *truncated_data ) * encoder->frame_length );
#endif

        memcpy(tmp_int_scf_fx_exp_requantized, tmp_int_scf_fx_exp, sizeof(*tmp_int_scf_fx_exp) * encoder->bands_number);

        BASOP_sub_start( "Quant. 1::DetCurve" );
        IF( encoder->fs_idx == 6 && bits_per_sample == 24 )
        {
             process_deterministic_curve( deterministic_curve, gain32, gain_e, encoder->bands_offset[encoder->low_band_limit], tmp_int_scf_fx, tmp_int_scf_fx_exp, encoder->bands_offset, encoder->low_band_limit, scratch );
             lsb_mean_split(d_fx, encoder->low_band_limit, encoder->bands_number, encoder->bands_offset, deterministic_curve[ encoder->bands_offset[encoder->low_band_limit]-1 ],
                            ll_deltaCodesBits,
                            deterministic_curve );
        }
        ELSE
        {
            process_deterministic_curve( deterministic_curve, gain32, gain_e, encoder->frame_length, tmp_int_scf_fx, tmp_int_scf_fx_exp, encoder->bands_offset, encoder->bands_number, scratch );
        }
        BASOP_sub_end();
        BASOP_sub_start( "Quant. 1::LsbRemove" );
        process_lsb_remove(d_fx, deterministic_curve, encoder->frame_length, truncated_data);
        BASOP_sub_end();
        BASOP_sub_start( "Quant. 1::QuantSpec" );
        processQuantizeSpec_fx( truncated_data, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, encoder->ll_ari_bits,
                            h_EncSetup->targetBitsAri, &h_EncSetup->mem_specBits, &nBits, encoder->fs_idx, &lastnz,
                            codingdata, &lsbMode, -1, encoder->hrmode,1
                            #ifdef RATE_FLAG_TUNING
                            ,encoder->lossless
                            #endif
#ifdef LL_INCL_HPVC
                            , &(h_EncSetup->hpvcEncCfg), hpvc_L_bits_tcxQ9
#endif
                            );
        BASOP_sub_end();
#ifndef LL_INCL_HPVC
        truncated_data = (Word32*) lc3_scratch_pop( scratch, truncated_data );
#endif
    }
    ELSE
    {
        BASOP_sub_start( "Quant. 1::QuantSpec" );
    #  ifdef ENABLE_HR_MODE
        processQuantizeSpec_fx( d_fx, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, h_EncSetup->targetBitsQuant,
                            h_EncSetup->targetBitsAri, &h_EncSetup->mem_specBits, &nBits, encoder->fs_idx, &lastnz,
                            codingdata, &lsbMode, -1,
                            encoder->hrmode,
                            0
                            #ifdef RATE_FLAG_TUNING
                            ,encoder->lossless
                            #endif
#ifdef LL_INCL_HPVC
                            , &(h_EncSetup->hpvcEncCfg), NULL
#endif
                            );
    #  else
        processQuantizeSpec_fx( d_fx, d_fx_exp, gain, gain_e, q_d_fx16, encoder->yLen, h_EncSetup->targetBitsQuant,
                            h_EncSetup->targetBitsAri, &h_EncSetup->mem_specBits, &nBits, encoder->fs_idx, &lastnz,
                            codingdata, &lsbMode, -1 );
    #  endif
        BASOP_sub_end();
}

#else /* CR14_A_ADD_LOSSLESS_MODE */
    BASOP_sub_start("Quant. 1");
#ifdef ENABLE_HR_MODE
    processQuantizeSpec_fx(d_fx, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, h_EncSetup->targetBitsQuant,
                           h_EncSetup->targetBitsAri, &h_EncSetup->mem_specBits, &nBits, encoder->fs_idx, &lastnz,
                           codingdata, &lsbMode, -1, encoder->hrmode);

#else
    processQuantizeSpec_fx(d_fx, d_fx_exp, gain, gain_e, q_d_fx16, encoder->yLen, h_EncSetup->targetBitsQuant,
                           h_EncSetup->targetBitsAri, &h_EncSetup->mem_specBits, &nBits, encoder->fs_idx, &lastnz,
                           codingdata, &lsbMode, -1);
#endif /* ENABLE_HR_MODE */
#endif /* CR14_A_ADD_LOSSLESS_MODE */
    BASOP_sub_end();
      
//    for (int k = 0; k < encoder->yLen; k++)
//    {
//        printf("q_d_fx24[%d] = %d\n", k, q_d_fx24[k]);
//    }
//      
//    printf("quantizedGainMin = %d\n", quantizedGainMin); 
//    printf("h_EncSetup->quantizedGainOff = %d\n", h_EncSetup->quantizedGainOff); 
//    printf("h_EncSetup->targetBitsQuant = %d\n", h_EncSetup->targetBitsQuant); 
//    printf("h_EncSetup->mem_specBits = %d\n", h_EncSetup->mem_specBits); 
//    printf("encoder->hrmode = %d\n", encoder->hrmode); 

#  ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag == 0)
    {
        processAdjustGlobalGain_fx( &quantizedGain, quantizedGainMin, h_EncSetup->quantizedGainOff, &gain32, &gain_e,
                                    h_EncSetup->targetBitsQuant, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx, encoder->hrmode, encoder->frame_dms, encoder->lossless );
    }
    ELSE 
    {
        IF( h_EncSetup->mem_specBits > encoder->ll_ari_bits && encoder->ll_cbr)
        {
            IF( encoder->fs_idx == 6 && bits_per_sample == 24 )
            {
#ifdef CR15_C_VARIOUS_LOSSLESS_FIXES
                processAdjustGlobalGain_fx( &quantizedGain, quantizedGainMin, encoder->quantizedGainOff_ll, &gain32, &gain_e,
                                            L_min(ll_ari_bits_lb, h_EncSetup->targetBitsAri), h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx, encoder->hrmode, encoder->frame_dms, encoder->lossless );
#else
                processAdjustGlobalGain_fx( &quantizedGain, quantizedGainMin, encoder->quantizedGainOff_ll, &gain32, &gain_e,
                                            ll_ari_bits_lb, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx, encoder->hrmode, encoder->frame_dms, encoder->lossless );
#endif
            }
            ELSE
            {
                processAdjustGlobalGain_fx( &quantizedGain, quantizedGainMin, encoder->quantizedGainOff_ll, &gain32, &gain_e,
                                            encoder->ll_ari_bits, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx, encoder->hrmode, encoder->frame_dms, encoder->lossless );
        
            }    
        }                
        ELSE
        {
            gainChange = 0; move16();
        }
    }
    gain = round_fx( gain32 );
#else /* CR14_A_ADD_LOSSLESS_MODE */
    BASOP_sub_start("Adj. Global Gain");
#ifdef ENABLE_HR_MODE
    //gain32 = L_shl_pos((Word32)gain, 16);
    processAdjustGlobalGain_fx(&quantizedGain, quantizedGainMin, h_EncSetup->quantizedGainOff, &gain32, &gain_e,
                               h_EncSetup->targetBitsQuant, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx
                               , encoder->hrmode, encoder->frame_dms
                               );
    gain = round_fx(gain32);
#else
    processAdjustGlobalGain_fx(&quantizedGain, quantizedGainMin, h_EncSetup->quantizedGainOff, &gain, &gain_e,
                               h_EncSetup->targetBitsQuant, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx);
#endif /* ENABLE_HR_MODE */
#endif /* CR14_A_ADD_LOSSLESS_MODE */
    BASOP_sub_end();


    BASOP_sub_start("Quant. 2");
    IF (sub(gainChange, 1) == 0)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        IF (ll_adap_flag)
        {
            Word16 mode;
#ifdef LL_INCL_HPVC
            /* Reuse outer-scope buffer so contents persist for HPVC analysis. */
            Word32 *truncated_data = hpvc_truncated_data;
#else
            Word32* truncated_data = (Word32*) lc3_scratch_push( scratch, sizeof( *truncated_data ) * encoder->frame_length );
#endif
            BASOP_sub_start( "Quant. 2::DetCurve" );
            IF( encoder->fs_idx == 6 && bits_per_sample == 24 )
            {
                process_deterministic_curve( deterministic_curve, gain32, gain_e, encoder->bands_offset[encoder->low_band_limit], tmp_int_scf_fx, tmp_int_scf_fx_exp_requantized, encoder->bands_offset, encoder->low_band_limit, scratch );
                lsb_mean_split(d_fx, encoder->low_band_limit, encoder->bands_number, encoder->bands_offset, deterministic_curve[ encoder->bands_offset[encoder->low_band_limit]-1 ],
                            ll_deltaCodesBits,
                            deterministic_curve );
            }
            ELSE
            {
                process_deterministic_curve( deterministic_curve, gain32, gain_e, encoder->frame_length, tmp_int_scf_fx, tmp_int_scf_fx_exp_requantized, encoder->bands_offset, encoder->bands_number, scratch );
            }
            BASOP_sub_end();
            BASOP_sub_start( "Quant. 2::LsbRemove" );
            process_lsb_remove(d_fx, deterministic_curve, encoder->frame_length, truncated_data);
            BASOP_sub_end();

            IF( nBits > h_EncSetup->targetBitsQuant )
            {
                mode = 0;
            }
            ELSE
            {
                mode = -1;
            }
            BASOP_sub_start( "Quant. 2::QuantSpec" );
            processQuantizeSpec_fx( truncated_data, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, encoder->ll_ari_bits,
                h_EncSetup->targetBitsAri, &nBits, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
                mode, encoder->hrmode,1
                #ifdef RATE_FLAG_TUNING
                ,encoder->lossless
                #endif
#ifdef LL_INCL_HPVC
                , &(h_EncSetup->hpvcEncCfg), hpvc_L_bits_tcxQ9
#endif
                );
            BASOP_sub_end();
#ifndef LL_INCL_HPVC
            truncated_data = (Word32*) lc3_scratch_pop( scratch, truncated_data );
#endif
        }
        ELSE
        {
            BASOP_sub_start( "Quant. 2::QuantSpec" );
            processQuantizeSpec_fx( d_fx, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, h_EncSetup->targetBitsQuant,
            h_EncSetup->targetBitsAri, NULL, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
            0, encoder->hrmode,0
            #ifdef RATE_FLAG_TUNING
            ,encoder->lossless
            #endif
#ifdef LL_INCL_HPVC
            , &(h_EncSetup->hpvcEncCfg), NULL
#endif
            );
            BASOP_sub_end();
        }
      
#else /* CR14_A_ADD_LOSSLESS_MODE */
      
#ifdef ENABLE_HR_MODE
        processQuantizeSpec_fx(d_fx, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, h_EncSetup->targetBitsQuant,
                               h_EncSetup->targetBitsAri, NULL, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
                               0, encoder->hrmode);
#else
        processQuantizeSpec_fx(d_fx, d_fx_exp, gain, gain_e, q_d_fx16, encoder->yLen, h_EncSetup->targetBitsQuant,
                               h_EncSetup->targetBitsAri, NULL, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
                               0);
#endif /* ENABLE_HR_MODE */
#endif /* CR14_A_ADD_LOSSLESS_MODE */
    }
    BASOP_sub_end();

#ifdef LL_INCL_HPVC
    h_EncSetup->hpvcEncCfg.mode = -1;
    move16();
#ifdef LL_HPVC_FORCE_LEGACY_TCX_ENC
    IF (0)
#else
    IF (h_EncSetup->hpvcEncCfg.active_flag != 0 && encoder->lossless != 0 && ll_adap_flag != 0)
#endif
    {
        Word16 Nqm_ana;
        Word16 L1_signal[MAX_LEN / LL_HPVC_N_SIGNAL];   /* L1 norms per Nsignal block */

        h_EncSetup->hpvcEncCfg.Tx_dec = &(Tx_dec[0]);
        h_EncSetup->hpvcEncCfg.Tx_splitRule = &(Tx_splitRule[0]);

        h_EncSetup->hpvcEncCfg.mode = 0;
        move16();

        Nqm_ana = shr_pos(sub(lastnz, h_EncSetup->hpvcEncCfg.startCoef), N_SIGNAL_LOG);

        if (Nqm_ana <= 0)
        {
            /* all-zero tail, disable further HPVC analysis and encoding */
            h_EncSetup->hpvcEncCfg.mode = -1;
            move16();
        }

        IF (h_EncSetup->hpvcEncCfg.mode >= 0)
        {
            /* Initial OL HPVC analysis: decides mode (0 or 1) and startCoef. */
            Word16 local_mode = hpvc_analyze_start_fx(hpvc_truncated_data, lastnz, L1_signal, &h_EncSetup->hpvcEncCfg);
            UNUSED(local_mode);
        }

        IF (h_EncSetup->hpvcEncCfg.mode >= 0)
        {
            Word32 L_bitsLongQ9[MAX_LEN / LL_HPVC_N_SIGNAL];
            Word16 nb_hpvc_trees;
            Word32 *Qm_ptr;
            Word16 *L1_sig_ptr;
            Word16 startBlockOffset, allTCX_flag, coef0;
            UNUSED(nb_hpvc_trees);

            Qm_ptr = &(hpvc_truncated_data[h_EncSetup->hpvcEncCfg.startCoef]);
            Nqm_ana = s_max(0, sub(lastnz, h_EncSetup->hpvcEncCfg.startCoef));
            Nqm_ana = shr_pos(Nqm_ana, N_SIGNAL_LOG);

#ifdef LL_HPVC_ALIGN_STARTCOEFF_TO_LASTNZ
            coef0 = h_EncSetup->hpvcEncCfg.startCoefList[0];  /* lastnz-adjusted */
#else
            coef0 = h_EncSetup->hpvcEncCfg.startCoefListNom[0];
            move16();
#endif
            startBlockOffset = shr_pos(sub(h_EncSetup->hpvcEncCfg.startCoef, coef0), N_SIGNAL_LOG);
            L1_sig_ptr = &(L1_signal[startBlockOffset]);

            /* L1-norm based per-block pre-decision TCX vs HPVC. */
            hpvc_segment_tcx_or_hpvc(L1_sig_ptr, Nqm_ana, Tx_dec);

            nb_hpvc_trees = hpvc_collect_hpvc_bitrates(
                L1_sig_ptr,
                Qm_ptr,
                Tx_dec,
                Tx_splitRule,
                Nqm_ana,
                &(L_bitsLongQ9[0]));

            {
                Word32 L_totalBits_fulltcxQ9;
                Word32 L_totalBits_mixed_onlyQ9;
                Word32 *L_bits_tcxQ9_ptr = &(hpvc_L_bits_tcxQ9[startBlockOffset]);

                allTCX_flag = hpvc_segment_accept_reject(L_bits_tcxQ9_ptr, L_bitsLongQ9,
                    L1_sig_ptr, Nqm_ana, Tx_dec,
                    &L_totalBits_fulltcxQ9, &L_totalBits_mixed_onlyQ9
#ifdef HPVC_APPLY_TREE_LIMIT
                    , h_EncSetup->hpvcEncCfg.nomTreeLim
#endif
#ifdef HPVC_MAXTREE_LIMIT
                    , h_EncSetup->hpvcEncCfg.maxTreeLim
#endif
                );

                if (allTCX_flag != 0)
                {
#ifdef LL_HPVC_FORCE_NP_TCX_ENC
                    h_EncSetup->hpvcEncCfg.mode = h_EncSetup->hpvcEncCfg.mode;
#else
                    /* Disable HPVC due to total bit-rate cost. */
                    h_EncSetup->hpvcEncCfg.mode = -1;
                    move16();
#endif
                }
            }

            /* Wire tree config storage from outer-scope buffer to keep it valid through processAriEncoder_fx. */
            {
                Word16 nTrees = 0;
                h_EncSetup->hpvcEncCfg.HpvcTreeEnumCfgPtr = &(hpvc_trees[0]);

                IF (h_EncSetup->hpvcEncCfg.mode >= 0)
                {
                    nTrees = hpvc_enumerate_trees(
                        Qm_ptr,
                        Tx_dec,
                        Tx_splitRule,
                        L1_sig_ptr,
                        Nqm_ana,
                        &(hpvc_trees[0]));
                    UNUSED(nTrees);
                }
            }
        }
    }
#endif

//    for (int k = 0; k < encoder->yLen; k++)
//    {
//        printf("q_d_fx24[%d] = %d\n", k, q_d_fx24[k]);
//    }

#ifdef CR14_A_ADD_LOSSLESS_MODE
        max_resBits_len = MAX_RESBITS_LEN;
        if (encoder->lossless) {
            max_resBits_len = L_add (L_mult0(encoder->frame_length,add(bitsPerSample,1)),1);
        }
#endif

    BASOP_sub_start("Res. Cod.");
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag)
    {
        numResBits = 0;
        basop_memset( resBits, 0, sizeof( *resBits ) * max_resBits_len );
    } ELSE IF( lsbMode == 0 )
#else
    IF (lsbMode == 0)
#endif
    {
        processResidualCoding_fx(d_fx_exp, d_fx,
#ifdef ENABLE_HR_MODE
                                 q_d_fx24,
#else
                                 q_d_fx16,
#endif
#ifdef ENABLE_HR_MODE
                                 gain32,
#else
                                 gain,
#endif
                                 gain_e, encoder->yLen, h_EncSetup->targetBitsQuant, nBits, resBits, &numResBits
#ifdef ENABLE_HR_MODE
                                 , encoder->hrmode
#endif
#if defined (CR9_C_ADD_1p25MS)
                                 , encoder->frame_dms
#endif
        );
    }
    ELSE
    {
        numResBits = 0;
        move16();
    }
    BASOP_sub_end();

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_adap_flag == 1)
    {
        goto skip_nf;
    }
#endif

    BASOP_sub_start("Noise fac");
    IF (h_EncSetup->lfe == 0)
    {
        processNoiseFactor_fx(&fac_ns_idx, d_fx_exp, d_fx,
#ifdef ENABLE_HR_MODE
                              q_d_fx24,
#else
                              q_d_fx16,
#endif
                              gain, gain_e, BW_cutoff_idx, encoder->frame_dms, h_EncSetup->targetBytes, scratch
#ifdef ENABLE_HR_MODE
                              , encoder->hrmode
#endif
                             );
    }
    ELSE
    {
        fac_ns_idx = 7;
        move16();
    }
    BASOP_sub_end();
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
skip_nf:
#endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (encoder->lossless)
    {
        IF (encoder->ll_tns == 0)
        {
            tns_numfilters = 0;
            tns_order[0] = -1;
            tns_order[1] = -1;
        }
        
        IF (encoder->ll_cbr == 0 || ll_adap_flag)
        {
            IF (encoder->ll_tns == 0)
            {
                BW_cutoff_idx = -1;
            }
            
            fac_ns_idx = -1;
        }
    }
#endif
      
//    printf("h_EncSetup->targetBitsAri = %d\n", h_EncSetup->targetBitsAri);
//    printf("h_EncSetup->targetBytes = %d\n", h_EncSetup->targetBytes);
//    printf("encoder->yLen = %d\n", encoder->yLen);
//    printf("encoder->BW_cutoff_bits = %d\n", encoder->BW_cutoff_bits);
//    printf("tns_numfilters = %d\n", tns_numfilters);
//    printf("lsbMode = %d\n", lsbMode);
//    printf("lastnz = %d\n", lastnz);
//    printf("fac_ns_idx = %d\n", fac_ns_idx);
//    printf("quantizedGain = %d\n", quantizedGain);
//      
//    printf("BW_cutoff_idx = %d\n", BW_cutoff_idx);
//    printf("encoder->frame_dms = %d\n", encoder->frame_dms);
//    printf("encoder->lossless = %d\n", encoder->lossless);
//    printf("BW_cutoff_idx = %d\n", BW_cutoff_idx);
//    printf("encoder->ll_cbr = %d\n", encoder->ll_cbr);
//    printf("ll_adap_flag = %d\n", ll_adap_flag);
//    printf("encoder->ll_tns = %d\n", encoder->ll_tns);
//    printf("off_idx = %d\n", off_idx);
//    printf("encoder->ll_offQuant = %d\n", encoder->ll_offQuant);
//    printf("encoder->ll_tns_remove = %d\n", encoder->ll_tns_remove);
//    printf("tns_lsb_num_remove = %d\n", tns_lsb_num_remove);
//    printf("encoder->wavFormat = %d\n", encoder->wavFormat);
//      
//    for (int k = 0; k < tns_numfilters; k++)
//    {
//        printf("tns_order[%d] = %d\n", k, tns_order[k]);
//    }
//      
//    for (int k = 0; k < 7; k++)
//    {
//        printf("L_scf_idx[%d] = %d\n", k, L_scf_idx[k]);
//    }
//
//    for (int k = 0; k < 3; k++)
//    {
//        printf("ltpf_idx[%d] = %d\n", k, ltpf_idx[k]);
//    }
//      
//    for (int k = 0; k < encoder->yLen; k++)
//    {
//        printf("q_d_fx24[%d] = %d\n", k, q_d_fx24[k]);
//    }
//      
//    for (int k = 0; k < numResBits; k++)
//    {
//        printf("resBits[%d] = %d\n", k, resBits[k]);
//    }

    BASOP_sub_start("Entropy cod");
    processEncoderEntropy(bytes, &bp_side, &mask_side, h_EncSetup->targetBitsAri, h_EncSetup->targetBytes,
                          encoder->yLen, encoder->BW_cutoff_bits, tns_numfilters, lsbMode, lastnz, tns_order,
                          fac_ns_idx, quantizedGain, BW_cutoff_idx, ltpf_idx, L_scf_idx, bfi_ext, encoder->fs_idx
#ifdef CR9_C_ADD_1p25MS
#ifdef FIX_TX_RX_STRUCT_STEREO
                          , encoder->frame_dms, &h_EncSetup->Tx_ltpf
#else
                          , encoder->frame_dms, &encoder->Tx_ltpf
#endif
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            ,  encoder->lossless,
                               ll_adap_flag,
                               encoder->b_relative,
                               encoder->ll_tns,
                               off_idx,
                               encoder->ll_offQuant,
                               encoder->ll_tns_remove,
                               tns_lsb_num_remove,
                               ll_deltaCodesBits
                            ,  bits_per_sample
                            , 0
                            , (UWord16) scaleSignal
#ifdef LL_INCL_HPVC
                            , &(h_EncSetup->hpvcEncCfg)
#endif
#endif
                         );
    BASOP_sub_end();
      
//    for (int k = 0; k < h_EncSetup->targetBytes; k++)
//    {
//        printf("bytes[%d] = %d\n", k, bytes[k]);
//    }
//      
//    printf("lsbMode = %d\n", lsbMode);
//    printf("h_EncSetup->enable_lpc_weighting = %d\n", h_EncSetup->enable_lpc_weighting);
//    printf("ll_side[0] = %d\n", ll_side[0]);
//    printf("ll_side[1] = %d\n", ll_side[1]);
//    printf("max_resBits_len = %d\n", max_resBits_len);
//    printf("tns_lsb_num_remove = %d\n", tns_lsb_num_remove);
//    printf("encoder->ll_cbr = %d\n", encoder->ll_cbr);
//    printf("ll_adap_flag = %d\n", ll_adap_flag);
//      
//    for (int k = 0; k < encoder->yLen; k++)
//    {
//        printf("deterministic_curve[%d] = %d\n", k, deterministic_curve[k]);
//    }
//      
//    for (int k = 0; k < encoder->yLen; k++)
//    {
//        printf("d_fx_orig[%d] = %d\n", k, d_fx_orig[k]);
//    }
      
    if ( scratch->max_scratch_calculation_only )
    {
        lastnz = encoder->yLen;
    }

    BASOP_sub_start("Ari cod");
    processAriEncoder_fx(bytes, bp_side, mask_side, h_EncSetup->targetBitsAri,
#ifdef ENABLE_HR_MODE
                         q_d_fx24,
#else
                         q_d_fx16,
#endif
                         tns_order, tns_numfilters, indexes,
                         lastnz,
                         codingdata, resBits, numResBits, lsbMode,
                         h_EncSetup->enable_lpc_weighting,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                          ll_side, deterministic_curve, d_fx, max_resBits_len
                          , tns_lsb_num_remove, encoder->ll_cbr, d_fx_orig, tns_lsb_remove
                          , ll_adap_flag
#ifdef LL_INCL_HPVC
                          , &(h_EncSetup->hpvcEncCfg)
#endif
                          , encoder->b_relative, encoder->bands_offset, encoder->bands_number, encoder->yLen,
#endif
                         scratch);
    BASOP_sub_end();
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
IF ( encoder->padding && ll_adap_flag && encoder->ll_cbr )
{
    IF( !scratch->max_scratch_calculation_only && ll_side[1] > ll_side[0] )
    {
        // number of bytes at buffer end (side info)
        Word16 lenBytesEnd = (h_EncSetup->targetBitsAri >> 3) - ll_side[1];
        UWord8* tmp_buf = (UWord8*) lc3_scratch_push( scratch, sizeof( *tmp_buf ) * lenBytesEnd );
        // total number of bytes without gap between fw and bw buffer parts
        encoder->ll_totalBytes = ll_side[0] + 1 + lenBytesEnd;
        // copy end (=bw) part to tmp_buf
        basop_memcpy( tmp_buf, bytes + ll_side[1], lenBytesEnd );
        // append it to first (=fw) part of buff
        basop_memcpy( bytes + ll_side[0] + 1, tmp_buf, lenBytesEnd );
        // fill with 0x7 from index ll_totalBytes to targetBitsAri/8 (bits=>bytes)
        basop_memset( bytes + encoder->ll_totalBytes, 0x7, (h_EncSetup->targetBitsAri >> 3) - encoder->ll_totalBytes );
        // total number of bytes without gap between fw and bw buffer parts but with end padding
        encoder->ll_totalBytes = (h_EncSetup->targetBitsAri >> 3);
        // Padding => encoder->ll_totalBytes is GROSS number of data bytes = full frame length
        tmp_buf = (UWord8*) lc3_scratch_pop( scratch, tmp_buf );
    }
    ELSE IF ( encoder->lossless )
    {
        encoder->ll_totalBytes = (h_EncSetup->targetBitsAri >> 3);
    }
}
ELSE
{
    IF( !scratch->max_scratch_calculation_only && ll_adap_flag && ll_side[1] > ll_side[0] )
    {
        // number of bytes at buffer end (side info)
        Word16 lenBytesEnd = (h_EncSetup->targetBitsAri >> 3) - ll_side[1];
        UWord8* tmp_buf = (UWord8*) lc3_scratch_push( scratch, sizeof( *tmp_buf ) * lenBytesEnd );
         // total number of bytes without gap between fw and bw buffer parts
        encoder->ll_totalBytes = ll_side[0] + 1 + lenBytesEnd;
        // copy end (=bw) part to tmp_buf
        basop_memcpy( tmp_buf, bytes + ll_side[1], lenBytesEnd );
        // append it to first (=fw) part of buff
        basop_memcpy( bytes + ll_side[0] + 1, tmp_buf, lenBytesEnd );
        // fill with zeros from index ll_totalBytes to targetBitsAri/8 (bits=>bytes)
        basop_memset( bytes + encoder->ll_totalBytes, 0, (h_EncSetup->targetBitsAri >> 3) - encoder->ll_totalBytes );
        /*basop_memcpy( bitBuffState_bw.ptr + ll_side[0] + 1, bitBuffState_bw.ptr + ll_side[1], lenBytesEnd );
        basop_memset( bitBuffState_bw.ptr + encoder->ll_totalBytes, 0, (h_EncSetup->targetBitsAri >> 3) - encoder->ll_totalBytes );*/

        // No padding => encoder->ll_totalBytes is NET number of data bytes
        if( encoder->ll_carryOver )
            encoder->ll_carryOverBytes = (h_EncSetup->targetBitsAri >> 3) - encoder->ll_totalBytes;
        tmp_buf = (UWord8*) lc3_scratch_pop( scratch, tmp_buf );
    }
    ELSE IF ( encoder->lossless )
    {
        encoder->ll_totalBytes = (h_EncSetup->targetBitsAri >> 3);
    }
}
#endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
    /* End lossless/lossy counting here. */
    /* In fact there is some more below but this does not yield a relevant contribution */
    BASOP_sub_sub_end(); // BASOP_sub_sub_start( "Enc(lossless)"/lossy );

    IF( encoder->lossless )
    {
        IF( bits_per_sample == 16 )
        {
            fallback_bit_planes = 17;
        }
        ELSE IF( bits_per_sample == 24 )
        {
            fallback_bit_planes = 25;
        }
#ifdef CR15_A_LOSSLESS_1p25MS
        Word16 fallback_threshold_bytes = (fallback_bit_planes * encoder->frame_length) / 8 + 1;
        fallback_num_bytes = (fallback_bit_planes * encoder->frame_length + 7) / 8 + 1;

        #ifdef CR14_A_ADD_LOSSLESS_MODE
        fallback_num_bytes += 1;
        #endif 

        IF( encoder->ll_totalBytes >= fallback_threshold_bytes )
#else
        fallback_num_bytes = (fallback_bit_planes * encoder->frame_length) / 8 + 1;

        IF( encoder->ll_totalBytes >= fallback_num_bytes )
#endif
        {
            ll_adap_flag = 1;
            BASOP_sub_start( "EncEntropyFallback" );
#ifdef CR15_A_LOSSLESS_1p25MS
            h_EncSetup->Tx_ltpf = 0;
#endif
            processEncoderEntropy( bytes, &bp_side, &mask_side, fallback_num_bytes*8, fallback_num_bytes,
                encoder->yLen, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, bfi_ext, encoder->fs_idx,
                encoder->frame_dms, NULL,
                encoder->lossless, ll_adap_flag, encoder->b_relative, 0, 0, 0, 0, 0 ,
                NULL,
                bits_per_sample,
                1
#ifdef CR14_A_ADD_LOSSLESS_MODE
                , (UWord16) scaleSignal
#ifdef LL_INCL_HPVC
                , &(h_EncSetup->hpvcEncCfg)
#endif
#endif
            );

            fallback_encoder( bytes, &bp_side, &mask_side, d_tda_fx, encoder->frame_length, fallback_bit_planes );

            encoder->ll_totalBytes = fallback_num_bytes;
            BASOP_sub_end(); // BASOP_sub_start( "EncEntropyFallback" );
        }
    }
#endif

    BASOP_sub_start("Reorder Bitstream Enc");
    test();
    IF (encoder->combined_channel_coding == 0 && h_EncSetup->n_pc > 0)
    {
        BASOP_sub_start("Reorder Ari dec");

#ifdef ENABLE_HR_MODE
#    ifdef ENABLE_HR_MODE
        Word32* xbuf = (Word32*) lc3_scratch_push( scratch, encoder->frame_length * sizeof( *xbuf ) );
#    else
        Word16* xbuf = (Word16*) lc3_scratch_push( scratch, encoder->frame_length * sizeof( *xbuf ) );
#    endif

        processAriDecoder_fx(bytes, &bp_side, &mask_side, h_EncSetup->total_bits, encoder->yLen, encoder->fs_idx,
                             h_EncSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                             fac_ns_idx, quantizedGain, encoder->frame_dms, h_EncSetup->n_pc, 0,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                              shr_pos( h_EncSetup->total_bits, 3 ), 1, &gain, &b_left, (Word32*)&gain,
#else
                              shr_pos( h_EncSetup->total_bits, 3 ), 1, &gain, &b_left, &gain,
#endif
                             xbuf,
                             &gain, resBits, indexes, &gain,
                             scratch
                             , encoder->hrmode
#    ifdef CR14_A_ADD_LOSSLESS_MODE
                              , deterministic_curve, ll_adap_flag
                              , encoder->lossless
#    endif
#ifdef LL_INCL_HPVC
                              , NULL
#endif
        );
  
#    ifdef ENABLE_HR_MODE
        xbuf = (Word32*) lc3_scratch_pop( scratch, xbuf );
#    else
        xbuf = (Word16*) lc3_scratch_pop( scratch, xbuf );
#    endif
#else

        processAriDecoder_fx(bytes, &bp_side, &mask_side, h_EncSetup->total_bits, encoder->yLen, encoder->fs_idx,
                             h_EncSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                             fac_ns_idx, quantizedGain, encoder->frame_dms, h_EncSetup->n_pc, 0,
                             shr_pos(h_EncSetup->total_bits, 3), 1, &gain, &b_left, &gain,
                             codingdata,
                             &gain, resBits, indexes, &gain,
                             scratch
#ifdef LL_INCL_HPVC
                             , NULL
#endif
        );
#endif

        BASOP_sub_end(); /* Ari dec */
        processReorderBitstream_fx(bytes, h_EncSetup->n_pccw, h_EncSetup->n_pc, b_left, scratch);
    }
    BASOP_sub_end();

    /* end q_d_fx16 */
    BASOP_sub_end();

    int_scf_fx_exp = (Word16*) lc3_scratch_pop( scratch, int_scf_fx_exp );
#  ifdef ENABLE_HR_MODE
    int_scf_fx = (Word32*) lc3_scratch_pop( scratch, int_scf_fx );
#  else
    int_scf_fx = (Word16*) lc3_scratch_pop( scratch, int_scf_fx );
#  endif
    scf = (Word16*) lc3_scratch_pop( scratch, scf );
#  ifdef ENABLE_HR_MODE
    scf_q = (Word32*) lc3_scratch_pop( scratch, scf_q );
#  else
    scf_q = (Word16*) lc3_scratch_pop( scratch, scf_q );
#  endif
    s_12k8 = (Word16*) lc3_scratch_pop( scratch, s_12k8 );
    resBits = (UWord8*) lc3_scratch_pop( scratch, resBits );
#  ifdef ENABLE_HR_MODE
    s_in_scaled = (Word32*) lc3_scratch_pop( scratch, s_in_scaled );
#  else
    s_in_scaled = (Word16*) lc3_scratch_pop( scratch, s_in_scaled );
#  endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    tns_lsb_remove      = (UWord8*) lc3_scratch_pop( scratch, tns_lsb_remove );
    d_fx_orig           = (Word32*) lc3_scratch_pop( scratch, d_fx_orig );
    deterministic_curve = (UWord8*) lc3_scratch_pop( scratch, deterministic_curve);
#endif
    codingdata = (Word16*) lc3_scratch_pop( scratch, codingdata );
#  ifdef ENABLE_HR_MODE
    q_d_fx24 = (Word32*) lc3_scratch_pop( scratch, q_d_fx24 );
#  else
    q_d_fx16 = (Word16*) lc3_scratch_pop( scratch, q_d_fx16 );
#  endif


#ifdef CR14_A_ADD_LOSSLESS_MODE
    d_tda_fx = (Word32*) lc3_scratch_pop( scratch, d_tda_fx );
#endif 

    indexes = (Word16*) lc3_scratch_pop( scratch, indexes );
    L_scf_idx = (Word32*) lc3_scratch_pop( scratch, L_scf_idx );
    d_fx = (Word32*) lc3_scratch_pop( scratch, d_fx );

    Dyn_Mem_Deluxe_Out();
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
void Enc_LC3PLUS(LC3PLUS_Enc *encoder, void **input, int bits_per_sample, UWord8 *output, lc3_scratch_t scratch, Word16 bfi_ext, Word32* nBytes)
#else
int Enc_LC3PLUS(LC3PLUS_Enc *encoder, void **input, int bits_per_sample, UWord8 *output, lc3_scratch_t scratch, Word16 bfi_ext)
#endif
{
    int ch = 0, output_size = 0;
    int input_size = 0;
    int totalBytes = (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
    int output_size2;

    UWord8 *lc3buf = output;
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    encoder->ll_carryOverBytes = 0;
#endif

    for (ch = 0; ch < encoder->channels; ch++)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if( ch > 0 )
        {
            update_enc_payload_sizes(encoder, ch,  encoder->totalBytes, encoder->ll_carryOverBytes);
        }
#endif
        Enc_LC3PLUS_Channel(encoder, ch, bits_per_sample, input[ch], lc3buf, scratch, bfi_ext);
        if (encoder->epmode && encoder->combined_channel_coding == 0)
        {
            output_size2 = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));
            BASOP_sub_start("fec_enc");

            fec_encoder(encoder->epmode, encoder->epmr, lc3buf, encoder->channel_setup[ch]->targetBytes, output_size2,
                        encoder->channel_setup[ch]->n_pccw, scratch);

            BASOP_sub_end();

            lc3buf += output_size2;
            output_size += output_size2;
#ifdef CR14_A_ADD_LOSSLESS_MODE
            nBytes[ch] = output_size2;
#endif
        }
        else
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            if ( encoder->lossless )
            {
                lc3buf += encoder->ll_totalBytes;
                nBytes[ch] = encoder->ll_totalBytes;
            }
            else
#endif
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
              if ( encoder->lossless )
              {
                  lc3buf += encoder->ll_totalBytes;
                  nBytes[ch] = encoder->ll_totalBytes;
              }
              else
#endif
              {
                  lc3buf += encoder->channel_setup[ch]->targetBytes;
                  output_size += encoder->channel_setup[ch]->targetBytes;
#ifdef CR14_A_ADD_LOSSLESS_MODE
                  nBytes[ch] = encoder->channel_setup[ch]->targetBytes;
#endif
              }
            }
        }
      
        if ( scratch->max_scratch_calculation_only )
        {
            break; /* Only one channel for max-scratch calculation */
        }
    }

    if (encoder->epmode > 0 && encoder->combined_channel_coding)
    {
        input_size  = output_size;
        output_size = (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
        BASOP_sub_start("fec_enc");
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
        for (int ch = 0; ch < encoder->channels; ch++)
        {
            nBytes[ch] = L_shr_pos_pos(output_size, encoder->channels-1);
          
            if ( scratch->max_scratch_calculation_only )
            {
                break; /* Only one channel for max-scratch calculation */
            }
        }
#endif

        fec_encoder(encoder->epmode, encoder->epmr, output, input_size, output_size, encoder->channel_setup[0]->n_pccw,
                    scratch);

        BASOP_sub_end();
    }

#ifndef CR14_A_ADD_LOSSLESS_MODE
    return output_size;
#endif
}
