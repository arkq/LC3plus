/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "util.h"

static int Dec_LC3PLUS_Channel(LC3PLUS_Dec *decoder, int channel, int bits_per_sample, UWord8 *bs_in, void *s_out, Word16 bfi, lc3_scratch_t scratch)
{
    Word16 bfi_ext = bfi;
    Word16 scale;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 fill_bits = 0;
#else
    Word16 fill_bits = 0;
    Word32 offset = 0;
#endif
    Word16 nf_seed = 0, gg_idx = 0, fac_ns_idx = 0, q_fx_exp = 0;
    Word16 bp_side = 0, mask_side = 0;
    Word16 tns_numfilters = 0, lsbMode = 0, lastnz = 0, BW_cutoff_idx = 0, BW_cutoff_idx_nf = 0;
    Word16 zero_frame = 0;
#ifdef ENABLE_RFRAME
    Word16 rframe = 0;
#endif
    Word16 ltpf_idx[3] = {0};
    Word16 spec_inv_idx = 0;
    Counter i = 0;

    /* Buffers */
    Word16 *int_scf_fx_exp;
    UWord8 *resBitBuf;

    Word16 *  int_scf_fx;
#ifdef ENABLE_HR_MODE
    Word32 *sqQdec;
#else
    Word16 *  sqQdec;
#endif
    Word16 *  x_fx, *indexes;
    Word16    scf_q[M];
    Word32 *  L_scf_idx;
    Word32 *  q_d_fx;
    DecSetup *h_DecSetup = decoder->channel_setup[channel];

#ifdef CR9_C_ADD_1p25MS_LRSNS
    Word16 pitch_rx_fx;
    Word16 ltpf_rx_fx;
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
    Word32 scf_q_ip[M];
#   ifdef ENABLE_HR_MODE
    //Counter i;
    Word32* x_fx_ip;
    Word32 *int_scf_fx_ip;
#   endif
#else
#   ifdef ENABLE_HR_MODE
    Word32 *x_fx_ip;
    Word32 *int_scf_fx_ip;
    Word32 scf_q_ip[M];
#   endif
#endif /*  CR9_C_ADD_1p25MS_LRSNS  */

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Counter i;
        Word16  scale;
        Word32  offset;
        Word16  fill_bits;
        Word16  nf_seed, gg_idx, fac_ns_idx, q_fx_exp;
        Word16  bp_side, mask_side;
        Word16  tns_numfilters, lsbMode, lastnz, BW_cutoff_idx, BW_cutoff_idx_nf;
        Word16  zero_frame;
        Word16  ltpf_idx[3];
#ifdef ENABLE_RFRAME
        Word16 rframe;
#endif
        Word16 spec_inv_idx;

        /* Buffers */
        Word16 *int_scf_fx_exp;
        UWord8 *resBitBuf;
#ifdef ENABLE_HR_MODE
        Word32 *sqQdec;
#else
        Word16 *sqQdec;
#endif
        Word16 *int_scf_fx, *x_fx, *indexes;
        Word32 *L_scf_idx;
        Word32 *q_d_fx;
        Word16 scf_q[M];
#ifdef ENABLE_HR_MODE
        Word32 scf_q_ip[M];
#endif
    };
    Dyn_Mem_In("Dec_LC3_Channel", sizeof(struct _dynmem));
#endif

    Word16 tns_order[TNS_NUMFILTERS_MAX] = {0};

#ifdef DISABLE_PLC
    UNUSED(decoder->plcMeth);
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32* q_d_res;
    UWord8* deterministic_curve;
    Word16 ll_adap_flag = 0, off_idx = 0, tns_lsb_num_remove = 0;
    Word16 b_relative = 0;
    Word16 fallback = 0;
    Word16 scaleSignal = 0;
    Word16 fallback_bit_planes = 0;
    Word16 input_tda = 0;
    Word32* q_res_r;
    Word32* q_res;
    Word8 ll_dec_rounding = 1;
    Word8 ll_adap_flag_2nd = 0;
    UNUSED(ll_adap_flag_2nd);
    Word32* int_scf_ll;
    Word16* int_scf_exp_ll;
    UWord8* residualDataLossless;
    Word16 res_bit_pos = 0;
    Word16 index_b = -1;
    Word16 index_x = -1;
    UWord8* tns_lsb_add;

    UWord8 deltaCodedBits[HIGH_BANDS_NUMBER];
    basop_memset(deltaCodedBits, 0, HIGH_BANDS_NUMBER);
    UWord8* eff_det_curve;

#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    UWord8* resBitBufLossless;
#else
    UWord8 resBitBufLossless[48000] = {0};
#endif

    q_d_fx = (Word32*) lc3_scratch_push( scratch, sizeof( *q_d_fx ) * decoder->frame_length );

#  ifdef ENABLE_HR_MODE
    /* allocate memory for residual bits */
    IF ( decoder->hrmode )
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        IF ( decoder->lossless )
        {
            Word16 max_resBits_len = L_shr(L_add (L_mult0(decoder->frame_length,add(bits_per_sample,1)),1),3); /* ((bps+1) * N)/8 */
            resBitBuf = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBuf ) * max_resBits_len );
            residualDataLossless = (UWord8*) lc3_scratch_push( scratch, sizeof( *residualDataLossless ) * L_shl(max_resBits_len , 3) );
            resBitBufLossless = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );

            basop_memset( resBitBuf, 0, sizeof( *resBitBuf ) * max_resBits_len );
            basop_memset( residualDataLossless, 0, sizeof( *residualDataLossless ) * L_shl(max_resBits_len , 3) );
            basop_memset( resBitBufLossless, 0, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );
        } ELSE
#endif
        {
            Word16 max_resBits_len = L_shr(L_add (L_mult0(decoder->frame_length,add(bits_per_sample,1)),1),3); /* ((bps+1) * N)/8 */
            resBitBuf = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBuf ) * MAX_RESBITS_LEN );
#ifdef CR14_A_ADD_LOSSLESS_MODE
            resBitBufLossless = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );
#endif
            basop_memset( resBitBuf, 0, sizeof( *resBitBuf ) * MAX_RESBITS_LEN );
#ifdef CR14_A_ADD_LOSSLESS_MODE
            basop_memset( resBitBufLossless, 0, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );
#endif
        }
    }
    ELSE
#  endif
    {
        Word16 maxResBits = decoder->frame_length;
#  ifdef ENABLE_12p5_DMS_MODE
        IF( decoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS ) { maxResBits = i_mult( maxResBits, 3 ); }
#  endif

        resBitBuf = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBuf ) * MAX_RESBITS_LEN );
        basop_memset( resBitBuf, 0, sizeof( *resBitBuf ) * MAX_RESBITS_LEN );

#ifdef CR14_A_ADD_LOSSLESS_MODE
        Word16 max_resBits_len = L_shr(L_add (L_mult0(decoder->frame_length,add(bits_per_sample,1)),1),3); /* ((bps+1) * N)/8 */
        resBitBufLossless = (UWord8*) lc3_scratch_push( scratch, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );
        basop_memset( resBitBufLossless, 0, sizeof( *resBitBufLossless ) * L_shl(max_resBits_len , 3) );
#endif
    }

    indexes = (Word16*) lc3_scratch_push( scratch, sizeof( *indexes ) * TNS_NUMFILTERS_MAX * MAXLAG );
    memset( indexes, 0, sizeof( *indexes ) * TNS_NUMFILTERS_MAX * MAXLAG );
    L_scf_idx = (Word32*) lc3_scratch_push( scratch, sizeof( *L_scf_idx ) * SCF_MAX_PARAM );
#  ifdef ENABLE_HR_MODE
    sqQdec = (Word32*) lc3_scratch_push( scratch, sizeof( *sqQdec ) * decoder->frame_length );
#  else
    sqQdec = (Word16*) lc3_scratch_push( scratch, sizeof( *sqQdec ) * decoder->frame_length );
#  endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    deterministic_curve = (UWord8*) lc3_scratch_push( scratch, sizeof( *deterministic_curve ) * decoder->frame_length );
#endif
    int_scf_fx_exp = (Word16*) lc3_scratch_push( scratch, sizeof( *int_scf_fx_exp ) * MAX_BANDS_NUMBER );
    int_scf_fx = (Word16*) lc3_scratch_push( scratch, sizeof( *int_scf_fx ) * MAX_BANDS_NUMBER );
    x_fx = (Word16*) lc3_scratch_push( scratch, sizeof( *x_fx ) * ( decoder->frame_length + decoder->stDec_ola_mem_fx_len ) );
#  ifdef ENABLE_HR_MODE
    x_fx_ip = (Word32*) lc3_scratch_push( scratch, sizeof( *x_fx_ip ) * ( decoder->frame_length + decoder->stDec_ola_mem_fx_len ) );
    int_scf_fx_ip = (Word32*) lc3_scratch_push( scratch, sizeof( *int_scf_fx_ip ) * MAX_BANDS_NUMBER );



#endif

#ifdef DISABLE_PLC
    memset(q_d_fx, 0, decoder->frame_length * sizeof(*q_d_fx));
#endif

    BASOP_sub_start("Decoder");

#ifdef ENABLE_RFRAME
    IF (sub(bfi, 3) == 0)
    {
        bfi = 2;
        move16();
        rframe = 1;
        move16();
    }
#endif

    if (bfi != 1)
    {
        BASOP_sub_sub_start("Dec(bfi=0)");
    }
    else
    {
        BASOP_sub_sub_start("Dec(bfi=1)");
    }

    BASOP_sub_start("Entropy dec");
#ifdef NEW_SIGNALLING_SCHEME_1p25
    h_DecSetup->ltpfinfo_frame_cntr_fx = add_sat(h_DecSetup->ltpfinfo_frame_cntr_fx, 1);
    /*ltpfinfo_frame_cntr_fx increased always,  also for bfi=1  */  /* set or reset inside dec_entropy_fx() */
#endif

    if ( scratch->max_scratch_calculation_only )
    {
        bfi = 1;
    }

    IF (sub(bfi, 1) != 0)
    {

//    printf("h_DecSetup->total_bits = %d\n", h_DecSetup->total_bits);
//    printf("decoder->BW_cutoff_bits = %d\n", decoder->BW_cutoff_bits);
//    printf("decoder->lossless = %d\n", decoder->lossless);
//    printf("decoder->ll_tns = %d\n", decoder->ll_tns);
//    printf("h_DecSetup->ll_offQuant = %d\n", h_DecSetup->ll_offQuant);
//    printf("decoder->ll_tns_remove = %d\n", decoder->ll_tns_remove);

        processDecoderEntropy_fx(bs_in, &bp_side, &mask_side, h_DecSetup->total_bits, decoder->yLen, decoder->fs_idx,
                                 decoder->BW_cutoff_bits, &tns_numfilters, &lsbMode, &lastnz, &bfi, tns_order,
                                 &fac_ns_idx, &gg_idx, &BW_cutoff_idx, ltpf_idx, L_scf_idx, decoder->frame_dms
#ifdef CR9_C_ADD_1p25MS
#   ifdef FIX_TX_RX_STRUCT_STEREO
                                 , h_DecSetup->ltpf_rx_status, &h_DecSetup->ltpf_mem_continuation
#   else
                                 , decoder->ltpf_rx_status, &decoder->ltpf_mem_continuation
#    endif
#    ifdef NEW_SIGNALLING_SCHEME_1p25
                                 ,
                                 &h_DecSetup->ltpfinfo_frame_cntr_fx
#    endif
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                  ,  decoder->lossless, decoder->ll_tns, h_DecSetup->ll_offQuant, decoder->ll_tns_remove
                                  , &ll_adap_flag, &b_relative, &off_idx, &tns_lsb_num_remove, deltaCodedBits
                                  , bits_per_sample
                                  , &fallback
                                  , &scaleSignal
#ifdef LL_INCL_HPVC
                                  , &(h_DecSetup->hpvcDecCfg)
#endif
#endif
                                );

        BW_cutoff_idx_nf = BW_cutoff_idx;
        move16();
    }
    BASOP_sub_end(); /* Entropy dec */

#ifdef CR14_A_ADD_LOSSLESS_MODE
    decoder->ll_adap_flag = ll_adap_flag;
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF( fallback )
    {
        IF( bits_per_sample == 16 )
        {
            fallback_bit_planes = 17;
        }
        ELSE IF( bits_per_sample == 24 )
        {
            fallback_bit_planes = 25;
        }

        fallback_decoder( bs_in, &bp_side, &mask_side, q_d_fx, decoder->frame_length, fallback_bit_planes );
        input_tda = 1;

        goto jump_to_imdct;
    }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    if ( scratch->max_scratch_calculation_only && decoder->lossless )
    {
        ll_adap_flag = 1;
        basop_memset( L_scf_idx, 0, sizeof( L_scf_idx[0] ) * SCF_MAX_PARAM );
    }

    // Count WMOPS for lossless vs. lossy code path
    if ( 1 == ll_adap_flag )
    {
        BASOP_sub_sub_start( "Dec(lossless)" );
    }
    else
    {
        BASOP_sub_sub_start( "Dec(lossy)" );
    }

    IF( ll_adap_flag )
    {
        Word16 gg_idx_off_ll = h_DecSetup->quantizedGainOff_ll;

        IF (h_DecSetup->ll_offQuant)
        {
            Word32 tmp;
            IF( decoder->fs_idx == 4 && bits_per_sample == 24 )
            {
                /*gg_idx_off_ll = round( -off_idx * 60.f / 7 - 130 );*/
                tmp = Mpy_32_32( L_shl( (Word32)i_mult( -off_idx, 60 ), 1 ), 306783378 );
                gg_idx_off_ll = sub( (Word16) L_shr( L_add( tmp, 1 ), 1 ), 130 );
            }
            ELSE IF( decoder->fs_idx >= 5 && bits_per_sample == 24 )
            {
                //gg_idx_off_ll = round( -off_idx * 100.f / 7 - 135 );
                tmp = Mpy_32_32( L_shl( (Word32)i_mult( -off_idx, 100 ), 1), 306783378 );
                gg_idx_off_ll = sub( (Word16) L_shr(L_add(tmp, 1),1), 135 );
            }

            IF (bits_per_sample == 24)
            {
                gg_idx_off_ll = gg_idx_off_ll + 48;
            }
        }


        Word32 gg; Word16 gg_e;

        int_scf_ll = (Word32*) lc3_scratch_push( scratch, sizeof( *int_scf_ll ) * MAX_BANDS_NUMBER );
        int_scf_exp_ll = (Word16*) lc3_scratch_push( scratch, sizeof( *int_scf_exp_ll ) * MAX_BANDS_NUMBER );

        processCalculateGlobalGain_fx( &gg, &gg_e, gg_idx, gg_idx_off_ll );

#if defined(CR15_A_LOSSLESS_1p25MS) && defined(CR9_C_ADD_1p25MS_LRSNS)
        IF( sub( decoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS ) == 0 )
        {
            pitch_rx_fx = ltpf_idx[0]; move16();
#  ifdef LRSNS_CBC_NO_LTPF_DEPENDENCY
            ltpf_rx_fx = 0; move16();
#  else
            ltpf_rx_fx = ltpf_idx[1]; move16();
#  endif
            snsQuantScfDecLR_fx( L_scf_idx, scf_q_ip, scf_q, pitch_rx_fx, ltpf_rx_fx, scratch );
#  ifdef ENABLE_HR_MODE
            downshift_w32_arr( scf_q_ip, scf_q, 26 - 11, M );
#  endif
        }
        ELSE
#endif
        processSnsQuantizeScfDecoder_fx( L_scf_idx, scf_q_ip, scratch );

        IF( decoder->fs_idx == 6 && bits_per_sample == 24 )
        {

            processSnsInterpolateScf_fx( scf_q_ip, int_scf_ll, int_scf_exp_ll, 1, decoder->low_band_limit, scratch );

            process_deterministic_curve( deterministic_curve, gg, gg_e, decoder->bands_offset[decoder->low_band_limit],
                          int_scf_ll, int_scf_exp_ll, decoder->bands_offset, decoder->low_band_limit ,scratch );

            Word16 startband = decoder->low_band_limit; move16();
            Word16 stopband = decoder->bands_number; move16();
            UWord8 lastseg = deterministic_curve[decoder->bands_offset[decoder->low_band_limit]-1];
            Word16 counter = 0; move16();

            FOR( Word16 i = startband; i < stopband && counter < HIGH_BANDS_NUMBER; i++ )
            {
                FOR( Word16 j = decoder->bands_offset[i]; j < decoder->bands_offset[i+1]; j++)
                {
                    deterministic_curve[j] =  lastseg - deltaCodedBits[counter]; move16();
                }
                lastseg -= deltaCodedBits[counter++]; move16();
            }
        }
        ELSE
        {
            processSnsInterpolateScf_fx( scf_q_ip, int_scf_ll, int_scf_exp_ll, 1, decoder->bands_number, scratch );

            process_deterministic_curve( deterministic_curve, gg, gg_e, decoder->frame_length, int_scf_ll, int_scf_exp_ll, decoder->bands_offset, decoder->bands_number, scratch );
        }


        int_scf_exp_ll = (Word16*) lc3_scratch_pop( scratch,  int_scf_exp_ll);
        int_scf_ll = (Word32*) lc3_scratch_pop( scratch, int_scf_ll );

    } ELSE {
        basop_memset(deterministic_curve, 0, sizeof(*deterministic_curve) * decoder->frame_length);
    }
#endif

    IF( scratch->max_scratch_calculation_only )
    {
        bfi = 0;
        move16();
        if ( bfi_ext )
        {
            bfi = 1;
            move16();
        }

        memset( L_scf_idx, 0, sizeof( *L_scf_idx ) * SCF_MAX_PARAM );
        memset( sqQdec, 0, sizeof( *sqQdec ) * decoder->frame_length );
    }

    BASOP_sub_start("Ari dec");
#ifdef LL_INCL_HPVC
    HpvcTreeEnumCfg decTrees[4 * HPVC_NOMTREE_COUNT_FB];
#endif
    IF (sub(bfi, 1) != 0)
    {
#ifdef LL_INCL_HPVC
        IF (h_DecSetup->hpvcDecCfg.active_flag != 0 && decoder->lossless != 0 && ll_adap_flag != 0)
        {
            h_DecSetup->hpvcDecCfg.HpvcTreeEnumCfgPtr = &(decTrees[0]);

#ifndef LL_HPVC_GLOBAL_FRAC
            {
                Word16 Nqm_ana;
                Nqm_ana = shr_pos(sub(lastnz, h_DecSetup->hpvcDecCfg.startCoef), N_SIGNAL_LOG);
                if (Nqm_ana <= 0)
                {
                    /* all-zero tail, disable further mixed TCX+HPVC decoding */
                    ASSERT((scratch->max_scratch_calculation_only != 0 || (h_DecSetup->hpvcDecCfg.mode < 0)) && "global signaling and lastnz(rx) mismatch");
                }
            }
#endif
        }
#endif



        processAriDecoder_fx(bs_in, &bp_side, &mask_side, h_DecSetup->total_bits, decoder->yLen, decoder->fs_idx,
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &bfi, tns_order,
                             fac_ns_idx, gg_idx, decoder->frame_dms,
                             decoder->n_pc, decoder->be_bp_left, decoder->be_bp_right, 0, &spec_inv_idx, &scale,
                             &fill_bits, sqQdec, &nf_seed, resBitBufLossless, indexes, &zero_frame, scratch
#ifdef ENABLE_HR_MODE
                             , decoder->hrmode
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                              , deterministic_curve, ll_adap_flag
                              ,  decoder->lossless
#endif
#ifdef LL_INCL_HPVC
                              , &(h_DecSetup->hpvcDecCfg)
#endif
        );



        if ( scratch->max_scratch_calculation_only )
        {
            spec_inv_idx = decoder->yLen;  // required to trigger noise filling, which necessary to run the imdct properly
            move16();
        }

#ifdef ENABLE_RFRAME
        test();test();
        IF (sub(rframe, 1) == 0 && zero_frame == 0 && sub(bfi, 1) != 0)
        {
            bfi = 2;
            move16();
            Word16 max_bw_stopband = BW_cutoff_bin_all[BW_cutoff_idx];
            SWITCH (decoder->frame_dms)
            {
#ifdef CR9_C_ADD_1p25MS
              case LC3PLUS_FRAME_DURATION_1p25MS:
                  max_bw_stopband  = shr_pos(max_bw_stopband, 3);
                  BREAK;
#endif
              case LC3PLUS_FRAME_DURATION_2p5MS:
                  max_bw_stopband  = shr_pos(max_bw_stopband, 2);
                  BREAK;
              case LC3PLUS_FRAME_DURATION_5MS:
                  max_bw_stopband  = shr_pos(max_bw_stopband, 1);
                  BREAK;
              case LC3PLUS_FRAME_DURATION_7p5MS:
                  max_bw_stopband = add(shr_pos(max_bw_stopband, 2), add(shr_pos(max_bw_stopband, 2), shr_pos(max_bw_stopband, 2)));
                  BREAK;
              case LC3PLUS_FRAME_DURATION_10MS:
                  BREAK;
              case LC3PLUS_FRAME_DURATION_UNDEFINED:
                  assert(0);
            }

            spec_inv_idx = s_max(lastnz, max_bw_stopband);
            move16();
        }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
        IF( bfi == 0 && (decoder->lossless == 0 || ll_adap_flag == 0 ) )
#else
        IF( bfi == 0 )
#endif
        {
            processAriDecoderScaling_fx(sqQdec, decoder->yLen, q_d_fx, &q_fx_exp);
        }
    }
    BASOP_sub_end(); /* Ari dec */

    #ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (scaleSignal && !ll_adap_flag)
    {
        q_fx_exp = q_fx_exp + scaleSignal;
    }
    #endif

    IF( scratch->max_scratch_calculation_only )
    {
        bfi = 0;
        move16();
        if ( bfi_ext )
        {
            bfi = 1;
            move16();
        }
    }

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (sub(bfi, 1) != 0)
    {
        q_d_res = (Word32*) lc3_scratch_push( scratch,  sizeof(*q_d_res) * decoder->frame_length );
        IF (ll_adap_flag)
        {


            eff_det_curve = (UWord8*) lc3_scratch_push( scratch,  sizeof(*eff_det_curve) * decoder->frame_length );
            q_res_r = (Word32*) lc3_scratch_push( scratch,  sizeof(*q_res_r) * decoder->frame_length );
            q_res = (Word32*) lc3_scratch_push( scratch,  sizeof(*q_res) * decoder->frame_length );

            FOR (i = 0; i < fill_bits; i++)
            {
                IF (resBitBufLossless[i >> RESBITS_PACK_SHIFT] & (1 << (i & RESBITS_PACK_MASK)))
                {
                    residualDataLossless[i] = 1;
                }
                ELSE
                {
                    residualDataLossless[i] = 0;
                }
            }

            compute_resbits_priority( sqQdec, deterministic_curve, decoder->bands_offset, decoder->bands_number, decoder->frame_length, fill_bits, b_relative, eff_det_curve, scratch );

            res_bit_pos = residual_decoder_lossless( sqQdec, q_res, residualDataLossless, deterministic_curve, eff_det_curve, decoder->frame_length, fill_bits, &index_b, &index_x, scratch );
            UNUSED(res_bit_pos);

            IF ((tns_order[0] + tns_order[1]) <= 0)
            {
                processProperRounding_fx( deterministic_curve, index_b, index_x, q_res, decoder->frame_length, q_res_r, eff_det_curve, scratch );

                IF( ll_dec_rounding )
                {
                    basop_memcpy( q_res, q_res_r, sizeof( *q_res_r ) * decoder->frame_length );
                }
            }

            basop_memcpy(q_d_res, q_res, sizeof(*q_res) * decoder->frame_length);


            q_res  = (Word32*) lc3_scratch_pop( scratch,  q_res);
            q_res_r = (Word32*) lc3_scratch_pop( scratch,  q_res_r );
            eff_det_curve = (UWord8*) lc3_scratch_pop( scratch,  eff_det_curve );

        }
    }
#endif



    BASOP_sub_start("SnsQuantScfDec");

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (decoder->lossless == 0 || ll_adap_flag == 0)
#endif
    {
        IF (sub(bfi, 1) != 0)
    #ifdef  CR9_C_ADD_1p25MS_LRSNS
        {
            IF(sub(decoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)
            {
                pitch_rx_fx = ltpf_idx[0]; move16();

    #ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
                ltpf_rx_fx = 0;           move16(); /* CB_C with binary means ,  not dependent on LTPF activation */
    #else
                ltpf_rx_fx  = ltpf_idx[1]; move16();/*  CB_C, with ternary means  dependent on LTPF activation */
    #endif
                snsQuantScfDecLR_fx(L_scf_idx, scf_q_ip, scf_q, pitch_rx_fx, ltpf_rx_fx, scratch); /*  9,12,29,30,  bits decoding and 2 pitch info bits  */
    #ifdef ENABLE_HR_MODE
                downshift_w32_arr(scf_q_ip /* Q26 */, scf_q/* Q11 */, 26 - 11, M);  /* W16Q11 version required for PLC  */
    #endif
            }
            ELSE
    #endif   /* CR9_C_ADD_1p25MS_LRSNS */
            {
    #ifdef ENABLE_HR_MODE
            processSnsQuantizeScfDecoder_fx(L_scf_idx, scf_q_ip, scratch);
            downshift_w32_arr(scf_q_ip, scf_q, 15, M); /* required for PLC */
    #else
            processSnsQuantizeScfDecoder_fx(L_scf_idx, scf_q, scratch);
    #endif
            }
    #ifdef  CR9_C_ADD_1p25MS_LRSNS
        }
    #endif
    }
#ifdef CR14_A_ADD_LOSSLESS_MODE
    ELSE
    {
        IF( sub( bfi, 1 ) != 0 )
        {
            downshift_w32_arr( scf_q_ip, scf_q, 15, M ); /* required for PLC */
        }
    }
#endif
BASOP_sub_end();

    BASOP_sub_start("PLC::ComputeStabFac");
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (decoder->lossless == 0)
#endif
    {
        if (h_DecSetup->plcAd)
        {
            processPLCcomputeStabFac_main(scf_q, h_DecSetup->plcAd->old_scf_q, h_DecSetup->plcAd->old_old_scf_q, bfi,
                                          h_DecSetup->prev_bfi, h_DecSetup->prev_prev_bfi, &h_DecSetup->plcAd->stab_fac);
        }
    }
    BASOP_sub_end();

    BASOP_sub_start("Partial Concealment");
    IF (sub(bfi, 1) != 0)
    {
        scale = 32767;
        move16();

        IF (h_DecSetup->plcAd)
        {
            scale = h_DecSetup->plcAd->stab_fac;
        }

        processPCmain_fx(rframe, &bfi, decoder->yLen, decoder->frame_dms, h_DecSetup->q_old_res_fx,
                         &h_DecSetup->q_old_res_fx_exp, sqQdec, h_DecSetup->q_old_d_fx, spec_inv_idx, ltpf_idx[0],
                         scale, q_d_fx, &q_fx_exp, gg_idx, h_DecSetup->quantizedGainOff, &h_DecSetup->prev_gg,
                         &h_DecSetup->prev_gg_e, &BW_cutoff_idx_nf, &h_DecSetup->prev_BW_cutoff_idx_nf, fac_ns_idx,
                         &h_DecSetup->prev_fac_ns_fx, &h_DecSetup->pc_nbLostFramesInRow);
    }
    BASOP_sub_end();

#ifdef FIX_PLC_CONFORM_ISSUES
#ifdef CR9_C_ADD_1p25MS
    IF( sub( bfi, 1 ) == 0 )
    {
#ifdef FIX_TX_RX_STRUCT_STEREO
        h_DecSetup->ltpf_rx_status[0] = 0;
        h_DecSetup->ltpf_rx_status[1] = 0;
#else
        decoder->ltpf_rx_status[0] = 0;
        decoder->ltpf_rx_status[1] = 0;
#endif
    }
#endif
#endif

    IF( scratch->max_scratch_calculation_only )
    {
        bfi = 0;
        move16();
        if ( bfi_ext )
        {
            bfi = 1;
            move16();
        }

        zero_frame = 0;
    }

    IF (sub(bfi, 1) != 0)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
       IF (decoder->lossless == 0 || ll_adap_flag == 0)
#endif
       {
        BASOP_sub_start("Residual dec");
        processResidualDecoding_fx(q_d_fx, q_fx_exp, decoder->yLen, resBitBufLossless, fill_bits
#ifdef ENABLE_HR_MODE
                                   , decoder->hrmode
#endif
#if defined (CR9_C_ADD_1p25MS)
                                   , decoder->frame_dms
#endif
        );
        BASOP_sub_end();
#ifdef CR14_A_ADD_LOSSLESS_MODE
           basop_memcpy(q_d_res, q_d_fx, sizeof(*q_d_fx) * decoder->frame_length);
#endif
        }

#ifdef CR14_A_ADD_LOSSLESS_MODE
       IF (decoder->lossless == 0 || ll_adap_flag == 0)
#endif
       {
        BASOP_sub_start("Noisefill");
#ifdef CR9_C_ADD_1p25MS
        IF (zero_frame == 0)
#else
        IF (zero_frame == 0)
#endif
        {
            processNoiseFilling_fx(q_d_fx, nf_seed, q_fx_exp, fac_ns_idx, BW_cutoff_idx_nf, decoder->frame_dms,
                                   h_DecSetup->prev_fac_ns_fx, spec_inv_idx, scratch
#ifdef ENABLE_HR_MODE
                                   , decoder->hrmode
#endif
            );
        }
        BASOP_sub_end();
       }

#ifdef CR14_A_ADD_LOSSLESS_MODE
       IF (decoder->lossless == 0 || ll_adap_flag == 0)
#endif
       {
        BASOP_sub_start("applyGlobalGain");
        processApplyGlobalGain_fx(q_d_fx, &q_fx_exp, decoder->yLen, gg_idx, h_DecSetup->quantizedGainOff);
        BASOP_sub_end();
       }



#ifdef CR14_A_ADD_LOSSLESS_MODE
       IF( ll_adap_flag )
       {
           basop_memcpy( q_d_fx, q_d_res, sizeof( *q_d_fx ) * decoder->frame_length );
       }
       q_d_res = (Word32*) lc3_scratch_pop( scratch, q_d_res);




       IF (decoder->lossless == 0 || decoder->ll_tns)
#endif
       {
#ifdef CR9_C_ADD_1p25MS
        if (tns_numfilters > 0) {
#endif



        BASOP_sub_start("Tns_dec");
        processTnsDecoder_fx(indexes, q_d_fx, decoder->yLen, tns_order, &q_fx_exp, BW_cutoff_idx, decoder->frame_dms,
                             scratch
#ifdef ENABLE_HR_MODE
                             , decoder->hrmode
#endif
#  ifdef CR14_A_ADD_LOSSLESS_MODE
                             , ll_adap_flag
#  endif
        );




#  ifdef CR14_A_ADD_LOSSLESS_MODE

            tns_lsb_add = (UWord8*) lc3_scratch_push( scratch, sizeof(*tns_lsb_add) * decoder->frame_length );

            IF (tns_lsb_num_remove > 0)
            {
                FOR (i = 0; i < decoder->yLen; i++) {
                    tns_lsb_add[i] = tns_lsb_num_remove;
                }
                process_lsb_add(q_d_fx, tns_lsb_add, decoder->yLen, q_d_fx);
                res_bit_pos = residual_decoder_lossless( q_d_fx, q_d_fx, residualDataLossless + res_bit_pos, tns_lsb_add, tns_lsb_add, decoder->frame_length, fill_bits, &index_b, &index_x, scratch );
                UNUSED(res_bit_pos);
            }

            tns_lsb_add = (UWord8*) lc3_scratch_pop( scratch, tns_lsb_add );
#  endif

        BASOP_sub_end();
#ifdef CR9_C_ADD_1p25MS
        }
#endif
       }



#ifdef CR14_A_ADD_LOSSLESS_MODE
       IF (decoder->lossless == 0 || ll_adap_flag == 0)
#endif
       {
#ifdef ENABLE_HR_MODE
        BASOP_sub_start("SnsInterpScfDec");
        processSnsInterpolateScf_fx(scf_q_ip, int_scf_fx_ip, int_scf_fx_exp, 0, decoder->bands_number, scratch);

        BASOP_sub_end();

        BASOP_sub_start("Mdct shaping_dec");
        processScfScaling(int_scf_fx_exp, decoder->bands_number, &q_fx_exp);

        processMdctShaping_fx(q_d_fx, int_scf_fx_ip, int_scf_fx_exp, decoder->bands_offset, decoder->bands_number);
        BASOP_sub_end();
#else
        BASOP_sub_start("SnsInterpScfDec");
        processSnsInterpolateScf_fx(scf_q, int_scf_fx, int_scf_fx_exp, 0, decoder->bands_number, scratch);
        BASOP_sub_end();

        BASOP_sub_start("Mdct shaping_dec");
        processScfScaling(int_scf_fx_exp, decoder->bands_number, &q_fx_exp);
        processMdctShaping_fx(q_d_fx, int_scf_fx, int_scf_fx_exp, decoder->bands_offset, decoder->bands_number);
        BASOP_sub_end();
        /* end int_scf_fx */
#endif /* ENABLE_HR_MODE */
       }
    }



    /* x_fx_ip will be used to store h_DecSetup->stDec_ola_mem_fx returned by PLCmain_fx*/
    /* This will be upshifted to 32 bit overlap buffer outside of the PLCmain function */
#ifdef CR14_A_ADD_LOSSLESS_MODE
jump_to_imdct:
    UNUSED(bfi);
#endif
#ifdef ENABLE_HR_MODE
    Word16 *plc_ola_mem = (Word16 *)x_fx_ip;
    IF(sub(bfi, 1) == 0)
    {
        FOR(i = 0; i < decoder->stDec_ola_mem_fx_len; i++)
        {
            plc_ola_mem[i] =  round_fx(h_DecSetup->stDec_ola_mem_fx[i]);
        }
    }
#endif

    IF( scratch->max_scratch_calculation_only )
    {
        h_DecSetup->ltpf_mem_pitch_int = MAX_PITCH_FS( decoder->fs );

        bfi = 0;
        move16();
        if ( bfi_ext )
        {
            bfi = 1;
            move16();
        }
    }

    BASOP_sub_start("PLC::Main");
    processPLCmain_fx(decoder->plcMeth, &h_DecSetup->concealMethod, &h_DecSetup->nbLostFramesInRow, bfi,
                      h_DecSetup->prev_bfi, decoder->frame_length, decoder->la_zeroes, decoder->W_fx, x_fx,
#ifdef ENABLE_HR_MODE
                      plc_ola_mem,
#else
                      h_DecSetup->stDec_ola_mem_fx,
#endif
                      &h_DecSetup->stDec_ola_mem_fx_exp, h_DecSetup->q_old_d_fx,
                      &h_DecSetup->q_old_fx_exp, q_d_fx, &q_fx_exp, decoder->yLen, decoder->fs_idx,
                      decoder->bands_offset, decoder->bands_number, &h_DecSetup->plc_damping, h_DecSetup->ltpf_mem_pitch_int,
                      h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ns_cum_alpha, &h_DecSetup->ns_seed, h_DecSetup->plcAd,
                      decoder->frame_dms, scratch, &h_DecSetup->pc_nbLostFramesInRow
#ifdef ENABLE_HR_MODE
                      , decoder->hrmode
#endif
                      , h_DecSetup->rel_pitch_change
                      , decoder->alpha_type_2_table
#ifdef CR14_A_ADD_LOSSLESS_MODE
                      , ll_adap_flag
                      , bits_per_sample
                      , decoder->lossless
#endif
                      );
    BASOP_sub_end();

#ifdef ENABLE_HR_MODE
    IF(sub(bfi, 1) == 0)
    {
#ifdef CR15_A_LOSSLESS_1p25MS
        IF( decoder->lossless && sub( decoder->frame_dms, LC3PLUS_FRAME_DURATION_1p25MS ) == 0 )
        {
            FOR(i = 0; i < decoder->stDec_ola_mem_fx_len; i++)
            {
                IF( sub( plc_ola_mem[i], round_fx( h_DecSetup->stDec_ola_mem_fx[i] ) ) != 0 )
                {
                    h_DecSetup->stDec_ola_mem_fx[i] = L_deposit_h(plc_ola_mem[i]);
                    move32();
                }
            }
        }
        ELSE
#endif
        {
            FOR(i = 0; i < decoder->stDec_ola_mem_fx_len; i++)
            {
                h_DecSetup->stDec_ola_mem_fx[i] = L_deposit_h(plc_ola_mem[i]);
            }
        }
    }
#endif

    BASOP_sub_start("PLC/PC::DampingScrambling");
    if (h_DecSetup->plcAd)
    {
        processPLCDampingScrambling_main_fx(
            bfi, h_DecSetup->concealMethod, h_DecSetup->nbLostFramesInRow, &h_DecSetup->plcAd->cum_fflcAtten,
            h_DecSetup->pc_nbLostFramesInRow, &h_DecSetup->ns_seed, &h_DecSetup->pc_seed,
            h_DecSetup->ltpf_mem_pitch_int, ltpf_idx[0], q_d_fx, &q_fx_exp, h_DecSetup->q_old_d_fx,
            &h_DecSetup->q_old_fx_exp, decoder->yLen, h_DecSetup->plcAd->stab_fac, decoder->frame_dms,
            &h_DecSetup->plcAd->cum_fading_slow, &h_DecSetup->plcAd->cum_fading_fast, spec_inv_idx
            , h_DecSetup->plcAd->plc_fadeout_type
            #ifdef CR14_A_ADD_LOSSLESS_MODE
            , ll_adap_flag,
            bits_per_sample
            #endif
        );
    }
    BASOP_sub_end();

    BASOP_sub_start("Imdct");

#ifdef CR14_A_ADD_LOSSLESS_MODE
#ifdef CR15_A_LOSSLESS_1p25MS
        if (decoder->lossless) {
#else
        if (decoder->lossless && decoder->frame_dms != LC3PLUS_FRAME_DURATION_1p25MS) {
#endif
            test(); test(); test();
            IF( sub( bfi, 1 ) != 0 || sub( h_DecSetup->concealMethod, LC3_CON_TEC_NS_STD ) == 0 || sub( h_DecSetup->concealMethod, LC3_CON_TEC_NS_ADV ) == 0 || sub( h_DecSetup->concealMethod, LC3_CON_TEC_FREQ_MUTING ) == 0 )
            {
            IF ( ll_adap_flag == 1 && h_DecSetup->ll_adap_prev == 0 )
            {
                Word16 shift_ola = sub( 31, h_DecSetup->stDec_ola_mem_fx_exp );
                IF( bits_per_sample == 24 )
                {
                    shift_ola = sub( shift_ola, 8 - scaleSignal );
                }
            
                shift_ola = sub( shift_ola, 1);

                if ( decoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && decoder->lossless )
                {
                    Word16 hr = getScaleFactor32_0( h_DecSetup->stDec_ola_mem_fx,
                                                   decoder->stDec_ola_mem_fx_len );
                    Word16 min_safe = sub( 4, hr );
                    if ( sub( shift_ola, min_safe ) < 0 )
                    {
                        shift_ola = min_safe;
                    }
                }

                Word16 i;
                FOR( i = 0; i < decoder->frame_length>>1; i++ )
                {
                    h_DecSetup->stDec_ola_mem_fx[i] = L_shr( L_add( L_shr( h_DecSetup->stDec_ola_mem_fx[i], shift_ola ),1),1);        move32();
                }
                /* TODO why is there a diff of 2 ???*/
                shift_ola = sub( shift_ola, 2);

                if ( decoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && decoder->lossless )
                {
                    Word16 hr = getScaleFactor32_0( &h_DecSetup->stDec_ola_mem_fx[i],
                                                   sub( decoder->stDec_ola_mem_fx_len, i ) );
                    Word16 min_safe = sub( 4, hr );
                    if ( sub( shift_ola, min_safe ) < 0 )
                    {
                        shift_ola = min_safe;
                    }
                }

                FOR(; i < decoder->stDec_ola_mem_fx_len; i++ )
                {
                    h_DecSetup->stDec_ola_mem_fx[i] = L_shr( L_add( L_shr( h_DecSetup->stDec_ola_mem_fx[i], shift_ola ),1),1);        move32();
                }
            }
            IF (decoder->lossless && ll_adap_flag)
            {
                /*set exponent of overlap memory for possible transitions to lossy path*/
                IF( bits_per_sample == 24 )
                {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                    IF( h_DecSetup->ScaleSignal_Memory  > scaleSignal)
                    {
                            h_DecSetup->stDec_ola_mem_fx_exp = 23 + h_DecSetup->ScaleSignal_Memory  ;
                    }
                    ELSE
                    {
                            h_DecSetup->stDec_ola_mem_fx_exp = 23 + scaleSignal ;
                    }                
#endif 
                }
                ELSE
                {
                    h_DecSetup->stDec_ola_mem_fx_exp = 31;
                }

            }

#ifdef CR14_A_ADD_LOSSLESS_MODE
            IF ( ll_adap_flag == 1 && h_DecSetup->ll_adap_prev == 1 && (h_DecSetup->ScaleSignal_Memory - scaleSignal > 0) )
            {
                Word16 ii;
                FOR( ii = 0; ii < decoder->stDec_ola_mem_fx_len; ii++ )
                {
                    h_DecSetup->stDec_ola_mem_fx[ii] = L_shl(h_DecSetup->stDec_ola_mem_fx[ii], h_DecSetup->ScaleSignal_Memory - scaleSignal);
                }
            }
            IF ( ll_adap_flag == 1 && h_DecSetup->ll_adap_prev == 1 && (scaleSignal - h_DecSetup->ScaleSignal_Memory) > 0 )
            {
                Word16 ii;
                FOR( ii = 0; ii < decoder->stDec_ola_mem_fx_len; ii++ )
                {
                    h_DecSetup->stDec_ola_mem_fx[ii] = L_shr(h_DecSetup->stDec_ola_mem_fx[ii], scaleSignal - h_DecSetup->ScaleSignal_Memory);
                }
            }
#endif

            const int32_t *lift[4];
            getLiftingCoeffs2_fx(decoder->frame_dms, decoder->fs, lift);
            if (ll_adap_flag) {
                invIntMdct2_fx(q_d_fx, x_fx_ip, decoder->frame_length, lift, decoder->la_zeroes, h_DecSetup->stDec_ola_mem_fx, &h_DecSetup->stDec_ola_mem_fx[decoder->frame_length>>1],
                input_tda,
                scratch );

            } else {
                /* TODO why is there a diff of 2 ???*/
                IF ( h_DecSetup->ll_adap_prev == 1 ) {
                    FOR(Word16 i=decoder->frame_length>>1; i < decoder->stDec_ola_mem_fx_len; i++ )
                    {
                        h_DecSetup->stDec_ola_mem_fx[i] = L_shr( h_DecSetup->stDec_ola_mem_fx[i], 2 );        move32();
                    }
                }
                invIntMdct2lossy_fx(q_d_fx, x_fx_ip, decoder->frame_length, &q_fx_exp, &h_DecSetup->stDec_ola_mem_fx_exp, lift, decoder->la_zeroes, h_DecSetup->stDec_ola_mem_fx, &h_DecSetup->stDec_ola_mem_fx[decoder->frame_length>>1],
                input_tda,
                scratch, decoder->frame_dms );

                    IF( sub( bfi, 1 ) == 0)
                    {
                        if ( h_DecSetup->stDec_ola_mem_fx_exp < 0 )
                        {
                            h_DecSetup->stDec_ola_mem_fx_exp = 0;
                        }
                    }
            }
            }
        } ELSE {
#endif
    ProcessingIMDCT(q_d_fx, &q_fx_exp, decoder->W_fx, h_DecSetup->stDec_ola_mem_fx, &h_DecSetup->stDec_ola_mem_fx_exp,
#ifdef ENABLE_HR_MODE
                    x_fx_ip,
#else
                    x_fx,
#endif
                    decoder->W_size, decoder->frame_length, decoder->stDec_ola_mem_fx_len, decoder->frame_dms,
                    h_DecSetup->concealMethod, bfi, h_DecSetup->prev_bfi, h_DecSetup->nbLostFramesInRow,
                    h_DecSetup->plcAd,
                    scratch
#ifdef ENABLE_HR_MODE
                    , decoder->hrmode
#endif
    );
#ifdef CR14_A_ADD_LOSSLESS_MODE
        }
#endif

#ifdef ENABLE_HR_MODE
        IF(sub(bfi, 1) != 0 || sub(h_DecSetup->concealMethod, LC3_CON_TEC_NS_STD) == 0 || sub(h_DecSetup->concealMethod, LC3_CON_TEC_NS_ADV) == 0 || sub(h_DecSetup->concealMethod, LC3_CON_TEC_FREQ_MUTING) == 0)
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Word16 headroom = 1;
            Word16 y_s_o = getScaleFactor32_0(x_fx_ip, decoder->frame_length) - headroom;

#ifdef CR15_A_LOSSLESS_1p25MS
            if (ll_adap_flag || decoder->lossless)
#else
            if ((ll_adap_flag || decoder->lossless) && decoder->frame_dms != LC3PLUS_FRAME_DURATION_1p25MS)
#endif
            {
                q_fx_exp = sub(q_fx_exp,y_s_o);

                if(ll_adap_flag)
                {
                    q_fx_exp = 31 - y_s_o ;

                    IF( bits_per_sample == 24 )
                    {
                        /*additional shift by 8 for 24-bit lossless as lossy path in previous frame(s) was operating on 16-bit representation*/
                        q_fx_exp = sub( q_fx_exp, 8 - scaleSignal);
                    }

                }

                FOR(int i=0;i<decoder->frame_length;i++)
                {
                    x_fx_ip[i] = L_shl( x_fx_ip[i], y_s_o );
                }
            }

            round_w32tow16_arr( x_fx_ip, x_fx, decoder->frame_length );

            if(ll_adap_flag)
            {
                FOR(int i=0;i<decoder->frame_length;i++)
                {
                    x_fx_ip[i] = L_shr( x_fx_ip[i], y_s_o );
                }
            }

#else
            round_w32tow16_arr(x_fx_ip, x_fx, decoder->frame_length);
#endif
        }
        ELSE
        {
            FOR(i = 0; i < decoder->frame_length; i++)
            {
                x_fx_ip[i] = L_deposit_h(x_fx[i]);
            }
        }
#endif /* ENABLE_HR_MODE */



    BASOP_sub_end();

    BASOP_sub_start("PLC::Update");

    processPLCupdate_fx(h_DecSetup->plcAd, x_fx, q_fx_exp, h_DecSetup->concealMethod, decoder->frame_length,
                        decoder->fs_idx, &h_DecSetup->nbLostFramesInRow, &h_DecSetup->prev_prev_bfi,
                        &h_DecSetup->prev_bfi, bfi, scf_q, &h_DecSetup->ns_cum_alpha
#ifdef ENABLE_HR_MODE
                        , decoder->hrmode
#endif
                        );

    #ifdef CR14_A_ADD_LOSSLESS_MODE
    if(ll_adap_flag)
    {
        q_fx_exp = 0;
    }
    #endif

    BASOP_sub_end();

    if ( scratch->max_scratch_calculation_only )
    {
        h_DecSetup->ltpf_mem_active = 1;
        h_DecSetup->ltpf_mem_scale_fac_idx = 0;
    }

    BASOP_sub_start("LtpfDec");
    process_ltpf_decoder_fx(&q_fx_exp, decoder->frame_length, decoder->ltpf_mem_x_len, decoder->fs_idx,
                            decoder->ltpf_mem_y_len, &h_DecSetup->ltpf_mem_e, x_fx, h_DecSetup->ltpf_mem_x, x_fx,
                            h_DecSetup->ltpf_mem_y, ltpf_idx[0], ltpf_idx[1], ltpf_idx[2],
                            &h_DecSetup->ltpf_mem_pitch_int, &h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ltpf_mem_gain,
                            &h_DecSetup->ltpf_mem_active, h_DecSetup->ltpf_scale_fac_idx, bfi,
                            h_DecSetup->concealMethod,
                            h_DecSetup->plc_damping, &h_DecSetup->ltpf_mem_scale_fac_idx,
                            &h_DecSetup->rel_pitch_change, decoder->hrmode, decoder->frame_dms,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            decoder->lossless,
#endif
                            scratch
#ifdef CR9_C_ADD_1p25MS
#ifdef FIX_TX_RX_STRUCT_STEREO
                            ,&h_DecSetup->ltpf_mem_continuation, &h_DecSetup->ltpf_mem_pitch_int_prev,
                            &h_DecSetup->ltpf_mem_pitch_fr_prev, &h_DecSetup->ltpf_mem_beta_idx_prev, &h_DecSetup->ltpf_mem_gain_prev,
                            &h_DecSetup->ltpf_mem_active_prev, &h_DecSetup->ltpf_pitch_stability_counter
#else
                            , &decoder->ltpf_mem_continuation, &decoder->ltpf_mem_pitch_int_prev,
                            &decoder->ltpf_mem_pitch_fr_prev, &decoder->ltpf_mem_beta_idx_prev, &decoder->ltpf_mem_gain_prev,
                            &decoder->ltpf_mem_active_prev, &decoder->ltpf_pitch_stability_counter
#endif
#endif
                           );
    BASOP_sub_end();

#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_1p25MS_HR
    IF (!(decoder->hrmode) || (decoder->hrmode && decoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && h_DecSetup->ltpf_scale_fac_idx != -1))
#else
    IF (!(decoder->hrmode))
#endif
    {
        FOR (i = 0; i < decoder->frame_length; i++)
        {
            x_fx_ip[i] = L_deposit_h(x_fx[i]);
        }
    }
#endif

    BASOP_sub_start("Output scaling");

#ifdef CR14_A_ADD_LOSSLESS_MODE
    if ((scaleSignal || h_DecSetup->ScaleSignal_Memory) && ll_adap_flag)
    {
        rescale_signal_decoder(x_fx_ip, scaleSignal, h_DecSetup->ScaleSignal_Memory, decoder->frame_length);
    }
    h_DecSetup->ScaleSignal_Memory = scaleSignal;
#endif

#ifndef CR14_A_ADD_LOSSLESS_MODE
    {
        scale  = sub(sub(31 + 16, bits_per_sample), q_fx_exp);
        offset = L_shr_sat(32768, sub(16, scale));
        IF (bits_per_sample == 16)
        {
            scale = sub(15, q_fx_exp);
            FOR (i = 0; i < decoder->frame_length; i++)
            {
#ifdef ENABLE_HR_MODE
                ((Word16 *)s_out)[i] = round_fx_sat(L_shr_sat(x_fx_ip[i], scale));
#else
                ((Word16 *)s_out)[i] = round_fx_sat(L_shr_sat(L_deposit_h(x_fx[i]), scale));
#endif
                move16();
            }
        }
        ELSE
        {
            FOR (i = 0; i < decoder->frame_length; i++)
            {
#ifdef ENABLE_HR_MODE
                ((Word32 *)s_out)[i] = L_shr_sat(L_add_sat(x_fx_ip[i], offset), scale);
#else
                ((Word32 *)s_out)[i] = L_shr_sat(L_add_sat(L_deposit_h(x_fx[i]), offset), scale);
#endif
                move32();
            }
        }
    }
#else /* CR14_A_ADD_LOSSLESS_MODE */
#  ifdef ENABLE_HR_MODE
        format_out_pcm( bits_per_sample, x_fx_ip, q_fx_exp, s_out, decoder->frame_length
                       , decoder->lossless, ll_adap_flag);
#  else
        format_out_pcm( bits_per_sample, x_fx, q_fx_exp, s_out, decoder->frame_length );
#  endif
#endif /* CR14_A_ADD_LOSSLESS_MODE */
    BASOP_sub_end(); /* Output scaling */

    BASOP_sub_sub_end();

#ifdef CR14_A_ADD_LOSSLESS_MODE
    BASOP_sub_sub_end();
#endif

    BASOP_sub_end(); /* Decoder */

#ifdef CR14_A_ADD_LOSSLESS_MODE
    h_DecSetup->ll_adap_prev = ll_adap_flag;
#endif

#ifdef ENABLE_HR_MODE
    int_scf_fx_ip = (Word32*) lc3_scratch_pop( scratch, int_scf_fx_ip );
    x_fx_ip = (Word32*) lc3_scratch_pop( scratch, x_fx_ip );
#endif
    x_fx = (Word16*) lc3_scratch_pop( scratch, x_fx );
    int_scf_fx = (Word16*) lc3_scratch_pop( scratch, int_scf_fx );
    int_scf_fx_exp = (Word16*) lc3_scratch_pop( scratch, int_scf_fx_exp );
#  ifdef CR14_A_ADD_LOSSLESS_MODE
    deterministic_curve = (UWord8*) lc3_scratch_pop( scratch, deterministic_curve );
#endif
#  ifdef ENABLE_HR_MODE
    sqQdec = (Word32*) lc3_scratch_pop( scratch, sqQdec );
#  else
    sqQdec = (Word16*) lc3_scratch_pop( scratch, sqQdec );
#  endif
    L_scf_idx = (Word32*) lc3_scratch_pop( scratch, L_scf_idx );
    indexes = (Word16*) lc3_scratch_pop( scratch, indexes );
#ifdef CR14_A_ADD_LOSSLESS_MODE
    resBitBufLossless = (UWord8*) lc3_scratch_pop( scratch, resBitBufLossless );
#endif
#  ifdef CR14_A_ADD_LOSSLESS_MODE
#  ifdef ENABLE_HR_MODE
    IF ( decoder->hrmode && decoder->lossless )
    {
        residualDataLossless = (UWord8*) lc3_scratch_pop( scratch, residualDataLossless );
    }
#  endif
#  endif
    resBitBuf = (UWord8*) lc3_scratch_pop( scratch, resBitBuf );
    q_d_fx = (Word32*) lc3_scratch_pop( scratch, q_d_fx );

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return bfi;
}

/* num_bytes = 0 -> bad frame */
LC3PLUS_Error Dec_LC3PLUS(LC3PLUS_Dec *decoder, UWord8 *input,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                          int *num_bytes,
#else
                          int num_bytes,
#endif
                          void **output, int bits_per_sample, lc3_scratch_t scratch,
                  int bfi_ext)
{
    int       ch = 0, bfi = bfi_ext;
    LC3PLUS_Error err = LC3PLUS_OK;
    int       fec_num_bytes;
    int       lc3_num_bytes;
    int       lc3_channel_num_bytes;
    int       channel_bfi, out_bfi;
    Word16    channel_epmr;

    if (bfi == 0)
    {
        bfi = !num_bytes;
    }

    if (decoder->ep_enabled)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        Word32 num_bytes_total = 0;
        for ( ch = 0; ch < decoder->channels; ch++ )
        {
            num_bytes_total += num_bytes[ch];

            if ( scratch->max_scratch_calculation_only )
            {
                break;
            }
        }
        decoder->combined_channel_coding = decoder->channels > 1 && num_bytes_total <= 160;
#else
        decoder->combined_channel_coding = decoder->channels > 1 && num_bytes <= 160;
#endif

        if (decoder->combined_channel_coding)
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            fec_num_bytes = num_bytes_total;
#else
            fec_num_bytes = num_bytes;
#endif

            BASOP_sub_start("fec_dec");

            decoder->error_report =
                fec_decoder(input, fec_num_bytes, &lc3_num_bytes, &decoder->epmr, decoder->combined_channel_coding,
                            &decoder->n_pccw, &bfi, &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc,
                            &decoder->m_fec, scratch);

            BASOP_sub_end();

            for (ch = 0; ch < decoder->channels; ch++)
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                Word32 tmp = lc3_num_bytes >> (decoder->channels-1);
                lc3_channel_num_bytes = tmp + ( ch < ( lc3_num_bytes - (tmp<<(decoder->channels-1)) ));
#else
                lc3_channel_num_bytes = lc3_num_bytes / decoder->channels + (ch < (lc3_num_bytes % decoder->channels));
#endif

                if (bfi != 1 && lc3_channel_num_bytes != decoder->channel_setup[ch]->last_size)
                {
                    err = update_dec_bitrate(decoder, ch, lc3_channel_num_bytes);

                    if (err)
                    {
                        bfi = 1;
                    }
                    else
                    {
                        decoder->channel_setup[ch]->last_size = lc3_channel_num_bytes;
                    }
                }

                if ( scratch->max_scratch_calculation_only )
                {
                    bfi = bfi_ext;
                }

                bfi = Dec_LC3PLUS_Channel(decoder, ch, bits_per_sample, input, output[ch], bfi, scratch);
                if (input != NULL)
                {
                    input += decoder->channel_setup[ch]->targetBytes;
                }

                if ( scratch->max_scratch_calculation_only )
                {
                    break;
                }
            }
        }
        else
        {
            decoder->epmr = 12;
            out_bfi       = 0;
            decoder->error_report = 0;

            for (ch = 0; ch < decoder->channels; ch++)
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                Word32 tmp = num_bytes_total >> (decoder->channels-1);
                fec_num_bytes = tmp + ( ch < ( num_bytes_total - (tmp<<(decoder->channels-1)) ));
#else
                fec_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));
#endif

                BASOP_sub_start("fec_dec");

                channel_bfi = bfi;

                Word32 chan_error_report =
                    fec_decoder(input, fec_num_bytes, &lc3_num_bytes, &channel_epmr, decoder->combined_channel_coding,
                                &decoder->n_pccw, &channel_bfi, &decoder->be_bp_left, &decoder->be_bp_right,
                                &decoder->n_pc, &decoder->m_fec, scratch);

                if (chan_error_report < 0 || decoder->error_report < 0) {
                    decoder->error_report = -1; move16();
                } else {
                    decoder->error_report = add(decoder->error_report, chan_error_report);
                }

                BASOP_sub_end();

                decoder->epmr = MIN(decoder->epmr, channel_epmr);


#ifdef ENABLE_PADDING
                if (channel_bfi != 1)
                {
                    Word16 padding_len, np_zero;

                    if (paddingDec_fx(input, shl(lc3_num_bytes, 3), decoder->yLen, decoder->BW_cutoff_bits,
                                      decoder->ep_enabled, &padding_len, &np_zero))
                    {
                        channel_bfi = 1;
                    }

                    if (input != NULL)
                    {
                        input     = input + np_zero;
                    }
                    decoder->n_pc = s_max(decoder->n_pc - (2 * np_zero), 0);

                    if (channel_bfi == 2)
                    {
                        if (decoder->be_bp_right < (8 * np_zero))
                        {
                            channel_bfi          = 0;
                            decoder->be_bp_left  = -1;
                            decoder->be_bp_right = -1;
                        }
                        else
                        {
                            decoder->be_bp_right = decoder->be_bp_right - (8 * np_zero);
                            decoder->be_bp_left  = s_max(decoder->be_bp_left - (8 * np_zero), 0);
                        }
                    }

                    lc3_num_bytes = lc3_num_bytes - padding_len;
                }
#endif

                if (channel_bfi != 1 && lc3_num_bytes != decoder->channel_setup[ch]->last_size)
                {
                    err = update_dec_bitrate(decoder, ch, lc3_num_bytes);
                    if (err)
                    {
                        channel_bfi = 1;
                    }
                    else
                    {
                        decoder->channel_setup[ch]->last_size = lc3_num_bytes;
                    }
                }

                if ( scratch->max_scratch_calculation_only )
                {
                    channel_bfi = bfi_ext;
                }

                channel_bfi = Dec_LC3PLUS_Channel(decoder, ch, bits_per_sample, input, output[ch], channel_bfi, scratch);

                out_bfi |= channel_bfi;
                if (input != NULL)
                {
                    input += fec_num_bytes;
                }

                if ( scratch->max_scratch_calculation_only )
                {
                    break;
                }
            }

            bfi = out_bfi & 1;
        }
    }
    else
    {
        for (ch = 0; ch < decoder->channels; ch++)
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            lc3_num_bytes = num_bytes[ch];
#else
            lc3_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));
#endif

#if defined(CR14_A_ADD_LOSSLESS_MODE)
            IF (decoder->lossless == 0)
#endif
#ifdef ENABLE_PADDING
            if (bfi != 1)
            {
                Word16 padding_len, np_zero;

                if (paddingDec_fx(input, shl(lc3_num_bytes, 3), decoder->yLen, decoder->BW_cutoff_bits,
                                  decoder->ep_enabled, &padding_len, &np_zero))
                {
                    bfi = 1;
                }

                lc3_num_bytes = lc3_num_bytes - padding_len;
                if (lc3_num_bytes < 20 || lc3_num_bytes > LC3PLUS_MAX_BYTES)
                {
                    bfi = 1; /* mark frame as broken if frame sizeif below the minimum of 20 bytes */
                }
            }
#endif

            if (bfi != 1 && lc3_num_bytes != decoder->channel_setup[ch]->last_size)
            {
                err = update_dec_bitrate(decoder, ch, lc3_num_bytes);
                if (err)
                {
                    bfi = 1;
                }
                else
                {
                    decoder->channel_setup[ch]->last_size = lc3_num_bytes;
                }
            }

            if ( scratch->max_scratch_calculation_only )
            {
                bfi = bfi_ext;
            }

            bfi = Dec_LC3PLUS_Channel(decoder, ch, bits_per_sample, input, output[ch], bfi, scratch);
            if (input != NULL)
            {
                input += decoder->channel_setup[ch]->targetBytes;
            }

            if ( scratch->max_scratch_calculation_only )
            {
                break;
            }
        }
    }

    return bfi == 1 ? LC3PLUS_DECODE_ERROR : LC3PLUS_OK;
}
