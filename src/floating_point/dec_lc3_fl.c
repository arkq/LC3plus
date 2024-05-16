/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static int Dec_LC3PLUS_Channel_fl(LC3PLUS_Dec* decoder, int channel, uint8_t* bs_in, void* s_out, int bps, int bfi_ext)
{
    DecSetup* h_DecSetup;
    LC3_INT       mask_side = 0, bp_side = 0, bfi = 0, gg_idx = 0, fac_ns_idx = 0, tns_numfilters = 0, bw_cutoff_idx = 0,
        lastnz = 0, lsbMode = 0, nf_seed = 0, zero_frame = 0, residualPresent = 0, nbits_residual = 0, bitsRead = 0,
        i = 0, tns_order[2] = {0}, sqQdec[MAX_LEN] = {0};
    LC3_INT b_left;
    LC3_FLOAT stab_fac = 0;

    h_DecSetup = decoder->channel_setup[channel];
    
    memset(h_DecSetup->tns_idx, 0, sizeof(*h_DecSetup->tns_idx) * TNS_NUMFILTERS_MAX * MAXLAG);

    bfi = bfi_ext;

    decoder->rframe = 0;
    if (bfi == 3)
    {
        bfi = 2;
        decoder->rframe = 1;
    }

    /* Entropy decoding */
    if (bfi != 1) {
        processDecoderEntropy_fl(bs_in, h_DecSetup->targetBytes, &mask_side, &bp_side, decoder->yLen, decoder->fs_idx,
                                 decoder->BW_cutoff_bits, &bfi, &gg_idx, h_DecSetup->scf_idx, &fac_ns_idx,
                                 &tns_numfilters, tns_order, h_DecSetup->ltpf_param, &bw_cutoff_idx, &lastnz, &lsbMode, decoder->frame_dms
        );
        h_DecSetup->BW_cutoff_idx_nf = bw_cutoff_idx;
    }
    
    /* Arithmetic decoding */
    if (bfi != 1) {
        processAriDecoder_fl(bs_in, bp_side, mask_side, decoder->yLen, decoder->fs_idx,
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &bfi, tns_order, fac_ns_idx, gg_idx, h_DecSetup->resBits,
                             sqQdec, &nf_seed, h_DecSetup->tns_idx, &zero_frame, h_DecSetup->targetBytes, &nbits_residual, &residualPresent, decoder->frame_dms,
                             decoder->n_pc, decoder->be_bp_left, decoder->be_bp_right, 0, &b_left, &h_DecSetup->spec_inv_idx,
                             decoder->hrmode
        );
        
        if (decoder->rframe == 1 && zero_frame == 0 && bfi != 1)
        {
            LC3_INT32 max_bw_stopband = BW_cutoff_bin_all[bw_cutoff_idx];
            bfi = 2;
            switch (decoder->frame_dms)
            {
            case 25:
                max_bw_stopband  = max_bw_stopband >> 2;
                break;
            case 50:
                max_bw_stopband  = max_bw_stopband >> 1;
                break;
            case 75:
                max_bw_stopband = 3 * (max_bw_stopband >> 2);
                break;
            case 100:
                break;
            }
            
            h_DecSetup->spec_inv_idx = MAX(lastnz, max_bw_stopband);
        }

        /* Cast from int to float */
        for (i = 0; i < decoder->yLen; i++) {
            h_DecSetup->sqQdec_fl[i] = (LC3_FLOAT)sqQdec[i];
        }
    }
    
    if (bfi != 1)
    {
        /* SNS Quantize Decoder */
        process_snsQuantizesScf_Dec(h_DecSetup->scf_idx, h_DecSetup->scf_q);
    }
    if (h_DecSetup->PlcAdvSetup)
    {
        processPlcComputeStabFacMain_fl(h_DecSetup->scf_q, h_DecSetup->PlcAdvSetup->scf_q_old, h_DecSetup->PlcAdvSetup->scf_q_old_old, bfi, h_DecSetup->PlcSetup.prevBfi,
                                        h_DecSetup->PlcSetup.prevprevBfi, &h_DecSetup->PlcAdvSetup->stabFac);
    }
    
    if ( bfi != 1 )
    {
        stab_fac = 1;
        if (h_DecSetup->PlcAdvSetup)
        {
            stab_fac = h_DecSetup->PlcAdvSetup->stabFac;
        }

        /* Partial Concealment */
        processPcMain_fl(&bfi, decoder, h_DecSetup->sqQdec_fl, h_DecSetup, h_DecSetup->ltpf_param[0], stab_fac, gg_idx, h_DecSetup->quantizedGainOff,
             fac_ns_idx, &h_DecSetup->statePC, h_DecSetup->spec_inv_idx, decoder->yLen);
    }

    /* Decoding only if no bit error detected */
    if (bfi != 1) {
        /* Residual decoding */
        if (residualPresent) {
            processResidualDecoding_fl(&bitsRead, h_DecSetup->sqQdec_fl, decoder->yLen, h_DecSetup->resBits,
                                       nbits_residual
                                       , decoder->hrmode
            );
        }
        
        
        /* Noise filling */
        if (zero_frame == 0) {
            processNoiseFilling_fl(h_DecSetup->sqQdec_fl, nf_seed, fac_ns_idx, decoder->cutoffBins[h_DecSetup->BW_cutoff_idx_nf], decoder->frame_dms, h_DecSetup->prev_fac_ns, h_DecSetup->spec_inv_idx);
        }
        
        /* Application of global gain */
        processApplyGlobalGain_fl(h_DecSetup->sqQdec_fl, decoder->yLen, gg_idx, h_DecSetup->quantizedGainOff);

        /* TNS decoder */
        processTnsDecoder_fl(h_DecSetup->sqQdec_fl, h_DecSetup->tns_idx, tns_order, tns_numfilters,
                             decoder->cutoffBins[bw_cutoff_idx], h_DecSetup->N_red_tns, h_DecSetup->fs_red_tns);

        /* SNS interpolation */
        processSnsInterpolateScf_fl(h_DecSetup->scf_q, 0, decoder->bands_number, h_DecSetup->int_scf);

        /* MDCT shaping */
        processMdctShaping_fl(h_DecSetup->sqQdec_fl, h_DecSetup->int_scf, decoder->bands_offset, decoder->bands_number);
    }
    
    /* PLC */
    processPlcMain_fl(h_DecSetup->sqQdec_fl, h_DecSetup->x_fl, decoder, h_DecSetup, bfi, h_DecSetup->PlcAdvSetup, &h_DecSetup->PlcSetup,
              decoder->plcMeth, h_DecSetup->ltpf_mem_pitch, h_DecSetup->ltpf_mem_pitch_fr, decoder->tilt, decoder->bands_offset,
              decoder->bands_number, decoder->bands_offsetPLC, decoder->n_bandsPLC, decoder->hrmode, &h_DecSetup->statePC);

    processPlcDampingScramblingMain_fl(&h_DecSetup->PlcNsSetup.seed,
                                       &h_DecSetup->statePC.seed, h_DecSetup->statePC.ns_nbLostCmpt_pc,
                                       h_DecSetup->PlcSetup.nbLostCmpt, &h_DecSetup->PlcAdvSetup->stabFac,
                                       &h_DecSetup->PlcAdvSetup->cum_fading_slow, &h_DecSetup->PlcAdvSetup->cum_fading_fast,
                                       h_DecSetup->PlcSetup.q_d_prev, h_DecSetup->sqQdec_fl, h_DecSetup->spec_inv_idx, decoder->yLen, bfi,
                                       decoder->frame_dms, h_DecSetup->concealMethod, h_DecSetup->ltpf_mem_pitch, h_DecSetup->ltpf_param[0],
                                       &h_DecSetup->PlcAdvSetup->cum_fflcAtten
                                       , h_DecSetup->PlcAdvSetup->plc_fadeout_type
                                      );
    
    /* IMDCT */
    if (h_DecSetup->concealMethod == 4 || bfi != 1 )
    {
        ProcessingIMDCT_fl(h_DecSetup->sqQdec_fl, decoder->frame_length, decoder->imdct_win, decoder->imdct_winLen, decoder->imdct_laZeros,
                           h_DecSetup->imdct_mem, h_DecSetup->x_fl, &h_DecSetup->dct4structImdct);
    }

    processPlcUpdate_fl(h_DecSetup->PlcAdvSetup
                        , decoder->frame_length, h_DecSetup->x_fl, h_DecSetup->scf_q,
                        &h_DecSetup->PlcSetup.nbLostCmpt, &h_DecSetup->PlcNsSetup.cum_alpha, bfi, &h_DecSetup->PlcSetup.prevBfi, &h_DecSetup->PlcSetup.prevprevBfi);
    
    /* LTPF decoder */
    process_ltpf_decoder_fl(h_DecSetup->x_fl, decoder->frame_length, h_DecSetup->x_fl, decoder->fs,
                            h_DecSetup->ltpf_mem_x, h_DecSetup->ltpf_mem_y, &h_DecSetup->ltpf_mem_pitch,
                            &h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ltpf_mem_gain, &h_DecSetup->ltpf_mem_beta_idx,
                            bfi, h_DecSetup->ltpf_param, h_DecSetup->ltpf_param_mem, h_DecSetup->ltpf_conf_beta_idx,
                            h_DecSetup->ltpf_conf_beta, h_DecSetup->concealMethod, h_DecSetup->alpha
                            , &h_DecSetup->ltpf_mem_active
                            , &h_DecSetup->rel_pitch_change, decoder->hrmode, decoder->frame_dms
                           );

    {
        /* Round, scale and copy output to output buffer */
        if (bps == 16) {
            for (i = 0; i < decoder->frame_length; i++) {
                LC3_FLOAT tmp        = round(h_DecSetup->x_fl[i]);
                ((int16_t*)s_out)[i] = (int16_t)fmaxf(fminf(tmp, 32767), -32768);
            }
        } else {
            for (i = 0; i < decoder->frame_length; i++) {
                LC3_FLOAT tmp = round(LC3_CONST_POW_2_23 * LC3_CONST_POW_2_M15 * h_DecSetup->x_fl[i]);
                ((int32_t*)s_out)[i] = (int32_t)fmaxf(fminf(tmp, LC3_CONST_POW_2_23_RED), LC3_CONST_POW_2_23_NEG);
            }
        }
    }
    return bfi;
}

LC3PLUS_Error Dec_LC3PLUS_fl(LC3PLUS_Dec* decoder, uint8_t* input, LC3_INT32 num_bytes, void** output, LC3_INT32 bps, LC3_INT32 bfi_ext)
{
        LC3_INT32       ch, bfi, lc3_num_bytes;
        LC3PLUS_Error err;
        LC3_INT32       fec_num_bytes;
        LC3_INT32       lc3_channel_num_bytes;
        LC3_INT32       channel_bfi, out_bfi;
        LC3PLUS_EpModeRequest channel_epmr;
    
    bfi = bfi_ext;
    lc3_num_bytes = 0;
    err = LC3PLUS_OK;
    
    if (bfi == 0)
    {
        bfi = !num_bytes;
    }
    
    if (decoder->ep_enabled)
    {
        decoder->combined_channel_coding = decoder->channels > 1 && num_bytes <= 160;

        if (decoder->combined_channel_coding)
        {
            fec_num_bytes = num_bytes;

            decoder->error_report =
                fec_decoder(input, fec_num_bytes, &lc3_num_bytes, (LC3PLUS_EpModeRequest*)&decoder->epmr, decoder->combined_channel_coding,
                            &decoder->n_pccw, &bfi, &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc, &decoder->m_fec);

            for (ch = 0; ch < decoder->channels; ch++)
            {
                lc3_channel_num_bytes = lc3_num_bytes / decoder->channels + (ch < (lc3_num_bytes % decoder->channels));


                if (bfi != 1 && lc3_channel_num_bytes != decoder->channel_setup[ch]->last_size)
                {
                    err = update_dec_bitrate(decoder, ch, lc3_channel_num_bytes);
                    if (err)
                    {
                        bfi = 1;
                        decoder->last_error = err;
                    }
                    else
                    {
                        decoder->channel_setup[ch]->last_size = lc3_channel_num_bytes;                        
                    }
                }

                bfi = Dec_LC3PLUS_Channel_fl(decoder, ch, input, output[ch], bps, bfi);
                if (input != NULL)
                {
                    input += decoder->channel_setup[ch]->targetBytes;
                }
            }
        }
        else
        {
            decoder->epmr = LC3PLUS_EPMR_HIGH_NC;
            out_bfi       = 0;

            for (ch = 0; ch < decoder->channels; ch++)
            {
                fec_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));

                channel_bfi = bfi;              

                decoder->error_report = fec_decoder(input, fec_num_bytes, &lc3_num_bytes, &channel_epmr,
                                                    decoder->combined_channel_coding, &decoder->n_pccw, &channel_bfi,
                                                    &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc, &decoder->m_fec);

                decoder->epmr = MIN((LC3PLUS_EpModeRequest) decoder->epmr, channel_epmr);

    
#ifdef ENABLE_PADDING
                if (channel_bfi != 1)
                {
                    LC3_INT32 padding_len = 0, np_zero = 0;

                    if (paddingDec_fl(input, (lc3_num_bytes << 3), decoder->yLen, decoder->BW_cutoff_bits, decoder->ep_enabled, &padding_len, &np_zero))
                    {
                        channel_bfi = 1;
                    }

                    if (input != NULL)
                    {
                        input = input + np_zero;
                    }
                    
                    decoder->n_pc = MAX(decoder->n_pc - (2 * np_zero), 0);
                    
                    if (channel_bfi == 2)
                    {
                        if (decoder->be_bp_right < (8 * np_zero))
                        {
                            channel_bfi = 0;
                            decoder->be_bp_left = -1;
                            decoder->be_bp_right = -1;
                        }
                        else
                        {
                            decoder->be_bp_right = decoder->be_bp_right - (8 * np_zero);
                            decoder->be_bp_left  = MAX(decoder->be_bp_left - (8 * np_zero), 0);
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
                        decoder->last_error = err;
                    }
                    else
                    {
                        decoder->channel_setup[ch]->last_size = lc3_num_bytes;
                    }
                }

                channel_bfi = Dec_LC3PLUS_Channel_fl(decoder, ch, input, output[ch], bps, channel_bfi);

                out_bfi |= channel_bfi;
                if (input != NULL)
                {
                    input += fec_num_bytes;
                }
            }

            bfi = out_bfi & 1;
        }
    }
    else
    {
        for (ch = 0; ch < decoder->channels; ch++)
        {
            lc3_num_bytes = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));

#ifdef ENABLE_PADDING
            if (bfi != 1)
            {
                LC3_INT32 padding_len = 0, np_zero = 0;

                if (paddingDec_fl(input, (lc3_num_bytes << 3), decoder->yLen, decoder->BW_cutoff_bits, decoder->ep_enabled, &padding_len, &np_zero))
                {
                    bfi = 1;
                    decoder->last_error = LC3PLUS_PADDING_ERROR;
                }
                
                lc3_num_bytes = lc3_num_bytes - padding_len;
                if (lc3_num_bytes < 20 || lc3_num_bytes > LC3PLUS_MAX_BYTES)
                {
                    bfi = 1;    /* mark frame as broken if frame size is below the minimum of 20 bytes or above the maximum of LC3PLUS_MAX_BYTES */
                    decoder->last_error = FRAMESIZE_ERROR;
                }
            }
#endif 
            
            if (bfi != 1 && lc3_num_bytes != decoder->channel_setup[ch]->last_size)
            {
                err = update_dec_bitrate(decoder, ch, lc3_num_bytes);
                if (err)
                {
                    bfi = 1;
                    decoder->last_error = err;
                }
                else
                {
                    decoder->channel_setup[ch]->last_size = lc3_num_bytes;
                }
            }

            bfi = Dec_LC3PLUS_Channel_fl(decoder, ch, input, output[ch], bps, bfi);
            if (input != NULL)
            {
                input += decoder->channel_setup[ch]->targetBytes;
            }
        }
    }

    if ((decoder->last_error == LC3PLUS_OK) && bfi)
    {
        decoder->last_error = LC3PLUS_DECODE_ERROR;
    }
    return bfi == 1 ? LC3PLUS_DECODE_ERROR : LC3PLUS_OK;
}
