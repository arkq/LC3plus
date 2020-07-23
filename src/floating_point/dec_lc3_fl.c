/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

static int Dec_LC3_Channel_fl(LC3_Dec* decoder, int channel, uint8_t* bs_in, void* s_out, int bps, int bfi_ext)
{
    DecSetup* h_DecSetup;
    LC3_INT       mask_side = 0, bp_side = 0, bfi = 0, gg_idx = 0, fac_ns_idx = 0, tns_numfilters = 0, bw_cutoff_idx = 0,
        lastnz = 0, lsbMode = 0, nf_seed = 0, zero_frame = 0, residualPresent = 0, nbits_residual = 0, bitsRead = 0,
        i = 0, tns_order[2] = {0}, sqQdec[MAX_LEN] = {0};
    

    h_DecSetup = decoder->channel_setup[channel];

    bfi = bfi_ext;


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
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &bfi, tns_order, fac_ns_idx, gg_idx, h_DecSetup->resBits, sqQdec, &nf_seed, h_DecSetup->tns_idx, &zero_frame, h_DecSetup->targetBytes, &nbits_residual, &residualPresent, decoder->frame_dms
                             ,
                             decoder->hrmode
        );

        
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

    

    /* Decoding only if no bit error detected */
    if (bfi != 1) {
        /* Residual decoding */
        if (residualPresent) {
            processResidualDecoding_fl(&bitsRead, h_DecSetup->sqQdec_fl, decoder->yLen, h_DecSetup->resBits,
                                       nbits_residual
									   , decoder->hrmode
            );
        }
        
   h_DecSetup->spec_inv_idx = decoder->yLen;     
        
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
    processPlcMain_fl(h_DecSetup->sqQdec_fl, h_DecSetup->x_fl, decoder, h_DecSetup, bfi
                      , &h_DecSetup->PlcSetup,
                      &h_DecSetup->PlcNsSetup, decoder->plcMeth, h_DecSetup->ltpf_mem_pitch, h_DecSetup->ltpf_mem_pitch_fr, decoder->tilt,
                      decoder->bands_offset, decoder->bands_number, decoder->bands_offsetPLC, decoder->n_bandsPLC
                      );
    
    /* IMDCT */
    if (h_DecSetup->concealMethod == 0 || h_DecSetup->concealMethod == 4 || h_DecSetup->concealMethod == 5 || bfi != 1 )
    {
        ProcessingIMDCT_fl(h_DecSetup->sqQdec_fl, decoder->frame_length, decoder->imdct_win, decoder->imdct_winLen, decoder->imdct_laZeros,
                           h_DecSetup->imdct_mem, h_DecSetup->x_fl, &h_DecSetup->dct4structImdct);
    }

    processPlcUpdate_fl(decoder->plcMeth
                        , decoder->frame_length, h_DecSetup->x_fl, h_DecSetup->scf_q,
                        &h_DecSetup->PlcSetup.nbLostCmpt, &h_DecSetup->PlcNsSetup.cum_alpha, bfi, &h_DecSetup->PlcSetup.prevBfi, &h_DecSetup->PlcSetup.prevprevBfi);
    
    /* LTPF decoder */
    process_ltpf_decoder_fl(h_DecSetup->x_fl, decoder->frame_length, h_DecSetup->x_fl, decoder->fs,
                            h_DecSetup->ltpf_mem_x, h_DecSetup->ltpf_mem_y, &h_DecSetup->ltpf_mem_pitch,
                            &h_DecSetup->ltpf_mem_pitch_fr, &h_DecSetup->ltpf_mem_gain, &h_DecSetup->ltpf_mem_beta_idx,
                            bfi, h_DecSetup->ltpf_param, h_DecSetup->ltpf_param_mem, h_DecSetup->ltpf_conf_beta_idx,
                            h_DecSetup->ltpf_conf_beta, h_DecSetup->concealMethod
                            , h_DecSetup->alpha
                           );

    {
        /* Round, scale and copy output to output buffer */
        if (bps == 16) {
            for (i = 0; i < decoder->frame_length; i++) {
                LC3_FLOAT tmp            = round(LC3_POW(2, 16 - 1) * (h_DecSetup->x_fl[i] * LC3_POW(2, -15)));
                ((int16_t*)s_out)[i] = (int16_t)fmaxf(fminf(tmp, 32767), -32768);
            }
        } else {
            for (i = 0; i < decoder->frame_length; i++) {
                ((int32_t*)s_out)[i] = (int32_t)round(LC3_POW(2, bps - 1) * (h_DecSetup->x_fl[i] * LC3_POW(2, -15)));
            }
        }
    }
    return bfi;
}


/* num_bytes = 0 -> bad frame */
LC3_Error Dec_LC3_fl(LC3_Dec* decoder, uint8_t* input, int num_bytes, void** output, int bps, int bfi_ext)
{
    int ch = 0, bfi = bfi_ext;
    LC3_Error err = LC3_OK;
    int num_bytes2;

    if (bfi == 0)
    {
        bfi = !num_bytes;
    }

    for (ch = 0; ch < decoder->channels; ch++)
    {
        num_bytes2 = num_bytes / decoder->channels + (ch < (num_bytes % decoder->channels));
        if (bfi != 1 && num_bytes2 != decoder->channel_setup[ch]->last_size)
        {
            err = update_dec_bitrate(decoder, ch, num_bytes2);
            if (err)
                return err;
            decoder->channel_setup[ch]->last_size = num_bytes2;
        }
        
        bfi = Dec_LC3_Channel_fl(decoder, ch, input, output[ch], bps, bfi);
        input += decoder->channel_setup[ch]->targetBytes;
    }

    return bfi == 1 ? LC3_DECODE_ERROR : LC3_OK;
}
