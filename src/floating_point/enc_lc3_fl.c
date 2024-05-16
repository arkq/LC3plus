/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                            
#include "functions.h"

static void Enc_LC3PLUS_Channel_fl(LC3PLUS_Enc* encoder, int channel, int32_t* s_in, uint8_t* bytes, int bps
, LC3_INT32 bfi_ext
)
{
    EncSetup* h_EncSetup;

    LC3_INT s_12k8_len = 0, T0_out = 0, ltpfBits = 0, BW_cutoff_idx = 0, tns_numfilters = 0, quantizedGain = 0,
        quantizedGainMin = 0, nbits = 0, nbits2 = 0, lastnz = 0, lsbMode = 0, gainChange = 0, bp_side = 0,
        mask_side = 0, fac_ns_idx = 0, numResBits = 0, tns_order[2] = {0}, i = 0;
    LC3_FLOAT normcorr = 0, gain = 0;
    

    LC3_FLOAT d_fl[MAX_LEN]                        = {0};
    LC3_INT   q_d[MAX_LEN]                         = {0};
    LC3_INT   indexes[TNS_NUMFILTERS_MAX * MAXLAG] = {0};

    h_EncSetup = encoder->channel_setup[channel];
    memset(bytes, 0, sizeof(uint8_t) * h_EncSetup->targetBytes);

    if (bps == 24) {
        for (i = 0; i < encoder->frame_length; i++) {
            int32_t tmp = ((int32_t*)s_in)[i];

            if (tmp >= 0)
            {
                tmp = tmp & 0x007fffff;
            }
            else
            {
                tmp = tmp | (int32_t)0xff800000;
            }

            h_EncSetup->s_in_scaled[i] = ((LC3_FLOAT) tmp / (float) LC3_POW(2, 8));
        }
    } else if (bps == 16) {
        for (i = 0; i < encoder->frame_length; i++) {
            h_EncSetup->s_in_scaled[i] = (LC3_FLOAT)((int16_t*)s_in)[i];
        }
    }

    /* MDCT */
    processMdct_fl(h_EncSetup->s_in_scaled, d_fl, &h_EncSetup->mdctStruct);

    /* 12.8 kHz resampler */
    process_resamp12k8_fl(h_EncSetup->s_in_scaled, encoder->frame_length, h_EncSetup->r12k8_mem_in,
                          encoder->r12k8_mem_in_len, h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out,
                          encoder->r12k8_mem_out_len, h_EncSetup->s_12k8, &s_12k8_len, encoder->fs_idx,
                          encoder->frame_dms, encoder->fs);

    /* Pitch estimation */
    processOlpa_fl(h_EncSetup->s_12k8, h_EncSetup->olpa_mem_s12k8, h_EncSetup->olpa_mem_s6k4,
                   &h_EncSetup->olpa_mem_pitch, 
                   &h_EncSetup->pitch_flag, 
                   &T0_out, &normcorr, s_12k8_len, encoder->frame_dms);

    /* LTPF encoder */
    process_ltpf_coder_fl(h_EncSetup->s_12k8, s_12k8_len + 1, h_EncSetup->ltpf_enable, T0_out, normcorr,
                          encoder->frame_dms, h_EncSetup->ltpf_mem_in, encoder->ltpf_mem_in_len,
                          &h_EncSetup->ltpf_mem_normcorr, &h_EncSetup->ltpf_mem_ltpf_on,
                          &h_EncSetup->ltpf_mem_pitch, h_EncSetup->ltpf_param, &h_EncSetup->ltpf_mem_mem_normcorr,
                          &ltpfBits
                          , encoder->hrmode
);

    /* Attack detector */
    attack_detector_fl(h_EncSetup->s_in_scaled, encoder->frame_length, encoder->fs, &h_EncSetup->attdec_position,
                       &h_EncSetup->attdec_acc_energy, &h_EncSetup->attdec_detected, h_EncSetup->attdec_filter_mem,
                       h_EncSetup->attack_handling, encoder->attdec_nblocks, encoder->attdec_hangover_thresh);

    /* Per-band energy */
    processPerBandEnergy_fl(encoder->bands_number, encoder->bands_offset, encoder->hrmode, encoder->frame_dms, h_EncSetup->ener, d_fl);
    /* Near Nyquist detector */
    processNearNyquistdetector_fl(&encoder->near_nyquist_flag, encoder->fs_idx, encoder->near_nyquist_index, encoder->bands_number, h_EncSetup->ener                              
                                 , encoder->frame_dms, encoder->hrmode );                            
    /* Disable LTPF if nyquist detector triggers or -lfe mode is active*/
    if (encoder->near_nyquist_flag != 0 || h_EncSetup->lfe == 1)
    {
        h_EncSetup->ltpf_mem_ltpf_on = 0;
        h_EncSetup->ltpf_param[1] = 0;
    }

    /* Bandwidth cut-off detection */
    if (h_EncSetup->lfe == 0) {
    /* No BW Cutoff for 8 kHz and 96 kHz. No detection if bandwidth controller is active */
    if (encoder->fs_idx > 0 && encoder->hrmode == 0 && encoder->bw_ctrl_active == 0) {
        processDetectCutoffWarped_fl(h_EncSetup->ener, encoder->fs_idx, encoder->frame_dms, &BW_cutoff_idx);
    } else {
        BW_cutoff_idx = encoder->fs_idx;
    }
    } else {
        BW_cutoff_idx = 0;
    }

    processSnsComputeScf_fl(h_EncSetup->ener, encoder->bands_number, h_EncSetup->scf,
                            h_EncSetup->attdec_detected, encoder->sns_damping, encoder->attdec_damping, encoder->fs_idx);

    /* SNS Quantizer */
    process_snsQuantizesScf_Enc(h_EncSetup->scf, h_EncSetup->L_scf_idx, h_EncSetup->scf_q, h_EncSetup->dct2StructSNS);

    /* SNS Interpolation */
    processSnsInterpolateScf_fl(h_EncSetup->scf_q, 1, encoder->bands_number, h_EncSetup->int_scf);

    /* MDCT shaping */
    processMdctShaping_fl(d_fl, h_EncSetup->int_scf, encoder->bands_offset, encoder->bands_number);
        
    /* Bandwidth controller */
    if (encoder->bandwidth < encoder->fs / 2) {
        process_cutoff_bandwidth(d_fl, encoder->yLen, encoder->bw_ctrl_cutoff_bin);
        BW_cutoff_idx = MIN(BW_cutoff_idx, encoder->bw_index);
    }
        
    /* TNS encoder */
    if (h_EncSetup->lfe == 0)
    {
    processTnsCoder_fl(d_fl, BW_cutoff_idx, encoder->cutoffBins[BW_cutoff_idx], encoder->fs, encoder->frame_length,
                       encoder->frame_dms, h_EncSetup->total_bits, tns_order, indexes, &tns_numfilters,
                       &(h_EncSetup->tns_bits)
                       , encoder->near_nyquist_flag
    );
    }
    else
    {
        tns_numfilters = 1;
        tns_order[0] = 0;
        h_EncSetup->tns_bits = tns_numfilters;
    }
    /* Global Gain Estimation */
    h_EncSetup->targetBitsQuant = h_EncSetup->targetBitsInit - (h_EncSetup->tns_bits + ltpfBits);

    if (h_EncSetup->targetBitsQuant < 0 && ltpfBits > 1)
    {
        /* Disable LTPF */
        h_EncSetup->ltpf_mem_ltpf_on = 0;
        h_EncSetup->ltpf_param[1]    = 0;
        ltpfBits                     = 1;
        h_EncSetup->targetBitsQuant  = h_EncSetup->targetBitsInit - (h_EncSetup->tns_bits + ltpfBits);
    }

    processEstimateGlobalGain_fl(d_fl, encoder->yLen, h_EncSetup->targetBitsQuant, &gain, &quantizedGain,
                                 &quantizedGainMin, h_EncSetup->quantizedGainOff, &h_EncSetup->targetBitsOff,
                                 &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits
                                 , encoder->hrmode, h_EncSetup->regBits, encoder->frame_ms

    );

    /* 1. Quantization */
    processQuantizeSpec_fl(d_fl, gain, q_d, encoder->yLen, h_EncSetup->total_bits, &nbits, &nbits2, encoder->fs,
                           &lastnz, h_EncSetup->codingdata, &lsbMode, -1, h_EncSetup->targetBitsQuant, encoder->hrmode
    );

    h_EncSetup->mem_specBits = nbits;

    /* Global Gain Adjustment */
    processAdjustGlobalGain_fl(&quantizedGain, quantizedGainMin, h_EncSetup->quantizedGainOff, &gain,
                               h_EncSetup->targetBitsQuant, h_EncSetup->mem_specBits, &gainChange, encoder->fs_idx
                               , encoder->hrmode, encoder->frame_dms
                               );

    /* 2. Quantization */
    if (gainChange) {
        processQuantizeSpec_fl(d_fl, gain, q_d, encoder->yLen, h_EncSetup->total_bits, &nbits, &nbits2, encoder->fs,
                               &lastnz,
                               h_EncSetup->codingdata, 
                               &lsbMode, 0, h_EncSetup->targetBitsQuant
        , encoder->hrmode
        );
    }

    /* Noise factor */
    if (h_EncSetup->lfe == 0)
    {    
    processNoiseFactor_fl(&fac_ns_idx, d_fl, q_d, gain, encoder->cutoffBins[BW_cutoff_idx], encoder->frame_dms,
                          h_EncSetup->targetBytes
    );
    }
    else
    {
        fac_ns_idx = 7;
    }
    /* Residual Coding */
    if (lsbMode == 0) {
        processResidualCoding_fl(d_fl, q_d, gain, encoder->yLen, h_EncSetup->targetBitsQuant, nbits2,
                                 h_EncSetup->resBits, &numResBits
                                 , encoder->hrmode
        );
    } else {
        numResBits = 0;
    }

    /* Entropy encoding */
    processEncoderEntropy_fl(bytes, &bp_side, &mask_side, h_EncSetup->targetBytes, encoder->BW_cutoff_bits,
                             BW_cutoff_idx, lastnz, encoder->yLen, lsbMode, quantizedGain, tns_numfilters, tns_order,
                             h_EncSetup->ltpf_param, h_EncSetup->L_scf_idx, fac_ns_idx
                              , bfi_ext, encoder->fs_idx
    );

    /* Artithmetic encoding */
    processAriEncoder_fl(bytes, bp_side, mask_side, q_d, tns_order, tns_numfilters, indexes, lastnz,
                         h_EncSetup->codingdata, 
                         h_EncSetup->resBits, numResBits, lsbMode, h_EncSetup->targetBitsAri,
                         h_EncSetup->enable_lpc_weighting);
    
    if (encoder->combined_channel_coding == 0 && h_EncSetup->n_pc > 0)
    {
        LC3_INT32 xbuf[MAX_LEN] = {0}, nf_seed, tns_idx[M], zero_frame, nbits_residual, residualPresent, b_left, spec_inv_idx;
   
        memset(h_EncSetup->resBits, 0, sizeof(*(h_EncSetup->resBits)) * MAX_RESBITS_LEN);
        
        processAriDecoder_fl(bytes, bp_side, mask_side, encoder->yLen, encoder->fs_idx, h_EncSetup->enable_lpc_weighting,
                 tns_numfilters, lsbMode, lastnz, &bfi_ext, tns_order, fac_ns_idx, quantizedGain,
                 h_EncSetup->resBits, xbuf, &nf_seed, tns_idx, &zero_frame, h_EncSetup->targetBytes, &nbits_residual,
                 &residualPresent, encoder->frame_dms, h_EncSetup->n_pc, 0, h_EncSetup->total_bits >> 3, 1, &b_left,
                 &spec_inv_idx, encoder->hrmode);
        
        processReorderBitstream_fl(bytes, h_EncSetup->n_pccw, h_EncSetup->n_pc, b_left, h_EncSetup->targetBytes);
    }
        
}

int Enc_LC3PLUS_fl(LC3PLUS_Enc* encoder, void** input, uint8_t* output, int bps
, LC3_INT32 bfi_ext
)
{
    int      ch = 0, output_size = 0;
    uint8_t* lc3buf = output;

    LC3_INT32 totalBytes;
    LC3_INT32 output_size2, input_size;
    
    totalBytes = encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);

    for (ch = 0; ch < encoder->channels; ch++)
    {
        Enc_LC3PLUS_Channel_fl(encoder, ch, input[ch], lc3buf, bps, bfi_ext);
        if (encoder->epmode && encoder->combined_channel_coding == 0)
        {
            output_size2 = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));

            fec_encoder(encoder->epmode, encoder->epmr, lc3buf, encoder->channel_setup[ch]->targetBytes, output_size2,
                        encoder->channel_setup[ch]->n_pccw);

            lc3buf += output_size2;
            output_size += output_size2;
        }
        else
        {
            lc3buf += encoder->channel_setup[ch]->targetBytes;
            output_size += encoder->channel_setup[ch]->targetBytes;
        }
    }

    if (encoder->epmode > 0 && encoder->combined_channel_coding)
    {
        input_size  = output_size;
        output_size = encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);

        fec_encoder(encoder->epmode, encoder->epmr, output, input_size, output_size, encoder->channel_setup[0]->n_pccw);
    }
    
    return output_size;
}
