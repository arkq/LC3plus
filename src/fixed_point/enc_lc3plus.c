/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include <stdio.h>  // ToDo: probably to be removed

static void Enc_LC3PLUS_Channel(LC3PLUS_Enc *encoder, int channel, int bits_per_sample, Word32 *s_in, UWord8 *bytes,
                            Word8 *scratchBuffer, int bfi_ext)
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
                      Word16 nBits, numResBits;
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
                      Word8 * currentScratch;
                      Word16 hrmode;
                     );
#else
    Dyn_Mem_Deluxe_In(Word16 d_fx_exp; Word16 gain_e, gain, quantizedGain, quantizedGainMin; Word16 ener_fx_exp;
                      Word16 pitch, normcorr; Word16 ltpf_bits; Word16 tns_numfilters;
                      Word16 lsbMode, lastnz, BW_cutoff_idx; Word16 gainChange, fac_ns_idx; Word16 nBits, numResBits;
                      Word16 bp_side, mask_side; Word16 s_12k8_len; Word16 b_left;

                      Word32 * L_scf_idx; Word32 * d_fx, *ener_fx;
                      Word16 * s_12k8, *int_scf_fx_exp, *q_d_fx16, *int_scf_fx, tns_order[TNS_NUMFILTERS_MAX], *indexes;
                      Word16 * scf, *scf_q; Word16 * codingdata; Word16 * s_in_scaled; UWord8 * resBits;
                      Word8 * currentScratch; Word16 ltpf_idx[3]; EncSetup * h_EncSetup;);
#endif /* ENABLE_HR_MODE */
    
    h_EncSetup = encoder->channel_setup[channel];
    
#ifdef ENABLE_HR_MODE  
    hrmode = encoder->hrmode;
#endif

    BASOP_sub_start("Encoder");

    /* BUFFER INITIALISATION. Some buffers may overlap since they are not used in the whole encoding process */
    d_fx      = scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_LEN bytes */
    L_scf_idx = scratchAlign(
        d_fx, sizeof(*d_fx) * s_max(80, encoder->frame_length)); /* Size = 4 * SCF_MAX_PARAM -> aligned to 32 bytes */
    indexes  = scratchAlign(L_scf_idx,
                           sizeof(*L_scf_idx) * SCF_MAX_PARAM); /* Size = 2 * TNS_NUMFILTERS_MAX * MAXLAG = 32 bytes */

#ifdef ENABLE_HR_MODE
    q_d_fx24 = scratchAlign(indexes, sizeof(*indexes) * (TNS_NUMFILTERS_MAX * MAXLAG));
#else
    q_d_fx16 = scratchAlign(indexes, sizeof(*indexes) * (TNS_NUMFILTERS_MAX * MAXLAG)); /* Size = 2 * MAX_LEN bytes */
#endif /* ENABLE_HR_MODE */
    
#ifdef ENABLE_HR_MODE
    codingdata =
        scratchAlign(q_d_fx24, sizeof(*q_d_fx24) * s_max(80, encoder->frame_length)); /* Size = 3 * MAX_LEN bytes */
#else
    codingdata =
        scratchAlign(q_d_fx16, sizeof(*q_d_fx16) * s_max(80, encoder->frame_length)); /* Size = 3 * MAX_LEN bytes */
#endif
    
#ifdef ENABLE_HR_MODE
    ener_fx     = scratchAlign(q_d_fx24, 0); /* Size = 4 * MAX_BANDS_NUMBER = 256 bytes */
    s_in_scaled = scratchAlign(q_d_fx24, 0); /* Size = 2 * MAX_LEN bytes */
#else
    ener_fx     = scratchAlign(q_d_fx16, 0); /* Size = 4 * MAX_BANDS_NUMBER = 256 bytes */
    s_in_scaled = scratchAlign(q_d_fx16, 0); /* Size = 2 * MAX_LEN bytes */
#endif
    
#ifdef ENABLE_HR_MODE
    s_in_scaled_lp = (Word16 *)s_in_scaled;
#endif

#ifdef ENABLE_HR_MODE
    /* allocate memory for residual bits */
    if (encoder->hrmode)
    {
        resBits = scratchAlign(codingdata, sizeof(*codingdata) * (3 * s_max(80, encoder->frame_length) / 2));
    }
    else
#endif /* #ifdef ENABLE_HR_MODE */
    {
        resBits        = scratchAlign(codingdata,
                               sizeof(*codingdata) * (3 * s_max(80, encoder->frame_length) / 2)); /* Size = MAX_LEN bytes */
    }

    {
#ifdef ENABLE_HR_MODE
        currentScratch =
            scratchAlign(resBits, sizeof(*resBits) * MAX_RESBITS_LEN); /* Size = 4 * MAX_LEN */
#else
        currentScratch =
            scratchAlign(resBits, sizeof(*resBits) * s_max(80, encoder->frame_length)); /* Size = 4 * MAX_LEN */
#endif /* #ifdef ENABLE_HR_MODE */
    }

    s_12k8         = scratchAlign(
        s_in_scaled,
        sizeof(*s_in_scaled) *
            s_max(80, encoder->frame_length)); /* Size = 2 * (LEN_12K8 + 1) = 258 bytes -> aligned to 288 bytes */
    scf_q      = scratchAlign(ener_fx, sizeof(*ener_fx) * MAX_BANDS_NUMBER); /* Size = 2 * M */
    scf        = scratchAlign(scf_q, sizeof(*scf_q) * M);                    /* Size = 2 * M */
    int_scf_fx = scratchAlign(scf, sizeof(*scf) * M); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    int_scf_fx_exp =
        scratchAlign(int_scf_fx, sizeof(*int_scf_fx) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */

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

    BASOP_sub_start("Mdct");

    /* currentScratch Size = 4 * MAX_LEN */
    processMdct_fx(s_in_scaled, h_EncSetup->x_exp, encoder->frame_length,
#ifdef ENABLE_HR_MODE
                   hrmode, 
#endif
                   encoder->W_fx, encoder->W_size,
                   h_EncSetup->stEnc_mdct_mem, encoder->stEnc_mdct_mem_len, d_fx, &d_fx_exp, currentScratch);
    BASOP_sub_end();

    /* begin s_12k8 */
    BASOP_sub_start("Resamp12k8");
    /* currentScratch Size = 2.25 * MAX_LEN bytes */
#ifdef ENABLE_HR_MODE
    downshift_w32_arr(s_in_scaled, s_in_scaled_lp, 16, encoder->frame_length); /* s_in_scaled is no longer required */
    
    process_resamp12k8_fx(s_in_scaled_lp, encoder->frame_length, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len, s_12k8,
                          &s_12k8_len, encoder->fs_idx, encoder->frame_dms, currentScratch
                          , bits_per_sample
                         );
#else
    process_resamp12k8_fx(s_in_scaled, encoder->frame_length, h_EncSetup->r12k8_mem_in, encoder->r12k8_mem_in_len,
                          h_EncSetup->r12k8_mem_50, h_EncSetup->r12k8_mem_out, encoder->r12k8_mem_out_len, s_12k8,
                          &s_12k8_len, encoder->fs_idx, encoder->frame_dms, currentScratch
                          , bits_per_sample
                          );
#endif /* ENABLE_HR_MODE */
    BASOP_sub_end();

    BASOP_sub_start("Olpa");
    /* currentScratch Size = 392 bytes */
    process_olpa_fx(&h_EncSetup->olpa_mem_s6k4_exp, h_EncSetup->olpa_mem_s12k8, h_EncSetup->olpa_mem_s6k4, &pitch,
                    s_12k8, s_12k8_len, &normcorr, &h_EncSetup->olpa_mem_pitch, 
                    &h_EncSetup->pitch_flag,                  
                    h_EncSetup->resamp_exp, encoder->frame_dms, currentScratch);
    BASOP_sub_end();

    BASOP_sub_start("LtpfEnc");
    /* currentScratch Size = 512 bytes */
    process_ltpf_coder_fx(&ltpf_bits, pitch, h_EncSetup->ltpf_enable, &h_EncSetup->ltpf_mem_in_exp,
                          h_EncSetup->ltpf_mem_in, encoder->ltpf_mem_in_len, ltpf_idx, s_12k8, s_12k8_len,
                          &h_EncSetup->ltpf_mem_normcorr, &h_EncSetup->ltpf_mem_mem_normcorr, normcorr,
                          &h_EncSetup->ltpf_mem_ltpf_on, &h_EncSetup->ltpf_mem_pitch, h_EncSetup->resamp_exp,
                          encoder->frame_dms, currentScratch
                          , encoder->hrmode
);
    BASOP_sub_end();

    /* end s_12k8 */
    BASOP_sub_start("AttackDetector");
    /* currentScratch Size = ??? bytes */
#ifdef ENABLE_HR_MODE
    attack_detector_fx(encoder, h_EncSetup, s_in_scaled_lp, sub(h_EncSetup->x_exp, 15), currentScratch);
#else
    attack_detector_fx(encoder, h_EncSetup, s_in_scaled, sub(h_EncSetup->x_exp, 15), currentScratch);
#endif
    BASOP_sub_end();
    
    /* begin ener_fx */
    BASOP_sub_start("PerBandEnergy");
    /* currentScratch Size = 160 bytes */
    processPerBandEnergy_fx(ener_fx, &ener_fx_exp, d_fx, d_fx_exp, encoder->bands_offset, encoder->fs_idx,
                            encoder->bands_number, 0, encoder->frame_dms, currentScratch
#ifdef ENABLE_HR_MODE
                            , encoder->hrmode
#endif
    );
    BASOP_sub_end();

    BASOP_sub_start("Near Nyquist Detector");
        /* Near Nyquist Detector */
        processNearNyquistdetector_fx(&encoder->near_nyquist_flag, encoder->fs_idx, encoder->near_nyquist_index,
                                      encoder->bands_number, ener_fx, ener_fx_exp
#ifdef ENABLE_HR_MODE
                                  , encoder->frame_dms, encoder->hrmode );
#else
                                  );
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

    BASOP_sub_start("SnsCompScf");


    /* currentScratch Size = 512 bytes */
    processSnsComputeScf_fx(ener_fx, ener_fx_exp, encoder->fs_idx, encoder->bands_number, scf,
                            h_EncSetup->attdec_detected, encoder->attdec_damping, currentScratch, encoder->sns_damping
    );
    BASOP_sub_end();

    BASOP_sub_start("SnsQuantScfEnc");
    /* currentScratch Size = 500 bytes */
        
    processSnsQuantizeScfEncoder_fx(scf, L_scf_idx, scf_q, currentScratch);
    BASOP_sub_end();

    BASOP_sub_start("SnsInterpScfEnc");
    /* currentScratch Size = 128 bytes */
    processSnsInterpolateScf_fx(scf_q, int_scf_fx, int_scf_fx_exp, 1, encoder->bands_number, currentScratch);
    BASOP_sub_end();

    BASOP_sub_start("Mdct shaping_enc");
    processMdctShaping_fx(d_fx, int_scf_fx, int_scf_fx_exp, encoder->bands_offset, encoder->bands_number);

    BASOP_sub_end();
    /* end int_scf_fx_exp */
    BASOP_sub_start("BandwidthControl_enc");
    if (encoder->bandwidth < L_shr_pos(encoder->fs, 1))
    {
        process_cutoff_bandwidth(d_fx, encoder->yLen, encoder->bw_ctrl_cutoff_bin);
        BW_cutoff_idx = s_min(BW_cutoff_idx, encoder->bw_index);
    }
    BASOP_sub_end();
    BASOP_sub_start("Tns_enc");
    /* currentScratch Size = 2 * MAX_LEN + 220 */

    IF (h_EncSetup->lfe == 0)
    {
    processTnsCoder_fx(&(h_EncSetup->tns_bits), indexes, d_fx, BW_cutoff_idx, tns_order, &tns_numfilters,
                       h_EncSetup->enable_lpc_weighting, encoder->nSubdivisions, encoder->frame_dms,
                       encoder->frame_length, currentScratch
#ifdef ENABLE_HR_MODE
                       , encoder->hrmode
#endif
                       , encoder->near_nyquist_flag   
    );
        }
    ELSE
    {
        tns_numfilters = 1;
        move16();
        tns_order[0] = 0;
        move16();
        h_EncSetup->tns_bits = tns_numfilters;
        move16();
    }
    BASOP_sub_end();

    BASOP_sub_start("Est. Global Gain");
    /* currentScratch Size = 4 * MAX_LEN bytes */
    h_EncSetup->targetBitsQuant = sub(h_EncSetup->targetBitsInit, add(h_EncSetup->tns_bits, ltpf_bits));

    test();
    IF (h_EncSetup->targetBitsQuant < 0 && sub(ltpf_bits, 1) > 0)
    {
        /* Disable LTPF */
        h_EncSetup->ltpf_mem_ltpf_on = 0;  move16();
        ltpf_idx[1]                  = 0;  move16();
        ltpf_bits                    = 1;  move16();
        h_EncSetup->targetBitsQuant  = sub(h_EncSetup->targetBitsInit, add(h_EncSetup->tns_bits, ltpf_bits));
    }
        
#ifdef ENABLE_HR_MODE
    Word32 gain32;
#endif
        
    processEstimateGlobalGain_fx(d_fx, d_fx_exp, encoder->yLen, h_EncSetup->targetBitsQuant,
#ifdef ENABLE_HR_MODE
                                 &gain32,
#else
                                 &gain,
#endif
                                 &gain_e,
                                 &quantizedGain, &quantizedGainMin, h_EncSetup->quantizedGainOff,
                                 &h_EncSetup->targetBitsOff, &h_EncSetup->mem_targetBits, h_EncSetup->mem_specBits,
                                 currentScratch
#ifdef ENABLE_HR_MODE
                                 , encoder->hrmode, h_EncSetup->regBits, encoder->frame_dms
#endif
    );
    BASOP_sub_end();
    /* begin q_d_fx16 */
        
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
    BASOP_sub_end();

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
    BASOP_sub_end();

    BASOP_sub_start("Quant. 2");
    IF (sub(gainChange, 1) == 0)
    {
#ifdef ENABLE_HR_MODE
        processQuantizeSpec_fx(d_fx, d_fx_exp, gain32, gain_e, q_d_fx24, encoder->yLen, h_EncSetup->targetBitsQuant,
                               h_EncSetup->targetBitsAri, NULL, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
                               0, encoder->hrmode);
#else
        processQuantizeSpec_fx(d_fx, d_fx_exp, gain, gain_e, q_d_fx16, encoder->yLen, h_EncSetup->targetBitsQuant,
                               h_EncSetup->targetBitsAri, NULL, &nBits, encoder->fs_idx, &lastnz, codingdata, &lsbMode,
                               0);
#endif /* ENABLE_HR_MODE */
    }
    BASOP_sub_end();

    BASOP_sub_start("Res. Cod.");
    IF (lsbMode == 0)
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
        );
    }
    ELSE
    {
        numResBits = 0;
        move16();
    }
    BASOP_sub_end();

    BASOP_sub_start("Noise fac");
    /* currentScratch Size = 2 * MAX_LEN bytes */
    IF (h_EncSetup->lfe == 0)
    {
        processNoiseFactor_fx(&fac_ns_idx, d_fx_exp, d_fx,
#ifdef ENABLE_HR_MODE
                              q_d_fx24,
#else
                              q_d_fx16,
#endif
                              gain, gain_e, BW_cutoff_idx, encoder->frame_dms, h_EncSetup->targetBytes, currentScratch
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
    
    BASOP_sub_start("Entropy cod");
    processEncoderEntropy(bytes, &bp_side, &mask_side, h_EncSetup->targetBitsAri, h_EncSetup->targetBytes,
                          encoder->yLen, encoder->BW_cutoff_bits, tns_numfilters, lsbMode, lastnz, tns_order,
                          fac_ns_idx, quantizedGain, BW_cutoff_idx, ltpf_idx, L_scf_idx, bfi_ext, encoder->fs_idx);
    BASOP_sub_end();

    BASOP_sub_start("Ari cod");
    processAriEncoder_fx(bytes, bp_side, mask_side, h_EncSetup->targetBitsAri,
#ifdef ENABLE_HR_MODE
                         q_d_fx24,
#else
                         q_d_fx16,
#endif
                         tns_order, tns_numfilters, indexes, lastnz, codingdata, resBits, numResBits, lsbMode,
                         h_EncSetup->enable_lpc_weighting, currentScratch);
    BASOP_sub_end();

    BASOP_sub_start("Reorder Bitstream Enc");
    test();
    IF (encoder->combined_channel_coding == 0 && h_EncSetup->n_pc > 0)
    {
        BASOP_sub_start("Reorder Ari dec");
        
#ifdef ENABLE_HR_MODE
        Word32 *xbuf = (Word32 *) scratchAlign(scratchBuffer, 0);
        
        processAriDecoder_fx(bytes, &bp_side, &mask_side, h_EncSetup->total_bits, encoder->yLen, encoder->fs_idx,
                             h_EncSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                             fac_ns_idx, quantizedGain, encoder->frame_dms, h_EncSetup->n_pc, 0,
                             shr_pos(h_EncSetup->total_bits, 3), 1, &gain, &b_left, &gain,
                             xbuf,
                             &gain, resBits, indexes, &gain, 
                             currentScratch
                             , encoder->hrmode
        );
#else
        
        processAriDecoder_fx(bytes, &bp_side, &mask_side, h_EncSetup->total_bits, encoder->yLen, encoder->fs_idx,
                             h_EncSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                             fac_ns_idx, quantizedGain, encoder->frame_dms, h_EncSetup->n_pc, 0,
                             shr_pos(h_EncSetup->total_bits, 3), 1, &gain, &b_left, &gain,
                             codingdata,
                             &gain, resBits, indexes, &gain, 
                             currentScratch
        );
#endif
        
        BASOP_sub_end(); /* Ari dec */
        processReorderBitstream_fx(bytes, h_EncSetup->n_pccw, h_EncSetup->n_pc, b_left, currentScratch);
    }
    BASOP_sub_end();

    /* end q_d_fx16 */
    BASOP_sub_end();

    Dyn_Mem_Deluxe_Out();
}

int Enc_LC3PLUS(LC3PLUS_Enc *encoder, void **input, int bits_per_sample, UWord8 *output, void *scratch, Word16 bfi_ext)
{
    int ch = 0, output_size = 0;
    int input_size = 0;
    int totalBytes = (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
    int output_size2;

    UWord8 *lc3buf = output;

    for (ch = 0; ch < encoder->channels; ch++)
    {
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
        output_size = (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
        BASOP_sub_start("fec_enc");

        fec_encoder(encoder->epmode, encoder->epmr, output, input_size, output_size, encoder->channel_setup[0]->n_pccw,
                    scratch);

        BASOP_sub_end();
    }

    return output_size;
}

