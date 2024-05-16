/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "setup_enc_lc3plus.h"
#include "functions.h"
#include <stdio.h>

/* if encoder is null only size is reported */
int alloc_encoder(LC3PLUS_Enc* encoder, int channels)
{
    int    ch   = 0;
    size_t size = sizeof(LC3PLUS_Enc);

    for (ch = 0; ch < channels; ch++) {
        EncSetup* setup = balloc(encoder, &size, sizeof(EncSetup));
        if (encoder) {
            encoder->channel_setup[ch] = setup;
        }
    }

    return (int)size;
}

LC3PLUS_Error FillEncSetup(LC3PLUS_Enc* encoder, int samplerate, int channels
                            , int hrmode
                            , int32_t lfe_channel_array[]
)
{
    int ch = 0;
    memset(encoder, 0, lc3plus_enc_get_size(samplerate, channels));
    alloc_encoder(encoder, channels);

    encoder->fs     = CODEC_FS(samplerate);
    encoder->fs_in  = samplerate;
    encoder->fs_idx = FS2FS_IDX(encoder->fs);
    encoder->frame_dms = 100;

    if (encoder->fs_idx > 4) {
        encoder->fs_idx = 5;
    }

    encoder->hrmode = hrmode != 0;

    encoder->channels          = channels;
    encoder->frame_ms          = 10;
    encoder->envelope_bits     = 38;
    encoder->global_gain_bits  = 8;
    encoder->noise_fac_bits    = 3;
    encoder->BW_cutoff_bits    = BW_cutoff_bits_all[encoder->fs_idx];

    encoder->r12k8_mem_in_len  = 2 * 8 * encoder->fs / 12800;
    encoder->r12k8_mem_out_len = 24;

    if (lfe_channel_array != NULL)
    {
        for (ch = 0; ch < encoder->channels; ch++)
        {
            encoder->channel_setup[ch]->lfe = lfe_channel_array[ch] != 0;
        }
    }

    encoder->bw_ctrl_active   = 0;
    encoder->bandwidth        = encoder->fs / 2;
    encoder->bandwidth_preset = encoder->fs / 2;


    if (encoder->fs == 8000) {
        encoder->tilt = 14;
    } else if (encoder->fs == 16000) {
        encoder->tilt = 18;
    } else if (encoder->fs == 24000) {
        encoder->tilt = 22;
    } else if (encoder->fs == 32000) {
        encoder->tilt = 26;
    } else if (encoder->fs == 48000) {
        encoder->tilt = 30;
    }
    else if (encoder->fs == 96000) {
        encoder->tilt = 34;
    }

    set_enc_frame_params(encoder);
    return LC3PLUS_OK;
}

/* set frame config params */
void set_enc_frame_params(LC3PLUS_Enc* encoder)
{
    int       ch = 0;
    EncSetup* setup;

    encoder->frame_length       = ceil(encoder->fs * 10 / 1000); /* fs * 0.01*2^6 */
    if (encoder->hrmode == 1)
    {
        encoder->yLen = encoder->frame_length;
    }
    else
    {
        encoder->yLen = MIN(MAX_BW, encoder->frame_length);
        encoder->sns_damping = 0.85;
    }
    
    encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
    encoder->bands_number       = 64;
    encoder->nSubdivisions      = 3;
    encoder->near_nyquist_index = encoder->bands_number - 2;
    encoder->near_nyquist_flag  = 0;
    encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN;
    
    if (encoder->fs_idx == 5)
    {
        encoder->hrmode = 1;
    }

    if (encoder->hrmode)
    {
        encoder->BW_cutoff_bits = 0;
    }
    else
    {
        encoder->BW_cutoff_bits = BW_cutoff_bits_all[encoder->fs_idx];
    }

    if (encoder->frame_ms == 10) {
        encoder->la_zeroes = MDCT_la_zeroes[encoder->fs_idx];
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND[encoder->fs_idx];
        }
        encoder->cutoffBins   = BW_cutoff_bin_all;
        
        encoder->attdec_nblocks         = 4;
        encoder->attdec_damping         = 0.5;
        encoder->attdec_hangover_thresh = 2;
    }
    else if (encoder->frame_ms == 7.5) {
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_7_5ms_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_7_5ms[encoder->fs_idx];
        }
        encoder->la_zeroes    = MDCT_la_zeroes_7_5ms[encoder->fs_idx];
        encoder->cutoffBins   = BW_cutoff_bin_all_7_5ms;
        encoder->attdec_nblocks         = 3;
        encoder->attdec_damping         = 0.3;
        encoder->attdec_hangover_thresh = 1;
        
        encoder->frame_length = (encoder->frame_length >> 2) * 3;
        encoder->yLen = (encoder->yLen >> 2) * 3;
        
        encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
        if (encoder->hrmode)
        {
            encoder->bands_number       = bands_number_7_5ms_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_number       = bands_number_7_5ms[encoder->fs_idx];
        }
        encoder->nSubdivisions      = 3;
        encoder->near_nyquist_index = encoder->bands_number - 4;
        encoder->r12k8_mem_out_len = ceil(2.0 * ((LC3_FLOAT) encoder->frame_length / 2.0 - (LC3_FLOAT) encoder->la_zeroes) * 12800.0 / (LC3_FLOAT) encoder->fs - 8.0);
    }
    else if (encoder->frame_ms == 5) {
        encoder->frame_length = encoder->frame_length >> 1;
        encoder->yLen /= 2;
        encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
        encoder->bands_number       = bands_number_5ms[encoder->fs_idx];
        encoder->nSubdivisions      = 2;
        encoder->near_nyquist_index = encoder->bands_number - 3;
        encoder->la_zeroes = MDCT_la_zeroes_5ms[encoder->fs_idx];
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_5ms_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_5ms[encoder->fs_idx];
        }
        encoder->cutoffBins   = BW_cutoff_bin_all_5ms;
    }
    else if (encoder->frame_ms == 2.5) {
        encoder->la_zeroes = MDCT_la_zeroes_2_5ms[encoder->fs_idx];
        if (encoder->hrmode)
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms[encoder->fs_idx];
        }
        encoder->cutoffBins   = BW_cutoff_bin_all_2_5ms;
        encoder->frame_length = encoder->frame_length >> 2;
        encoder->yLen /= 4;
        encoder->stEnc_mdct_mem_len = encoder->frame_length - encoder->la_zeroes;
        if (encoder->hrmode)
        {
            encoder->bands_number       = bands_number_2_5ms_HR[encoder->fs_idx];
        }
        else
        {
            encoder->bands_number       = bands_number_2_5ms[encoder->fs_idx];
        }

        encoder->nSubdivisions      = 2;
        encoder->near_nyquist_index = encoder->bands_number - 2;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + (LEN_12K8 >> 2);
    }

    for (ch = 0; ch < encoder->channels; ch++) {
        setup = encoder->channel_setup[ch];

        setup->olpa_mem_pitch = 17;
        setup->pitch_flag = 0;
        if (setup->mdctStruct.mem != NULL) {
            mdct_free(&setup->mdctStruct);
            mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms, encoder->fs_idx, encoder->hrmode);

            dct2_free(&setup->dct2StructSNS);
            dct2_init(&setup->dct2StructSNS, M);
        }
        else
        {
            mdct_init(&setup->mdctStruct, encoder->frame_length, encoder->frame_dms, encoder->fs_idx, encoder->hrmode);
            dct2_init(&setup->dct2StructSNS, M);
        }
    }
}

/* change encoder bitrate */
LC3PLUS_Error update_enc_bitrate(LC3PLUS_Enc* encoder, int bitrate)
{
    int ch = 0, bitsTmp = 0, minBR = 0, maxBR = 0, totalBytes = 0;
    LC3_INT channel_bytes = 0, max_bytes = 0;

#ifdef ENABLE_HR_MODE_FL
    if (encoder->hrmode)
    {
        switch (encoder->frame_dms)
        {
        case 25:
            maxBR = 672000;
            if (encoder->fs == 48000) {minBR = MIN_BR_25MS_48KHZ_HR;}
            else if (encoder->fs == 96000) {minBR = MIN_BR_25MS_96KHZ_HR;}
            else { return LC3PLUS_HRMODE_ERROR;}
            break;
        case 50:
            maxBR = 600000;
            if (encoder->fs == 48000) {minBR = MIN_BR_50MS_48KHZ_HR;}
            else if (encoder->fs == 96000) {minBR = MIN_BR_50MS_96KHZ_HR;}
            else { return LC3PLUS_HRMODE_ERROR;}
            break;
        case  75:
            maxBR = 500000;
            if      (encoder->fs == 48000) {minBR = MIN_BR_075DMS_48KHZ_HR;}
            else if (encoder->fs == 96000) {minBR = MIN_BR_075DMS_96KHZ_HR;}
            else                           {return LC3PLUS_HRMODE_ERROR;}
            break;      
        case 100:
            maxBR = 500000;
            if (encoder->fs == 48000) {minBR = MIN_BR_100MS_48KHZ_HR;}
            else if (encoder->fs == 96000) {minBR = MIN_BR_100MS_96KHZ_HR;}
            else { return LC3PLUS_HRMODE_ERROR;}
            break;
        default:
            return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
#endif /* ENABLE_HR_MODE_FL */
    {
        minBR = (MIN_NBYTES << 3);
        maxBR = MAX_BR;

        switch (encoder->frame_dms)
        {
        case  25:
            minBR = MIN_BR_025DMS;
            maxBR = MAX_BR;
            break;
        case  50:
            minBR = MIN_BR_050DMS;
            maxBR = MAX_BR;
            /* have additional limitations for 5.0ms */
            switch (encoder->fs_in)
            {
            case  8000:  maxBR = MAX_BR_050DMS_NB;   break;
            default:                                 break;
            }
            break;
        case  75:
            minBR = MIN_BR_075DMS;
            maxBR = MAX_BR_075DMS;
            /* have additional limitations for 7.5ms */
            switch (encoder->fs_in)
            {
            case  8000:  maxBR = MAX_BR_075DMS_NB  ; break;
            case 16000:  maxBR = MAX_BR_075DMS_WB  ; break;
            case 24000:  maxBR = MAX_BR_075DMS_SSWB; break;
            default:                                 break;
            }
            break;
        case 100: 
            /* have additional limitations for 10ms */
            minBR = MIN_BR_100DMS;
            maxBR = MAX_BR;
            switch (encoder->fs_in)
            {
            case  8000:  maxBR = MAX_BR_100DMS_NB  ; break;
            case 16000:  maxBR = MAX_BR_100DMS_WB  ; break;
            case 24000:  maxBR = MAX_BR_100DMS_SSWB; break;
            default:     maxBR = MAX_BR;             break;
            }
            break;
        default: return LC3PLUS_FRAMEMS_ERROR;
        }
        maxBR *= (encoder->fs_in == 44100 ? 441. / 480 : 1);
    }
    minBR *= encoder->channels;
    maxBR *= encoder->channels;
    
    encoder->combined_channel_coding = 0;
    
    if (encoder->channels > 1 && encoder->epmode)
    {
        if (encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in) <= 160)
        {
            encoder->combined_channel_coding = 1;
        }
    }

    if (encoder->epmode > 0)
    {
        max_bytes = bitrate * encoder->frame_length / (8 * encoder->fs_in * encoder->channels);
        if (max_bytes < FEC_SLOT_BYTES_MIN || max_bytes > FEC_SLOT_BYTES_MAX)
        {
            encoder->lc3_br_set = 0;
            return LC3PLUS_BITRATE_ERROR;
        }
    }

    if (encoder->combined_channel_coding)
    {
        totalBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding,
                                          bitrate * encoder->frame_length / (8 * encoder->fs_in));

        encoder->channel_setup[0]->n_pccw =
            fec_get_n_pccw(bitrate * encoder->frame_length / (8 * encoder->fs_in), encoder->epmode,
                           encoder->combined_channel_coding);

        encoder->channel_setup[0]->n_pc = fec_get_n_pc(encoder->epmode, encoder->channel_setup[0]->n_pccw,
                                                       bitrate * encoder->frame_length / (8 * encoder->fs_in));
    }
    else
    {
        totalBytes = bitrate * encoder->frame_length / (8 * encoder->fs_in);
    }
    
    if (encoder->frame_dms <= 50)
    {
        encoder->tnsMaxOrder = 4;
    } else {
        encoder->tnsMaxOrder = 8;
    }

    if (bitrate < minBR || bitrate > maxBR) {
        return LC3PLUS_BITRATE_ERROR;
    }
    
    encoder->lc3_br_set = 1;
    for (ch = 0; ch < encoder->channels; ch++) {

        EncSetup* setup = encoder->channel_setup[ch];
        
        setup->targetBytes = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));
        channel_bytes = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));

        if (encoder->combined_channel_coding)
        {
            setup->targetBytes = channel_bytes;
        }
        else
        {
            setup->targetBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding, channel_bytes);
            setup->n_pccw = fec_get_n_pccw(channel_bytes, encoder->epmode, encoder->combined_channel_coding);
            setup->n_pc = fec_get_n_pc(encoder->epmode, setup->n_pccw, channel_bytes);
        }
        // reduce bandwith to 12kHz if bitrate is low
        if (encoder->frame_dms == 100 &&
            ((setup->targetBytes < 40 && encoder->fs == 48000) ||
             (setup->targetBytes < 36 && encoder->fs == 32000)))
        {
            encoder->bandwidth = MIN(12000, encoder->bandwidth_preset);
        }
        else
        {
            /* channel with highest index has lowest bitrate.
            For a second channel with lower targetBytes, bandwidth is overwritten */
            encoder->bandwidth = encoder->bandwidth_preset;
        }
        encoder->bw_ctrl_cutoff_bin = encoder->bandwidth * encoder->frame_dms / 5000;
        encoder->bw_index           = (encoder->bandwidth / 4000) - 1;
        setup->total_bits     = setup->targetBytes << 3;
        setup->targetBitsInit = setup->total_bits - encoder->envelope_bits - encoder->global_gain_bits -
                                encoder->noise_fac_bits - encoder->BW_cutoff_bits -
                                ceil(LC3_LOGTWO(encoder->frame_length / 2)) - 2 - 1;

        if (setup->total_bits > 1280) {
            setup->targetBitsInit = setup->targetBitsInit - 1;
        }
        if (setup->total_bits > 2560) {
            setup->targetBitsInit = setup->targetBitsInit - 1;
        }

        if (encoder->hrmode)
        {
            setup->targetBitsInit -= 1;
        }

        setup->targetBitsAri        = setup->total_bits;
        setup->enable_lpc_weighting = setup->total_bits < 480;

        if (encoder->frame_ms == 7.5) {
            setup->enable_lpc_weighting = setup->total_bits < 360;
        }
        if (encoder->frame_ms == 5) {
            setup->enable_lpc_weighting = setup->total_bits < 240;
        }
        if (encoder->frame_ms == 2.5) {
            setup->enable_lpc_weighting = setup->total_bits < 120;
        }

        setup->quantizedGainOff =
            -(MIN(115, setup->total_bits / (10 * (encoder->fs_idx + 1))) + 105 + 5 * (encoder->fs_idx + 1));

        if (encoder->hrmode && encoder->fs_idx == 5)
        {
            setup->quantizedGainOff = MAX(setup->quantizedGainOff, -181);
        }

        if (encoder->frame_ms == 10 && ((encoder->fs_in >= 44100 && setup->targetBytes >= 100) ||
                                        (encoder->fs_in == 32000 && setup->targetBytes >= 81)) && setup->targetBytes < 340 && encoder->hrmode == 0) {
            setup->attack_handling = 1;

        }     
        else if (encoder->frame_dms == 75 && ((encoder->fs_in >= 44100 && setup->targetBytes >= 75) ||
                (encoder->fs_in == 32000 && setup->targetBytes >= 61)) && setup->targetBytes < 150 && encoder->hrmode == 0)
        {
            setup->attack_handling = 1;
        }
        else
        {
            /* reset for bitrate switching */
            setup->attack_handling = 0;

            setup->attdec_filter_mem[0] = 0;
            setup->attdec_filter_mem[1] = 0;

            setup->attdec_detected   = 0;
            setup->attdec_position   = 0;
            setup->attdec_acc_energy = 0;
        }

        bitsTmp = setup->total_bits;
        if (encoder->frame_ms == 2.5) {
            bitsTmp = bitsTmp * 4.0 * (1.0 - 0.4);
        }
        if (encoder->frame_ms == 5) {
            bitsTmp = bitsTmp * 2 - 160;
        }
        if (encoder->frame_ms == 7.5) {
            bitsTmp = round(bitsTmp * 10 / 7.5);
        }

        if (bitsTmp < 400 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 480 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 560 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else if (bitsTmp < 640 + (encoder->fs_idx - 1) * 80) {
            setup->ltpf_enable = 1;
        } else {
            setup->ltpf_enable = 0;
        }
        if (encoder->hrmode) {
            setup->ltpf_enable = 0;
        }

        /* turn down SNS shaping for higher rates */
        if (encoder->hrmode == 0) {
            encoder->sns_damping = 0.85;
        } else {
            encoder->sns_damping = 0.6;
            if (encoder->fs_idx >= 4) {
                if (encoder->frame_ms == 10)
                {
                    if (setup->total_bits > 4400) {
                        encoder->sns_damping = 6881.0/32768.0;
                    }
                }
                if (encoder->frame_ms == 7.5)
                {
                    if (setup->total_bits > 3*4400/4) {
                        encoder->sns_damping = 5898.0/32768.0;
                    }
                }
                if (encoder->frame_ms == 5)
                {
                    if (setup->total_bits > 4600/2) {
                        encoder->sns_damping = 4915.0/32768.0;
                    }
                }
                if (encoder->frame_ms == 2.5)
                {
                    if (setup->total_bits > 4600/4) {
                        encoder->sns_damping = 4915.0/32768.0;
                    }
                }
            }
        }

        if (encoder->hrmode && encoder->fs_idx >= 4)
        {
            int real_rate = setup->targetBytes * 8000 / encoder->frame_ms;
            setup->regBits = real_rate / 12500;

            if (encoder->fs_idx == 5)
            {
                if (encoder->frame_ms == 10)
                {
                    setup->regBits +=2;
                }
                if (encoder->frame_ms == 7.5)
                {
                    setup->regBits +=1;
                }
                if (encoder->frame_ms == 2.5)
                {
                    setup->regBits -= 6;
                }
            }
            else
            {
                if (encoder->frame_ms == 2.5)
                {
                    setup->regBits -= 6;
                }
                else if (encoder->frame_ms == 5)
                {
                    setup->regBits += 0;
                }
                if (encoder->frame_ms == 7.5)
                {
                    setup->regBits +=2;
                }
                if (encoder->frame_ms == 10)
                {
                    setup->regBits += 5;
                }
            }
            if (setup->regBits < 6)
            {
                setup->regBits = 6;
            }
            if (setup->regBits > 23)
            {
                setup->regBits = 23;
            }
        }
        else
        {
            setup->regBits = -1;
        }
    }

    encoder->bitrate = bitrate;

    return LC3PLUS_OK;
}

void update_enc_bandwidth(LC3PLUS_Enc* encoder, int bandwidth)
{
    int index = 0;

    if (bandwidth >= encoder->fs_in) {
        encoder->bandwidth = 0;
    }
    else
    {
        encoder->bandwidth = bandwidth;
        index              = FS2FS_IDX(bandwidth);
        encoder->bw_ctrl_cutoff_bin = encoder->cutoffBins[index];
    }
}
