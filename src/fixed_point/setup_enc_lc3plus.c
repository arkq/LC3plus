/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "setup_enc_lc3plus.h"

/* if encoder is null only size is reported */
int alloc_encoder(LC3PLUS_Enc *encoder, int samplerate, int channels)
{
    int    ch         = 0;
    size_t size       = sizeof(LC3PLUS_Enc);
    void * mdct_mem32 = NULL, *stEnc_mdct_mem = NULL;

    for (ch = 0; ch < channels; ch++)
    {
        EncSetup *setup = balloc(encoder, &size, sizeof(EncSetup));
        mdct_mem32      = balloc(encoder, &size, sizeof(*setup->mdct_mem32) * DYN_MAX_MDCT_LEN(samplerate));
        stEnc_mdct_mem  = balloc(encoder, &size, sizeof(*setup->stEnc_mdct_mem) * DYN_MAX_MDCT_LEN(samplerate));
        if (encoder)
        {
            encoder->channel_setup[ch] = setup;
            setup->mdct_mem32          = mdct_mem32;
            setup->stEnc_mdct_mem      = stEnc_mdct_mem;
        }
    }

    return (int)size;
}

LC3PLUS_Error FillEncSetup(LC3PLUS_Enc *encoder, int samplerate, int channels
#ifdef ENABLE_HR_MODE
                           , int hrmode
#endif
                        , int32_t lfe_channel_array[]   
                          )
{
    int ch = 0;

    memset(encoder, 0, lc3plus_enc_get_size(samplerate, channels));
    alloc_encoder(encoder, samplerate, channels);

    encoder->fs     = CODEC_FS(samplerate);
    encoder->fs_in  = samplerate;
    encoder->fs_idx = FS2FS_IDX(encoder->fs);
    
#  ifdef ENABLE_HR_MODE
    if (encoder->fs_idx > 4)
    {
        encoder->fs_idx = 5; 
    }
    encoder->hrmode = hrmode != 0;
#  endif
    
    encoder->channels         = channels;
    encoder->frame_dms        = LC3PLUS_FRAME_DURATION_10MS;
    encoder->envelope_bits    = 38;
    #ifdef CR9_C_ADD_1p25MS_LRSNS
    /*  "38" SNS bit constant  kept here,  final LR-SNS bitrate adjusted in enc_lc3() after  LR-SNSVQ */
#endif /* CR9_C_ADD_1p25MS_LRSNS */
    encoder->global_gain_bits = 8;
    encoder->noise_fac_bits   = 3;
    encoder->r12k8_mem_in_len  = extract_l(L_shr_pos(Mpy_32_16(encoder->fs, 20972), 9));
    encoder->r12k8_mem_out_len = 24;
    encoder->epmr = LC3PLUS_EPMR_ZERO;
    encoder->bw_ctrl_active   = 0;
    encoder->bandwidth        = L_shr_pos(encoder->fs, 1);
    encoder->bandwidth_preset = L_shr_pos(encoder->fs, 1);

    if (lfe_channel_array != NULL)
    {
        for (ch = 0; ch < encoder->channels; ch++)
        {
            encoder->channel_setup[ch]->lfe = lfe_channel_array[ch] != 0;
        }
    }

    for (ch = 0; ch < encoder->channels; ch++)
    {
        encoder->channel_setup[ch]->x_exp      = 15;
        encoder->channel_setup[ch]->resamp_exp = 17;
    }

    set_enc_frame_params(encoder);

    return lc3plus_enc_set_ep_mode(encoder, LC3PLUS_EP_OFF); /* also calls update_enc_bitrate */
}

/* set frame config params */
void set_enc_frame_params(LC3PLUS_Enc *encoder)
{
    Word16 tmp;
#ifdef CR9_C_ADD_1p25MS
#ifndef FIX_TX_RX_STRUCT_STEREO
    encoder->Tx_ltpf           = 0;
#endif
    encoder->LT_normcorr       = 0xFFFF >> 2;
#endif
    
    encoder->frame_length = extract_l(L_shr_pos(Mpy_32_16(encoder->fs, 20972), 6)); /* fs * 0.01*2^6 */
    
#  ifdef ENABLE_HR_MODE
    if (encoder->hrmode)
    {
        encoder->BW_cutoff_bits = 0;
    }
    else
#  endif
    {
        encoder->BW_cutoff_bits = BW_cutoff_bits_all[encoder->fs_idx];
    }

    SWITCH (encoder->frame_dms)
    {
#ifdef CR9_C_ADD_1p25MS
    case LC3PLUS_FRAME_DURATION_1p25MS:
        encoder->BW_cutoff_bits     = 0;
        encoder->frame_length       = shr_pos(encoder->frame_length, 3);
        encoder->la_zeroes          = 0;
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->nSubdivisions      = 0;
#ifdef FIX_LTPF_PITCH_MEM_LEN
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + 112;
#else
        encoder->ltpf_mem_in_len = LTPF_MEMIN_LEN + ( LEN_12K8 >> 1 ) + 16;
#endif
        encoder->r12k8_mem_out_len  = 8;
#ifndef CR9_C_ADD_1p25MS_NOISEFILLING
        encoder->noise_fac_bits     = 0;
#endif
#    ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            assert(0);
        }
        else
#    endif
        {
            encoder->yLen         = s_min(MAX_BW >> 3, encoder->frame_length);
            encoder->W_fx         = LowDelayShapes_n960_1_25ms[encoder->fs_idx];
            encoder->W_size       = LowDelayShapes_n960_len_1_25ms[encoder->fs_idx];
            encoder->bands_number = bands_number_1_25ms[encoder->fs_idx];
            encoder->bands_offset = bands_offset_1_25ms[encoder->fs_idx];
            encoder->near_nyquist_index  = encoder->bands_number - 2;
            encoder->near_nyquist_flag = 0;
        }
    BREAK;  
#endif
    case LC3PLUS_FRAME_DURATION_2p5MS:
        encoder->frame_length       = shr_pos(encoder->frame_length, 2);
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes_2_5ms[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->nSubdivisions      = 2;
#ifdef FIX_LTPF_PITCH_MEM_LEN
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + 96;
#else
        encoder->ltpf_mem_in_len = LTPF_MEMIN_LEN + ( LEN_12K8 >> 2 );
#endif
        
#    ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_number = bands_number_2_5ms_HR[encoder->fs_idx];
            encoder->bands_offset = bands_offset_2_5ms_HR[encoder->fs_idx - 4];
            encoder->W_fx         = LowDelayShapes_n960_HRA_2_5ms[encoder->fs_idx - 4];		
            encoder->W_size       = LowDelayShapes_n960_len_2_5ms[encoder->fs_idx];
            encoder->yLen         = encoder->frame_length;
        }
        else
#    endif
        {
            encoder->yLen         = s_min(MAX_BW >> 2, encoder->frame_length);
            encoder->W_fx         = LowDelayShapes_n960_2_5ms[encoder->fs_idx];
            encoder->W_size       = LowDelayShapes_n960_len_2_5ms[encoder->fs_idx];
            encoder->bands_number = bands_number_2_5ms[encoder->fs_idx];
            encoder->bands_offset = bands_offset_2_5ms[encoder->fs_idx];
            encoder->near_nyquist_index  = encoder->bands_number - 2;
            encoder->near_nyquist_flag = 0;
        }

        BREAK;
    case LC3PLUS_FRAME_DURATION_5MS:
        encoder->frame_length       = shr_pos(encoder->frame_length, 1);
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes_5ms[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->nSubdivisions      = 2;
#ifdef FIX_LTPF_PITCH_MEM_LEN
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN + 64;
#else
        encoder->ltpf_mem_in_len = LTPF_MEMIN_LEN + ( LEN_12K8 >> 1 );
#endif
        
#    ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_offset = bands_offset_5ms_HR[encoder->fs_idx - 4];
            encoder->W_fx         = LowDelayShapes_n960_HRA_5ms[encoder->fs_idx - 4];		
            encoder->W_size       = LowDelayShapes_n960_len_5ms[encoder->fs_idx];
            encoder->yLen         = encoder->frame_length;
            encoder->bands_number = bands_number_5ms[encoder->fs_idx];
        }
        else
#    endif
        {
            encoder->yLen         = s_min(MAX_BW >> 1, encoder->frame_length);
            encoder->W_fx         = LowDelayShapes_n960_5ms[encoder->fs_idx];
            encoder->W_size       = LowDelayShapes_n960_len_5ms[encoder->fs_idx];
            encoder->bands_number = bands_number_5ms[encoder->fs_idx];
            encoder->bands_offset = bands_offset_5ms[encoder->fs_idx];
            encoder->near_nyquist_index  = encoder->bands_number - 3;
            encoder->near_nyquist_flag = 0;
        }
        BREAK;
        
    case LC3PLUS_FRAME_DURATION_7p5MS:
        tmp                         = shr_pos(encoder->frame_length, 2);
        encoder->frame_length       = add(tmp, add(tmp, tmp));
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes_7_5ms[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->bands_number       = bands_number_7_5ms[encoder->fs_idx];
        encoder->nSubdivisions      = 3;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN;
        encoder->attdec_nblocks		= 3;
        encoder->attdec_damping		= 9830;
        encoder->attdec_hangover_thresh = 1;
        encoder->r12k8_mem_out_len  = 44;
        encoder->near_nyquist_index  = encoder->bands_number - 4;
        encoder->near_nyquist_flag = 0;
#    ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->yLen               = encoder->frame_length;
            encoder->W_fx               = LowDelayShapes_n960_HRA_7_5ms[encoder->fs_idx - 4];
            encoder->W_size             = LowDelayShapes_n960_len_7_5ms[encoder->fs_idx];
            encoder->bands_number       = bands_number_7_5ms[encoder->fs_idx];
            encoder->bands_offset       = bands_offset_7_5ms_HR[encoder->fs_idx - 4];
        }
        else
#    endif
        {
            encoder->yLen               = s_min((MAX_BW >> 2) * 3, encoder->frame_length);
            encoder->W_fx               = LowDelayShapes_n960_7_5ms[encoder->fs_idx];
            encoder->W_size             = LowDelayShapes_n960_len_7_5ms[encoder->fs_idx];
            encoder->bands_number       = bands_number_7_5ms[encoder->fs_idx];
            encoder->bands_offset       = bands_offset_7_5ms[encoder->fs_idx];
        }
		BREAK;
        
    case LC3PLUS_FRAME_DURATION_10MS:
        encoder->la_zeroes          = LowDelayShapes_n960_la_zeroes[encoder->fs_idx];
        encoder->stEnc_mdct_mem_len = sub(encoder->frame_length, encoder->la_zeroes);
        encoder->bands_number       = 64;
        encoder->nSubdivisions      = 3;
        encoder->ltpf_mem_in_len    = LTPF_MEMIN_LEN;
        encoder->attdec_nblocks         = 4;
        encoder->attdec_damping         = 16384;
        encoder->attdec_hangover_thresh = 2;
        encoder->near_nyquist_index  = encoder->bands_number - 2;
        encoder->near_nyquist_flag = 0;
        
#    ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            encoder->bands_offset = bands_offset_HR[encoder->fs_idx - 4];		
            encoder->W_fx         = LowDelayShapes_n960_HRA[encoder->fs_idx - 4];		
            encoder->W_size       = LowDelayShapes_n960_len[encoder->fs_idx];
            encoder->yLen         = encoder->frame_length;
        }
        else
#    endif
        {
            encoder->yLen         = s_min(MAX_BW, encoder->frame_length);
            encoder->W_fx         = LowDelayShapes_n960[encoder->fs_idx];
            encoder->W_size       = LowDelayShapes_n960_len[encoder->fs_idx];
            encoder->bands_offset = bands_offset[encoder->fs_idx];
        }
        BREAK;
    case LC3PLUS_FRAME_DURATION_UNDEFINED:
        assert(0);
    }
}

/* change encoder bitrate */
LC3PLUS_Error update_enc_bitrate(LC3PLUS_Enc *encoder, int bitrate)
{
    int ch = 0;
    int totalBytes = 0, maxBR = 0, minBR = 0, max_bytes = 0;
    int channel_bytes = 0;

#   ifdef CR12_D_FIX_BITRATE_LIMITS
    int fec_slot_bytes_min = 0, check_bytes = 0;
#   endif
#  ifdef ENABLE_HR_MODE
    if (encoder->hrmode)
    {
        SWITCH (encoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
        case LC3PLUS_FRAME_DURATION_1p25MS:
            assert(0);
            maxBR = 672000;
            if (encoder->fs == 48000)
            {
                minBR = MIN_BR_25MS_48KHZ_HR;
            }
            else if (encoder->fs == 96000)
            {
                minBR = MIN_BR_25MS_96KHZ_HR;
            }
            else
            {
                return LC3PLUS_HRMODE_ERROR;
            }
            break;
#endif
        case LC3PLUS_FRAME_DURATION_2p5MS:
            maxBR = 672000;
            if (encoder->fs == 48000)
            {
                minBR = MIN_BR_25MS_48KHZ_HR;
            }
            else if (encoder->fs == 96000)
            {
                minBR = MIN_BR_25MS_96KHZ_HR;
            }
            else
            {
                return LC3PLUS_HRMODE_ERROR;
            }
            break;
        case LC3PLUS_FRAME_DURATION_5MS:
            maxBR = 600000;
            if (encoder->fs == 48000)
            {
                minBR = MIN_BR_50MS_48KHZ_HR;
            }
            else if (encoder->fs == 96000)
            {
                minBR = MIN_BR_50MS_96KHZ_HR;
            }
            else
            {
                return LC3PLUS_HRMODE_ERROR;
            }
            break;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            maxBR = 500000;
            if      (encoder->fs == 48000) {minBR = MIN_BR_075DMS_48KHZ_HR;}
            else if (encoder->fs == 96000) {minBR = MIN_BR_075DMS_96KHZ_HR;}
            else                           {return LC3PLUS_HRMODE_ERROR;}
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS:
            maxBR = 500000;
            if (encoder->fs == 48000)
            {
                minBR = MIN_BR_100MS_48KHZ_HR;
            }
            else if (encoder->fs == 96000)
            {
                minBR = MIN_BR_100MS_96KHZ_HR;
            }
            else
            {
                return LC3PLUS_HRMODE_ERROR;
            }
            break;
        case LC3PLUS_FRAME_DURATION_UNDEFINED: return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
    {
        minBR = (MIN_NBYTES << 3);
        maxBR = MAX_BR;
            
        SWITCH (encoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
        case  LC3PLUS_FRAME_DURATION_1p25MS:
            minBR = MIN_BR_0125DMS;
            maxBR = MAX_BR_0125DMS;
            BREAK;
#endif
        case  LC3PLUS_FRAME_DURATION_2p5MS:
            minBR = MIN_BR_025DMS;
            maxBR = MAX_BR;
            BREAK;
        case  LC3PLUS_FRAME_DURATION_5MS:
            minBR = MIN_BR_050DMS;
            maxBR = MAX_BR;
            /* have additional limitations for 5.0ms */
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_050DMS_NB;   BREAK;
#        endif
            default:                                 BREAK;
            }
            BREAK;
        case  LC3PLUS_FRAME_DURATION_7p5MS:
            minBR = MIN_BR_075DMS;
            maxBR = MAX_BR_075DMS; // special value for maxBR @ 7.5ms
            /* have additional limitations for 7.5ms */
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_075DMS_NB  ; BREAK;
#        endif
#        ifdef SUBSET_WB
            case 16000:  maxBR = MAX_BR_075DMS_WB  ; BREAK;
#        endif
#        ifdef SUBSET_SSWB
            case 24000:  maxBR = MAX_BR_075DMS_SSWB; BREAK;
#        endif
            default:                                 BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS: 
            /* have additional limitations for 10ms */
            minBR = MIN_BR_100DMS;
            maxBR = MAX_BR;
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_100DMS_NB  ; BREAK;
#        endif
#        ifdef SUBSET_WB
            case 16000:  maxBR = MAX_BR_100DMS_WB  ; BREAK;
#        endif
#        ifdef SUBSET_SSWB
            case 24000:  maxBR = MAX_BR_100DMS_SSWB; BREAK;
#        endif
            default:     maxBR = MAX_BR;             BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_UNDEFINED: return LC3PLUS_FRAMEMS_ERROR;
        }
 
        /* 441/480 in Q31 and 1000/75 in Q23 */
        if (encoder->fs_in == 44100)
        {
            minBR = Mpy_32_32(minBR, 1973000602);
            maxBR = Mpy_32_32(maxBR, 1973000602);
        }
    }
#else /* ENABLE_HR_MODE */
        minBR = (MIN_NBYTES << 3);
        maxBR = MAX_BR;
            
        SWITCH (encoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
        case  LC3PLUS_FRAME_DURATION_1p25MS:
            minBR = MIN_BR_0125DMS;
            maxBR = MAX_BR_0125DMS;
            BREAK;
#endif
        case  LC3PLUS_FRAME_DURATION_2p5MS:
            minBR = MIN_BR_025DMS;
            maxBR = MAX_BR;
            BREAK;
        case  LC3PLUS_FRAME_DURATION_5MS:
            minBR = MIN_BR_050DMS;
            maxBR = MAX_BR;
            /* have additional limitations for 5.0ms */
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_050DMS_NB;   BREAK;
#        endif
            default:                                 BREAK;
            }
            BREAK;
        case  LC3PLUS_FRAME_DURATION_7p5MS:
            minBR = MIN_BR_075DMS;
            maxBR = MAX_BR_075DMS; // special value for maxBR @ 7.5ms
            /* have additional limitations for 7.5ms */
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_075DMS_NB  ; BREAK;
#        endif
#        ifdef SUBSET_WB
            case 16000:  maxBR = MAX_BR_075DMS_WB  ; BREAK;
#        endif
#        ifdef SUBSET_SSWB
            case 24000:  maxBR = MAX_BR_075DMS_SSWB; BREAK;
#        endif
            default:                                 BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS: 
            /* have additional limitations for 10ms */
            minBR = MIN_BR_100DMS;
            maxBR = MAX_BR;
            SWITCH (encoder->fs_in)
            {
#        ifdef SUBSET_NB
            case  8000:  maxBR = MAX_BR_100DMS_NB  ; BREAK;
#        endif
#        ifdef SUBSET_WB
            case 16000:  maxBR = MAX_BR_100DMS_WB  ; BREAK;
#        endif
#        ifdef SUBSET_SSWB
            case 24000:  maxBR = MAX_BR_100DMS_SSWB; BREAK;
#        endif
            default:     maxBR = MAX_BR;             BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_UNDEFINED: return LC3PLUS_FRAMEMS_ERROR;
        }
 
        /* 441/480 in Q31 and 1000/75 in Q23 */
        if (encoder->fs_in == 44100)
        {
            minBR = Mpy_32_32(minBR, 1973000602);
            maxBR = Mpy_32_32(maxBR, 1973000602);
        }
#endif /* ENABLE_HR_MODE */

    minBR *= encoder->channels;
    maxBR *= encoder->channels;

    if (bitrate < minBR || bitrate > maxBR)
    {
        return LC3PLUS_BITRATE_ERROR;
    }

    encoder->bitrate    = bitrate;
    encoder->lc3_br_set = 1;

    /* move stuff to encoder->channel_setup */

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
#ifdef CR12_D_FIX_BITRATE_LIMITS
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode){
            SWITCH( encoder->frame_dms )
            {
            case LC3PLUS_FRAME_DURATION_2p5MS:
                IF( encoder->fs_in == 48000){
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_025DMS_48KHZ_HR;
                } ELSE {
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_025DMS_96KHZ_HR;
                }
                BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:
                IF( encoder->fs_in == 48000){
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_050DMS_48KHZ_HR;
                } ELSE {
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_050DMS_96KHZ_HR;
                }
                BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:
                IF( encoder->fs_in == 48000){
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_075DMS_48KHZ_HR;
                } ELSE {
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_075DMS_96KHZ_HR;
                }
                BREAK;
            case LC3PLUS_FRAME_DURATION_10MS:
                IF( encoder->fs_in == 48000){
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_100DMS_48KHZ_HR;
                } ELSE {
                    fec_slot_bytes_min = FEC_SLOT_BYTES_MIN_100DMS_96KHZ_HR;
                }
                BREAK;
            default:
                return LC3PLUS_FRAMEMS_ERROR;
            }
        }
        else 
#endif        
        {
            fec_slot_bytes_min = FEC_SLOT_BYTES_MIN;
        }

        check_bytes = bitrate * encoder->frame_length / ( 8 * encoder->fs_in * encoder->channels );
        maxBR = FEC_SLOT_BYTES_MAX * ( 8 * encoder->fs_in * encoder->channels ) / encoder->frame_length;
        if ( check_bytes < fec_slot_bytes_min || bitrate > maxBR )
#else        
        max_bytes = bitrate * encoder->frame_length / (8 * encoder->fs_in * encoder->channels);
        if (max_bytes < FEC_SLOT_BYTES_MIN || max_bytes > FEC_SLOT_BYTES_MAX)
#endif /* CR12_D_FIX_BITRATE_LIMITS */
        {
            return LC3PLUS_BITRATE_ERROR;
        }
    }

    if (encoder->combined_channel_coding)
    {
        totalBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding,
                                       bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in));

        encoder->channel_setup[0]->n_pccw =
            fec_get_n_pccw(bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in), encoder->epmode,
                           encoder->combined_channel_coding);

        encoder->channel_setup[0]->n_pc = fec_get_n_pc(encoder->epmode, encoder->channel_setup[0]->n_pccw,
                                                       bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in));
    }
    else
    {
        totalBytes = bitrate * (Word32)encoder->frame_length / (8 * encoder->fs_in);
    }

    for (ch = 0; ch < encoder->channels; ch++)
    {
        EncSetup *setup = encoder->channel_setup[ch];
        channel_bytes = totalBytes / encoder->channels + (ch < (totalBytes % encoder->channels));

        if (encoder->combined_channel_coding)
        {
            setup->targetBytes = channel_bytes;
        }
        else
        {
            setup->targetBytes = fec_get_data_size(encoder->epmode, encoder->combined_channel_coding, channel_bytes);
            setup->n_pccw      = fec_get_n_pccw(channel_bytes, encoder->epmode, encoder->combined_channel_coding);
            setup->n_pc        = fec_get_n_pc(encoder->epmode, setup->n_pccw, channel_bytes);
        }
        /* reduce bandwith to 12kHz if bitrate is low */
        if (sub(encoder->frame_dms, LC3PLUS_FRAME_DURATION_10MS) == 0 &&
            ((sub(setup->targetBytes, 40) < 0 && L_sub(encoder->fs, 48000) == 0) ||
             (sub(setup->targetBytes, 36) < 0 && L_sub(encoder->fs, 32000) == 0)))
        {
            encoder->bandwidth = L_min(12000, encoder->bandwidth_preset);
        }
        else
        {
            /* channel with highest index has lowest bitrate.
            For a second channel with lower targetBytes, bandwidth is overwritten */
            encoder->bandwidth = encoder->bandwidth_preset;
        }

        {
            Word16 tmp = 0;

            SWITCH(encoder->frame_dms)
            {
#ifdef CR9_C_ADD_1p25MS
                case  LC3PLUS_FRAME_DURATION_1p25MS: tmp = 1; BREAK;
#endif
                case  LC3PLUS_FRAME_DURATION_2p5MS:  tmp = 2; BREAK;
                case  LC3PLUS_FRAME_DURATION_5MS:    tmp = 4; BREAK;
                case  LC3PLUS_FRAME_DURATION_7p5MS:  tmp = 6; BREAK;
                case  LC3PLUS_FRAME_DURATION_10MS:   tmp = 8; BREAK;
                case  LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
            }
            
            encoder->bw_ctrl_cutoff_bin = L_mult0(Mpy_32_32(encoder->bandwidth, 10737420>>1), tmp); /* bandwidth * frame_dms / 5000 */
        }

        encoder->bw_index           = sub( Mpy_32_32(encoder->bandwidth, 536871), 1);                     /* (bandwidth / 4000 ) - 1      */
        SWITCH (encoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
            case LC3PLUS_FRAME_DURATION_1p25MS: max_bytes = MAX_NBYTES_025; break;
#endif
            case LC3PLUS_FRAME_DURATION_2p5MS: max_bytes = MAX_NBYTES_025; break;
            case LC3PLUS_FRAME_DURATION_5MS: max_bytes = MAX_NBYTES_050; break;
            case LC3PLUS_FRAME_DURATION_7p5MS: max_bytes = MAX_NBYTES_075; BREAK;
            case LC3PLUS_FRAME_DURATION_10MS: max_bytes = MAX_NBYTES_100; break;
            case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
        }
#ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            max_bytes = MAX_NBYTES_RED_HR;
        }
#endif
        if (setup->targetBytes < MIN_NBYTES || setup->targetBytes > max_bytes)
        {
            return LC3PLUS_BITRATE_ERROR;
        }

        setup->total_bits = shl(setup->targetBytes, 3);
        setup->targetBitsInit =
            sub(setup->total_bits,
                add(encoder->envelope_bits,
                    add(encoder->global_gain_bits, add(encoder->noise_fac_bits, encoder->BW_cutoff_bits))));
        setup->targetBitsInit = sub(setup->targetBitsInit, getLastNzBits_fx (encoder->frame_length) + 3);
        if (setup->total_bits > 1280)
        {
            setup->targetBitsInit = sub(setup->targetBitsInit, 1);
        }
        if (setup->total_bits > 2560)
        {
            setup->targetBitsInit = sub(setup->targetBitsInit, 1);
        }
        
#  ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            setup->targetBitsInit -= 1;
        }
#  endif

        setup->targetBitsAri = setup->total_bits;

        SWITCH (encoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
            case LC3PLUS_FRAME_DURATION_1p25MS:
                /* 13763 = 3.36 * 2^12 */
                setup->ltpf_enable =
                    sub(extract_l(L_shr(L_mult0(13763, setup->total_bits), 12)), add(560, i_mult(80, encoder->fs_idx))) < 0;
                setup->enable_lpc_weighting = 0;
                BREAK;
#endif
            case LC3PLUS_FRAME_DURATION_2p5MS:
                /* 9830 = 2.4 * 2^12 */
                setup->ltpf_enable =
                    sub(extract_l(L_shr(L_mult0(9830, setup->total_bits), 12)), add(560, i_mult(80, encoder->fs_idx))) < 0;
                setup->enable_lpc_weighting = 0;
                BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:
                setup->ltpf_enable = sub(sub(i_mult(setup->total_bits, 2), 160), add(560, i_mult(80, encoder->fs_idx))) < 0;
                setup->enable_lpc_weighting = setup->total_bits < 240;
                BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:
                setup->ltpf_enable = sub(L_shr(L_mult0(10923, setup->total_bits), 13), add(560, i_mult(80, encoder->fs_idx))) < 0;
                setup->enable_lpc_weighting = setup->total_bits < 360;
                BREAK;
            case LC3PLUS_FRAME_DURATION_10MS:
                setup->enable_lpc_weighting = setup->total_bits < 480;
                setup->ltpf_enable          = sub(setup->total_bits, add(560, i_mult(80, encoder->fs_idx))) < 0;
                BREAK;
            case LC3PLUS_FRAME_DURATION_UNDEFINED:
                assert(0);
        }

#ifdef FIX_BOTH_1p25_WB_GLOBGAINOFFSET_NONBE  
        IF(encoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
        {
            setup->quantizedGainOff = calc_GGainOffset_1p25_fx(setup->total_bits, encoder->fs_idx); /* enc/dec common function */
        }
        ELSE
        {
            setup->quantizedGainOff =
                -(s_min(115, setup->total_bits / (10 * (encoder->fs_idx + 1))) + 105 + 5 * (encoder->fs_idx + 1));
        }
#else 

        setup->quantizedGainOff =
            -(s_min(115, setup->total_bits / (10 * (encoder->fs_idx + 1))) + 105 + 5 * (encoder->fs_idx + 1));
#endif

#    ifdef ENABLE_HR_MODE
    if (encoder->hrmode && encoder->fs_idx == 5)
    {
        setup->quantizedGainOff = MAX(setup->quantizedGainOff, -181);
    }
#    endif

#  ifdef ENABLE_HR_MODE
        if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS &&
            ((encoder->fs_in >= 44100 && setup->targetBytes >= 100) ||
             (encoder->fs_in == 32000 && setup->targetBytes >= 81))
            && setup->targetBytes < 340
            && (encoder->hrmode == 0))
#  else
        if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS &&
            ((encoder->fs_in >= 44100 && setup->targetBytes >= 100) ||
             (encoder->fs_in == 32000 && setup->targetBytes >= 81))
            && setup->targetBytes < 340
        )
#  endif
        {
            setup->attack_handling = 1;
        }
        else if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS && ((encoder->fs_in >= 44100 && setup->targetBytes >= 75) ||
                (encoder->fs_in == 32000 && setup->targetBytes >= 61)) && setup->targetBytes < 150
#  ifdef ENABLE_HR_MODE
                && encoder->hrmode == 0
#  endif
                 )
        {
            setup->attack_handling = 1;
        }
        else
        {
            /* reset attack detector for bitrate switching */
            setup->attack_handling      = 0;
            setup->attdec_filter_mem[0] = 0;
            setup->attdec_filter_mem[1] = 0;
            setup->attdec_detected      = 0;
            setup->attdec_position      = 0;
            setup->attdec_acc_energy    = 0;
            setup->attdec_scaling       = 0;
        }
        
#  ifdef ENABLE_HR_MODE
        if (encoder->hrmode)
        {
            setup->ltpf_enable = 0;
        }
#  endif
        encoder->sns_damping = SNS_DAMPING;

#  ifdef ENABLE_HR_MODE
        IF (encoder->hrmode)
        {
            encoder->sns_damping = SNS_DAMPING_HRMODE;
            IF (encoder->fs_idx >= 4) 
            {
                IF ((encoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS) & (setup->total_bits > 4400))
                {
                    encoder->sns_damping = SNS_DAMPING_HRMODE_UB_10MS;
                }
                IF ((encoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS) & (setup->total_bits > 3300))
                {
                    encoder->sns_damping = SNS_DAMPING_HRMODE_UB_7_5MS;
                }
                IF ((encoder->frame_dms == LC3PLUS_FRAME_DURATION_5MS) & (setup->total_bits > 2300))
                {
                    encoder->sns_damping = SNS_DAMPING_HRMODE_UB_5MS;
                }
                IF ((encoder->frame_dms == LC3PLUS_FRAME_DURATION_2p5MS) & (setup->total_bits > 1150))
                {
                    encoder->sns_damping = SNS_DAMPING_HRMODE_UB_2_5MS;
                }
            }
        }

        if (encoder->hrmode && encoder->fs_idx >= 4)
        {
            int real_rate  = setup->targetBytes * 8 * 10000 / (encoder->frame_dms * 1.25 * 10);
            setup->regBits = real_rate / 12500;

            if (encoder->fs_idx == 5)
            {
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS)
                {
                    setup->regBits += 2;
                }
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS)
                {
                    setup->regBits +=1;
                }
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_2p5MS)
                {
                    setup->regBits -= 6;
                }
            }
            else
            {
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_2p5MS)
                {
                    setup->regBits -= 6;
                }
                else if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_5MS)
                {
                    setup->regBits += 0;
                }
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS)
                {
                    setup->regBits +=2;
                }
                if (encoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS)
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
#  endif
    }

    encoder->bitrate = bitrate;

    return LC3PLUS_OK;
}

