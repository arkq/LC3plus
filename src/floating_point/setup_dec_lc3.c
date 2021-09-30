/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "setup_dec_lc3.h"
#include "functions.h"
#include <stdio.h>
#include <assert.h>

/* if decoder is null only size is reported */
#  include "fft/iis_fft.h"

int alloc_decoder(LC3PLUS_Dec* decoder, int samplerate, int channels)
{
    int    ch   = 0;
    size_t size = sizeof(LC3PLUS_Dec);
   size_t frame_len = DYN_MAX_LEN_EXT(samplerate);

   void *PlcAdvSetup = NULL;
   LC3_FLOAT *pcmbufHist, *harmonicBuf;
   LC3_FLOAT *PhECU_oold_grp_shape, *PhECU_old_grp_shape;
   LC3_FLOAT *PhECU_xfp;
   Complex   *PhECU_X_sav_m;
   LC3_INT32   *PhECU_plocs;
   LC3_FLOAT *PhECU_f0est, *PhECU_mag_chg_1st, *PhECU_Xavg;
   LC3_FLOAT *sine_table1_phecu, *sine_table2_phecu;
   HANDLE_IIS_FFT handle_fft_phaseecu;
   HANDLE_IIS_FFT handle_ifft_phaseecu;

    for (ch = 0; ch < channels; ch++) {
        DecSetup* setup = balloc(decoder, &size, sizeof(DecSetup));

         size_t max_pitch = ceilf(228.0 * CODEC_FS(samplerate) / 12800.0);
         size_t pcm_plc_len = max_pitch + frame_len;
         pcmbufHist = balloc(decoder, &size, sizeof(LC3_FLOAT) * pcm_plc_len);
         harmonicBuf = balloc(decoder, &size, sizeof(LC3_FLOAT) * max_pitch);
         PlcAdvSetup = balloc(decoder, &size, sizeof(*setup->PlcAdvSetup));
         PhECU_oold_grp_shape = balloc(decoder, &size, sizeof(LC3_FLOAT) *MAX_LGW);          /* BASOP Word16  PhECU_oold_grp_shape_fx[MAX_LGW]; */
         PhECU_old_grp_shape = balloc(decoder, &size, sizeof(LC3_FLOAT) *MAX_LGW);          /* BASOP Word16  PhECU_old_grp_shape_fx[MAX_LGW] ; */
         PhECU_xfp = balloc(decoder, &size, sizeof(LC3_FLOAT) *(frame_len * 16 / 10));
         PhECU_X_sav_m = balloc(decoder, &size, sizeof(Complex) *(((frame_len * 16 / 10) / 2) + 1));/*MAX_PLC_LMSPEC*/
         PhECU_plocs = balloc(decoder, &size, sizeof(LC3_INT32) * (((frame_len * 16 / 10) / 4) + 1 + 1));   /* BASOP Word16 *PhECU_plocs;    */
          
         handle_fft_phaseecu = balloc(decoder, &size, sizeof(IIS_FFT) * 1);
         handle_ifft_phaseecu = balloc(decoder, &size, sizeof(IIS_FFT) * 1);
         PhECU_f0est = balloc(decoder, &size, sizeof(LC3_FLOAT) * (((frame_len * 16 / 10) / 4) + 1));        /*BASOP Word32 *PhECU_f0est;*/
         PhECU_mag_chg_1st = balloc(decoder, &size, sizeof(LC3_FLOAT) *MAX_LGW);                /*  BASOP Word16  PhECU_mag_chg_1st[MAX_LGW];*/
         PhECU_Xavg = balloc(decoder, &size, sizeof(LC3_FLOAT) * MAX_LGW);               /*  BASOP Word16  PhECU_Xavg[MAX_LGW] ; */
          
         sine_table1_phecu = balloc(decoder, &size, sizeof(LC3_FLOAT) * (((CODEC_FS(samplerate) * 16) / 1000) / 2 + 1));
         sine_table2_phecu = balloc(decoder, &size, sizeof(LC3_FLOAT) * (((CODEC_FS(samplerate) * 16) / 1000) / 2 + 1));

         LC3_FLOAT *q_old_res = balloc(decoder, &size, sizeof(LC3_FLOAT) * frame_len);

        if (decoder) {
            decoder->channel_setup[ch] = setup;

            setup->PlcAdvSetup = PlcAdvSetup;

            setup->PlcAdvSetup->pcmbufHist = pcmbufHist;
            setup->PlcAdvSetup->PlcTdcSetup.harmonicBuf = harmonicBuf;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_oold_grp_shape = PhECU_oold_grp_shape;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_old_grp_shape = PhECU_old_grp_shape;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp = PhECU_xfp;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_X_sav_m = PhECU_X_sav_m;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_plocs = PhECU_plocs;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_f0est = PhECU_f0est;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_mag_chg_1st = PhECU_mag_chg_1st;
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Xavg = PhECU_Xavg;
            setup->PlcAdvSetup->PlcPhEcuSetup.handle_fft_phaseecu = handle_fft_phaseecu;
            setup->PlcAdvSetup->PlcPhEcuSetup.handle_ifft_phaseecu = handle_ifft_phaseecu;
             
            setup->PlcAdvSetup->PlcPhEcuSetup.handle_fft_phaseecu->sine_table = sine_table1_phecu;
            setup->PlcAdvSetup->PlcPhEcuSetup.handle_ifft_phaseecu->sine_table = sine_table2_phecu;
             
            setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot = (CODEC_FS(samplerate) * 16) / 1000;
            real_fft_init(&(setup->PlcAdvSetup->PlcPhEcuSetup.PhEcu_Fft), setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot, &(setup->PlcAdvSetup->PlcPhEcuSetup.handle_fft_phaseecu));
            real_ifft_init(&(setup->PlcAdvSetup->PlcPhEcuSetup.PhEcu_Ifft), setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot, &(setup->PlcAdvSetup->PlcPhEcuSetup.handle_ifft_phaseecu));
            setup->statePC.q_old_res = q_old_res;
        }
    }

    return (int)size;
}

LC3PLUS_Error FillDecSetup(LC3PLUS_Dec* decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode)
{
    memset(decoder, 0, lc3plus_dec_get_size(samplerate, channels));
    alloc_decoder(decoder, samplerate, channels);

    decoder->fs     = CODEC_FS(samplerate);
    decoder->fs_out = samplerate;
    decoder->fs_idx = FS2FS_IDX(decoder->fs);
    decoder->plcMeth = plc_mode;
    
    
    if (decoder->fs_idx > 4) {
        decoder->fs_idx = 5;
    }
    decoder->channels       = channels;
    decoder->frame_ms       = 10;
    decoder->frame_dms      = 100;
    decoder->BW_cutoff_bits = BW_cutoff_bits_all[decoder->fs_idx];
    
    if (decoder->fs == 8000) {
        decoder->tilt = 14;
    } else if (decoder->fs == 16000) {
        decoder->tilt = 18;
    } else if (decoder->fs == 24000) {
        decoder->tilt = 22;
    } else if (decoder->fs == 32000) {
        decoder->tilt = 26;
    } else if (decoder->fs == 48000) {
        decoder->tilt = 30;
    }
    else if (decoder->fs == 96000) {
        decoder->tilt = 34;
    }

    set_dec_frame_params(decoder);

    lc3plus_dec_set_ep_enabled(decoder, 0);
    
    return LC3PLUS_OK;
}

/* set frame config params */
void set_dec_frame_params(LC3PLUS_Dec* decoder)
{
    int ch = 0;
    
    if (decoder->fs_idx == 5)
    {
        decoder->hrmode = 1;
    }

    decoder->frame_length = ceil(decoder->fs * 10 / 1000); /* fs * 0.01*2^6 */
    if (decoder->hrmode == 1)
    {
        decoder->yLen = decoder->frame_length;
    }
    else
    {
        decoder->yLen = MIN(MAX_BW, decoder->frame_length);
    }

    decoder->bands_number = 64;
    if (decoder->frame_ms == 2.5) 
    {
        decoder->frame_length = decoder->frame_length >> 2;
        decoder->yLen /= 4;
        if (decoder->hrmode)
        {
            decoder->bands_number = bands_number_2_5ms_HR[decoder->fs_idx];
        } 
        else
        {
            decoder->bands_number = bands_number_2_5ms[decoder->fs_idx];
        }
    }
    if (decoder->frame_ms == 5) 
    {
        decoder->frame_length = decoder->frame_length >> 1;
        decoder->yLen /= 2;
        decoder->bands_number = bands_number_5ms[decoder->fs_idx];
    }

    if (decoder->hrmode)
    {
        decoder->BW_cutoff_bits    = 0;
    }
    else
    {
        decoder->BW_cutoff_bits    = BW_cutoff_bits_all[decoder->fs_idx];
    }

    if (decoder->frame_ms == 10) 
    {
        if (decoder->hrmode)
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND_HR[decoder->fs_idx];
        }
        else
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND[decoder->fs_idx];
        }
        decoder->cutoffBins   = BW_cutoff_bin_all;
    }
    else if (decoder->frame_ms == 2.5) 
    {
        if (decoder->hrmode)
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms_HR[decoder->fs_idx];
        }
        else
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND_2_5ms[decoder->fs_idx];
        }
        decoder->cutoffBins   = BW_cutoff_bin_all_2_5ms;
    }
    else if (decoder->frame_ms == 5) 
    {
        if (decoder->hrmode)
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND_5ms_HR[decoder->fs_idx];
        }
        else
        {
            decoder->bands_offset = ACC_COEFF_PER_BAND_5ms[decoder->fs_idx];
        }
        decoder->cutoffBins   = BW_cutoff_bin_all_5ms;
    }
    
    decoder->n_bandsPLC = MIN(decoder->frame_length, 80);
    
    if (decoder->frame_ms == 10)
    {
        decoder->bands_offsetPLC = ACC_COEFF_PER_BAND_PLC[decoder->fs_idx];
    }
    else if (decoder->frame_ms == 5)
    {
        decoder->bands_offsetPLC = ACC_COEFF_PER_BAND_PLC_5ms[decoder->fs_idx];

        if (decoder->fs == 24000)
        {
            decoder->n_bandsPLC = 40;
        }
    }
    else if (decoder->frame_ms == 2.5)
    {
        decoder->bands_offsetPLC = ACC_COEFF_PER_BAND_PLC_2_5ms[decoder->fs_idx];

        if (decoder->fs == 48000)
        {
            decoder->n_bandsPLC = 60;
        }
    }
    assert(decoder->bands_offsetPLC);

    if (decoder->frame_ms == 10) {
        decoder->imdct_win     = MDCT_WINS_10ms[decoder->hrmode][decoder->fs_idx];
        decoder->imdct_laZeros = MDCT_la_zeroes[decoder->fs_idx];
        decoder->imdct_winLen  = MDCT_WINDOWS_LENGTHS_10ms[decoder->fs_idx];
    }
    else if (decoder->frame_ms == 2.5) {
        decoder->imdct_win     = MDCT_WINS_2_5ms[decoder->hrmode][decoder->fs_idx];
        decoder->imdct_laZeros = MDCT_la_zeroes_2_5ms[decoder->fs_idx];
        decoder->imdct_winLen  = MDCT_WINDOWS_LENGTHS_2_5ms[decoder->fs_idx];
    }
    else if (decoder->frame_ms == 5) {
        decoder->imdct_win     = MDCT_WINS_5ms[decoder->hrmode][decoder->fs_idx];
        decoder->imdct_laZeros = MDCT_la_zeroes_5ms[decoder->fs_idx];
        decoder->imdct_winLen  = MDCT_WINDOWS_LENGTHS_5ms[decoder->fs_idx];
    }

    decoder->la_zeroes = decoder->imdct_laZeros;

    decoder->imdct_memLen = decoder->frame_length - decoder->imdct_laZeros;

    for (ch = 0; ch < decoder->channels; ch++) {
        DecSetup* setup = decoder->channel_setup[ch];
        
        setup->ltpf_mem_beta_idx = -1;
        
        setup->statePC.seed = 24607;

        if (decoder) {
            /* Init DCT4 structs */
            if (setup->dct4structImdct.length != 0) {
                dct4_free(&setup->dct4structImdct);
                dct4_init(&setup->dct4structImdct, decoder->frame_length);
            } else {
                dct4_init(&setup->dct4structImdct, decoder->frame_length);
            }
            
            setup->PlcNsSetup.cum_alpha = 1;
            setup->PlcNsSetup.seed = 24607;
            setup->alpha = 1;
            if (setup->PlcAdvSetup)
            {
            LC3_INT32 pitch_max = 0, pitch_ana_len = 0, tdc_synt_len = 0;
            pitch_max     = ceil(228.0 * (LC3_FLOAT) decoder->fs / 12800.0);
                pitch_ana_len = pitch_max + decoder->frame_length * (LC3_FLOAT) 100 / decoder->frame_dms;
                tdc_synt_len  = 16 + 1 + pitch_max  + ceil(decoder->frame_length / 2);
                setup->PlcAdvSetup->max_len_pcm_plc = MAX(pitch_ana_len, tdc_synt_len);
                setup->PlcAdvSetup->PlcTdcSetup.preemphFac = plc_preemph_fac[decoder->fs_idx];
                setup->PlcAdvSetup->PlcTdcSetup.seed = 24607;
                setup->PlcAdvSetup->PlcTdcSetup.lpcorder = 16;

                if (decoder->fs_idx == 0 && decoder->frame_dms == 25)
                {
                    setup->PlcAdvSetup->PlcTdcSetup.lpcorder = 8;
                }

                setup->PlcAdvSetup->stabFac = 1;
                setup->PlcAdvSetup->cum_fading_fast = 1;
                setup->PlcAdvSetup->cum_fading_slow = 1;
                setup->PlcAdvSetup->cum_fflcAtten = 1;

                if (decoder->fs_idx <= 4 && decoder->frame_dms == 100)
                {
                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot = (decoder->fs * 16) / 1000;  /* 16 ms of samples at fs*/

                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_f0hzLtpBin = 0;
                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_norm_corr = 0;

                   set_vec(PHECU_GRP_SHAPE_INIT, setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_oold_grp_shape, MAX_LGW);
                   set_vec(PHECU_GRP_SHAPE_INIT, setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_old_grp_shape, MAX_LGW);

                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_L_oold_xfp_w_E = (LC3_FLOAT)PHECU_LTOT_MIN_MAN * LC3_POW(2.0, PHECU_LTOT_MIN_EXP - 31);
                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_L_old_xfp_w_E = (LC3_FLOAT)PHECU_LTOT_MIN_MAN * LC3_POW(2.0, PHECU_LTOT_MIN_EXP - 31);

                   /* CFL uses separate buffers for pcmHist,  xfp and Xsav and q_d , BASOP uses an optimized joint buffer*/
                   zero_float(setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp, setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot);
                   zero_cmplx(setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_X_sav_m, (setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot/2 + 1));

                   set_vec(POS_ONE_Q15, setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_mag_chg_1st,  MAX_LGW);

                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_beta_mute = (16384.0/32768.0);
                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_seed = 21845;

                   assert(decoder->frame_dms == 100);
                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_LDWIN_OLAP =  (decoder->frame_length / 4 );   /* 2.5 ms for regular 10 ms MDCT */

                   setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_t_adv = (
                      decoder->frame_length
                      + setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot
                      + setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_LDWIN_OLAP )/ 2;
                }
                setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_short_flag_prev  = 0;    /* fullband transient  */
                setup->PlcAdvSetup->PlcPhEcuSetup.PhECU_num_plocs = 0;
            }
        }
    }
}

LC3PLUS_Error update_dec_bitrate(LC3PLUS_Dec* decoder, int ch, int nBytes)
{
    int totalBits = 0, bitsTmp = 0, channel_bytes = 0, maxBytes = 0, minBytes = 0;

    if (decoder->hrmode)
    {
        switch (decoder->frame_dms)
        {
        case 25:
            maxBytes = 210;
            minBytes = MIN_NBYTES;
            break;
        case 50:
            maxBytes = 375;
            minBytes = MIN_NBYTES;
            break;
        case 100:
            maxBytes = 625;
            minBytes = MIN_NBYTES;
            break;
        default:
            return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
    {
        minBytes = MIN_NBYTES;
        maxBytes = MAX_NBYTES_100; // for backward compatibility, MAX_NBYTES_100 is used for all frame lengths
    }

    channel_bytes = nBytes;

        DecSetup* setup = decoder->channel_setup[ch];

        if (channel_bytes < minBytes || channel_bytes > maxBytes)
        {
            return LC3PLUS_NUMBYTES_ERROR;
        }
    
        setup->targetBytes          = channel_bytes;
        setup->total_bits           = setup->targetBytes << 3;
        setup->enable_lpc_weighting = (setup->total_bits < 480);
        setup->quantizedGainOff =
            -(MIN(115, setup->total_bits / (10 * (decoder->fs_idx + 1))) + 105 + 5 * (decoder->fs_idx + 1));

        if (decoder->hrmode && decoder->fs_idx == 5)
        {
            setup->quantizedGainOff = MAX(setup->quantizedGainOff, -181);
        }

        totalBits = setup->total_bits;
        
        if (decoder->frame_ms == 2.5) {
            setup->enable_lpc_weighting = setup->total_bits < 120;
            totalBits                   = setup->total_bits * 4.0 * (1.0 - 0.4);
        }
        if (decoder->frame_ms == 5) {
            setup->enable_lpc_weighting = (setup->total_bits < 240);
            totalBits                   = setup->total_bits * 2 - 160;
        }
    
        if (decoder->frame_length > 40 * ((LC3_FLOAT) (decoder->frame_dms) / 10.0)) {
            setup->N_red_tns  = 40 * ((LC3_FLOAT) (decoder->frame_dms) / 10.0);
            setup->fs_red_tns = 40000;
        } else {
            setup->N_red_tns = decoder->frame_length;
            setup->fs_red_tns = decoder->fs;
        }

        bitsTmp = totalBits;
        
        if (bitsTmp < 400 + (decoder->fs_idx - 1) * 80) {
            setup->ltpf_conf_beta     = 0.4;
            setup->ltpf_conf_beta_idx = 0;
        } else if (bitsTmp < 480 + (decoder->fs_idx - 1) * 80) {
            setup->ltpf_conf_beta     = 0.35;
            setup->ltpf_conf_beta_idx = 1;
        } else if (bitsTmp < 560 + (decoder->fs_idx - 1) * 80) {
            setup->ltpf_conf_beta     = 0.3;
            setup->ltpf_conf_beta_idx = 2;
        } else if (bitsTmp < 640 + (decoder->fs_idx - 1) * 80) {
            setup->ltpf_conf_beta     = 0.25;
            setup->ltpf_conf_beta_idx = 3;
        } else {
            setup->ltpf_conf_beta     = 0;
            setup->ltpf_conf_beta_idx = -1;
        }

        /* No LTPF in hrmode */
        if (decoder->hrmode == 1) {
            setup->ltpf_conf_beta     = 0;
            setup->ltpf_conf_beta_idx = -1;
        }
    
    return LC3PLUS_OK;
}
