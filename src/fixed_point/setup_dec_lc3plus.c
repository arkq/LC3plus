/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "setup_dec_lc3plus.h"

/* if decoder is null only size is reported */
/*assume 10ms  for memory allocation for now */
int alloc_decoder(LC3PLUS_Dec *decoder, int samplerate, int channels)
{
    int    ch         = 0;
    size_t size       = sizeof(LC3PLUS_Dec);
    void * ltpf_mem_x = NULL, *ltpf_mem_y = NULL, *stDec_ola_mem_fx = NULL;
    int    max_len = DYN_MAX_LEN_EXT(samplerate); /*NB  sing 80 as minimum value changes BE for NB */
    Word32 * plc_longterm_advc_tdc = NULL, *plc_longterm_advc_ns = NULL;
    Word16 longterm_analysis_counter_max = 0, longterm_analysis_counter_max_bytebuffer = 0;

    void *q_old_res_fx = NULL;
    Word16 *sharedBuf  = NULL;
#if defined(ENABLE_PC)
    void *q_old_d_fx = NULL;
#endif
    void *plcAd = NULL, *PhECU_f0est = NULL, *PhECU_plocs = NULL;
    
    for (ch = 0; ch < channels; ch++)
    {
        DecSetup *setup = balloc(decoder, &size, sizeof(DecSetup));
#ifdef ENABLE_HR_MODE
        int q_old_len = MIN(CODEC_FS(samplerate) * 100 / 10000, MAX_BW_HR);
#else
        int q_old_len = MIN(CODEC_FS(samplerate) * 100 / 10000, MAX_BW   );
#endif
        int x_old_len = DYN_MAX_LEN_PCM_PLC(samplerate); /* max(pitchmax + frame ms,  M+1 + pitchmax + frame/2) */
        int fs_idx    = (samplerate / 10000);            /* floor= integer truncation is needed  here */

        ltpf_mem_x = balloc(decoder, &size, sizeof(*setup->ltpf_mem_x) * (max_len + max_len / 40 - 2));
        ltpf_mem_y =
            balloc(decoder, &size,
                   sizeof(*setup->ltpf_mem_y) * (max_len + CEILING(MAX_PITCH_12K8 * max_len, 128) + (max_len / 80)));

        stDec_ola_mem_fx = balloc(decoder, &size, sizeof(*setup->stDec_ola_mem_fx) * DYN_MAX_MDCT_LEN(samplerate));

#if defined(ENABLE_PC)
        q_old_d_fx = balloc(decoder, &size, sizeof(*setup->q_old_d_fx) * q_old_len);
#endif
        sharedBuf  = balloc(decoder, &size, sizeof(*setup->q_old_d_fx) * q_old_len + sizeof(*setup->plcAd->x_old_tot_fx) * x_old_len);

        /* To save static RAM, a large buffer sharedBuf is used to share
         * space for q_old_d_fx, x_old_tot_fx and X_sav_fx (note:
         * PHASE_ECU operates only at 10ms and standard precision)
         *
         * If partial Concealment is active, q_old_d_fx is saved in a separate buffer
         *
         * The following graph provides an overview about the memory
         * sharing for 10ms and normal resolution:
         */

    /*    Total buffer     | <---     sharedBuf                                                    --->  |   */
    /* ConcealMethod 3or4  | <---q_old_d_fx---> | <-------------------- x_old_tot_fx ------------------> |   */
    /* ConcealMethod 3     | <---q_old_d_fx---> |      |<------- M+1+pitch_max+ (frame_len/2)----------> |   */ /* first BFI-frame only? */
    /* Meth 2,prevBfi=0    |                                |               <-------   16 ms xfp ------->|   */
    /* Meth 2              | ^ X_sav_fx ptr                 | <-------- maintained old_tot_fx ---------> |   */ /* pitchmax samples maintained */
    /* Meth 2              | <----- X_sav 16 ms ----------> |               ^xfp_fx ptr                  |   */
    /* Meth2  ,prevBfi=1   |                                |                <-12.25(forTDC),3.75(PLC2)->|   */ /*the part used by PhECU  is 3.75ms */
    /* ConcealMethod 2                          |<- <10ms ->|                                                */
    int max_plocs = DYN_MAX_PLOCS(samplerate);
    plcAd         = balloc(decoder, &size, sizeof(*setup->plcAd));
    PhECU_f0est   = balloc(decoder, &size, sizeof(*setup->plcAd->PhECU_f0est) * max_plocs);
    PhECU_plocs   = balloc(decoder, &size, sizeof(*setup->plcAd->PhECU_plocs) * max_plocs);
#ifdef ENABLE_HR_MODE
    q_old_res_fx = balloc(decoder, &size, sizeof(*setup->q_old_res_fx) * max_len);
#else
    q_old_res_fx = balloc(decoder, &size, sizeof(*setup->q_old_res_fx) * MIN(max_len, MAX_BW));
#endif
        
/*        longterm_analysis_counter_max = PLC_LONGTERM_ANALYSIS_MS * (100.0f / 25.0f);*/
/*        assert(longterm_analysis_counter_max == (PLC_LONGTERM_ANALYSIS_MS * (100 /25)));  */ /* test integer division for compatibility */  
/*        longterm_analysis_counter_max_bytebuffer = floor(longterm_analysis_counter_max / 30.0);  */
/*        assert(longterm_analysis_counter_max_bytebuffer == (longterm_analysis_counter_max / 30)); */ /* test integer division for compatibility */ 

/*        longterm_analysis_counter_max = (longterm_analysis_counter_max_bytebuffer+1) * 30;  */
        
        longterm_analysis_counter_max = plc_fadeout_param_maxlen[0];
        longterm_analysis_counter_max_bytebuffer = plc_fadeout_param_maxbytes[0];
        
        plc_longterm_advc_tdc = balloc(decoder, &size, sizeof(Word32) * longterm_analysis_counter_max_bytebuffer);
        plc_longterm_advc_ns  = balloc(decoder, &size, sizeof(Word32) * longterm_analysis_counter_max_bytebuffer);
        
        if (decoder)
        {
            decoder->channel_setup[ch] = setup;
            setup->ltpf_mem_x          = ltpf_mem_x;
            setup->ltpf_mem_y          = ltpf_mem_y;
            setup->stDec_ola_mem_fx    = stDec_ola_mem_fx;
            setup->q_old_res_fx        = q_old_res_fx;
#if defined(ENABLE_PC)
            setup->q_old_d_fx          = q_old_d_fx;
#else
            setup->q_old_d_fx          = sharedBuf;
#endif
            setup->plcAd = plcAd;
        }
        if (decoder && plcAd)
        {
#ifdef ENABLE_HR_MODE
            if (fs_idx > 4)
            {
                fs_idx = 4;
            }
#endif

            setup->plcAd->x_old_tot_fx   = &sharedBuf[q_old_len];
            setup->plcAd->PhECU_f0est    = PhECU_f0est;
            setup->plcAd->PhECU_xfp_fx   =  &(setup->plcAd->x_old_tot_fx[x_old_len-LprotSzPtr[fs_idx]]); /* point to the last 16ms of the  x_old_tot_fx      */
            setup->plcAd->PhECU_X_sav_fx =  sharedBuf; /* reuse  of lprot(=num_FsByResQ0[fs_idx]) values from this point fwd, i.e beyond the end of q_old_fx  */
            setup->plcAd->PhECU_plocs    = PhECU_plocs;
            
            setup->plcAd->longterm_analysis_counter_max = longterm_analysis_counter_max;
            setup->plcAd->longterm_analysis_counter_max_bytebuffer = longterm_analysis_counter_max_bytebuffer;
                
            setup->plcAd->plc_longterm_advc_tdc = plc_longterm_advc_tdc;
            setup->plcAd->plc_longterm_advc_ns  = plc_longterm_advc_ns;
        }
    }

    return (int)size;
}

LC3PLUS_Error FillDecSetup(LC3PLUS_Dec *decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode
#ifdef ENABLE_HR_MODE
                           , int hrmode
#endif
                          )
{
    int ch = 0;

    memset(decoder, 0, lc3plus_dec_get_size(samplerate, channels, plc_mode));
    alloc_decoder(decoder, samplerate, channels);
    
#ifdef ENABLE_HR_MODE
    decoder->hrmode = hrmode != 0;
#endif

    decoder->fs     = CODEC_FS(samplerate);
    decoder->fs_out = samplerate;
    decoder->fs_idx = FS2FS_IDX(decoder->fs);
    decoder->channels  = channels;
    decoder->frame_dms = 100;
    decoder->plcMeth   = plc_mode;
    {
        decoder->BW_cutoff_bits = BW_cutoff_bits_all[decoder->fs_idx];
    }
    decoder->ltpf_mem_x_len = extract_l(L_shr_pos(Mpy_32_16(L_max(16000, decoder->fs), 16778), 11)) - 2;
    decoder->ltpf_mem_y_len = extract_l(L_shr_pos(Mpy_32_16(decoder->fs, 18678) - 1, 5)) + 1 +
                              extract_l(L_shr_pos(Mpy_32_16(L_max(16000, decoder->fs), 16778), 12));

    set_dec_frame_params(decoder);

    for (ch = 0; ch < decoder->channels; ch++)
    {
        DecSetup *setup               = decoder->channel_setup[ch];
        setup->plc_damping            = 32767;
        setup->ltpf_mem_scale_fac_idx = -1;
        move16();

        setup->pc_seed = 24607;
        setup->ns_seed      = 24607;
        setup->ns_cum_alpha = 32767;

    int i = 0;
    /* 0 = 0kHz  1= 8kHz, 2= 16 kHz   3= 24 , 4 = 32  5=40 6=48kHz */
#if (PHECU_XFP_LA == 0)
    Word16 oneMsTab_LA[5] = {0, 0, 0, 0, 0};
#else
#  if (PHECU_XFP_LA == 4)
    Word16 oneMsTab_LA[7] = {0 /*unused*/, 2, 4, 6, 8, 10 /*unused*/, 12};
#  else
    Word16 oneMsTab_LA[7] = {0 /*unused*/, 8, 16, 24, 32, 40 /*unused*/, 48};
#  endif
#endif
    Word16 oneMsTab[5] = {8, 16, 24, 32, 48};

    setup->plcAd->stab_fac        = 32767;
    setup->plcAd->tdc_seed        = 24607;
    setup->plcAd->tdc_preemph_fac = plc_preemph_fac[decoder->fs_idx];
    setup->plcAd->tdc_lpc_order   = 16;
    setup->plcAd->PhECU_fs_idx_fx =
        mult(decoder->frame_length,
         (Word16)(32768.0 / 99.0)); /* truncation needed , i.e no rounding can be applied here */
    /* idx=frame/80,  0=8kHZ, 1=16kHz, 2=24 kHz, 3=32 kHz 5=*, 4=48 */

    setup->plcAd->max_len_pcm_plc     = DYN_MAX_LEN_PCM_PLC(decoder->fs);
    
    if ((decoder->hrmode == 0) && (samplerate <= 48000))
    {
        setup->plcAd->PhECU_frame_ms = (Word16)(
            decoder->frame_dms * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
        setup->plcAd->PhECU_seed_fx = 21845;
        setup->plcAd->PhECU_LprotOrg_fx =
            shl_pos(oneMsTab[setup->plcAd->PhECU_fs_idx_fx], 4); /* 16 *1ms = 1.6 *framelength */
        setup->plcAd->PhECU_Lprot_fx      = setup->plcAd->PhECU_LprotOrg_fx;
        setup->plcAd->PhECU_LA            = oneMsTab_LA[setup->plcAd->PhECU_fs_idx_fx];
        setup->plcAd->PhECU_whr_tot_taper = sub(setup->plcAd->PhECU_Lprot_fx, decoder->frame_length); /* 3+3 ms */
        setup->plcAd->PhECU_whr_tot_flat  = decoder->frame_length;                                    /* 10 ms */
        setup->plcAd->PhECU_LDWIN_OLAP    = shr_pos(decoder->frame_length, 2);                        /* 2.5 ms */
        setup->plcAd->max_lprot           = DYN_MAX_LPROT(decoder->fs);
        setup->plcAd->max_plocs           = DYN_MAX_PLOCS(decoder->fs);
        setup->plcAd->PhECU_margin_xfp    = 0;
        setup->plcAd->PhECU_L_oold_xfp_w_E_fx = LTOT_MIN_MAN;
        setup->plcAd->PhECU_L_old_xfp_w_E_fx  = LTOT_MIN_MAN;

        setup->plcAd->PhECU_oold_xfp_w_E_exp_fx = UNINIT_OR_UNSAFE_OOLD_SENTINEL;
        setup->plcAd->PhECU_old_xfp_w_E_exp_fx  = LTOT_INIT_FLAG;

        setup->plcAd->PhECU_oold_Ltot_exp_fx = UNINIT_OR_UNSAFE_OOLD_SENTINEL;
        setup->plcAd->PhECU_old_Ltot_exp_fx  = LTOT_INIT_FLAG;


        for (i = 0; i < MAX_LGW; i++)
        {
            setup->plcAd->PhECU_oold_grp_shape_fx[i] =
            GRP_SHAPE_INIT; /* negative value will be replaced in very first calculation */
            setup->plcAd->PhECU_old_grp_shape_fx[i] = GRP_SHAPE_INIT;
        }
        /* t_adv=ctrl.FRAME/2 + ctrl.PhECU.LprotOrg/2 - ctrl.PhECU.LA + LDWIN_OLAP/2; */
        i = add(add(decoder->frame_length, setup->plcAd->PhECU_LprotOrg_fx), setup->plcAd->PhECU_LDWIN_OLAP);
        setup->plcAd->PhECU_t_adv = sub(shr_pos(i, 1), setup->plcAd->PhECU_LA);

        for (i = 0; i < MAX_LGW; i++)
        {
            setup->plcAd->PhECU_mag_chg_1st[i] = 32767;
        }
        setup->plcAd->PhECU_beta_mute = 16384;
		setup->plcAd->PhECU_nonpure_tone_flag = -1;
        }
    }

    lc3plus_dec_set_ep_enabled(decoder, 0);
    return LC3PLUS_OK;
}

/* set frame config params */
void set_dec_frame_params(LC3PLUS_Dec *decoder)
{
    Word16 tmp = 0;
    Word16 n;
    
    decoder->frame_length = extract_l(L_shr_pos(Mpy_32_16(decoder->fs, 20972), 6)); /* fs * 0.01*2^6 */
    
#ifdef ENABLE_HR_MODE
    if (decoder->hrmode)
    {
        decoder->BW_cutoff_bits = 0;
    }
#endif

    SWITCH (decoder->frame_dms)
    {
    case 25:
        decoder->frame_length         = shr_pos(decoder->frame_length, 2);
        decoder->la_zeroes    = LowDelayShapes_n960_la_zeroes_2_5ms[decoder->fs_idx];
        decoder->stDec_ola_mem_fx_len = sub(decoder->frame_length, decoder->la_zeroes);
        
#ifdef ENABLE_HR_MODE
        if (decoder->hrmode)
        {
            decoder->bands_number = bands_number_2_5ms_HR[decoder->fs_idx];
            decoder->bands_offset = bands_offset_2_5ms_HR[decoder->fs_idx - 4];
            decoder->W_fx         = LowDelayShapes_n960_HRA_2_5ms[decoder->fs_idx - 4];
            decoder->W_size       = LowDelayShapes_n960_len_2_5ms[decoder->fs_idx];
            decoder->yLen         = decoder->frame_length;            
        }
        else
#endif
        {
            decoder->bands_number = bands_number_2_5ms[decoder->fs_idx];
            decoder->bands_offset = bands_offset_2_5ms[decoder->fs_idx];
            decoder->W_fx         = LowDelayShapes_n960_2_5ms[decoder->fs_idx];
            decoder->W_size       = LowDelayShapes_n960_len_2_5ms[decoder->fs_idx];
            decoder->yLen         = s_min(MAX_BW >> 2, decoder->frame_length);            
        }
        
        if (decoder->fs_idx == 0)
        {
            int ch;
            for (ch = 0; ch < decoder->channels; ch++)
            {
                DecSetup *setup = decoder->channel_setup[ch];
                if (setup->plcAd != NULL)
                {
                    setup->plcAd->tdc_lpc_order = 8;
                }
            }
        }
        BREAK;
    case 50:
        decoder->frame_length         = shr_pos(decoder->frame_length, 1);
        decoder->la_zeroes    = LowDelayShapes_n960_la_zeroes_5ms[decoder->fs_idx];
        decoder->stDec_ola_mem_fx_len = sub(decoder->frame_length, decoder->la_zeroes);
        
#ifdef ENABLE_HR_MODE
        if (decoder->hrmode)
        {
            decoder->bands_offset = bands_offset_5ms_HR[decoder->fs_idx - 4];
            decoder->W_fx         = LowDelayShapes_n960_HRA_5ms[decoder->fs_idx - 4];
            decoder->W_size       = LowDelayShapes_n960_len_5ms[decoder->fs_idx];
            decoder->yLen         = decoder->frame_length;
            decoder->bands_number = bands_number_5ms[decoder->fs_idx];
        }
        else
#endif
        {
            decoder->bands_offset = bands_offset_5ms[decoder->fs_idx];
            decoder->W_fx         = LowDelayShapes_n960_5ms[decoder->fs_idx];
            decoder->W_size       = LowDelayShapes_n960_len_5ms[decoder->fs_idx];
            decoder->yLen         = s_min(MAX_BW >> 1, decoder->frame_length);
            decoder->bands_number         = bands_number_5ms[decoder->fs_idx];
        }
        BREAK;
        
    case 75:
        tmp                           = shr_pos(decoder->frame_length, 2);
        decoder->frame_length         = add(tmp, add(tmp, tmp));
        decoder->la_zeroes            = LowDelayShapes_n960_la_zeroes_7_5ms[decoder->fs_idx];
        decoder->stDec_ola_mem_fx_len = sub(decoder->frame_length, decoder->la_zeroes);
#    ifdef ENABLE_HR_MODE
        if (decoder->hrmode)
        {
            decoder->bands_offset = bands_offset_7_5ms_HR[decoder->fs_idx - 4];
            decoder->W_fx         = LowDelayShapes_n960_HRA_7_5ms[decoder->fs_idx - 4];
            decoder->W_size       = LowDelayShapes_n960_len_7_5ms[decoder->fs_idx];
            decoder->yLen         = decoder->frame_length;
            decoder->bands_number = bands_number_7_5ms[decoder->fs_idx];
        }
        else
#    endif
        {
            decoder->bands_offset = bands_offset_7_5ms[decoder->fs_idx];
            decoder->W_fx         = LowDelayShapes_n960_7_5ms[decoder->fs_idx];
            decoder->W_size       = LowDelayShapes_n960_len_7_5ms[decoder->fs_idx];
            decoder->yLen         = s_min((MAX_BW >> 2) * 3, decoder->frame_length);
            decoder->bands_number = bands_number_7_5ms[decoder->fs_idx];
        }
        BREAK;
        
    case 100:
        decoder->la_zeroes            = LowDelayShapes_n960_la_zeroes[decoder->fs_idx];
        decoder->stDec_ola_mem_fx_len = sub(decoder->frame_length, decoder->la_zeroes);
        decoder->bands_number         = 64;
        
#ifdef ENABLE_HR_MODE
        if (decoder->hrmode)
        {
            decoder->bands_offset = bands_offset_HR[decoder->fs_idx - 4];
            decoder->W_fx         = LowDelayShapes_n960_HRA[decoder->fs_idx - 4];
            decoder->W_size       = LowDelayShapes_n960_len[decoder->fs_idx];
            decoder->yLen         = decoder->frame_length;
        }
        else
#endif
        {
            decoder->W_fx         = LowDelayShapes_n960[decoder->fs_idx];
            decoder->W_size       = LowDelayShapes_n960_len[decoder->fs_idx];
            decoder->bands_offset = bands_offset[decoder->fs_idx];
            decoder->yLen         = s_min(MAX_BW, decoder->frame_length);
        }
        BREAK;
    }

    {
        int ch;
        for (ch = 0; ch < decoder->channels; ch++)
        {
            DecSetup *setup = decoder->channel_setup[ch];
            if (setup->plcAd != NULL)
            { /*only set  if plcAd was actually allocated */
                setup->plcAd->longterm_analysis_counter_max = plc_fadeout_param_maxlen[(decoder->frame_dms / 25) - 1];
                setup->plcAd->longterm_analysis_counter_max_bytebuffer = plc_fadeout_param_maxbytes[(decoder->frame_dms / 25) - 1];
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    decoder->frame_dms *
                    0.1); /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
            }
        }
    }
    FOR (n=0; n < PLC_FADEOUT_TYPE_1_IN_MS*10/decoder->frame_dms;n++){
        decoder->alpha_type_2_table[n] = type_2_fadeout_fx(n, decoder->frame_dms);
    }
}

LC3PLUS_Error update_dec_bitrate(LC3PLUS_Dec *decoder, int ch, Word16 nBytes)
{
    int tmp = 0, totalBits = 0;
    int channel_bytes = 0;
    int minBytes, maxBytes;

    DecSetup *setup = decoder->channel_setup[ch];
    channel_bytes   = nBytes;
    
#ifdef ENABLE_HR_MODE
    if (decoder->hrmode)
    {
        SWITCH (decoder->frame_dms)
        {
        case 25:
            maxBytes = 210;
            minBytes = MIN_NBYTES;
            BREAK;
        case 50:
            maxBytes = 375;
            minBytes = MIN_NBYTES;
            BREAK;
        case 75:
            maxBytes = 625;
            minBytes = MIN_NBYTES;
            BREAK;
        case 100:
            maxBytes = 625;
            minBytes = MIN_NBYTES;
            BREAK;
        default: return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
#endif
    {
        minBytes = MIN_NBYTES;
        maxBytes = MAX_NBYTES_100;
    }
    
    if (channel_bytes < minBytes || channel_bytes > maxBytes)
    {
        return LC3PLUS_NUMBYTES_ERROR;
    }

    setup->targetBytes = channel_bytes;
    move16();
    setup->total_bits           = shl(setup->targetBytes, 3);
    setup->enable_lpc_weighting = (setup->total_bits < 480);
    setup->quantizedGainOff =
        -(s_min(115, setup->total_bits / (10 * (decoder->fs_idx + 1))) + 105 + 5 * (decoder->fs_idx + 1));
    tmp       = i_mult(80, decoder->fs_idx);
    
#ifdef ENABLE_HR_MODE
    if (decoder->hrmode && decoder->fs_idx == 5)
    {
        setup->quantizedGainOff = MAX(setup->quantizedGainOff, -181);
    }
#endif
    
    totalBits = setup->total_bits;


    SWITCH (decoder->frame_dms)
    {
    case 25:
        setup->enable_lpc_weighting = 0;
        /* total_bits * 2.4 */
        totalBits = extract_l(L_shr(L_mult0(19661, setup->total_bits), 13));
        BREAK;
    case 50:
        setup->enable_lpc_weighting = setup->total_bits < 240;
        totalBits                   = sub(i_mult(setup->total_bits, 2), 160);
        BREAK;
    case 75:
        setup->enable_lpc_weighting = setup->total_bits < 360;
        totalBits = L_shr(L_mult0(10923, setup->total_bits), 13);
        BREAK;
    case 100:
        BREAK;
    }


    if (sub(totalBits, add(320, tmp)) < 0)
    {
        setup->ltpf_scale_fac_idx = 0;
        move16();
    }
    else if (sub(totalBits, add(400, tmp)) < 0)
    {
        setup->ltpf_scale_fac_idx = 1;
        move16();
    }
    else if (sub(totalBits, add(480, tmp)) < 0)
    {
        setup->ltpf_scale_fac_idx = 2;
        move16();
    }
    else if (sub(totalBits, add(560, tmp)) < 0)
    {
        setup->ltpf_scale_fac_idx = 3;
        move16();
    }
    else
    {
        setup->ltpf_scale_fac_idx = -1;
        move16();
    }
    
#ifdef ENABLE_HR_MODE
    /* No LTPF in hrmode */
    if (decoder->hrmode)
    {
        setup->ltpf_scale_fac_idx = -1;
        move16();
    }
#endif

    return LC3PLUS_OK;
}

