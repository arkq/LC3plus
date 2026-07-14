/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
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
#ifdef CR14_A_ADD_LOSSLESS_MODE
                           , int wavFormat
#endif
                          )
{
    int ch = 0;

    memset(decoder, 0, lc3plus_dec_get_size(samplerate, channels, plc_mode));
    alloc_decoder(decoder, samplerate, channels);
    
#ifdef CR14_A_ADD_LOSSLESS_MODE
    decoder->wavFormat = wavFormat;
    
    IF (hrmode == 2)
    {
        decoder->lossless = 1;
        decoder->hrmode = 1;
    
        decoder->ll_tns = 1;
        decoder->ll_tns_remove = 1;
    }
#endif
    
#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (!decoder->lossless)
    {
        decoder->hrmode = hrmode != 0;
    }
#else
    decoder->hrmode = hrmode != 0;
#endif
#endif

    decoder->fs     = CODEC_FS(samplerate);
    decoder->fs_out = samplerate;
    decoder->fs_idx = FS2FS_IDX(decoder->fs);
    decoder->channels  = channels;
    decoder->frame_dms = LC3PLUS_FRAME_DURATION_10MS;
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
        
#  ifdef CR9_C_ADD_1p25MS
#    ifdef CR9_C_ADD_1p25MS_LRSNS
        setup->ltpfinfo_frame_cntr_fx = -32768; move16();
#    endif 
#  endif

        #  ifdef CR9_C_ADD_1p25MS
        #    if 0
                /*  @FhG  @alex  ToAddhere channel dependent ::  setup->ltpf_mem_continuation = 0;  */
        #    endif 
        #    ifdef CR9_C_ADD_1p25MS_LRSNS
                setup->ltpfinfo_frame_cntr_fx = -32768; move16();
        #    endif 
        #  endif

        setup->pc_seed = 24607;
        setup->ns_seed      = 24607;
        setup->ns_cum_alpha = 32767;

#  ifdef CR9_C_ADD_1p25MS
#    ifdef NEW_SIGNALLING_SCHEME_1p25
        setup->ltpfinfo_frame_cntr_fx = -32768; move16();
#    endif 
#  endif

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
#ifdef CR14_B_REMOVE_FLOAT_IN_BASOP_CODE
        SWITCH (decoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
            case LC3PLUS_FRAME_DURATION_1p25MS:
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    LC3PLUS_FRAME_DURATION_1p25MS * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
            BREAK;
#endif
            case LC3PLUS_FRAME_DURATION_2p5MS:
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    LC3PLUS_FRAME_DURATION_2p5MS * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
            BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    LC3PLUS_FRAME_DURATION_5MS * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
            BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    LC3PLUS_FRAME_DURATION_7p5MS * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
            BREAK;
            case LC3PLUS_FRAME_DURATION_10MS:
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    LC3PLUS_FRAME_DURATION_10MS * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
            BREAK;
            case LC3PLUS_FRAME_DURATION_UNDEFINED:
                setup->plcAd->PhECU_frame_ms = 0;
            assert(0);
        }
#else
        setup->plcAd->PhECU_frame_ms = (Word16)(
            decoder->frame_dms * 1.25 * 10 * 0.1); /* needed in PLCUpdate and PLC main functions, adjusted in first frame */
#endif
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
#ifdef CR14_B_REMOVE_FLOAT_IN_BASOP_CODE
    Word16 n_limit;
#endif

#ifndef FIX_TX_RX_STRUCT_STEREO
#ifdef CR9_C_ADD_1p25MS
    decoder->ltpf_rx_status[0] = 0;
    decoder->ltpf_rx_status[1] = 0;
    decoder->ltpf_mem_continuation = 0;
    decoder->ltpf_mem_active_prev = 0;
#endif
#endif
    
    decoder->frame_length = extract_l(L_shr_pos(Mpy_32_16(decoder->fs, 20972), 6)); /* fs * 0.01*2^6 */
    
#ifdef ENABLE_HR_MODE
    if (decoder->hrmode)
    {
        decoder->BW_cutoff_bits = 0;
    }
#endif

    SWITCH (decoder->frame_dms)
    {
#ifdef CR9_C_ADD_1p25MS
    case LC3PLUS_FRAME_DURATION_1p25MS:
        decoder->frame_length         = shr_pos(decoder->frame_length, 3);
        decoder->la_zeroes            = 0;
        decoder->stDec_ola_mem_fx_len = sub(decoder->frame_length, decoder->la_zeroes);
        decoder->BW_cutoff_bits       = 0;
#    ifdef ENABLE_HR_MODE
        if (decoder->hrmode)
        {
#ifdef CR14_A_ADD_1p25MS_HR
            decoder->bands_number = bands_number_1_25ms_HR[decoder->fs_idx];
            decoder->bands_offset = bands_offset_1_25ms_HR[decoder->fs_idx];
            decoder->W_fx = LowDelayShapes_n960_1_25ms_HR[decoder->fs_idx];
            decoder->W_size = LowDelayShapes_n960_len_1_25ms[decoder->fs_idx];
            decoder->yLen = decoder->frame_length;
#ifdef CR15_A_LOSSLESS_1p25MS
            decoder->low_band_limit = 40;
#endif
#else
            assert(0);
#endif
        }
        else
#    endif
        {
            decoder->bands_number = bands_number_1_25ms[decoder->fs_idx];
            decoder->bands_offset = bands_offset_1_25ms[decoder->fs_idx];
            decoder->W_fx         = LowDelayShapes_n960_1_25ms[decoder->fs_idx];
            decoder->W_size       = LowDelayShapes_n960_len_1_25ms[decoder->fs_idx];
            decoder->yLen         = s_min(MAX_BW >> 3, decoder->frame_length);            
        }
#ifdef FIX_PLC_CONFORM_ISSUES
        if (decoder->fs_idx == 0 || decoder->frame_length <= 20)
#else
        if (decoder->fs_idx == 0)
#endif
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
#endif
    case LC3PLUS_FRAME_DURATION_2p5MS:
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
            #ifdef CR14_A_ADD_LOSSLESS_MODE
            decoder->low_band_limit = 40;  
            #endif          
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
        
#ifdef FIX_PLC_CONFORM_ISSUES
        if (decoder->fs_idx == 0 || decoder->frame_length <= 20)
#else
        if (decoder->fs_idx == 0)
#endif
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
    case LC3PLUS_FRAME_DURATION_5MS:
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
            #ifdef CR14_A_ADD_LOSSLESS_MODE
            decoder->low_band_limit = 47;
            #endif 
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
        
    case LC3PLUS_FRAME_DURATION_7p5MS:
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
            #ifdef CR14_A_ADD_LOSSLESS_MODE
            decoder->low_band_limit = 52;
            #endif

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
        
    case LC3PLUS_FRAME_DURATION_10MS:
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
            #ifdef CR14_A_ADD_LOSSLESS_MODE
            decoder->low_band_limit = 52;
            #endif

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
    case LC3PLUS_FRAME_DURATION_UNDEFINED:
        assert(0);
    }

    {
        int ch;
        for (ch = 0; ch < decoder->channels; ch++)
        {
            DecSetup *setup = decoder->channel_setup[ch];
            if (setup->plcAd != NULL)
            { /*only set  if plcAd was actually allocated */
#ifdef CR14_B_REMOVE_FLOAT_IN_BASOP_CODE
                int idx;

                SWITCH (decoder->frame_dms)
                {
#ifdef CR9_C_ADD_1p25MS
                    case LC3PLUS_FRAME_DURATION_1p25MS:
                        idx = (LC3PLUS_FRAME_DURATION_1p25MS * 1.25 * 10 / 25) - 1;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = (Word16)(LC3PLUS_FRAME_DURATION_1p25MS * 1.25 * 10 * 0.1);
                        n_limit = PLC_FADEOUT_TYPE_1_IN_MS*10/(LC3PLUS_FRAME_DURATION_1p25MS*1.25*10);
                    BREAK;
#endif
                    case LC3PLUS_FRAME_DURATION_2p5MS:
                        idx = (LC3PLUS_FRAME_DURATION_2p5MS * 1.25 * 10 / 25) - 1;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = (Word16)(LC3PLUS_FRAME_DURATION_2p5MS * 1.25 * 10 * 0.1);
                        n_limit = PLC_FADEOUT_TYPE_1_IN_MS*10/(LC3PLUS_FRAME_DURATION_2p5MS*1.25*10);
                    BREAK;
                    case LC3PLUS_FRAME_DURATION_5MS:
                        idx = (LC3PLUS_FRAME_DURATION_5MS * 1.25 * 10 / 25) - 1;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = (Word16)(LC3PLUS_FRAME_DURATION_5MS * 1.25 * 10 * 0.1);
                        n_limit = PLC_FADEOUT_TYPE_1_IN_MS*10/(LC3PLUS_FRAME_DURATION_5MS*1.25*10);
                    BREAK;
                    case LC3PLUS_FRAME_DURATION_7p5MS:
                        idx = (LC3PLUS_FRAME_DURATION_7p5MS * 1.25 * 10 / 25) - 1;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = (Word16)(LC3PLUS_FRAME_DURATION_7p5MS * 1.25 * 10 * 0.1);
                        n_limit = PLC_FADEOUT_TYPE_1_IN_MS*10/(LC3PLUS_FRAME_DURATION_7p5MS*1.25*10);
                    BREAK;
                    case LC3PLUS_FRAME_DURATION_10MS:
                        idx = (LC3PLUS_FRAME_DURATION_10MS * 1.25 * 10 / 25) - 1;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = (Word16)(LC3PLUS_FRAME_DURATION_10MS * 1.25 * 10 * 0.1);
                        n_limit = PLC_FADEOUT_TYPE_1_IN_MS*10/(LC3PLUS_FRAME_DURATION_10MS*1.25*10);
                    BREAK;
                    case LC3PLUS_FRAME_DURATION_UNDEFINED:
                        idx = 0;
                        /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
                        setup->plcAd->PhECU_frame_ms = 0;
                        n_limit = 0;
                        assert(0);
                }
#else // #ifdef CR14_B_REMOVE_FLOAT_IN_BASOP_CODE
                int idx = (decoder->frame_dms * 1.25 * 10 / 25) - 1;
                setup->plcAd->PhECU_frame_ms = (Word16)(
                    decoder->frame_dms * 1.25 * 10 *
                    0.1); /* needed in processPLCupdate_fx(),  now set properly set in first frame /second time  */
#endif
                setup->plcAd->longterm_analysis_counter_max = plc_fadeout_param_maxlen[idx];
                setup->plcAd->longterm_analysis_counter_max_bytebuffer = plc_fadeout_param_maxbytes[idx];
            }
        }
    }
#ifdef CR14_B_REMOVE_FLOAT_IN_BASOP_CODE
    FOR (n=0; n < n_limit; n++){
#else
    FOR (n=0; n < PLC_FADEOUT_TYPE_1_IN_MS*10/(decoder->frame_dms*1.25*10);n++){
#endif
        decoder->alpha_type_2_table[n] = type_2_fadeout_fx(n, decoder->frame_dms);
    }
}

LC3PLUS_Error update_dec_bitrate(LC3PLUS_Dec *decoder, int ch, Word16 nBytes)
{
    int tmp = 0, totalBits = 0;
    int channel_bytes = 0;
    int minBytes = 0, maxBytes = 0;

    DecSetup *setup = decoder->channel_setup[ch];
    channel_bytes   = nBytes;
    
#ifdef ENABLE_HR_MODE
    if (decoder->hrmode)
    {
        SWITCH (decoder->frame_dms)
        {
#ifdef CR9_C_ADD_1p25MS
            case LC3PLUS_FRAME_DURATION_1p25MS:
                maxBytes = 210;
                minBytes = MIN_NBYTES;
                BREAK;
#endif
            case LC3PLUS_FRAME_DURATION_2p5MS:
                maxBytes = 210;
                minBytes = MIN_NBYTES;
                BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:
                maxBytes = 375;
                minBytes = MIN_NBYTES;
                BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:
                maxBytes = 625;
                minBytes = MIN_NBYTES;
                BREAK;
            case LC3PLUS_FRAME_DURATION_10MS:
                maxBytes = 625;
                minBytes = MIN_NBYTES;
                BREAK;
            case LC3PLUS_FRAME_DURATION_UNDEFINED: return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
#endif
    {
        minBytes = MIN_NBYTES;
        maxBytes = MAX_NBYTES_100;
    }
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (decoder->lossless)
    {
        maxBytes = 7000;
#ifdef CR15_A_LOSSLESS_1p25MS
        minBytes = 6;
#else
        minBytes = 8;
#endif
    }
#endif
    
    if (channel_bytes < minBytes || channel_bytes > maxBytes)
    {
        return LC3PLUS_NUMBYTES_ERROR;
    }

    setup->targetBytes = channel_bytes;
    move16();
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    setup->total_bits           = L_shl(setup->targetBytes, 3);
#else
    setup->total_bits           = shl(setup->targetBytes, 3);
#endif
    setup->enable_lpc_weighting = (setup->total_bits < 480);


#ifdef FIX_BOTH_1p25_WB_GLOBGAINOFFSET_NONBE  
    IF(decoder->frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        setup->quantizedGainOff = calc_GGainOffset_1p25_fx(setup->total_bits, decoder->fs_idx);  /* enc/dec common function */
    }
    ELSE
    {
        setup->quantizedGainOff =
            -(s_min(115, setup->total_bits / (10 * (decoder->fs_idx + 1))) + 105 + 5 * (decoder->fs_idx + 1));
    }
#else 
    setup->quantizedGainOff =
        -(s_min(115, setup->total_bits / (10 * (decoder->fs_idx + 1))) + 105 + 5 * (decoder->fs_idx + 1));
#endif
    
#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
    if ( (decoder->hrmode && decoder->fs_idx == 5) || decoder->lossless )
#else
    if (decoder->hrmode && decoder->fs_idx == 5)
#endif
    {
        setup->quantizedGainOff = MAX(setup->quantizedGainOff, -181);
    }
#endif
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (decoder->fs_idx == 4)
    {
        IF (decoder->wavFormat == 16)
        {
            setup->ll_ari_bits = 1800;
        } ELSE {
            setup->ll_ari_bits = 1750;
            setup->ll_offQuant = 0;
        }
    } ELSE {
        IF (decoder->wavFormat == 16)
        {
            setup->ll_ari_bits = 4100;
        } ELSE {
            setup->ll_ari_bits = 3500;
            setup->ll_offQuant = 1;
        }
    }
    
    setup->ll_ari_bits = (Word16)( (setup->ll_ari_bits*decoder->frame_dms) >> 3 );

    {
        Word32 fact[7] = {214748365, 107374183, 71582789, 53687092, 42949673, 35791395, 30678338}; // ceil(2^31 / [10 20 30 40 50 60 70])
        setup->quantizedGainOff_ll = -( s_min( 115, Mpy_32_32(setup->ll_ari_bits, fact[decoder->fs_idx])) + 105 + 5 * ( decoder->fs_idx + 1 ) );
    }
    IF (decoder->wavFormat == 24)
    {
        setup->quantizedGainOff_ll = setup->quantizedGainOff_ll + 48;
    }

#ifdef LL_INCL_HPVC
    setup->hpvcDecCfg.active_flag = 0;
    if ((decoder->fs >= LL_HPVC_LOWEST_FS) &&
        (decoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS || decoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS))
    {
        setup->hpvcDecCfg.active_flag = 1;
        setup->hpvcDecCfg.startCoefListNom[0] = LL_HPVC_BW_STARTA_FB10MS;       /* lim0 */
        setup->hpvcDecCfg.startCoefListNom[1] = LL_HPVC_BW_STARTB_FB10MS;       /* lim1 */

#ifdef HPVC_CORRECT7p5_START
        if (decoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS)
        {
            Word16 tot_blocks_fb = (480 - LL_HPVC_BW_STARTA_FB10MS) >> N_SIGNAL_LOG;
            Word16 top_blocks_fb = (480 - LL_HPVC_BW_STARTB_FB10MS) >> N_SIGNAL_LOG;
            Word16 tot_blocks = (3 * tot_blocks_fb * (decoder->fs / 48000)) >> 2;  /* reduce by 75% for 7.5 */
            Word16 top_blocks = (3 * top_blocks_fb * (decoder->fs / 48000)) >> 2;

#ifdef HPVC_CORRECT7p5_START_PLUS
            if (decoder->fs >= 96000)
            {
                tot_blocks += 0; /* allow slighly more  than 75% of 10ms at  7p5ms */
                tot_blocks = MIN(tot_blocks, tot_blocks_fb * (decoder->fs / 48000));  /* never more blocks than 10ms */
            }
#endif
            ASSERT(tot_blocks > top_blocks);
            ASSERT(top_blocks >= 1);

            /* adjust so the HPVC Nsignal Phase is maintained in Nominal 7.5ms Coeffs by adding nb_ blocks from the top */
            setup->hpvcDecCfg.startCoefListNom[0] = decoder->frame_length - (tot_blocks * LL_HPVC_N_SIGNAL);
            setup->hpvcDecCfg.startCoefListNom[1] = decoder->frame_length - (top_blocks * LL_HPVC_N_SIGNAL);

            ASSERT(setup->hpvcDecCfg.startCoefListNom[1] >= (LL_HPVC_N_SIGNAL + setup->hpvcDecCfg.startCoefListNom[0]));
        }
        else
        { /* 10ms */
            setup->hpvcDecCfg.startCoefListNom[0] *= (decoder->fs / 48000);
            setup->hpvcDecCfg.startCoefListNom[1] *= (decoder->fs / 48000);
        }
        setup->hpvcDecCfg.startCoef = setup->hpvcDecCfg.startCoefListNom[0];  /* lim0 */
#else
        if (decoder->frame_dms == LC3PLUS_FRAME_DURATION_7p5MS)
        {
            Word16 tot_blocks_fb = (480 - LL_HPVC_BW_STARTA_FB10MS) >> N_SIGNAL_LOG;
            Word16 top_blocks_fb = (480 - LL_HPVC_BW_STARTB_FB10MS) >> N_SIGNAL_LOG;

            setup->hpvcDecCfg.startCoefListNom[0] = 320 - (tot_blocks_fb * LL_HPVC_N_SIGNAL);
            setup->hpvcDecCfg.startCoefListNom[1] = 320 - (top_blocks_fb * LL_HPVC_N_SIGNAL);
        }
        if (decoder->fs == 96000)
        {
            setup->hpvcDecCfg.startCoefListNom[0] *= 2;
            setup->hpvcDecCfg.startCoefListNom[1] *= 2;
        }
        if (decoder->fs == 192000)
        {
            setup->hpvcDecCfg.startCoefListNom[0] *= 4;
            setup->hpvcDecCfg.startCoefListNom[1] *= 4;
        }
        setup->hpvcDecCfg.startCoef = setup->hpvcDecCfg.startCoefListNom[0];   /* received over bitstream */
#endif

        setup->hpvcDecCfg.nomTreeLim = HPVC_NOMTREE_COUNT_FB;   /* soft limit enc side */
#ifdef HPVC_MAXTREE_LIMIT
        setup->hpvcDecCfg.maxTreeLim = HPVC_MAXTREE_LIMIT_FB;   /* hard limit */
        if (decoder->fs == 96000)
        {
            setup->hpvcDecCfg.maxTreeLim = (HPVC_MAXTREE_LIMIT_FB - 1);
        }
        if (decoder->fs == 192000)
        {
            setup->hpvcDecCfg.maxTreeLim = (HPVC_MAXTREE_LIMIT_FB - 4); /* hard limit */ /*3 trees*/
        }
        if (decoder->frame_dms == LC3PLUS_FRAME_DURATION_10MS)
        {
            setup->hpvcDecCfg.nomTreeLim += 1;
            setup->hpvcDecCfg.maxTreeLim += 1;
        }
#endif
    }
    setup->hpvcDecCfg.mode = -1;   /* default to TCX-only path */
#endif

#endif

    tmp       = i_mult(80, decoder->fs_idx);

    SWITCH (decoder->frame_dms)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE  
        UWord16 low; // dummy only - not used
#endif 
#ifdef CR9_C_ADD_1p25MS
        case LC3PLUS_FRAME_DURATION_1p25MS:
            setup->enable_lpc_weighting = 0;
            /* total_bits * 3.36 */
#ifdef CR14_A_ADD_LOSSLESS_MODE            
            Mpy_32_16_ss(L_shl(setup->total_bits,2), 27525, &totalBits, &low);
#else
            totalBits = extract_l(L_shr(L_mult0(27526, setup->total_bits), 13));
#endif
            BREAK;
#endif
        case LC3PLUS_FRAME_DURATION_2p5MS:
            setup->enable_lpc_weighting = 0;
            /* total_bits * 2.4 */
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Mpy_32_16_ss(L_shl(setup->total_bits,2), 19661, &totalBits, &low);
#else
            totalBits = extract_l(L_shr(L_mult0(19661, setup->total_bits), 13));
#endif
            BREAK;
        case LC3PLUS_FRAME_DURATION_5MS:
            setup->enable_lpc_weighting = setup->total_bits < 240;
            /* total_bits * 2 */
#ifdef CR14_A_ADD_LOSSLESS_MODE
            totalBits = L_sub( L_shl( setup->total_bits, 1 ), 160 );
#else
            totalBits = sub(i_mult(setup->total_bits, 2), 160);
#endif
            BREAK;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            setup->enable_lpc_weighting = setup->total_bits < 360;
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Mpy_32_16_ss(L_shl(setup->total_bits,2), 10923, &totalBits, &low);
#else
            totalBits = L_shr(L_mult0(10923, setup->total_bits), 13);
#endif
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS:
            totalBits = setup->total_bits;
            BREAK;
        case LC3PLUS_FRAME_DURATION_UNDEFINED:
            totalBits = 0;
            assert(0);
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
#ifdef CR14_A_ADD_1p25MS_HR
#ifdef CR15_A_LOSSLESS_1p25MS
        if (decoder->frame_dms != LC3PLUS_FRAME_DURATION_1p25MS || decoder->lossless)
#else
        if (decoder->frame_dms != LC3PLUS_FRAME_DURATION_1p25MS)
#endif
        {
            setup->ltpf_scale_fac_idx = -1;
        }
#else
        setup->ltpf_scale_fac_idx = -1;
#endif
        move16();
    }
#endif

    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_get_decoder_min_max_bytes( LC3PLUS_Dec* decoder, Word32* min_bytes, Word32* max_bytes )
{
    Word32 min_bytes_loc = 0, max_bytes_loc = 0;
    UNUSED( decoder );
    
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (decoder->lossless)
    {
        *max_bytes = 7000;
        *min_bytes = 8;
        
        return LC3PLUS_OK;
    }
#endif

#  ifdef ENABLE_HR_MODE
    if ( decoder->hrmode )
    {
        SWITCH( decoder->frame_dms )
        {
        case LC3PLUS_FRAME_DURATION_2p5MS:
            max_bytes_loc = 210;
            min_bytes_loc = MIN_NBYTES;
            BREAK;
        case LC3PLUS_FRAME_DURATION_5MS:
            max_bytes_loc = 375;
            min_bytes_loc = MIN_NBYTES;
            BREAK;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            max_bytes_loc = 625;
            min_bytes_loc = MIN_NBYTES;
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS:
            max_bytes_loc = 625;
            min_bytes_loc = MIN_NBYTES;
            BREAK;
        default:
            return LC3PLUS_HRMODE_ERROR;
        }
    }
    else
#  endif
    {
        min_bytes_loc = MIN_NBYTES;
        max_bytes_loc = MAX_NBYTES_100;  // for backward compatibility, MAX_NBYTES_100 is used for all frame lengths
    }

    *min_bytes = min_bytes_loc;
    *max_bytes = max_bytes_loc;

    return LC3PLUS_OK;
}
