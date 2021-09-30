/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void processPlcMain_fl(LC3_FLOAT *q_d_fl_c, LC3_FLOAT *syntM_fl_c, LC3PLUS_Dec* decoder, DecSetup* h_DecSetup, LC3_INT bfi,
               PlcAdvSetup *PlcAdvSetup, PlcSetup *PlcSetup, LC3_INT plcMeth, LC3_INT ltpf_pitch_int, LC3_INT ltpf_pitch_fr,
               LC3_INT tilt, const LC3_INT *bands_offset, LC3_INT bands_number, const LC3_INT *bands_offsetPLC,
               LC3_INT n_bandsPLC, LC3_INT16 hrmode, pcState *statePC
)
{
    LC3_FLOAT r[MAX_BANDS_NUMBER_PLC], A[M + 1], synth[MAX_LEN + MDCT_MEM_LEN_MAX], energies[MAX_BANDS_NUMBER_PLC];
    LC3_INT32 pitch_classifier;
    LC3_FLOAT xcorr;
    LC3_INT32 yLen;
    LC3_INT16 prev_bfi_plc2;
    LC3_FLOAT  phEcu_env_stab_local[1];
    LC3_FLOAT  phEcu_pfind_sens[1];

    prev_bfi_plc2 = 1;
    if (PlcSetup->nbLostCmpt == 0)
    {
        prev_bfi_plc2 = 0;
    }
    assert((h_DecSetup->PlcSetup.prevBfi == 1) == (prev_bfi_plc2 == 1));

    if (bfi == 1 && PlcAdvSetup)
    {
        /* FFLC increases the PFLC counter */
        statePC->ns_nbLostCmpt_pc = statePC->ns_nbLostCmpt_pc + 1;
    }

    pitch_classifier = ltpf_pitch_int;
#ifdef NONBE_PLC_CLASSIFER_LAG_FIX
    if (ltpf_pitch_fr > 2)
    {
        pitch_classifier++;
    }
#endif

    processPlcClassify_fl(plcMeth, &h_DecSetup->concealMethod, &PlcSetup->nbLostCmpt, bfi, &xcorr,
                          decoder->frame_length, decoder->frame_dms, pitch_classifier, decoder->fs,
                          bands_offset, bands_number, tilt, PlcAdvSetup, hrmode
        );

    if (bfi == 1)
    {
        switch (h_DecSetup->concealMethod)
        {
        case 2:
        {
            assert(decoder->fs_idx == floor(decoder->fs / 10000));
            // phaseECU supports only 10ms framing
            assert(PlcSetup->nbLostCmpt != 0 || decoder->frame_dms == 100);
            
            if (decoder->frame_dms != 100)
            {
                // muting, if frame size changed during phaseECU concealment
                memset(q_d_fl_c, 0, sizeof(LC3_FLOAT) * decoder->frame_length);
                h_DecSetup->alpha = 0;
                break;
            }

            /*  call phaseEcu  */
            LC3_FLOAT pitch_fl_c = (LC3_FLOAT)ltpf_pitch_int + (LC3_FLOAT)ltpf_pitch_fr / 4.0; /* use non-rounded pitch indeces  */


            if (prev_bfi_plc2 == 0)
            {
                /* convert fractional pitch lag info at current fs to a normalized fractional bin-frequency   */
                PlcAdvSetup->PlcPhEcuSetup.PhECU_f0hzLtpBin = plc_phEcuSetF0Hz(decoder->fs, &pitch_fl_c);
                /* several buffers  used in Cflt , a  copy  pcmbufHist,  right  before calling PhEcu in bad frames    */
                assert(bfi == 1);
                move_float(PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp,
                           &(PlcAdvSetup->pcmbufHist[PlcAdvSetup->max_len_pcm_plc - PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot]),
                           PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot);

                /* a first bfi frame:: calc windowed 16 ms energy twice in a 26 ms buffer separated by 10 ms*/
                {
                    const LC3_FLOAT *w, *prev_xfp;
                    LC3_INT32 i, oold_start;

                    oold_start = PlcAdvSetup->max_len_pcm_plc - (decoder->frame_length + PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot);
                    assert(oold_start > 0);
                    w = PhECU_whr16ms_wins[decoder->fs_idx]; /* hammrect table */
                    prev_xfp = &(PlcAdvSetup->pcmbufHist[oold_start + 0]);

                    PlcAdvSetup->PlcPhEcuSetup.PhECU_L_oold_xfp_w_E = 0;
                    PlcAdvSetup->PlcPhEcuSetup.PhECU_L_old_xfp_w_E = 0;
                    for (i = 0; i < PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot; i++)
                    {
                        PlcAdvSetup->PlcPhEcuSetup.PhECU_L_oold_xfp_w_E += sqrf(prev_xfp[i] * w[i]);
                        PlcAdvSetup->PlcPhEcuSetup.PhECU_L_old_xfp_w_E  += sqrf(PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp[i] * w[i]);
                    }

                }

            } /* (prev_bfi_plc2 == 0)*/
            else 
            {
                /* overwrite  last 3.75 ms of  xfp with most recent pcmbufHist   tail ,  right  before calling PhEcu in bursts of bad frames   */           
                LC3_INT32  lenCopyOla =  decoder->la_zeroes;  /*copy_part + ola_part = 3.75 ms for 10 ms frame*/

                assert(bfi == 1 && prev_bfi_plc2);
                move_float(&(PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp[PlcAdvSetup->PlcPhEcuSetup.PhECU_Lprot-lenCopyOla]),
                           &(PlcAdvSetup->pcmbufHist[PlcAdvSetup->max_len_pcm_plc - lenCopyOla]), lenCopyOla);
                   
            }

            {
                LC3_FLOAT x_tda[MAX_LEN]; /* 960/2 */
                PlcAdvSetup->PlcPhEcuSetup.PhECU_norm_corr = xcorr;
                phEcu_env_stab_local[0] = (LC3_FLOAT)PHECU_ENV_STAB_LOCAL;
                phEcu_pfind_sens[0]     = (LC3_FLOAT)PHECU_PFIND_SENS;

                plc_phEcu_hq_ecu(&(PlcAdvSetup->PlcPhEcuSetup.PhECU_f0hzLtpBin),
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_norm_corr),
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_xfp,
                                 prev_bfi_plc2,
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_short_flag_prev),
                                 decoder->fs,
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_time_offs),
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_X_sav_m, /* Complex */
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_num_plocs),
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_plocs,
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_f0est,
                                 MDCT_WINS_10ms[hrmode][decoder->fs_idx],

                                 phEcu_env_stab_local,
                                 PHECU_DELTA_CORR,
                                 phEcu_pfind_sens, 
                                 PHECU_LA,
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_t_adv,
                                 PhECU_whr16ms_wins[decoder->fs_idx],
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_oold_grp_shape,
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_L_oold_xfp_w_E),
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_old_grp_shape,
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_L_old_xfp_w_E),
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhECU_beta_mute),
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_mag_chg_1st,
                                 PlcAdvSetup->PlcPhEcuSetup.PhECU_Xavg,
                                 decoder->la_zeroes,
                                 x_tda, /* time domain aliased output */
                                 NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL
                                 ,
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhEcu_Fft),
                                 &(PlcAdvSetup->PlcPhEcuSetup.PhEcu_Ifft)
                 );


                ProcessingITDA_WIN_OLA_fl(x_tda, decoder->frame_length, decoder->imdct_win, decoder->imdct_winLen, decoder->imdct_laZeros,
                                          h_DecSetup->imdct_mem, synth);
                move_float(syntM_fl_c, synth, decoder->frame_length);

                 

            }
        }
        break;
        case 3:
            if (PlcSetup->nbLostCmpt == 1)
            {
                PlcAdvSetup->PlcTdcSetup.fract = ltpf_pitch_fr;
            }

            processPerBandEnergy_fl(n_bandsPLC, bands_offsetPLC, hrmode, decoder->frame_dms, energies, PlcSetup->q_d_prev);
            processTdcPreemphasis_fl(energies, &PlcAdvSetup->PlcTdcSetup.preemphFac, n_bandsPLC);
            processTdcInverseOdft_fl(energies, n_bandsPLC, r, PlcAdvSetup->PlcTdcSetup.lpcorder);
            processTdcLpcEstimation_fl(r, decoder->fs_idx, PlcAdvSetup->PlcTdcSetup.lpcorder + 1, A, decoder->frame_dms);
            processTdcApply_fl(ltpf_pitch_int, &PlcAdvSetup->PlcTdcSetup.preemphFac, A, PlcAdvSetup->PlcTdcSetup.lpcorder, PlcAdvSetup->pcmbufHist, PlcAdvSetup->max_len_pcm_plc, decoder->frame_length,
                               decoder->frame_dms, decoder->fs, PlcSetup->nbLostCmpt, decoder->frame_length - decoder->la_zeroes, &PlcAdvSetup->stabFac, PlcAdvSetup->PlcTdcSetup.harmonicBuf,
                               PlcAdvSetup->PlcTdcSetup.synthHist, &PlcAdvSetup->PlcTdcSetup.fract, &PlcAdvSetup->PlcTdcSetup.seed, &PlcAdvSetup->PlcTdcSetup.gain_c,
                               &h_DecSetup->alpha, synth);

            processTdcTdac_fl(synth, decoder->imdct_win, decoder->frame_length, decoder->la_zeroes, h_DecSetup->imdct_mem);
            memmove(syntM_fl_c, synth, sizeof(LC3_FLOAT) * decoder->frame_length);
            break;
        case 4:
            processNoiseSubstitution_fl(q_d_fl_c, PlcSetup->q_d_prev, decoder->yLen);
            break;
        default:
            assert("Invalid PLC method!");
        }
    }

    if (bfi == 0)
    {
        processPlcUpdateSpec_fl(PlcSetup->q_d_prev, q_d_fl_c, decoder->yLen);
    }

    yLen = MIN(decoder->frame_length, MAX_PLC_LMEM);
    if (PlcAdvSetup != NULL && (decoder->frame_dms == 100) && (hrmode == 0))
        {
            /*  BASOP processPLCspec2shape_fx(prev_bfi, bfi, q_old_d_fx, yLen, plcAd->PhECU_oold_grp_shape_fx, plcAd->PhECU_old_grp_shape_fx);*/
            plc_phEcu_processPLCspec2shape(prev_bfi_plc2, bfi, q_d_fl_c, yLen,
                                           PlcAdvSetup->PlcPhEcuSetup.PhECU_oold_grp_shape, PlcAdvSetup->PlcPhEcuSetup.PhECU_old_grp_shape);
        }
}

