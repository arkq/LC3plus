/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"

void plc_phEcu_hq_ecu(
    LC3_FLOAT *f0binPtr, LC3_FLOAT *f0ltpGainPtr, LC3_FLOAT *xfp,
    LC3_INT16 prev_bfi, LC3_INT32 *short_flag_prev, LC3_INT32 fs,
    LC3_INT32 *time_offs, Complex *X_sav_m, LC3_INT32 *n_plocs, LC3_INT32 *plocs, LC3_FLOAT *f0est,
    const LC3_FLOAT *mdctWin, LC3_FLOAT *env_stabPtr, LC3_INT32 delta_corr, LC3_FLOAT *pfind_sensPtr,
    LC3_INT32 PhECU_LA, LC3_INT32 t_adv, const LC3_FLOAT *winWhr,
    LC3_FLOAT *oold_grp_shape, LC3_FLOAT *oold_EwPtr, LC3_FLOAT *old_grp_shape, LC3_FLOAT *old_EwPtr,
    LC3_FLOAT *st_beta_mute, LC3_FLOAT *st_mag_chg_1st, LC3_FLOAT *st_Xavg, LC3_INT32 LA_ZEROS, LC3_FLOAT *x_tda,
    LC3_FLOAT *xsubst_dbg, Complex *X_out_m_dbg,
    LC3_INT32 *seed_dbg, LC3_FLOAT *mag_chg_dbg, LC3_INT32 *tr_dec_dbg, LC3_FLOAT *gpc_dbg, LC3_FLOAT *X_i_new_re_dbg, LC3_FLOAT *X_i_new_im_dbg,
    LC3_FLOAT *corr_phase_dbg,
    Fft *PhEcu_Fft, Fft *PhEcu_Ifft    
    , LC3_UINT8 plc_fadeout_type, LC3_INT16 *nonpure_tone_flag_ptr   /* nonpure tone flag */
    
    )
{
   LC3_INT32 i;
   LC3_INT32 fs_idx, L, Lprot, n_grp, Lecu, LXsav, Lxfp_inuse;
   LC3_FLOAT alpha[8];
   LC3_FLOAT beta[8];
   LC3_FLOAT mag_chg[8];
   LC3_FLOAT xfp_local_rnd[2*MAX_LEN];
   Complex   X_out_m[2*MAX_LEN];
   LC3_INT32 seed;
   LC3_INT32 burst_len;


    fs_idx = (LC3_INT32)floor(fs / 10000.0);
    L = (LC3_INT32)floor(0.01 * fs);
    Lprot = (LC3_INT32)(1.6 * L);
    n_grp = xavg_N_grp[fs_idx];
    Lecu = 2 * L;
    LXsav = Lprot / 2 + 1;   /* 48 kHz may be optimized , to save only up to 20 kHz as in BASOP */
    Lxfp_inuse = Lprot ;
    if (prev_bfi == 1){
       Lxfp_inuse = (LC3_INT32)(L*(3.75/10.0));
    }

      

    UNUSED(env_stabPtr);
    UNUSED(xsubst_dbg);
    UNUSED(X_out_m_dbg);
    UNUSED(seed_dbg);
    UNUSED(mag_chg_dbg);
    UNUSED(tr_dec_dbg);
    UNUSED(gpc_dbg);
    UNUSED(X_i_new_re_dbg);
    UNUSED(X_i_new_im_dbg);
    UNUSED(corr_phase_dbg);
    

    if (prev_bfi != 1)
    {
       for (i = (Lprot-Lxfp_inuse); i < Lprot; i++) {
          xfp_local_rnd[i] = xfp[i];
          /* hysteresis of low level input  aligns float fft analysis and peak  picking to BASOP performance for low level noisy signals  */
          if (xfp[i] >= -0.5 && xfp[i] <= 0.5) {
             xfp_local_rnd[i] = 0.0;
          }
       }
	    *nonpure_tone_flag_ptr = -1;  /* set  nonpure tone flag for new analysis */
       
        *time_offs = 0;
        burst_len = (*time_offs / L + 1);
        plc_phEcu_trans_burst_ana_sub(fs_idx, burst_len, n_grp, oold_grp_shape, oold_EwPtr , old_grp_shape,  old_EwPtr, st_beta_mute,
                                      st_mag_chg_1st, st_Xavg, alpha, beta, mag_chg, NULL, NULL
                                      , plc_fadeout_type                                               
              );
 
        plc_phEcu_spec_ana(xfp_local_rnd, Lprot, winWhr, pfind_sensPtr, plocs, n_plocs, f0est, X_sav_m, &LXsav, f0binPtr, f0ltpGainPtr, fs_idx, PhEcu_Fft);
    }
    else
    {
        *time_offs = *time_offs + L;
        *time_offs = imin(32767 ,*time_offs); /* limit  to Word16 range as in BASOP ~= 70 10ms frames@48kHz */ 
        burst_len = ((*time_offs / L) + 1);

        plc_phEcu_trans_burst_ana_sub(fs_idx, burst_len, n_grp, oold_grp_shape, oold_EwPtr, old_grp_shape, old_EwPtr, st_beta_mute,
                                      st_mag_chg_1st, st_Xavg, alpha, beta, mag_chg, NULL, NULL                                    
                                      , plc_fadeout_type                                          
                                      );

    }

    seed = *time_offs;
   
    if (*short_flag_prev != 0)
    {
        *n_plocs = 0;
    }

    move_cmplx( X_out_m, X_sav_m, LXsav);

    /* inplace X_out_m update */
    plc_phEcu_subst_spec(plocs, *n_plocs, f0est, *time_offs, X_out_m, LXsav, mag_chg, &seed, alpha, beta, st_Xavg, t_adv, Lprot, delta_corr,
                         plc_fadeout_type,
                         nonpure_tone_flag_ptr,  /*  nonpure_tone_flag   , a state  updated here  */
    
    
                         NULL, NULL, NULL);
  



        plc_phEcu_rec_frame(X_out_m, L, Lecu, winWhr, mdctWin, Lprot,
           xfp, /* last 3.75ms of non-rounded xfp used here  */
           *time_offs,
           x_tda /* output */, 
           NULL, NULL, NULL,
           LA_ZEROS, PhECU_LA, PhEcu_Ifft);
 
}

