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

void plc_phEcu_trans_burst_ana_sub(LC3_INT32 fs_idx, LC3_INT32 burst_len, LC3_INT32 n_grp, LC3_FLOAT *oold_spect_shape, 
                                   LC3_FLOAT *oold_EwPtr, LC3_FLOAT *old_spect_shape, 
                                   LC3_FLOAT *old_EwPtr,   LC3_FLOAT *stPhECU_beta_mute,
                                   LC3_FLOAT *stPhECU_mag_chg_1st, LC3_FLOAT *stPhECU_Xavg,  LC3_FLOAT *alpha, LC3_FLOAT *beta, LC3_FLOAT *mag_chg, LC3_INT32 *tr_dec_dbg, LC3_FLOAT *gpc_dbg
								   , LC3_UINT8 plc_fadeout_type
) 
{
   LC3_FLOAT gr_pow_left[MAX_LGW];
   LC3_FLOAT gr_pow_right[MAX_LGW];
   LC3_FLOAT trans[MAX_LGW];
   LC3_FLOAT grp_pow_change[MAX_LGW];
   LC3_FLOAT ph_dith[MAX_LGW];
   LC3_FLOAT att_val[MAX_LGW];
   LC3_INT32 tr_dec[MAX_LGW];

   LC3_INT32 attDegreeFrames;
   LC3_FLOAT thresh_dbg;
    
    UNUSED(tr_dec_dbg);
    UNUSED(gpc_dbg);

   if (burst_len <= 1)
    {
        plc_phEcu_tba_spect_Xavg(fs_idx, n_grp, oold_spect_shape, oold_EwPtr, old_spect_shape, old_EwPtr, gr_pow_left, gr_pow_right, stPhECU_Xavg);

        plc_phEcu_tba_per_band_gain(n_grp, gr_pow_left, gr_pow_right, trans, grp_pow_change);

    }

    plc_phEcu_tba_trans_dect_gains(burst_len, n_grp, grp_pow_change, stPhECU_beta_mute, stPhECU_mag_chg_1st, alpha, beta, mag_chg, ph_dith, tr_dec, att_val, &attDegreeFrames, &thresh_dbg
		                           , plc_fadeout_type
	);



    return;
}

