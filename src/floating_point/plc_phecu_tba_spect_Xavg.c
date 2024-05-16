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

void plc_phEcu_tba_spect_Xavg(LC3_INT32 fs_idx, LC3_INT32 n_grp, LC3_FLOAT *oold_spec_shape, 
   LC3_FLOAT *oold_EwPtr, LC3_FLOAT *old_spec_shape, 
   LC3_FLOAT *old_EwPtr,  LC3_FLOAT *gr_pow_left, LC3_FLOAT *gr_pow_right, LC3_FLOAT *Xavg) 
{
      LC3_INT32 i;
      LC3_FLOAT XavgEn[MAX_LGW];
      LC3_FLOAT xfp_w_scale, oold_Escale, old_Escale;

    /* 8k       16k    24k    32k      48k   */
    LC3_FLOAT flt_xfp_wE_MDCT2FFT_target[5] = { (LC3_FLOAT) 1.9906, (LC3_FLOAT) 4.0445, (LC3_FLOAT) 6.0980, (LC3_FLOAT) 8.1533, (LC3_FLOAT) 12.2603 };
    LC3_INT32   gw_0[10] = { 1, 3, 5, 9, 17, 33, 49, 65, 81, 97 }; 
   
    /* prepare scale factor */

    xfp_w_scale = LC3_ROUND(flt_xfp_wE_MDCT2FFT_target[fs_idx]/(LC3_FLOAT)16.0*(LC3_FLOAT) 32768.0) / (LC3_FLOAT) LC3_POW(2,11);

    /* prepare left and right subband energies */
    oold_Escale = (*oold_EwPtr) * xfp_w_scale;
    old_Escale  = (*old_EwPtr) * xfp_w_scale;
    for (i = 0;i < n_grp;i++) {
        gr_pow_left[i] = oold_spec_shape[i] * oold_Escale;
        gr_pow_right[i] = old_spec_shape[i] * old_Escale;

        XavgEn[i] = ((LC3_FLOAT) 0.5) * (gr_pow_left[i] + gr_pow_right[i]) / (gw_0[i + 1] - gw_0[i]);
        Xavg[i] = LC3_SQRT(XavgEn[i]);
    }

    return;
}

