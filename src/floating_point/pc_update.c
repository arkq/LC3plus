/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPcUpdate_fl(LC3_INT32 bfi, LC3_FLOAT *q_res, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_INT32 rframe, LC3_INT32 *BW_cutoff_idx_nf, LC3_INT32 *prev_BW_cutoff_idx_nf,
      LC3_INT32 fac_ns_idx, LC3_FLOAT *prev_fac_ns, LC3_FLOAT *fac, LC3_FLOAT *q_old_res, LC3_FLOAT *prev_gg, LC3_INT32 spec_inv_idx, LC3_INT32 yLen)
{
        LC3_FLOAT gg;
    
    gg = LC3_POW(10.0, ((LC3_FLOAT) (gg_idx + gg_idx_off)) / 28.0);
    *prev_gg = gg;
    move_float(q_old_res, q_res, yLen);
    
    if (rframe == 0)
    {
        *prev_BW_cutoff_idx_nf = *BW_cutoff_idx_nf;
        *prev_fac_ns = (8.0 - (LC3_FLOAT) fac_ns_idx) / 16.0;
    } else if ((bfi == 2) && (*BW_cutoff_idx_nf != *prev_BW_cutoff_idx_nf) && (spec_inv_idx < yLen))
    {
        *BW_cutoff_idx_nf = *prev_BW_cutoff_idx_nf;
        *prev_fac_ns = *prev_fac_ns * (*fac);
        *prev_fac_ns = MAX(*prev_fac_ns, (8.0 - 7.0) / 16.0);
        *prev_fac_ns = MIN(*prev_fac_ns, (8.0 - 0.0) / 16.0); 
    }
    
}

