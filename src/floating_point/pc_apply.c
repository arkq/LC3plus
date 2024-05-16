/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPcApply_fl(LC3_FLOAT *q_res, LC3_FLOAT *q_old_res, LC3_FLOAT *q_d_prev, LC3_INT32 spec_inv_idx, LC3_INT32 yLen, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_FLOAT *prev_gg, LC3_FLOAT *fac, LC3_INT32 *pc_nbLostCmpt)
{
        LC3_FLOAT gg, mean_nrg_low, mean_nrg_high, ener_prev, ener_curr, fac_local;
        LC3_INT32 i;
    
    ener_prev = 0; ener_curr = 0; mean_nrg_low = 0; mean_nrg_high = 0;
    
    *pc_nbLostCmpt = *pc_nbLostCmpt + 1;
    
    assert(spec_inv_idx > 1);
    
    gg = LC3_POW(10, (((LC3_FLOAT) (gg_idx + gg_idx_off)) / 28.0));

    for (i = 0; i < spec_inv_idx; i++)
    {
        mean_nrg_low += LC3_POW(q_d_prev[i], 2);
    }
    
    mean_nrg_low /= (LC3_FLOAT) spec_inv_idx;
    
    if (spec_inv_idx < yLen)
    {
        for (i = spec_inv_idx; i < yLen; i++)
        {
            mean_nrg_high += LC3_POW(q_d_prev[i], 2);
        }
    }
    
    mean_nrg_high /= (LC3_FLOAT) (yLen - spec_inv_idx);

    for (i = 0; i < spec_inv_idx; i++)
    {
        ener_prev += LC3_POW(q_old_res[i], 2);
        ener_curr += LC3_POW(q_res[i], 2);
    }

    *fac = 1;
    if (ener_prev > 0)
    {
        *fac = LC3_SQRT(ener_curr / ener_prev);
    }
    
    fac_local = *fac;
    if (mean_nrg_low <= mean_nrg_high || ener_prev * LC3_POW(*prev_gg, 2) <= ener_curr * LC3_POW(gg, 2))
    {
        fac_local = *prev_gg / gg;
    }
    
    for (i = spec_inv_idx; i < yLen; i++)
    {
        q_res[i] = q_old_res[i] * fac_local;
        
        if (LC3_FABS(q_res[i]) < (1 - 0.375))
        {
            q_res[i] = 0;
        }
    }
}

