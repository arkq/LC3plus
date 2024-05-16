/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPcMain_fl(LC3_INT32 *bfi, LC3PLUS_Dec* decoder, LC3_FLOAT *sqQdec, DecSetup* h_DecSetup, LC3_INT32 pitch_present, LC3_FLOAT stab_fac, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_INT32 fac_ns_idx, pcState *statePC, LC3_INT32 spec_inv_idx, LC3_INT32 yLen)
{
        LC3_FLOAT fac;
    
    /* PC Classifier */
    if (*bfi == 2)
    {
        processPcClassify_fl(pitch_present, decoder->frame_dms, h_DecSetup->PlcSetup.q_d_prev, statePC->q_old_res, decoder->yLen, h_DecSetup->spec_inv_idx, stab_fac, bfi);
    }

    /* PC Apply */
    if (*bfi == 2)
    {
        processPcApply_fl(sqQdec, statePC->q_old_res, h_DecSetup->PlcSetup.q_d_prev, h_DecSetup->spec_inv_idx, decoder->yLen, gg_idx, gg_idx_off, &statePC->prev_gg, &fac, &statePC->ns_nbLostCmpt_pc);
    }
    
    /* PC Update */
    if (*bfi != 1)
    {
        processPcUpdate_fl(*bfi, sqQdec, gg_idx, gg_idx_off, decoder->rframe, &h_DecSetup->BW_cutoff_idx_nf, &h_DecSetup->prev_BW_cutoff_idx_nf, fac_ns_idx, &h_DecSetup->prev_fac_ns, &fac, statePC->q_old_res, &statePC->prev_gg, spec_inv_idx, yLen);
    }
    
    /* Reset counter */
    if (*bfi == 0)
    {
        statePC->ns_nbLostCmpt_pc = 0;
    }
}

