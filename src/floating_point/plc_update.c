/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void processPlcUpdate_fl(LC3_INT plcMeth
, LC3_INT frame_length, LC3_FLOAT *syntM, LC3_FLOAT *scf_q, LC3_INT *nbLostCmpt, LC3_FLOAT *cum_alpha, LC3_INT bfi, LC3_INT *prevBfi, LC3_INT *prevprevBfi)
{
    UNUSED(plcMeth);
    UNUSED(frame_length);
    UNUSED(syntM);
    UNUSED(scf_q);
    
    
    if (bfi != 1)
    {   
        *nbLostCmpt = 0;
        *cum_alpha = 1;
        
    }
    
    *prevprevBfi = *prevBfi;
    *prevBfi = bfi;
}

void processPlcUpdateSpec_fl(LC3_FLOAT *q_d_prev, LC3_FLOAT *q_d_fl_c, LC3_INT yLen)
{
    move_float(q_d_prev, q_d_fl_c, yLen);
}
