/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void processPlcComputeStabFac_fl(LC3_FLOAT *scf_q, LC3_FLOAT *old_scf_q, LC3_INT32 prev_bfi, LC3_FLOAT *stab_fac);

void processPlcComputeStabFacMain_fl(LC3_FLOAT *scf_q, LC3_FLOAT *old_scf_q, LC3_FLOAT *old_old_scf_q, LC3_INT32 bfi, LC3_INT32 prev_bfi,
                                      LC3_INT32 prev_prev_bfi, LC3_FLOAT *stab_fac)
{
    if (bfi == 1)
    {
        if (prev_bfi != 1)
        {
            processPlcComputeStabFac_fl(old_scf_q, old_old_scf_q, prev_prev_bfi, stab_fac);
        }
    }
    else if (bfi == 2)
    {
        processPlcComputeStabFac_fl(scf_q, old_scf_q, prev_bfi, stab_fac);
    }
}

static void processPlcComputeStabFac_fl(LC3_FLOAT *scf_q, LC3_FLOAT *old_scf_q, LC3_INT32 prev_bfi, LC3_FLOAT *stab_fac)
{
        LC3_FLOAT tmp;
        LC3_INT32 i;
    
    tmp = 0;
    
    if (prev_bfi == 1)
    {
        *stab_fac = 0.8;
    }
    else
    {
        for (i = 0; i < M; i++)
        {
            tmp += (scf_q[i] - old_scf_q[i]) * (scf_q[i] - old_scf_q[i]);
        }

        *stab_fac = 1.25 - tmp / 25.0;

        if (*stab_fac > 1)
        {
            *stab_fac = 1;
        }

        if (*stab_fac < 0)
        {
            *stab_fac = 0;
        }
    }
}

