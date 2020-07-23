/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void processResidualDecoding_fl(LC3_INT* bitsRead, LC3_FLOAT x[], LC3_INT L_spec, uint8_t prm[], LC3_INT resQBits
                , LC3_INT hrmode
)
{
    LC3_INT k = 0, n = 0;
    LC3_FLOAT offset1 = 0, offset2 = 0;
    LC3_FLOAT offset = 0;
    LC3_INT nz_idx[MAX_LEN] = {0};
    LC3_INT N_nz = 0, idx = 0;

    LC3_INT iter = 0, iter_max = 1;

    if (hrmode)
    {
        iter_max = EXT_RES_ITER_MAX;
        offset = offset1 = offset2 = 0.25;
    }
    else
    {
        offset1 = 0.1875;
        offset2 = 0.3125;
    }

    if (hrmode)
    {
        /* enumerat non-zero coefficients */
        for (k = 0; k < L_spec; k ++)
        {
            if (x[k])
            {
                nz_idx[N_nz ++] = k;
            }
        }
        /* apply residual corrections */
        while (n < resQBits && iter < iter_max)
        {
            for (k = 0; k < N_nz; k ++)
            {
                idx = nz_idx[k];
                
                if ((prm[n >> 3] & 1 << (n & 7)) == 0)
                {
                    x[idx] -= offset;
                }
                else
                {
                    
                    x[idx] += offset;
                }
                if (++n >= resQBits)
                {
                    break;
                }
            }
            offset /= 2;
            iter ++;
        }
    }
    else
    {
        while (k < L_spec && n < resQBits) {
            if (x[k] != 0) {
                if ((prm[n >> 3] & 1 << (n & 7)) == 0)
                {
                    if (x[k] > 0) {
                        x[k] -= offset1;
                    } else {
                        x[k] -= offset2;
                    }
                } else {
                    if (x[k] > 0) {
                        x[k] += offset2;
                    } else {
                        x[k] += offset1;
                    }
                }
                n++;
            }
            
            k++;
        }
    }
    *bitsRead = n;
}
