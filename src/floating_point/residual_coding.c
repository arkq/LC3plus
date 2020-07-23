/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void processResidualCoding_fl(LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gain, LC3_INT L_spec, LC3_INT targetBits, LC3_INT nBits, uint8_t* resBits, LC3_INT* numResBits
                , LC3_INT hrmode
)
{
    LC3_INT n = 0, m = 0, k = 0;
    LC3_INT iter=0;
    LC3_FLOAT offset;
    LC3_INT iter_max = 1;
    LC3_INT nz_idx[MAX_LEN] = {0};
    LC3_INT N_nz = 0, idx = 0;


    memset(resBits, 0, MAX_RESBITS_LEN);

    m = targetBits - nBits + 4;
    if (hrmode)
    {
        m += 10;
    }

    assert(m <= MAX_RESBITS);

    offset = .25;
    if (hrmode)
    {
        iter_max = EXT_RES_ITER_MAX;

    }
    for (k = 0; k < L_spec; k ++)
    {
        if (xq[k])
        {
            nz_idx[N_nz ++] = k;
        }
    }
    while (iter < iter_max && n < m)
    {
        k = 0;
        while (k < N_nz && n < m)
        {
            idx = nz_idx[k];
            
            if (x[idx] >= (LC3_FLOAT)xq[idx] * gain)
            {
                resBits[n >> 3] |= 1 << (n & 7);
                x[idx] -= gain * offset;
            }
            else
            {
                resBits[n >> 3] &= ~(1 << (n & 7));
                x[idx] += gain * offset;
            }
            
            n++;
            
            k++;
        }
        iter ++;
        offset *= .5;
    }

    *numResBits = n;
}
