/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processReorderBitstream_fl(LC3_UINT8* bytes, LC3_INT32 n_pccw, LC3_INT32 n_pc, LC3_INT32 b_left, LC3_INT32 len)
{
    LC3_UINT8 bytes_local[MAX_NBYTES2];
    LC3_INT32 i, block_bytes;
    
    assert(b_left > 0);
    
    memcpy(bytes_local, bytes, len * sizeof(LC3_UINT8));
    
    if (n_pccw == 0)
    {
        return;
    }
    
    block_bytes = ceil((LC3_FLOAT) n_pc / 2.0);
    
    for (i = 0; i < block_bytes; i++)
    {
        bytes[i] = bytes_local[b_left + i];
    }
    
    for (i = 0; i < b_left; i++)
    {
        bytes[block_bytes + i] = bytes_local[i];
    }
}

