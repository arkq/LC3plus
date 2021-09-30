/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               
#include "functions.h"


/* return pointer to aligned base + base_size, *base_size += size + 8 bytes align */
void *balloc(void *base, size_t *base_size, size_t size)
{
    uintptr_t ptr = ((uintptr_t)base + *base_size + ALIGNMENT_BALLOC_RED) & (uintptr_t)~ALIGNMENT_BALLOC_RED;
    assert((uintptr_t)base % ALIGNMENT_BALLOC == 0); /* base must be 8-byte aligned */
    *base_size = (ptr - (uintptr_t)base) + size;
    //*base_size = (*base_size + size + ALIGNMENT_BALLOC_RED) & (size_t)~ALIGNMENT_BALLOC_RED;
    return (void *)ptr;
}


LC3_FLOAT array_max_abs(LC3_FLOAT *in, LC3_INT32 len)
{
        LC3_FLOAT max;
        LC3_INT32 i;
    
    max = LC3_FABS(in[0]);
    
    for (i = 0; i < len; i++)
    {
        if (LC3_FABS(in[i]) > LC3_FABS(max))
        {
            max = LC3_FABS(in[i]);
        }
    }
    
    return max;
}

