/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                          
#include "functions.h"

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

