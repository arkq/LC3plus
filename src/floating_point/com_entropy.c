/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

LC3_INT16 getLastNzBits (LC3_INT16 N)
{
    LC3_INT16 minBits;

    minBits = ceil(LC3_LOGTWO(N >> 1));
#ifdef CR9_C_ADD_1p25MS
    /* minimum of 2 spare bits */
    if (((1 << minBits) - (N >> 1)) < 2) 
    {
        minBits++;
    }
#endif

    return minBits;
}

