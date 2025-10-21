/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

Word16 getLastNzBits_fx (Word16 N)
{
    Dyn_Mem_Deluxe_In( Word16 minBits;);
    minBits = sub(14, norm_s(negate(N)));
#ifdef CR9_C_ADD_1p25MS
    /* minimum of 2 spare bits */
    if (sub(sub(shl(1,minBits),shr(N,1)),2) < 0) 
    {
        minBits = add(minBits,1);
    }
#endif
    Dyn_Mem_Deluxe_Out();
    return minBits;
}
