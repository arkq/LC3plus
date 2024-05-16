/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNoiseSubstitution_fl(LC3_FLOAT* spec, LC3_FLOAT* spec_prev, LC3_INT32 yLen)
{
    memmove(spec, spec_prev, sizeof(LC3_FLOAT) * yLen);
    
    spec[0] *= 0.2;
    spec[1] *= 0.5;
}

