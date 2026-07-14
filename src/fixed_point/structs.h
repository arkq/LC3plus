/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef STRUCTS_H
#define STRUCTS_H

#include "defines.h"

#ifdef CR14_A_ADD_LOSSLESS_MODE
typedef struct {
    const Word32 *CS1;  // Coefficient for the 1st lifting step
    const Word32 *CS2;  // Coefficient for the 2nd lifting step
    const Word32 *CS3;  // Coefficient for the 3rd lifting step
    const Word32 *CS4;  // Coefficient for the 4th lifting step
} Lifting_fx;
  
typedef struct
{
    Word16 start_freq[2];
    Word16 stop_freq;
} TnsStartStopFreqs;
#endif

#endif /* STRUCTS_H */
