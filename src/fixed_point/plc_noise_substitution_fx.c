/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "defines.h"

#include "functions.h"


void processPLCNoiseSubstitution_fx(Word32 spec[], Word16 spec_prev[], Word16 L_spec)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    FOR (i = 0; i < L_spec; i++)
    {
        spec[i] = L_deposit_h(spec_prev[i]);
    }

    /* High pass to prevent overflows */
    spec[0] = Mpy_32_16(spec[0], 6553 /* 0.2 Q15*/);  move32();
    spec[1] = Mpy_32_16(spec[1], 16384 /* 0.5 Q15*/); move32();

    Dyn_Mem_Deluxe_Out();
}



