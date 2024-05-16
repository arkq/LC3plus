/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"

LC3_FLOAT plc_phEcuSetF0Hz(LC3_INT32 fs, LC3_FLOAT * old_pitchPtr)
{
    LC3_FLOAT result;

    result = 0;
    if (*old_pitchPtr != 0)
    {
      result = LC3_ROUND(fs/(*old_pitchPtr)/PHECU_FRES * 128.0) / 128.0;
    }

    return result;
}

