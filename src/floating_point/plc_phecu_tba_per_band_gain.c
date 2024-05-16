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

void plc_phEcu_tba_per_band_gain(LC3_INT32 n_grp, LC3_FLOAT *gr_pow_left, LC3_FLOAT *gr_pow_right, LC3_FLOAT *trans, LC3_FLOAT *grp_pow_change) 
{
    LC3_INT32 i;

    /* per band gain */
    for (i = 0; i < n_grp; i++) {
        if (gr_pow_left[i] > 0)
        {
            trans[i] = gr_pow_right[i] / gr_pow_left[i];
        }
        else
        {
            /* handle division by zero case */
            if  (gr_pow_right[i] > 0)
            {
                trans[i] = 10.0;   /* positive/0  transient  */
            }
            else
            {
                trans[i] = 1.0;  /* 0/0  no transient , no power change */
            }
        }
        grp_pow_change[i] = (LC3_FLOAT) 10.0  * LC3_LOGTEN(trans[i]);

    }

    return;  
}
