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

void plc_phEcu_F0_refine_first( LC3_INT32 *plocs,            /* i/o  0 ... Lprot/2 +1*/
                                LC3_INT32 n_plocs,
                                LC3_FLOAT *f0est,        /* i/o  f0est */
                                const LC3_INT32 Xabs_len,
                                LC3_FLOAT *f0binPtr,  /* i    */
                                LC3_FLOAT *f0gainPtr, /* i    */
                                const LC3_INT32 nSubm
                                )
{
   LC3_FLOAT sens;
   LC3_INT32 i, j, high_idx, breakflag;
   LC3_FLOAT f0est_lim[MAX_PLC_NPLOCS];
   LC3_FLOAT f0bin;
   LC3_FLOAT f0gain;

   f0bin  = *f0binPtr;
   f0gain = *f0gainPtr;

    if (n_plocs > 0 && f0gain > 0.25) {
        
        sens = 0.5;
        if (f0gain < 0.75) {
            sens = 0.25;
        }
        
        high_idx = -1;
        for (i = 0; i < n_plocs; i++) {
            if (plocs[i] <= 25) { /* 25 ~= 1550 Hz */
                high_idx = MAX(high_idx, i);
            } else {
                /* Optimization, only works if plocs vector is sorted. Which it should be. */
                break;
            }
        }
        
        if (high_idx != -1) {
            high_idx++;
            move_float(f0est_lim, f0est, high_idx);
            
            breakflag = 0;
            for (i = 0; i < nSubm; i++) {
                for (j = 0; j < high_idx; j++) {
               if (LC3_FABS(f0est_lim[j] - (i+1) * f0bin) < sens) {
                        f0est[j] = (i+1)*f0bin;
                  plocs[j] = MIN(Xabs_len-1, MAX(1,(LC3_INT32) LC3_ROUND(f0est[j])));
                        breakflag = 1;
                        break;
                    }
                }
                if (breakflag) {
                    break;
                }
                sens *= 0.875;
            }
        }
    }

   return;
}

