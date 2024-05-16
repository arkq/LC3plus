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

void plc_phEcu_LF_peak_analysis(LC3_INT32 *plocs,        /* i/o  0 ... Lprot/2 +1*/
                                LC3_INT32 *n_plocs,      /* i/o   0.. MAX_PLOCS  */
                                LC3_FLOAT *f0est,        /* i/o  Q16*/
                                const LC3_FLOAT *Xabs,
                                LC3_FLOAT *f0binPtr,
                                LC3_FLOAT *f0gainPtr,
                                const LC3_INT32 nSubm
)
{
    LC3_INT32 i, j, fin, f_ind, prel_low, prel_high, start;
    LC3_FLOAT f0est_prel[3];
    LC3_INT32 plocs_prel[3];
    LC3_INT32 n_prel;
    LC3_FLOAT f0est_old[MAX_PLC_NPLOCS];
    LC3_INT32 plocs_old[MAX_PLC_NPLOCS];
    LC3_FLOAT peakLF_Xval, f;
    LC3_FLOAT f0bin  ;
    LC3_FLOAT f0gain ;
  
    f0bin  = *f0binPtr;
    f0gain = *f0gainPtr;

    if (*n_plocs > 0 && f0gain > 0.25 && f0bin < 2.75) {

        /* Assumes sorted plocs */
        if (plocs[0] < 3) {
            fin = MIN(3, *n_plocs);
            peakLF_Xval = Xabs[plocs[0]];
            for (i = 1; i < fin; i++) {
                peakLF_Xval = MAX(peakLF_Xval, Xabs[plocs[i]]);
            }

            n_prel = 0;
            for (i = 0; i < nSubm; i++) {
                f = (i+1)*f0bin;
            f_ind = (LC3_INT32)LC3_ROUND(f);
                if (f*PHECU_FRES <= 400 && Xabs[f_ind] > peakLF_Xval*0.375) {
                    f0est_prel[n_prel] = f;
                    plocs_prel[n_prel] = f_ind;
                    n_prel++;
                }
            }

            if (n_prel > 0) {
                prel_low = plocs_prel[0];
                prel_high = plocs_prel[n_prel-1];

                /*   initial assumption:: all original peaks (1 or 2 of them)  are positioned   below  prel_low  */   
                start =  (*n_plocs);  /* at this point  'start' is the  location_c where to add any harmonics peaks */      
                for (i = (*n_plocs)-1; i >= 0; i--) {
                    if (plocs[i] >= prel_low) {
                        start = i;
                    }
                }

                /* found position_c where to start adding */
                start =  (start-1 );        /*  one step lower,  now start is  of original LF peaks to keep  */  
                start =  MAX(start, -1);    /*  limit for loop   */ 

                if (prel_high < plocs[0]) {
                    fin = 0;
                } else {
                    fin = (*n_plocs)+1;
                    for (i = 0; i < *n_plocs; i++) {
                        if (plocs[i] <= prel_high) {
                            fin = i;
                        }
                    }
                    fin++;
                }

            move_int(plocs_old, plocs, *n_plocs);
            move_float(f0est_old, f0est, *n_plocs);

            j = (start+1);      /*   [0..(j-1)] of original LF peaks will be kept */  
            /* j now points to first location_c where to add  peaks */

            for (i = 0; i < n_prel; i++) {
               plocs[j] = plocs_prel[i];
               f0est[j] = f0est_prel[i];
               j++;
            }
            for (i = fin; i < *n_plocs; i++) {
               plocs[j] = plocs_old[i];
               f0est[j] = f0est_old[i];
               j++;
            }

                *n_plocs = j;

            }
        }
    }

    return;
}

