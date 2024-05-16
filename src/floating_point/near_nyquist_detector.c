/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNearNyquistdetector_fl(LC3_INT16* near_nyquist_flag, const LC3_INT fs_idx, const LC3_INT near_nyquist_index,
                                   const LC3_INT bands_number, const LC3_FLOAT* ener         , const LC3_INT16 frame_dms, const LC3_INT16 hrmode)                               
{
    *near_nyquist_flag = 0;
    if (hrmode == 0){
        if (fs_idx < 4)
        { 
            LC3_INT   i = 0;
            LC3_FLOAT ener_low = FLT_EPSILON, ener_high = 0;

            for (i=0; i<near_nyquist_index; i++)
            {
                ener_low += ener[i];
            }

            for (i=near_nyquist_index; i<bands_number; i++)
            {
                ener_high += ener[i];
            }

            if (ener_high > NN_thresh * ener_low){
                *near_nyquist_flag = 1;
            }
        }
    } 
    else // hrmode == 1 
    {
        // inverse spectral flatness = mean(energy) ./ 2^(mean(log2(energy)));
        LC3_INT32 td_thresh, i = 0;
        LC3_FLOAT mean_ener = 0, mean_ener_log2 = 0, inv_flatness = 0;

        switch (frame_dms)
        {
        case 25:
            td_thresh = TD_HR_thresh_2_5ms;
            break;
        case 50:
            td_thresh = TD_HR_thresh_5ms;
            break;
        case 75:
            td_thresh = TD_HR_thresh_7_5ms;
            break;
        default:                             /* 100 */
            td_thresh = TD_HR_thresh_10ms;
            break;
        }

        // calculate arithmetic mean
        for (i = 0; i < bands_number; i++)
        {
            mean_ener += ener[i];
        }
        mean_ener = mean_ener / bands_number;

        // calculate geometric mean
        for (i = 0; i < bands_number; i++)
        {    
            if (ener[i] != 0) {
                mean_ener_log2 += LC3_LOGTWO(ener[i]);
            }
        }
        mean_ener_log2 = mean_ener_log2 / bands_number;

        inv_flatness = mean_ener / LC3_POW(2,mean_ener_log2);

        if (inv_flatness > td_thresh) {
            *near_nyquist_flag = 1;
        } 
    }
}
