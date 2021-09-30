#include "functions.h"


void processNearNyquistdetector_fl(LC3_INT16* near_nyquist_flag, const LC3_INT fs_idx, const LC3_INT near_nyquist_index,
                                   const LC3_INT bands_number, const LC3_FLOAT* ener)
{
    *near_nyquist_flag = 0;
    
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
