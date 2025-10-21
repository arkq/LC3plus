/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                          
#include "functions.h"

#ifdef FIX_BOTH_1p25_WB_GLOBGAINOFFSET_NONBE  
/* tilt factor in gainOffset quantized and adjusted for low Fs  and  1p25ms framing */
LC3_INT16 calc_GGainOffset_1p25(LC3_INT16 total_bits, LC3_INT16 fs_idx)
{
    LC3_INT16 gain_off_tilt_1p25_Q19[6] = { 20480, 17408, 17476, 13107, 10486, 8738 }; /*  vector of 1p25 tilts for  NB to UB */
   /* Corresponding FLT =  gain_off_tilt_1p25 = {0.039062500000000   0.033203125000000   0.033333333333333   0.025000000000000   0.020000000000000   0.016666666666667 }*/


    LC3_INT16 tmp1 = (LC3_INT16)(  ( ((LC3_INT32)total_bits)*((LC3_INT32)gain_off_tilt_1p25_Q19[fs_idx]) ) >> (3 + 16)  ); /*no rounding on purpose */
    LC3_INT16 tmp2 = 105 + 5 * (fs_idx + 1);

    tmp2 = -(MIN(115, tmp1) + tmp2);

#ifdef  FIX_BOTH_1p25_WB_GLOBGAINOFFSET_LOWLIM_NONBE 
    if (fs_idx <= 1)
    {   /* only NB and WB additionally limited to -135 */
        tmp2 = MAX(tmp2, FIX_BOTH_1p25_WB_GLOBGAINOFFSET_LOWLIM_NONBE);
    }
#endif 

    return tmp2;
}
#endif 

LC3_FLOAT array_max_abs(LC3_FLOAT *in, LC3_INT32 len)
{
        LC3_FLOAT max;
        LC3_INT32 i;
    
    max = LC3_FABS(in[0]);
    
    for (i = 0; i < len; i++)
    {
        if (LC3_FABS(in[i]) > LC3_FABS(max))
        {
            max = LC3_FABS(in[i]);
        }
    }
    
    return max;
}

