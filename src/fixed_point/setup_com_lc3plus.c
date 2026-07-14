/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef FIX_BOTH_1p25_WB_GLOBGAINOFFSET_NONBE
/* tilt factor in gainOffset quantized and adjusted for low Fs  and  1p25ms framing */
#ifdef CR15_A_LOSSLESS_1p25MS
Word16 calc_GGainOffset_1p25_fx(Word32 total_bits, Word16 fs_idx)
#else
Word16 calc_GGainOffset_1p25_fx(Word16 total_bits, Word16 fs_idx)
#endif
{
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word16 gain_off_tilt_1p25_Q19[7] = { 20480, 17408, 17476, 13107, 10486, 8738, 7490 }; /*  vector for NB to 192k */
    /* Corresponding FLT =  gain_off_tilt_1p25 = {0.039062500000000   0.033203125000000   0.033333333333333   0.025000000000000   0.020000000000000   0.016666666666667   0.0142857142857143}*/
#else
    Word16 gain_off_tilt_1p25_Q19[6] = { 20480, 17408, 17476, 13107, 10486, 8738 }; /*  vector for NB to UB */
    /* Corresponding FLT =  gain_off_tilt_1p25 = {0.039062500000000   0.033203125000000   0.033333333333333   0.025000000000000   0.020000000000000   0.016666666666667 }*/
#endif

#ifdef CR15_A_LOSSLESS_1p25MS
    Word16 tmp1 = extract_h(Mpy_32_16(L_shl(total_bits, 12), gain_off_tilt_1p25_Q19[fs_idx]));  /* extract bits in Q0 , no rounding for interop*/
#else
    Word16 tmp1 = extract_h(L_shr_pos(L_mult0(total_bits, gain_off_tilt_1p25_Q19[fs_idx]), 3));  /* extract bits in Q0 , no rounding for interop*/
#endif
    Word16 tmp2 = extract_l(L_mac0(105L, 5, add(fs_idx, 1)));

    tmp2 = negate(add(s_min(115, tmp1), tmp2)); move16();

 
#ifdef  FIX_BOTH_1p25_WB_GLOBGAINOFFSET_LOWLIM_NONBE 
    if (sub(fs_idx, 1) <= 0)
    {   /* only NB and WB additionally limited to -135 */
        tmp2 = s_max(tmp2, FIX_BOTH_1p25_WB_GLOBGAINOFFSET_LOWLIM_NONBE);
    }
#endif 
 

    return tmp2;
}
#endif 


