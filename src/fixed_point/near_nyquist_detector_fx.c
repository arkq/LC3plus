/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNearNyquistdetector_fx(Word16 *near_nyquist_flag, const Word16 fs_idx, const Word16 near_nyquist_index,
                                   const Word16 bands_number, const Word32 *ener_fx, const Word16 ener_fx_exp
#ifdef ENABLE_HR_MODE
                                   , Word16 frame_dms, Word16 hrmode)
#else
                                   )
#endif
{
    *near_nyquist_flag = 0;
#ifdef ENABLE_HR_MODE
    IF (hrmode == 0){
#endif
        IF (sub(fs_idx, 4) < 0)
        {
            Dyn_Mem_Deluxe_In(
                Word16 i;
                Word16 nrg_above_thresh;
                Word16 ener_low_exp;
                Word16 ener_high_exp; 
                Word16 comp_energy_exp;
                Word32 comp_energy;
                Word32 ener_low; 
                Word32 ener_high;
            );

            ener_low = 0; move32();
            ener_low_exp = 0; move16();
            FOR (i = 0; i < near_nyquist_index; i++)
            {
                ener_low = BASOP_Util_Add_Mant32Exp(ener_fx[i], ener_fx_exp, ener_low, ener_low_exp, &ener_low_exp);
            }

            ener_high = 0; move32();
            ener_high_exp = 0; move16();
            FOR (i = near_nyquist_index; i < bands_number; i++)
            {
                ener_high = BASOP_Util_Add_Mant32Exp(ener_fx[i], ener_fx_exp, ener_high, ener_high_exp, &ener_high_exp);
            }
            
            comp_energy = Mpy_32_16(ener_low, NN_thresh); /* Mpy_32_16 -> 32Q15 */
            comp_energy_exp = add(add(ener_low_exp, NN_thresh_exp),15);

            nrg_above_thresh = BASOP_Util_Cmp_Mant32Exp(ener_high, ener_high_exp, comp_energy, comp_energy_exp); /* 1 if firstNumber > secondNumber */

            if (sub(nrg_above_thresh, 1) == 0)
            {
                *near_nyquist_flag = 1;
            }

            Dyn_Mem_Deluxe_Out();
        }
#ifdef ENABLE_HR_MODE
    } 
    ELSE // hrmode == 1 
    {
        // inverse spectral flatness = mean(energy) ./ 2^(mean(log2(energy)));
        Word32 td_thresh;

        SWITCH (frame_dms)
        {
        case 25:
            td_thresh = TD_HR_thresh_2_5ms;
            BREAK;
        case 50:
            td_thresh = TD_HR_thresh_5ms;
            BREAK;
        case 75:
            td_thresh = TD_HR_thresh_7_5ms;
            BREAK;
        default:                             /* 100 */
            td_thresh = TD_HR_thresh_10ms;
            BREAK;
        }

        Word16 mean_ener_exp = 0;

        Word32 sum_ener = 0; move32();
        Word16 sum_ener_exp = 0; move16();
        FOR (Word16 i = 0; i < bands_number; i++)
        {
            sum_ener = BASOP_Util_Add_Mant32Exp(ener_fx[i], ener_fx_exp, sum_ener, sum_ener_exp, &sum_ener_exp);
        }

        Word16 denom = sub(14,norm_s(bands_number));
        IF (sub(frame_dms, 50) == 0){
            denom = sub(15,norm_s(bands_number)); move16();
        }

        Word32 mean_ener = L_shr(sum_ener, denom); move32(); // = sum_ener / bands_number
        mean_ener_exp = sum_ener_exp; move16();

        Word32 sum_ener_log2 = 0;move32();
        Word16 sum_ener_log2_exp = 0;move16();
        Word32 mean_ener_log2 = 0;move32();
        
        FOR (Word16 i = 0; i < bands_number; i++)
        {    
            IF (ener_fx[i] != 0) {
                Word32 log2Value = L_add(BASOP_Util_Log2(ener_fx[i]), L_shl_pos(L_deposit_l(ener_fx_exp), 25));
                /* input argument is in Q7.25 , returns pow(2,(x/64) 
                    floatingpoint value log2_fl = log2Value/pow(2,31-6)  */
                sum_ener_log2 = BASOP_Util_Add_Mant32Exp(log2Value, 6, sum_ener_log2, sum_ener_log2_exp, &sum_ener_log2_exp); move32(); 
            }
        }
        mean_ener_log2 = L_shr(sum_ener_log2, denom); move32(); //mean_ener_log2 = sum_ener_log2 / bands_number
        Word16 mean_ener_log2_exp = sum_ener_log2_exp; 
        Word32 mean_ener_log2_fl = L_shr_pos(mean_ener_log2 ,s_min(31, sub(31,mean_ener_log2_exp))); //mean_ener_log2_fl = mean_ener_log2 / 2^(31-mean_ener_log2_exp)
        Word32 inv_flatness = 0;
        if (L_sub(norm_l(mean_ener),sub(sub(mean_ener_exp,31),mean_ener_log2_fl)) < 0 ) {
            inv_flatness = maxWord32;
        } 
        else {
            inv_flatness = L_shl(mean_ener, s_max(-31,sub(sub(mean_ener_exp,31),mean_ener_log2_fl)));
        }
        IF (L_sub(inv_flatness, td_thresh) > 0) {
            *near_nyquist_flag = 1; move16();
        } 
    }
    #endif // ENABLE_HR_MODE
}
