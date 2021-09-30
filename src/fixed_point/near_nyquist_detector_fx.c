
#include "functions.h"

void processNearNyquistdetector_fx(Word16 *near_nyquist_flag, const Word16 fs_idx, const Word16 near_nyquist_index,
                                   const Word16 bands_number, const Word32 *ener_fx, const Word16 ener_fx_exp)
{
    *near_nyquist_flag = 0;
    
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
}
