/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                              
#include "functions.h"

void processApplyGlobalGain_fx(Word32 x[], Word16 *x_e, Word16 xLen, Word16 global_gain_idx, Word16 global_gain_off)
{
    Counter i;
#ifdef ENABLE_HR_MODE
    Word32 global_gain;
#else
    Word16 global_gain;
#endif
    Word16 global_gain_e;
    Word32 tmp32;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processApplyGlobalGain_fx", sizeof(struct {
                   Counter i;
                   Word16  global_gain, global_gain_e;
                   Word32  tmp32;
               }));
#endif

#ifdef ENABLE_HR_MODE
    /* 1 / (28 * log 2) is 0x797D in Q18, L_shl_pos by 7 results in Q25 tmp32 */
    /* round(2^31 / (28 * log10(2))) = 254778081 */
    //tmp32 = L_shl_pos(Mpy_32_16(254778081, add(global_gain_idx, global_gain_off)), 9);
    Word32  mh;
    UWord16 ml;

    Mpy_32_16_ss(254778081, add(global_gain_idx, global_gain_off), &mh, &ml);
    tmp32 = L_shl_pos(mh, 9) | L_deposit_l((shr((Word16)ml, 7)) & 0x1ff);
    move16();
    /* Uses an argument in Q25 */
    global_gain = BASOP_Util_InvLog2(L_or(tmp32, (Word32)0xFE000000));
#else
    tmp32         = L_shl_pos(L_mult0(add(global_gain_idx, global_gain_off), 0x797D), 7);
    /* Uses an argument in Q25 */
    global_gain = round_fx(BASOP_Util_InvLog2(L_or(tmp32, (Word32)0xFE000000)));
#endif
    global_gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);

    FOR (i = 0; i < xLen; i++)
    {
#ifdef ENABLE_HR_MODE
        x[i] = Mpy_32_32(x[i], global_gain);
#else
        x[i] = Mpy_32_16(x[i], global_gain);
#endif
        move32();
    }

    *x_e = add(*x_e, global_gain_e); move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

