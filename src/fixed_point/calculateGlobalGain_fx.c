/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef CR14_A_ADD_LOSSLESS_MODE

void processCalculateGlobalGain_fx( Word32 *gg, Word16 *gg_e, Word16 global_gain_idx, Word16 global_gain_off )
{
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 tmp32;

    tmp32 = Mpy_32_16( 0x797CD707, L_shl_pos( add( global_gain_idx, global_gain_off ), 6 ) );
#else
    Word32 mh, tmp32;
    UWord16 ml;

    Mpy_32_16_ss( 254778081, add( global_gain_idx, global_gain_off ), &mh, &ml );
    tmp32 = L_shl_pos( mh, 9 ) | L_deposit_l( ( shr( (Word16) ml, 7 ) ) & 0x1ff );
#endif
    move16();
    /* Uses an argument in Q25 */
    *gg = BASOP_Util_InvLog2( L_or( tmp32, (Word32) 0xFE000000 ) );

    /* Mpy_32_16 result is in Q16 converted to Q18 using the upshift */
    *gg_e = add( extract_l( L_shr_pos( tmp32, 25 ) ), 1 );
}

#endif
