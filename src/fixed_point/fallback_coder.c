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

void fallback_encoder( UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word32 x[], Word16 L_spec, Word16 fallback_bit_planes )
{
    Dyn_Mem_Deluxe_In(
        Word16 i, j, sign;
        Word32 abs_value;
        UWord8 byte;
        Word16 byte_pointer;
        Word8 fallback_byte_planes;
    );
    
    byte_pointer = 0; move16();
    fallback_byte_planes = shr(sub(fallback_bit_planes, 1), 3);

    FOR ( i = 0; i < L_spec; i++ )
    {
        /* sign: write 1 for negative numbers, else 0 */
        if( x[i] < 0 )
        {
            sign = 1; move16();
        }
        else
        {
            sign = 0; move16();
        }
        
        write_bit_backward( bytes, bp_side, mask_side, sign );

        abs_value = L_abs( x[i] );
        FOR ( j = 0; j < fallback_byte_planes; j++ )
        {
            byte = L_and(L_shr_pos(abs_value, i_mult(8, j)), 0xFF);
            bytes[byte_pointer++] = byte; move16();
        }
    }
  
    Dyn_Mem_Deluxe_Out();
}

void fallback_decoder( UWord8* bytes, Word16* bp_side, Word16* mask_side, Word32 y[], Word16 L_spec, Word16 fallback_bit_planes )
{
    Dyn_Mem_Deluxe_In(
        Word16 i, j, sign;
        UWord8* ptr;
        Word16 byte_pointer;
        Word32 byte;
        Word8 fallback_byte_planes;
    );

    byte_pointer = 0; move16();
    ptr = bytes; move16();
    fallback_byte_planes = shr(sub(fallback_bit_planes, 1), 3);
    
    FOR ( i = 0; i < L_spec; i++ )
    {
        sign = read_bit( ptr, bp_side, mask_side );

        byte = 0; move16();
        FOR ( j = 0; j < fallback_byte_planes; j++ )
        {
            byte = L_or(L_shl_pos(ptr[byte_pointer++], i_mult(8, j)), byte);
        }
        
        IF (sub(sign, 1) == 0)
        {
            y[i] = L_negate(byte);
        } ELSE {
            y[i] = byte; move16();
        }
    }
  
    Dyn_Mem_Deluxe_Out();
}

#endif
