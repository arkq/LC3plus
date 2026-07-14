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
Word32 low_zero_search_lossless(Word32 x[],Word32 frame_length)
{
       /* Estimate bit usage */
        Word16 idx = 0;
        Word32 cnt = 0;

        FOR (Word16 k = 0; k < frame_length; k++)
        {
            IF (L_abs(x[k]) > 0)
            {
                IF( k > idx )
                {
                    idx = k;
                }
            }
        }

        IF (idx != 0)
        {        
            cnt = 0;
            FOR (Word16 k = 0; k < idx-1; k++)
            {
                IF (x[k] == 0)
                {
                    cnt++;
                }
            }     
        }   
        return cnt;
}
#endif 
