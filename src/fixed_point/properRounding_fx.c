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

static Word32 sign( Word32 x )
{
    IF ( x > 0 )
    {
        return 1;
    }
    IF ( x < 0 )
    {
        return -1;
    }

    return 0;
}

void processProperRounding_fx( UWord8* detcurve, Word16 b_index, Word16 x_index, Word32* q_res, Word16 framesize, Word32* q_res_r, UWord8* eff_det_curve, lc3_scratch_t scratch )
{
    Word16 k, missing;
    Word32 offset;

    basop_memcpy(q_res_r, q_res, sizeof(*q_res_r) * framesize);

    /* rounding for partially transmitted coefficients */
    IF ( eff_det_curve != NULL )
    {
        FOR ( k = 0; k < framesize; k++ )
        {
            missing = sub(detcurve[k], eff_det_curve[k]);
            IF ( missing > 0 && q_res[k] != 0 )
            {
                offset = L_sub(L_shl_pos(1, sub(missing, 1)), 1);
                IF ( q_res[k] > 0 )
                {
                    q_res_r[k] = L_add(q_res[k], offset);
                }
                ELSE
                {
                    q_res_r[k] = L_sub(q_res[k], offset);
                }
            }
        }
        return;
    }

    {
    Word32 sumlsb = 0;
    Word32* new_curve;
    Word32* ind;
    UWord8* b;
    UWord8 fill_b;

    IF (b_index <= 0)
    {
        return;
    }

    new_curve = (Word32*) lc3_scratch_push( scratch, sizeof(*new_curve) * framesize );
    ind       = (Word32*) lc3_scratch_push( scratch, sizeof(*ind) * framesize );
    b         = (UWord8*) lc3_scratch_push( scratch, sizeof(*b) * framesize );

    basop_memset(new_curve, 0, sizeof(*new_curve) * framesize);
    basop_memset(ind, 0, sizeof(*ind) * framesize);

    fill_b = 1 * b_index;
    basop_memset(b, fill_b, sizeof(*b) * framesize);

    IF (x_index > 0)
    {
        FOR (int k = x_index + 1; k < framesize; k++)
        {
            b[k] = b[k] + 1;
        }
    }

    FOR (int k = 0; k < framesize; k++)
    {
        new_curve[k] = s_min(detcurve[k], b[k]);
    }

    /* Processing */
    int counter = 0;
    FOR (int k = 0; k < framesize; k++)
    {
        IF ((q_res[k] != 0) && (new_curve[k] > 0))
        {
            ind[counter] = k;
            counter++;
        }
    }

    FOR (int k = 0; k < counter; k++)
    {
        q_res_r[ind[k]] = q_res[ind[k]] + (L_shl_pos(1, new_curve[ind[k]]) - 1) * sign(q_res[ind[k]]);
    }

    /* Sanity checks */
    FOR (int k = 0; k < framesize; k++)
    {
        sumlsb += L_and(q_res[k], L_shl_pos(1, new_curve[k]) - 1);
    }

    assert(sumlsb == 0);
    UNUSED(sumlsb);

    b         = (UWord8*) lc3_scratch_pop( scratch, b );
    ind       = (Word32*) lc3_scratch_pop( scratch, ind );
    new_curve = (Word32*) lc3_scratch_pop( scratch, new_curve );
    }
}

#endif
