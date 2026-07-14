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

/* emit the less significant bits of the integer values that the prioritization rule selected */
void residual_encoder_lossless( Word32 x[], UWord8 res_bits[], UWord8 det_curve[], UWord8 eff_det_curve[], Word16 L_spec, Word32 max_resBits_len, Word32 *bit_pos, lc3_scratch_t scratch)
{
    Dyn_Mem_Deluxe_In(
        Word16 i, j, d, t, j_stop;
        Word32 abs_x;
        );

    UNUSED(scratch);
    UNUSED(max_resBits_len);

    if ( scratch->max_scratch_calculation_only )
    {
        Dyn_Mem_Deluxe_Out();
        return;
    }

    FOR ( i = 0; i < L_spec; i++ )
    {
        d = det_curve[i];
        t = eff_det_curve[i];
        IF ( d == 0 || t == 0 )
        {
            CONTINUE;
        }
        abs_x = L_abs( x[i] );

        /* no MSBs above det_curve -> write sign bit */
        IF ( L_shr( abs_x, d ) == 0 )
        {
            res_bits[*bit_pos] = (x[i] < 0) ? 1 : 0;
            (*bit_pos)++;
        }

        /* write top t data bits out of d, MSB first */
        j_stop = sub(d, t);
        FOR ( j = sub(d, 1); j >= j_stop; j-- )
        {
            res_bits[*bit_pos] = s_and( L_shr( abs_x, j ), 1 );
            (*bit_pos)++;
        }
    }
    assert (*bit_pos <= max_resBits_len); /* check allocated resbits buffer size */
    Dyn_Mem_Deluxe_Out();
}

/* reconstruct the less significant bits that the prioritization rule included in the variable length information */
Word32 residual_decoder_lossless( Word32 x[], Word32 y[], UWord8 res_bits[], UWord8 det_curve[], UWord8 eff_det_curve[], Word16 L_spec, Word32 nbits_res, Word16 *index_b, Word16 *index_x, lc3_scratch_t scratch )
{
    Dyn_Mem_Deluxe_In(
        Word16 i, j, d, t, avail;
        Word32 bit_pos, abs_x, lsb, sgn;
        );

    UNUSED(scratch);

    *index_x = -1;
    *index_b = 0;

    IF ( nbits_res == 0 )
    {
        basop_memmove( y, x, sizeof(*x) * L_spec );
        Dyn_Mem_Deluxe_Out();
        return 0;
    }

    bit_pos = 0;
    FOR ( i = 0; i < L_spec; i++ )
    {
        d = det_curve[i];
        t = eff_det_curve[i];

        IF ( d == 0 || t == 0 )
        {
            y[i] = x[i];
            CONTINUE;
        }

        abs_x = L_abs( x[i] );

        /* no MSBs above det_curve -> read sign bit */
        IF ( L_shr( abs_x, d ) == 0 )
        {
            sgn = (res_bits[bit_pos] == 0) ? 1 : -1;
            bit_pos++;
        }
        ELSE
        {
            sgn = (x[i] > 0) ? 1 : -1;
        }

        /* read t data bits, MSB first */
        lsb = 0;
        avail = t;
        FOR ( j = sub(avail, 1); j >= 0; j-- )
        {
            lsb = L_add( L_shl_pos( lsb, 1 ), res_bits[bit_pos] );
            bit_pos++;
        }
        if ( avail < d )
        {
            lsb = L_shl_pos( lsb, sub(d, avail) );
        }

        y[i] = L_add( abs_x, lsb );
        if ( sgn < 0 )
        {
            y[i] = L_negate( y[i] );
        }
    }

    Dyn_Mem_Deluxe_Out();
    return bit_pos;
}

/* fine per-bin granularity: number of bins of the last block that still receive the next plane */
static Word16 compute_binIdxSplit(
    Word16 b,
    Word16 transmitted_b,
    Word16 d_band_b,
    Word32 x[],
    const Word16* bands_offset,
    Word32 budget
)
{
    Word16 start = bands_offset[b];
    Word16 stop  = bands_offset[b + 1];
    Word16 width = sub(stop, start);
    Word16 k, split_k;
    Word32 running, bin_cost;

    IF (transmitted_b > 0)
    {
        return extract_l(budget);
    }

    split_k = 0;
    running = 0;
    FOR (k = 0; k < width; k++)
    {
        bin_cost = 1;
        IF (L_shr(L_abs(x[start + k]), d_band_b) == 0)
        {
            bin_cost = L_add(bin_cost, 1);
        }
        IF (L_add(running, bin_cost) > budget)
        {
            BREAK;
        }
        running = L_add(running, bin_cost);
        split_k = add(k, 1);
    }
    return split_k;
}

/* prioritization rule for less significant bits; b_relative in the bitstream selects first (absolute numeric weight) or second (relative numeric weight) rule, frequency is the second-ranking criterion in both */
void compute_resbits_priority(
    Word32 x[],
    UWord8 det_curve[],
    const Word16* bands_offset,
    Word16 bands_number,
    Word16 L_spec,
    Word32 budget,
    Word16 b_relative,
    UWord8 eff_det_curve[],
    lc3_scratch_t scratch
)
{
    Word16 b, k, cost;
    Word16 d_band[MAX_BANDS_NUMBER];
    Word16 w_band[MAX_BANDS_NUMBER];
    Word16 s_band[MAX_BANDS_NUMBER];
    Word16 transmitted[MAX_BANDS_NUMBER];
    Word16 max_d, level, layer;
    Word16 split_band, binIdxSplit, done;
    Word32 total_cost;

    UNUSED(scratch);

    /* pre-compute per-band info */
    total_cost = 0;
    max_d = 0;
    FOR ( b = 0; b < bands_number; b++ )
    {
        Word16 start = bands_offset[b];
        Word16 stop  = bands_offset[b + 1];
        d_band[b] = det_curve[start];
        w_band[b] = sub(stop, start);
        s_band[b] = 0;

        IF ( d_band[b] > 0 )
        {
            FOR ( k = start; k < stop; k++ )
            {
                IF ( L_shr( L_abs( x[k] ), d_band[b] ) == 0 )
                {
                    s_band[b] = add(s_band[b], 1);
                }
            }
            total_cost = L_add(total_cost, (Word32)add(i_mult(d_band[b], w_band[b]), s_band[b]));
        }
        transmitted[b] = 0;
        if ( d_band[b] > max_d ) { max_d = d_band[b]; }
    }

    /* if everything fits, copy det_curve directly */
    IF ( total_cost <= budget )
    {
        basop_memcpy( eff_det_curve, det_curve, sizeof(*det_curve) * L_spec );
        return;
    }

    split_band  = -1;
    binIdxSplit = 0;
    done        = 0;

    IF ( b_relative == 0 )
    {
        /* first prioritization rule: absolute numeric weight as highest-ranking criterion */
        FOR ( level = max_d; level >= 1; level-- )
        {
            IF ( done ) { BREAK; }
            FOR ( b = 0; b < bands_number; b++ )
            {
                IF ( d_band[b] < level ) { CONTINUE; }
                IF ( sub(d_band[b], transmitted[b]) <= 0 ) { CONTINUE; }

                cost = w_band[b];
                if ( transmitted[b] == 0 )
                {
                    cost = add(cost, s_band[b]);
                }

                IF ( (Word32)cost > budget )
                {
                    /* last-block split: switch to per-bin granularity for this band */
                    binIdxSplit = compute_binIdxSplit(b, transmitted[b], d_band[b], x, bands_offset, budget);
                    if ( binIdxSplit > 0 )
                    {
                        split_band = b;
                    }
                    done = 1;
                    BREAK;
                }

                budget = L_sub(budget, (Word32)cost);
                transmitted[b] = add(transmitted[b], 1);
            }
        }
    }
    ELSE
    {
        /* second prioritization rule: relative numeric weight as highest-ranking criterion */
        FOR ( layer = 1; layer <= max_d; layer++ )
        {
            IF ( done ) { BREAK; }
            FOR ( b = 0; b < bands_number; b++ )
            {
                IF ( d_band[b] < layer ) { CONTINUE; }
                IF ( sub(d_band[b], transmitted[b]) <= 0 ) { CONTINUE; }

                cost = w_band[b];
                if ( transmitted[b] == 0 )
                {
                    cost = add(cost, s_band[b]);
                }

                IF ( (Word32)cost > budget )
                {
                    /* last-block split: switch to per-bin granularity for this band */
                    binIdxSplit = compute_binIdxSplit(b, transmitted[b], d_band[b], x, bands_offset, budget);
                    if ( binIdxSplit > 0 )
                    {
                        split_band = b;
                    }
                    done = 1;
                    BREAK;
                }

                budget = L_sub(budget, (Word32)cost);
                transmitted[b] = add(transmitted[b], 1);
            }
        }
    }

    /* expand to per-coefficient, split band gets two values */
    FOR ( b = 0; b < bands_number; b++ )
    {
        Word16 start = bands_offset[b];
        Word16 stop  = bands_offset[b + 1];
        IF ( b == split_band )
        {
            FOR ( k = start; k < add(start, binIdxSplit); k++ )
            {
                eff_det_curve[k] = (UWord8) add(transmitted[b], 1);
            }
            FOR ( k = add(start, binIdxSplit); k < stop; k++ )
            {
                eff_det_curve[k] = (UWord8) transmitted[b];
            }
        }
        ELSE
        {
            FOR ( k = start; k < stop; k++ )
            {
                eff_det_curve[k] = (UWord8) transmitted[b];
            }
        }
    }
    /* trailing coefficients */
    FOR ( k = bands_offset[bands_number]; k < L_spec; k++ )
    {
        eff_det_curve[k] = 0;
    }
}

#endif
