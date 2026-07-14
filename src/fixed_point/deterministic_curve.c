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

void process_deterministic_curve( UWord8* deterministic_curve, Word32 gg, Word16 gain_e, Word16 frame_length, Word32* int_scf_fx, Word16* int_scf_fx_exp, const Word16* bands_offset, Word16 bands_number, lc3_scratch_t scratch ) 
{
    Dyn_Mem_Deluxe_In(     
        Word32 *tmp;
        Word32 *invInput;
        Word32 gg_full_scale;
        Word16 max_scale_scf;
        Word16 max_scale_scf_tmp;
        Word16 shift;
        Word32 shift_round;
        Word32 tmp_mult;
        Word32 tmp_rounded;
    );
    
    tmp = (Word32*) lc3_scratch_push( scratch, sizeof( *tmp ) * frame_length );
    invInput = (Word32*) lc3_scratch_push( scratch, sizeof( *invInput ) * frame_length );
    
    /* Deterministic_Curve = log2(ceil(gg*mdct_shaping(ones(1,ctrl.FRAME), 1./q_spec_env_inv, ctrl.ACC_COEFF_PER_BAND))); */

    /* 1./int_scf_fx */
    FOR (int k = 0; k < bands_number; k++)
    {
        invInput[k] = invFixp(int_scf_fx[k], &int_scf_fx_exp[k]); move32();
    }
    
    max_scale_scf = 0; move16();
    /* Get maximum scale number for scf */
    FOR (int k = 0; k < bands_number; k++)
    {
        max_scale_scf = s_max(max_scale_scf, int_scf_fx_exp[k]);
    }
    
    /* Create array of ones in format of max_scale_scf */
    
    max_scale_scf_tmp = sub(31, max_scale_scf);
    tmp_mult =  L_shl_pos(1, max_scale_scf_tmp);
    FOR (int k = 0; k < frame_length; k++)
    {
        tmp[k] = tmp_mult; move32();
    }

    /* modifies tmp */
    processMdctShaping_fx( tmp, invInput, int_scf_fx_exp, bands_offset, bands_number );
    
    /* global gain is scaled with Q31 and additional Q-gain_e -> bring global gain to Q31 */
    IF (gain_e < 0)
    {
        gg_full_scale = L_shr_pos(gg, abs_s(gain_e));
        gain_e = 0; move16();
    } ELSE {
        gg_full_scale = gg; move32();
    }

    shift = sub(sub(31, gain_e), max_scale_scf);
    shift_round = L_sub(L_shl_pos( 1, shift ), 1);
    
    FOR (int k = 0; k < frame_length; k++)
    {
        /* Mpy_32_32: gg*mdct_shaping, mdct_shaping coeffs scaled to Q31 */
        /* L_shr_pos and + 1: ceil() */
        /* 31 - norm_l: log2(...) */

        /* multiplication with global gain */
        tmp_mult = Mpy_32_32( gg_full_scale, tmp[k] );

        /* ceil operation: add 1 if floor(tmp_mult) != tmp_mult */
        tmp_rounded = L_shr_pos( tmp_mult, shift );

        test();
        tmp_rounded = L_add (tmp_rounded, (tmp_mult & shift_round) != 0);
        tmp_rounded = L_sub(tmp_rounded, 1);

        /* ceil(log2(...))*/
        if ( tmp_rounded <= 0 ) 
        {
            tmp_rounded = -1; move16(); /* results in deterministic_curve[k] = 0; */
        }
        deterministic_curve[k] = sub(31, norm_l( tmp_rounded )); move16();
    }
  
    invInput = (Word32*) lc3_scratch_pop( scratch, invInput );
    tmp = (Word32*) lc3_scratch_pop( scratch, tmp );
  
    Dyn_Mem_Deluxe_Out();
}

void process_lsb_remove(Word32* q_d, UWord8* detcurve, Word16 frame_length, Word32* out)
{
    FOR (int k = 0; k < frame_length; k++)
    {
        IF (q_d[k] < 0)
        {
            out[k] = L_negate(L_shr_pos(L_abs(q_d[k]), detcurve[k]));       move32();
        } ELSE {
            out[k] = L_shr_pos(q_d[k], detcurve[k]);                        move32();
        }
    }
}

void process_lsb_add(Word32* q_d, UWord8* detcurve, Word16 frame_length, Word32* out)
{
    FOR (int k = 0; k < frame_length; k++)
    {
        out[k] = L_shl_pos(q_d[k], detcurve[k]);                            move32();
    }
}

static Word32 getBitWidth(Word32 value)
{
  IF (INT_MIN == value  || INT_MAX == value) {
    return sizeof(Word32) * 8 - 1;
  }
  else {
    UWord8 i = 0; 

    if (value < 0) {
      value = -value;   move32();
    }

    WHILE ((L_shr(value,i)) > 0) i++;

    return i;
  }
}

void lsb_mean_split(Word32* d_fx, Word16 startband, Word16 stopband, const Word16* band_offset,UWord8 lastSegementDetCurve, UWord8* deltaCodedBits, UWord8* determinstic_curve)
{

    Dyn_Mem_Deluxe_In(
    Word8 counter;
    UWord8 lastsegment;
    Word16 scale;
    );

    counter = 0; 
    lastsegment = 0; move16();

    //Fill detcurve with zeros in case of last det curve segment is zero
    IF( lastSegementDetCurve == 0 )
    {

        FOR (Word16 bin = band_offset[startband]; bin < band_offset[stopband]; bin++)
        {
            determinstic_curve[bin] = lastsegment; move16();
        }
        Dyn_Mem_Deluxe_Out();
        return;
    }


    FOR( Word16 i = startband ; i < stopband; i++ )
    {
        assert( counter < HIGH_BANDS_NUMBER);
        Word16 start = band_offset[i];          move16();
        Word16 stop  = band_offset[i + 1];      move16();
        Word16 width = stop - start;            move16();

        Word16 sumBitWidths = 0;                move16();
        Word16 avgBitWidths = 0;                move16();


        FOR (Word16 bin = start; bin < stop; bin++) {
            Word16 bitWidth = getBitWidth(d_fx[bin]);
            sumBitWidths = add(sumBitWidths,bitWidth);
        }

        /* average rounded to next lower integer */
        IF(sumBitWidths > 0 && (sumBitWidths > width) )
        {
            avgBitWidths = BASOP_Util_Divide1616_Scale(sumBitWidths, width, &scale); move16();
            avgBitWidths = shr(avgBitWidths, 15 - scale ); move16();
        }
        IF(avgBitWidths < 0)
        {
            avgBitWidths = 0;   move16();
        }
        
        IF (i == startband)
        {
            IF( avgBitWidths < lastSegementDetCurve)
            {
                deltaCodedBits[counter] = 1; move16();
                lastsegment = lastSegementDetCurve - 1; move16();
            }
            ELSE
            {
                lastsegment = lastSegementDetCurve; move16();
            }
        }
        ELSE
        {
            IF( (avgBitWidths < lastsegment) && (lastsegment > 0) )
            {
                    deltaCodedBits[counter] = 1; move16();
                    lastsegment--; move16();
            }            
        }

        FOR (int bin2 = start; bin2 < stop; bin2++)
        {
            determinstic_curve[bin2] = lastsegment; move16();
        }


        counter++; move16();

    }
    Dyn_Mem_Deluxe_Out();
    return;
}
#  endif
