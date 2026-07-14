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

Word32 estimate_active_bits_after_tns(Word32 x[], Word16 frame_length, Word32 predictionGain[], TnsStartStopFreqs* startstopfreqs)
{
    Word16 active_bits_before_1, active_bits_before_2; 
    Word16 est_active_bits_after_1, est_active_bits_after_2;
    Word16 first_segment, last_segment;
    Word32 all_est_active_bits_after;

    Word16 start_1, start_2, stop_2;
    start_1 = startstopfreqs->start_freq[0];
    start_2 = startstopfreqs->start_freq[1];
    stop_2 = startstopfreqs->stop_freq;

    //    function all_est_active_bits_after = estimate_active_bits_after_tns( x, predictionGain, ctrl)
    //        active_bits_1_before = obj.get_active_bits(x(ctrl.tns_conf.startfreq(1):ctrl.tns_conf.stopfreq(1)));
    //active_bits_before_1 = calculate_active_bits(x + start_1, start_2 - start_1);
    active_bits_before_1 = calculate_active_bits(x + start_1, sub(start_2, start_1));
    //        active_bits_2_before = obj.get_active_bits(x(ctrl.tns_conf.startfreq(2):ctrl.tns_conf.stopfreq(2)));
    //active_bits_before_2 = calculate_active_bits(x + start_2, stop_2 - start_2);
    active_bits_before_2 = calculate_active_bits(x + start_2, sub(stop_2, start_2));
    // set to bits before in case TNS filter is not applied
    est_active_bits_after_1 = active_bits_before_1; move16();
    est_active_bits_after_2 = active_bits_before_2; move16();

   IF ( L_sub( predictionGain[0], 98304 ) > 0 )
   {
        Word32 pG_dB_1, active_bit_ratio_1;
        pG_dB_1 = L_add(BASOP_Util_Log2(predictionGain[0]),(15<<25));
        active_bit_ratio_1 = L_add(Mpy_32_32(pG_dB_1, -156442592),32484046);
        est_active_bits_after_1 = extract_l(L_shr_pos_pos(Mpy_32_16(active_bit_ratio_1, active_bits_before_1), 10));
   }

   IF ( L_sub( predictionGain[1], 98304 ) > 0 )
   {
        Word32 pG_dB_2, active_bit_ratio_2;
        // log10(x) = log2(x) * 1/log2(10)
        // predictionGain is in Q15.16, BASOP_Util_Log2 expects Q0.31
        // => compensate by adding 15 (in Q6.25 => '15<<25') to result
        pG_dB_2 = L_add(BASOP_Util_Log2(predictionGain[1]),(15<<25)); // in: Q0.31, out: Q6.25
    //   active_bit_ratio(2) = -0.0356 * pG_dB(2) + 1.0014 + est_buffer;
        // round(2^31*(-0.0356*10*1/log2(10))) = -230138690
        // round(2^25*(1.0014)) = 33601408
        active_bit_ratio_2 = L_add(Mpy_32_32(pG_dB_2, -230138690),33601408); // Q6.25
    //   est_active_bits_after(2) = active_bits_before(2) .* active_bit_ratio(2); 
        // Mpy_32_16: Q6.25*Q15.0 => Q21.10, shr by 10 to get Q31.0
        est_active_bits_after_2 = extract_l(L_shr_pos_pos(Mpy_32_16(active_bit_ratio_2, active_bits_before_2), 10));
   }       

   //   first_segment = obj.get_active_bits(x(1:ctrl.tns_conf.startfreq(1)-1));
   first_segment = calculate_active_bits(x, start_1);
   //   last_segment = obj.get_active_bits(x(ctrl.tns_conf.stopfreq(2)+1:end));
   //last_segment = calculate_active_bits(x + stop_2, frame_length - stop_2);
   last_segment = calculate_active_bits(x + stop_2, sub(frame_length,stop_2));
   //   all_est_active_bits_after = first_segment + est_active_bits_after(1) + est_active_bits_after(2) + last_segment;
   //all_est_active_bits_after = first_segment + est_active_bits_after_1 + est_active_bits_after_2 + last_segment;
   all_est_active_bits_after = add(first_segment,add(est_active_bits_after_1,add(est_active_bits_after_2,last_segment)));
   //   all_est_active_bits_after = all_est_active_bits_after + 130;
   //all_est_active_bits_after = all_est_active_bits_after + 130;
   all_est_active_bits_after = add(all_est_active_bits_after, 130);
   return all_est_active_bits_after;
}

Word32 calculate_active_bits(Word32 x[], Word16 frame_length)
{
    Word32 zero_cnt_after_last_non_zero = 0; move32();
    Word32 active_bits=0; move32();

    FOR(Word32 i = 0; i < frame_length; i++)
    {
        IF(x[i]==0)
        {
            //active_bits += 1;
            active_bits = L_add(active_bits,1);
            //zero_cnt_after_last_non_zero += 1;
            zero_cnt_after_last_non_zero = L_add(zero_cnt_after_last_non_zero,1);
        }
        ELSE
        {
            //active_bits += 31 - norm_l(L_abs(x[i]));
            active_bits = L_add(active_bits,sub(31, norm_l(L_abs(x[i]))));
            zero_cnt_after_last_non_zero = 0; move32();
        }
    }

    active_bits = L_sub(active_bits,zero_cnt_after_last_non_zero);

    return active_bits;
}

Word32 estimate_bit_usage(Word16 fs_idx, Word32 active_bits, Word32 entropy_bits)
{
    Word32 bit_usage_estimate = 0; move32();
    assert(active_bits <= 960 * 24); // debugging
    Word32 active_bits_2 = L_mult0(active_bits,active_bits); // L_mult0: Q15.0*Q15.0=>Q30.0 - should work for 192kHz due to splitting the frame 
    assert( fs_idx >= 4 );
    IF (sub(fs_idx, 4) == 0)
    {
        // bit_usage_estimate = -3.152e-12 * powf(active_bits, 4) + 7.027e-08 * powf(active_bits, 3) - 4.951e-04 * powf(active_bits, 2) + 2.406 * active_bits - 88.05;
        // => Change to integer arithmetic:

        /* 
        bit_usage_estimate =-3.152e-12 * powf(active_bits, 4)    | = -( (active_bits^2 * sqrt(3.152e-12)) * (active_bits^2 * sqrt(3.152e-12)) )
                            + 7.027e-08 * powf(active_bits, 3)   |   +  (active_bits * 7.027e-8^(2/3)) * active_bits^2 * 7.027e-8^(1/3) )
                            - 4.951e-04 * powf(active_bits, 2)   |   -(  active_bits^2 * 4.951e-4  )
                            + 2.406 * active_bits                |   + 2.046 * active_bits
                            - 88.05                              |   - 88.05
        */
       

      
       // sqrt(3.152e-12) = 1.775387281693772e-06 | round(x*2^(31+19)) = 1998908375, using max range of Word32
       // active_bits_2: Q30.0
       // Mpy_32_32 (inner): Q30.0*Q0.50=Q50 => Q19
       // Mpy_32_32 (outer): Q19*Q19=Q38     => Q7  ==> sum all terms in Q7
       Word32 p4_1 = Mpy_32_32(Mpy_32_32(active_bits_2, 1998908375), 
                               Mpy_32_32(active_bits_2, 1998908375));
       
       // 7.027e-08^(2/3) = 1.702864018826761e-05 | round(x*2^31) = 36569
       // 7.027e-08^(1/3) = 0.004126577297018     | round(x*2^(31+7)) = 1134304930
       // max active bits for fs=48kHz is max. frame_size (48000/s * 0.01s = 480) * 31 = 14880 => 36569*14880, fits into Word32
       // active_bits_2: Q30.0
       // UL_Mpy_32_32: Q0.31*Q15.0=>Q31
       // Mpy_32_32 (inner): Q30.0*Q0.31   => Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+7) => Q7 
       Word32 p3_1 = Mpy_32_32(Mpy_32_32(active_bits_2, UL_Mpy_32_32(36569UL,L_deposit_l(active_bits))), 1134304930);
       
       // 4.951e-4<<(31+7) = 136092052
       // active_bits_2: Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+7) => Q7 
       Word32 p2_1 = Mpy_32_32(active_bits_2, 136092052);
       
       // round(2.406*2^29) = 1291711414; lsh(active_bits,2+7) to get to precision of summation (Q7)
       Word32 p1_1 = Mpy_32_32( L_shl(L_deposit_l(active_bits),2+7), 1291711414); 
       
       // round(88.05*2^7) = 11270
       Word32 p0_1 = 11270; // no move32(), compiler should replace with constant when used
       
       // -p4_1 + p3_1 - p2_1 + p1_1 - p0_1 = (p3_1+p1_1) - (p4_1+p2_1+p0_1)
       Word32 p = L_sub(L_add(p3_1,p1_1),(L_add(p4_1,L_add(p2_1,p0_1))));
       // active_bits=0 => p<0 for the first frame(s?)
       bit_usage_estimate = L_shr_pos(p,7);
    }
    ELSE IF (sub(fs_idx, 5) == 0)
    {
        // bit_usage_estimate = 7.412e-09 * powf(active_bits, 3) - 1.670e-04 * powf(active_bits, 2) + 2.182 * active_bits + 228.8;
        // => Change to integer arithmetic:
        /* 
        bit_usage_estimate =  7.412e-09 * powf(active_bits, 3)   |   (active_bits * 7.412e-09^(2/3)) * active_bits^2 * 7.412e-09^(1/3)
                            - 1.670e-04 * powf(active_bits, 2)   |   - 1.670e-04 * active_bits^2
                            + 2.182 * active_bits                |   + 2.182 * active_bits
                            + 228.8                              |   + 228.8
        */

       // 7.412e-9^(2/3) = 3.801517030488451e-06 | round(x*2^31) = 8164
       // 7.412e-9^(1/3) = 0.001949747940245     | round(x*2^(31+9)) = 2143770532
       // active_bits_2: Q30.0
       // L_mult0: Q0.31*Q15.0=>Q31
       // Mpy_32_32 (inner): Q30.0*Q0.31   => Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+9) => Q9     => sum all terms in Q9
       Word32 p3_1 = Mpy_32_32(Mpy_32_32(active_bits_2, L_mult0(8164,active_bits)), 2143770532);
       
       // 1.670e-4<<(31+9) = 183618442
       // active_bits_2: Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+9) => Q9
       Word32 p2_1 = Mpy_32_32(active_bits_2, 183618442);

       // round(2.182*2^29) = 1171452330; lsh(active_bits,2+9) to get to precision of summation (Q9)
       Word32 p1_1 = Mpy_32_32( L_shl(L_deposit_l(active_bits),2+9), 1171452330);

       // round(228.8*2^9) = 117146
       Word32 p0_1 = 117146; // no move32(), compiler should replace with constant when used

       Word32 p = L_sub(L_add(p3_1,L_add(p1_1,p0_1)), p2_1);
       bit_usage_estimate = L_shr_pos_pos(p,9);
    }
    ELSE IF (sub(fs_idx, 6) == 0  )
    {
        // bit_usage_estimate = 7.412e-09 * powf(active_bits, 3) - 1.670e-04 * powf(active_bits, 2) + 2.182 * active_bits + 228.8;
        // => Change to integer arithmetic:
        /* 
        bit_usage_estimate =  7.412e-09 * powf(active_bits, 3)   |   (active_bits * 7.412e-09^(2/3)) * active_bits^2 * 7.412e-09^(1/3)
                            - 1.670e-04 * powf(active_bits, 2)   |   - 1.670e-04 * active_bits^2
                            + 2.182 * active_bits                |   + 2.182 * active_bits
                            + 228.8                              |   + 228.8
        */

       // 7.412e-9^(2/3) = 3.801517030488451e-06 | round(x*2^31) = 8164
       // 7.412e-9^(1/3) = 0.001949747940245     | round(x*2^(31+9)) = 2143770532
       // active_bits_2: Q30.0
       // L_mult0: Q0.31*Q15.0=>Q31
       // Mpy_32_32 (inner): Q30.0*Q0.31   => Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+9) => Q9     => sum all terms in Q9
       Word32 p3_1 = Mpy_32_32(Mpy_32_32(active_bits_2, L_mult0(8164,active_bits)), 2143770532);
       
       // 1.670e-4<<(31+9) = 183618442
       // active_bits_2: Q30.0
       // Mpy_32_32 (outer): Q30.0*Q(31+9) => Q9
       Word32 p2_1 = Mpy_32_32(active_bits_2, 183618442);

       // round(2.182*2^29) = 1171452330; lsh(active_bits,2+9) to get to precision of summation (Q9)
       Word32 p1_1 = Mpy_32_32( L_shl(L_deposit_l(active_bits),2+9), 1171452330);

       // round(228.8*2^9) = 117146
       Word32 p0_1 = 117146; // no move32(), compiler should replace with constant when used

       Word32 p = L_sub(L_add(p3_1,L_add(p1_1,p0_1)), p2_1);
       bit_usage_estimate = L_shr_pos_pos(p,9) + L_shr_pos_pos(L_mult0(entropy_bits >> 1, 26214), 14) ;

    }
    return bit_usage_estimate;
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
/* Inverse of the encoder's L_shr(input, scaleSignal). For the per-frame /
 * switching variant the prev frame's scale and the lookahead range
 * (la_zeroes) and the prev-was-lossy flag would let us only rescale the
 * non-overlap part of the buffer, but in fixed-shift mode every frame uses
 * the same scaleSignal so a plain L_shl on the whole buffer is correct. */
void rescale_signal_decoder(Word32* x_fx_ip, Word16 scaleSignal, Word16 scaleSignalMemory, Word16 frame_length)
{
    UNUSED(scaleSignalMemory);
    IF(scaleSignal)
    {
        FOR ( int n = 0 ; n < frame_length ; n++)
        {
            x_fx_ip[n] = L_shl(x_fx_ip[n], scaleSignal);
        }
    }
}
#endif

Word16 get_ll_adap_flag( Word32 bit_usage_estimate, Word32 total_bits, Word16 bit_balance )
{
    Word16 ll_adap_flag;

    //IF( bit_usage_estimate <= (total_bits + bit_balance) )
    IF( L_sub(bit_usage_estimate, L_add(total_bits, L_deposit_l(bit_balance))) <= 0 )
    {
        ll_adap_flag = 1; move16();
    }
    ELSE
    {
        ll_adap_flag = 0; move16();
    }

    return ll_adap_flag;
}
#endif
