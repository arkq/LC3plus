/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void   Parcor2Index(const Word16 parCoeff[] /*Q15*/, Word16 index[], Word16 order);
static void   Index2Parcor(const Word16 index[], Word16 parCoeff[], Word16 order);
static Word32 FIRLattice(Word16 order, const Word16 *parCoeff /*Q15*/, Word32 *state, Word32 x /* Q0 */);

/*************************************************************************/

#ifdef CR14_A_ADD_LOSSLESS_MODE
void processAdaptiveTns_fx( Word16* bits, Word16 indexes[], Word32 x[], Word16 BW_cutoff_idx, Word16 order[], Word16* numfilters, Word16 enable_lpc_weighting, Word16 nSubdivisions, LC3PLUS_FrameDuration frame_dms, Word16 maxLen, lc3_scratch_t scratch
                         , Word16 hrmode, Word16 near_nyquist_flag, Word16 *ll_adap_flag, Word16 fs_idx, Word16 total_bits, Word16 ll_tns_lsb_num_remove_limit, Word16 *tns_lsb_num_remove)
{
    Dyn_Mem_Deluxe_In(     Word16 *RC[2];
    Word16 RC_1[MAXLAG];
    Word16 RC_2[MAXLAG];
    Word32 predictionGain[2];);

    /* TNS coder call with two execution modes. 1: only prediction gain and RC , 2: apply TNS filter,*/
    predictionGain[0] = 0;
    predictionGain[1] = 0;

    RC[0] = RC_1;
    RC[1] = RC_2;
    Word32* x_before_tns;
    Word32 active_bits;
    Word32 bit_usage_estimate;
    UWord8* tns_lsb_remove;
    *tns_lsb_num_remove = 0;
    TnsStartStopFreqs startstopfreqs;

    if(fs_idx==6)
    {
         maxLen = maxLen >> 1;
    }

    x_before_tns   = (Word32*) lc3_scratch_push( scratch, sizeof( *x_before_tns ) * maxLen );
    tns_lsb_remove = (UWord8*) lc3_scratch_push( scratch, sizeof( *tns_lsb_remove ) * maxLen );

    /* get predictionGain and RC coefficients from tns filter 1 and 2 */
    processTnsCoder_fx(bits, indexes, x, BW_cutoff_idx, order, numfilters, enable_lpc_weighting, nSubdivisions, frame_dms,
                    maxLen, scratch, hrmode, near_nyquist_flag, predictionGain, 1, RC, &startstopfreqs, 1);

    /* only apply TNS filter if prediction gain for at least one filter is > 1.5 */
    IF ( L_sub( predictionGain[0], 98304 ) > 0 || L_sub( predictionGain[1], 98304 ) > 0 ){

        /* copy values for later comparison */
        memcpy(x_before_tns, x, maxLen * sizeof(Word32));
        *ll_adap_flag = 1; move16();
        /* estimate active bits and bit usage after applied TNS filter */
        active_bits = estimate_active_bits_after_tns(x, maxLen, predictionGain, &startstopfreqs);
        bit_usage_estimate = L_add(estimate_bit_usage( fs_idx, active_bits, 0 ), (Word32)(( maxLen >> 2 ) + 10));

        IF (L_sub(bit_usage_estimate, (Word32)total_bits) < 0) {
            /* apply TNS filter directly */
            processTnsCoder_fx(bits, indexes, x, BW_cutoff_idx, order, numfilters, enable_lpc_weighting, nSubdivisions, frame_dms,
                        maxLen, scratch, hrmode, near_nyquist_flag, predictionGain, 2, RC, &startstopfreqs, 1);
        }
        ELSE {
            /* estimate the number of needed LSB remove bits to make TNS lossless */ 
            //*tns_lsb_num_remove = (UWord8) ((L_add(L_sub(bit_usage_estimate, (Word32)total_bits), (Word32)(maxLen - 1)))/(Word32)maxLen); /*480-1 = ceil*/
            Word32 tmp = L_sub(bit_usage_estimate, (Word32)total_bits);
            *tns_lsb_num_remove = 0;
            WHILE( tmp > 0 )
            {
                tmp = L_sub(tmp, maxLen);
                *tns_lsb_num_remove = add(*tns_lsb_num_remove, 1);
            }
            IF (*tns_lsb_num_remove > ll_tns_lsb_num_remove_limit)
            {
                /* apply TNS filter without LSB reduction */
                processTnsCoder_fx(bits, indexes, x, BW_cutoff_idx, order, numfilters, enable_lpc_weighting, nSubdivisions, frame_dms,
                            maxLen, scratch, hrmode, near_nyquist_flag, predictionGain, 2, RC, &startstopfreqs, 1);
                *tns_lsb_num_remove = 0; move16();
                *ll_adap_flag = 0; move16();
            }
            ELSE
            {
                /* remove LSBs and apply filter */
                basop_memset(tns_lsb_remove, *tns_lsb_num_remove, maxLen * sizeof(tns_lsb_remove[0]));
                process_lsb_remove(x, tns_lsb_remove, maxLen, x);
                processTnsCoder_fx(bits, indexes, x, BW_cutoff_idx, order, numfilters, enable_lpc_weighting, nSubdivisions, frame_dms,
                            maxLen, scratch, hrmode, near_nyquist_flag, predictionGain, 2, RC, &startstopfreqs, 1);
            }
        }

        Word32 energy_diff = 0; move32();

        Word16 scale_x = getScaleFactor32(x, maxLen);
        Word16 scale_x_before_tns = getScaleFactor32(x_before_tns, maxLen);
        
        Word16 scale_total = s_min(scale_x, scale_x_before_tns);

        FOR (Word16 i = 0; i < maxLen; i++) {
            Word16 tmp16 = extract_h(L_shl_pos(x[i], scale_total));
            Word32 e_x = L_mult0(tmp16, tmp16);
            
            tmp16 = extract_h(L_shl_pos(x_before_tns[i], scale_total));
            Word32 e_x_before_tns = L_mult0(tmp16, tmp16);
            
            energy_diff = L_add_sat(energy_diff, L_sub_sat(e_x, e_x_before_tns));
        }

        // Check if energy_after is greater than energy_before
        
        IF (energy_diff > 0) {
            memcpy(x, x_before_tns, maxLen * sizeof(Word32));
            *ll_adap_flag = 0; move16();
            *tns_lsb_num_remove = 0; move16();
            order[0] = 0; move16();
            order[1] = 0; move16();
            *bits = 2; move16();
        }
    }

    tns_lsb_remove = (UWord8*) lc3_scratch_pop( scratch, tns_lsb_remove );
    x_before_tns   = (Word32*) lc3_scratch_pop( scratch, x_before_tns );

    Dyn_Mem_Out();
}
#endif

void processTnsCoder_fx(Word16 *bits, Word16 indexes[], Word32 x[], Word16 BW_cutoff_idx, Word16 order[],
                        Word16 *numfilters, Word16 enable_lpc_weighting, Word16 nSubdivisions, LC3PLUS_FrameDuration frame_dms,
                        Word16 maxLen, lc3_scratch_t scratch
#ifdef ENABLE_HR_MODE
                        , Word16 hrmode
#endif
                       , Word16 near_nyquist_flag
#ifdef CR14_A_ADD_LOSSLESS_MODE
                        /* can take values 0: execute prediction gain and tns filter. 1: only prediciton gain. 2: only tns filter*/
                         , Word32 predictionGain[], Word16 execution_mode, Word16 ** RC, TnsStartStopFreqs* startstopfreqs, Word16 lossless
#endif
)
{
    Dyn_Mem_Deluxe_In(Word16 * tmpbuf; Word32 * rxx, epsP, *state, L_tmp, *A, alpha; 
                      Word16 n, n2, headroom, shift, tmp, shifts, facs, facs_e, stopfreq, xLen, maxOrder;
                      Word16 startfreq[TNS_NUMFILTERS_MAX]; const Word16 *subdiv_startfreq, *subdiv_stopfreq;
                      Counter i, j, iSubdivisions, lag;);
                      
#ifndef CR14_A_ADD_LOSSLESS_MODE
    Word32 predictionGain;
    Word16 *RC;
    Word16 inv;
#else
    Word16 inv;
#endif

    /* Buffer alignment */
    tmpbuf = (Word16*) lc3_scratch_push( scratch, sizeof( *tmpbuf ) * maxLen );
    rxx = (Word32*) lc3_scratch_push( scratch, sizeof( *rxx ) * ( MAXLAG + 1 ) );
    state = (Word32*) lc3_scratch_push( scratch, sizeof( *state ) * ( MAXLAG ) );
    A = (Word32*) lc3_scratch_push( scratch, sizeof( *A ) * ( MAXLAG + 1 ) );

#ifndef CR14_A_ADD_LOSSLESS_MODE
    RC = (Word16*) lc3_scratch_push( scratch, sizeof( *RC ) * MAXLAG );
#endif

    /* Init */
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (sub(execution_mode, 2) != 0)
#endif
    {
        *bits = 0;
        move16();
        maxOrder = MAXLAG;
        move16();
        *numfilters = 1;
        move16();
    }
    
    subdiv_startfreq = 0; 
    move16();
    subdiv_stopfreq  = 0; 
    move16();

#ifdef ENABLE_HR_MODE
    if (hrmode)
    {
        xLen = BW_cutoff_bin_all_HR[BW_cutoff_idx];
    }
    else
#endif
    {
        xLen = BW_cutoff_bin_all[BW_cutoff_idx];
    }
    move16();

    SWITCH (frame_dms)
    {
#ifdef CR9_C_ADD_1p25MS
    case LC3PLUS_FRAME_DURATION_1p25MS:
        *bits    = 0;
        order[0] = 0;
        order[1] = 0;
        *numfilters = 0;
        goto tns_exit;
#endif
    case LC3PLUS_FRAME_DURATION_2p5MS:
        startfreq[0] = 3;
        move16();
        
#ifdef ENABLE_HR_MODE
        if (hrmode)
        {
            subdiv_startfreq = tns_subdiv_startfreq_2_5ms_HR[BW_cutoff_idx - 4];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq_2_5ms_HR[BW_cutoff_idx - 4];
            move16();
        }
        else
#endif
        {
            subdiv_startfreq = tns_subdiv_startfreq_2_5ms[BW_cutoff_idx];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq_2_5ms[BW_cutoff_idx];
            move16();
        }
        xLen     = shr_pos(xLen, 2);
        maxOrder = 4;
        move16();
        BREAK;
    case LC3PLUS_FRAME_DURATION_5MS:
        startfreq[0] = 6;
        move16();
        
#ifdef ENABLE_HR_MODE
        if (hrmode)
        {
            subdiv_startfreq = tns_subdiv_startfreq_5ms_HR[BW_cutoff_idx - 4];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq_5ms_HR[BW_cutoff_idx - 4];
            move16();
        }
        else
#endif
        {
            subdiv_startfreq = tns_subdiv_startfreq_5ms[BW_cutoff_idx];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq_5ms[BW_cutoff_idx];
            move16();
        }
        xLen     = shr_pos(xLen, 1);
        maxOrder = 4;
        BREAK;
    case LC3PLUS_FRAME_DURATION_7p5MS:
        startfreq[0] = 9;
        move16();
        subdiv_startfreq = tns_subdiv_startfreq_7_5ms[BW_cutoff_idx];
        move16();
        subdiv_stopfreq = tns_subdiv_stopfreq_7_5ms[BW_cutoff_idx];
        move16();
        tmp      = shr_pos(xLen, 2);
        xLen     = add(tmp, add(tmp, tmp));
        maxOrder = 8;
        BREAK;
    case LC3PLUS_FRAME_DURATION_10MS:
        startfreq[0] = 12;
        move16();
        
#ifdef ENABLE_HR_MODE
        if (hrmode)
        {
            subdiv_startfreq = tns_subdiv_startfreq_HR[BW_cutoff_idx - 4];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq_HR[BW_cutoff_idx - 4];
            move16();
        }
        else
#endif
        {
            subdiv_startfreq = tns_subdiv_startfreq[BW_cutoff_idx];
            move16();
            subdiv_stopfreq = tns_subdiv_stopfreq[BW_cutoff_idx];
            move16();
        }
        BREAK;
    case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
    }

    IF (sub(BW_cutoff_idx, 3) >= 0 && frame_dms >= LC3PLUS_FRAME_DURATION_5MS)
    {
        *numfilters  = 2;
        startfreq[1] = shr_pos(xLen, 1);
    }
    
#ifdef CR14_A_ADD_LOSSLESS_MODE
    startstopfreqs->start_freq[0] = startfreq[0]; move16();
    startstopfreqs->stop_freq = startfreq[1]; move16();
    startstopfreqs->start_freq[1] = shr_pos( startfreq[1], 1 ); move16();
#endif 

    basop_memset(state, 0, MAXLAG * sizeof(*state));

#if defined(CR14_A_ADD_LOSSLESS_MODE) 
    if (lossless)
    {
        if ( near_nyquist_flag )
        {
            goto tns_exit;
        }
    }
#endif
    
    FOR (j = 0; j < *numfilters; j++)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if (sub(execution_mode, 2) != 0)
        {
#endif
        basop_memset(rxx, 0, (maxOrder + 1) * sizeof(*rxx));

        FOR (iSubdivisions = 0; iSubdivisions < nSubdivisions; iSubdivisions++)
        {
            n = sub(subdiv_stopfreq[nSubdivisions * j + iSubdivisions],
                    subdiv_startfreq[nSubdivisions * j + iSubdivisions]);

            /*norms[iFilter][iSubdivisions] = norm2FLOAT(pSpectrum+iStartLine, iEndLine-iStartLine);*/
            headroom = getScaleFactor32(x + subdiv_startfreq[nSubdivisions * j + iSubdivisions], n);

            /* Calculate norm of spectrum band */
            L_tmp = Norm32Norm(x + subdiv_startfreq[nSubdivisions * j + iSubdivisions], headroom, n, &shift);

            /* Rounding to avoid overflow when computing the autocorrelation below */
            tmp   = sub(norm_l(L_tmp), 1);
            L_tmp = L_shl(L_tmp, tmp);
            shift = sub(shift, tmp);
            L_tmp = L_add(L_tmp, 0x8000);
            L_tmp = L_and(L_tmp, 0x7FFF0000);

            IF (L_tmp == 0)
            {
                rxx[0] = 0x7FFFFFFF;
                move32();
                basop_memset(&rxx[1], 0, (maxOrder) * sizeof(*rxx));
                BREAK;
            }

            /* get pre-shift for autocorrelation */
            tmp    = sub(shift, norm_l(L_tmp)); /* exponent for normalized L_tmp */
            tmp    = shr_pos(sub(1, tmp), 1);   /* pre-shift to apply before autocorrelation */
            shifts = s_min(tmp, headroom);

            /* calc normalization factor */
            facs_e = shl_pos(sub(tmp, shifts), 1);

            SWITCH (frame_dms)
            {
#ifdef CR9_C_ADD_1p25MS
            case LC3PLUS_FRAME_DURATION_1p25MS: assert(0);
#endif
            case LC3PLUS_FRAME_DURATION_2p5MS: facs_e = add(facs_e, 1); BREAK;
            case LC3PLUS_FRAME_DURATION_5MS: facs_e = add(facs_e, 1); BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS: BREAK;
            case LC3PLUS_FRAME_DURATION_10MS: BREAK;
            case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
            }

            tmp   = sub(1, shl_pos(tmp, 1));       /* exponent of autocorrelation */
            L_tmp = L_shl(L_tmp, sub(shift, tmp)); /* shift L_tmp to that exponent */
            /* calc factor (with 2 bits headroom for sum of 3 subdivisions) */
            facs = div_s(0x2000, round_fx(L_tmp)); /* L_tmp is >= 0x2000000 */

            FOR (i = 0; i < n; i++)
            {
                tmpbuf[i] = round_fx_sat(L_shl_sat(x[subdiv_startfreq[nSubdivisions * j + iSubdivisions] + i], shifts));
                move16();
            }

            FOR (lag = 0; lag <= maxOrder; lag++)
            {
                n2 = sub(n, lag);
                L_tmp = L_deposit_l(0);
                FOR (i = 0; i < n2; i++)
                {
                    L_tmp = L_mac0(L_tmp, tmpbuf[i], tmpbuf[i + lag]);
                }
                if (lag != 0)
                    L_tmp = Mpy_32_32(L_tmp, tnsAcfWindow[lag - 1]);

                L_tmp = Mpy_32_16(L_tmp, facs);
                L_tmp = L_shl(L_tmp, facs_e);

                rxx[lag] = L_add(rxx[lag], L_tmp);
                move32();
            }
        }

        /* Levinson-Durbin */
#ifdef CR14_A_ADD_LOSSLESS_MODE
        processLevinson_fx( A, rxx, maxOrder, RC[j], &epsP, scratch );
#else
        processLevinson_fx(A, rxx, maxOrder, RC, &epsP, scratch);
#endif

        /* Prediction Gain */
        shift          = norm_l(epsP);
        inv            = div_s(16383, extract_h(L_shl_pos(epsP, shift)));
#ifdef CR14_A_ADD_LOSSLESS_MODE
        predictionGain[j] = Mpy_32_32( rxx[0], Mpy_32_16( L_sub( MAX_32, Mpy_32_16( L_shl( epsP, shift ), inv ) ), inv ) );
        IF( L_sub( predictionGain[j], L_shr_pos_pos( 0x30000000, shift ) ) > 0 && near_nyquist_flag == 0 )
        {
            /* If Prediction Gain is low */
            test();
            IF( enable_lpc_weighting != 0 && L_sub( predictionGain[j], L_shr_pos_pos( 0x40000000, shift ) ) < 0 )
            {
                /* LPC weighting */
                alpha = L_add( 0x6CCCCCCD,
                               Mpy_32_32( 0x13333333, L_shl_pos( L_sub( L_shl_pos( predictionGain[j], shift ), 0x30000000 ), 3 ) ) );
               L_tmp = alpha; move32();
#else  
        predictionGain = Mpy_32_32(rxx[0], Mpy_32_16(L_sub(MAX_32, Mpy_32_16(L_shl(epsP, shift), inv)), inv));

        IF (L_sub(predictionGain, L_shr_pos_pos(0x30000000, shift)) > 0 && near_nyquist_flag == 0)
        {
            /* If Prediction Gain is low */
            test();
            IF (enable_lpc_weighting != 0 && L_sub(predictionGain, L_shr_pos_pos(0x40000000, shift)) < 0)
            {
                /* LPC weighting */
                alpha = L_add(0x6CCCCCCD,
                              Mpy_32_32(0x13333333, L_shl_pos(L_sub(L_shl_pos(predictionGain, shift), 0x30000000), 3)));
                L_tmp = alpha;
#endif
                FOR (i = 1; i < maxOrder; i++)
                {
                    A[i] = Mpy_32_32(A[i], L_tmp);
                    move32();
                    L_tmp = Mpy_32_32(L_tmp, alpha);
                }
                A[maxOrder] = Mpy_32_32(A[maxOrder], L_tmp);
                move32();

                /* LPC -> RC */
#ifdef CR14_A_ADD_LOSSLESS_MODE
                lpc2rc( A, RC[j], maxOrder );
#else
                lpc2rc(A, RC, maxOrder);
#endif
            }

            /* Reflection Coefficients Quantization */
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Parcor2Index( RC[j], &indexes[MAXLAG * j], maxOrder );
#else
            Parcor2Index(RC, &indexes[MAXLAG * j], maxOrder);
#endif

            /* reduce filter order by truncating trailing zeros */
            i = sub(maxOrder, 1);
            WHILE ((i >= 0) && (indexes[MAXLAG * j + i] == INDEX_SHIFT))
            {
                i = sub(i, 1);
            }
            order[j] = add(i, 1);

            // Disable TNS if order is 0:
            IF (order[j] == 0) {
                // Jump to else statement
                goto tns_disabled;
            }
            /* Count bits */
            L_tmp = L_deposit_l(ac_tns_order_bits[enable_lpc_weighting][order[j] - 1]);
            FOR (i = 0; i < order[j]; i++)
            {
                L_tmp = L_add(L_tmp, L_deposit_l(ac_tns_coef_bits[i][indexes[MAXLAG * j + i]]));
            }
            
            *bits = add(*bits, add(2, extract_l(L_shr_pos(L_sub(L_tmp, 1), 11))));
            move16();

            /* Unquantize Reflection Coefficients */
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Index2Parcor( &indexes[MAXLAG * j], RC[j], order[j] );
#else
            Index2Parcor(&indexes[MAXLAG * j], RC, order[j]);
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
            IF (sub(execution_mode, 1) == 0)
            {
                /* write RC and order to output, don't apply filter */
            }
            ELSE
            { /* execution_mode != 1 : apply filter */
#endif 
            /* Stop frequency */
            stopfreq = xLen;
            move16();
            IF (sub(*numfilters, 2) == 0 && j == 0)
            {
                stopfreq = startfreq[1];
            }

            /* Filter */
            FOR (i = startfreq[j]; i < stopfreq; i++)
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                x[i] = FIRLattice( order[j], RC[j], state, x[i] );
#else
                x[i] = FIRLattice(order[j], RC, state, x[i]);
#endif
                move32();
            }
#ifdef CR14_A_ADD_LOSSLESS_MODE
            }
#endif 
        }
        ELSE
        {
tns_disabled:
            /* TNS disabled */
            *bits    = add(*bits, 1);
            order[j] = 0;
        }
#ifdef CR14_A_ADD_LOSSLESS_MODE
        }
        ELSE /* from here execution mode == 2 */
        {
            if (j == 0){
                basop_memset(state, 0, MAXLAG * sizeof(*state));
            }
            if (order[j] > 0){
                /* Stop frequency */
                stopfreq = xLen;
                move16();
                IF (sub(*numfilters, 2) == 0 && j == 0)
                {
                    stopfreq = startfreq[1];
                }
                /* Filter */
                FOR (i = startfreq[j]; i < stopfreq; i++)
                {
                    x[i] = FIRLattice(order[j], RC[j], state, x[i]);
                    move32();
                }
            }
        }
        IF (sub(execution_mode, 1) == 0) {
            predictionGain[j] = L_shr_pos_pos(predictionGain[j], 13-shift);
        }
#endif /* CR14_A_ADD_LOSSLESS_MODE_H */
    }

#ifdef CR9_C_ADD_1p25MS
tns_exit:
#endif

    Dyn_Mem_Deluxe_Out();

#ifndef CR14_A_ADD_LOSSLESS_MODE
    RC = (Word16*) lc3_scratch_pop( scratch, RC );
#endif
    A = (Word32*) lc3_scratch_pop( scratch, A );
    state = (Word32*) lc3_scratch_pop( scratch, state );
    rxx = (Word32*) lc3_scratch_pop( scratch, rxx );
    tmpbuf = (Word16*) lc3_scratch_pop( scratch, tmpbuf );
    

}

/*************************************************************************/

static void Parcor2Index(const Word16 parCoeff[] /*Q15*/, Word16 index[], Word16 order)
{
    Dyn_Mem_Deluxe_In(Counter i; Word16 iIndex; Word16 x;);

    FOR (i = 0; i < order; i++)
    {
        move16();
        move16();
        iIndex = 1;
        x      = parCoeff[i];

        WHILE ((iIndex < TNS_COEF_RES) && (x > tnsQuantThr[iIndex - 1]))
        {
            iIndex = add(iIndex, 1);
        }
        index[i] = sub(iIndex, 1);
        move16();
    }

    Dyn_Mem_Deluxe_Out();
}

static void Index2Parcor(const Word16 index[], Word16 parCoeff[], Word16 order)
{
    Counter i;
    FOR (i = 0; i < order; i++)
    {
        parCoeff[i] = tnsQuantPts[index[i]];
        move16();
    }
}

static Word32 FIRLattice(Word16 order, const Word16 *parCoeff /*Q15*/, Word32 *state, Word32 x /* Q0 */)
{
    Dyn_Mem_Deluxe_In(Counter i; Word32 tmpSave, tmp;);

    tmpSave = L_add(x, 0);

    FOR (i = 0; i < order - 1; i++)
    {
        tmp      = L_add(state[i], Mpy_32_16(x, parCoeff[i]));
        x        = L_add(x, Mpy_32_16(state[i], parCoeff[i])); /* exponent: 31+0 */
        state[i] = tmpSave;
        move32();
        tmpSave = L_add(tmp, 0);
    }

    /* last stage: only need half operations */
    x                = L_add(x, Mpy_32_16(state[order - 1], parCoeff[order - 1]));
    state[order - 1] = tmpSave;
    move32();
    Dyn_Mem_Deluxe_Out();
    return x;
}

