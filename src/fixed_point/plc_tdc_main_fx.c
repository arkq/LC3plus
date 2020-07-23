/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "defines.h"
#include "functions.h"


/*****************************************************************************/

static Word16 TDC_Dot_product(const Word16 x[], const Word16 y[], const Word16 lg);
static Word32 TDC_L_Dot_product(const Word16 x[], const Word16 y[], const Word16 lg);
static void   TDC_highPassFiltering_fx(const Word16 L_buffer, Word16 exc2[], const Word16 l_fir_fer,
                                       const Word16 *hp_filt);
static Word32 TDC_calcGainp(Word16 x[], Word16 y[], Word16 lg);
static void   TDC_calcGainc(Word16 *exc, Word16 Q_exc, Word16 old_fpitch, Word16 lg, Word16 frame_dms, Word16 lp_gainp, Word32 *lp_gainc);
static void   TDC_random_fx(Word16 *seed, Word16 lg, Word16 *y);
static Word16 TDC_preemph(Word16 *x, const Word16 fac, const Word16 lg);
static void   TDC_LPC_residu_fx(const Word16 *a, Word16 *x, Word16 *y, Word16 lg, Word16 m);
static void   TDC_deemph_fx(const Word16 *x, Word16 *y, const Word16 fac, const Word16 lg, const Word16 mem);
static void   TDC_LPC_synthesis_fx(const Word16 sh, const Word16 a[], const Word16 x[], Word16 y[], const Word16 lg,
                                   const Word16 m);
static void   TDC_normalize_energy_fx(Word16 *gain, Word16 *gain_exp, const Word16 *x, const Word16 lg);
 

/*****************************************************************************/


/*
 * processTimeDomainConcealment_Apply_fx
 *
 * Parameters:
 *   pitch_int                i  : integer pitch lag                 Q0
 *   preemphFac_fx            i  : preemphase factor                 Q15
 *   A_fx                     i  : lp filter coefficients            Qx
 *   pcmbufHist_fx            i  : pointer to input signal           Qq_fx_old_exp
 *   frame_length             i  : frame length                      Q0
 *   fs_idx                   i  : sample rate index                 Q0
 *   nbLostFramesInRow        i  : number of consecutive lost frames Q0
 *   overlap                  i  : overlap length                    Q0
 *   stabFac_fx               i  : stability factor                  Q15
 *   fract                    i/o: fraction of lag                   Q0
 *   seed_fx                  i/o: pointer to seed                   Q0
 *   gain_p_fx                i/o: pointer to gainp                  Q15
 *   gain_c_fx                i/o: pointer to gainc                  15Q16
 *   synth_fx                 o  : pointer to synthesized signal     Q_syn
 *   Q_syn                    o  : exponent for synthesized signal   Q0
 *   alpha                    o  : damping factor                    Q15
 *   scratchBuffer            i  : scratch buffer
 *
 * Function:
 *    Perform the time domain concealment.
 *
 * Returns:
 *    void
 */
void processTimeDomainConcealment_Apply_fx(const Word16 pitch_int, const Word16 preemphFac_fx, const Word16 *A_fx,
                                           const Word16 lpc_order, const Word16 *pcmbufHist_fx, const Word16 frame_length,
                                           const Word16 frame_dms, const Word16 fs_idx, const Word16 nbLostFramesInRow,
                                           const Word16 overlap, const Word16 stabFac_fx, Word16 *fract,
                                           Word16 *seed_fx, 
                                           Word32 *gain_c_fx, Word16 *synth_fx, Word16 *Q_syn, Word16 *alpha, Word16 max_len_pcm_plc,
                                           Word16 harmonicBuf_fx[MAX_PITCH], Word16 synthHist_fx[M], Word16 *const harmonicBuf_Q,
                                           Word8 *scratchBuffer)
{
    Counter       i;
    Word16        s, s1, c1, c2, len, cnt, g_fx, ilen, Tc, nextInc, beforeNextInc;
    Word32        tmp32, tmp32_2, gainc_tmp;
    Word16        gain_p_fx;
    Word32        gain_c_32_fx;
    Word16        gain_c_16_fx, gain_c_16_fx_exp, gain_inov_fx, gain_inov_fx_exp, ilen_exp;
    Word16        hpBlendFac;
    Word16        len_pi_lf_2, frame_length_2, step_fx, step_n_fx, gain_h_fx, nbLostCmpt_loc, mem_deemph;
    Word16 *      synth_mem_fx, *synth_tmp_fx, *exc2_fx, *exc_fx, *pt_exc, *pt1_exc, *x_pre_fx;
    Word16 *      harmonicBufPtr;
    Word16        Q_exc = 0, exp_scale;
    const Word16 *hp_filt_fx, *TDC_high_harm;
    Word16        alphaPrev_fx;
    Word16        throttle;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processTimeDomainConcealment_Apply_fx", sizeof(struct {
                   Counter i;
                   Word16 s, s1, c1, c2, len, cnt, g_fx, ilen, Tc, nextInc, beforeNextInc;
                   Word32 tmp32, tmp32_2, gainc_tmp;
                   Word16  gain_p_fx;
                   Word32  gain_c_32_fx;
                   Word16  gain_c_16_fx, gain_c_16_fx_exp, gain_inov_fx, gain_inov_fx_exp, ilen_exp;
                   Word16  hpBlendFac;
                   Word16  len_pi_lf_2, frame_length_2, step_fx, step_n_fx, gain_h_fx, nbLostCmpt_loc, mem_deemph;
                   Word16 *      synth_mem_fx, *synth_tmp_fx, *exc2_fx, *exc_fx, *pt_exc, *pt1_exc, *x_pre_fx;
                   Word16 *      harmonicBufPtr;
                   Word16        Q_exc, exp_scale;
                   const Word16 *hp_filt_fx, *TDC_high_harm;
                   Word16        alphaPrev_fx;
                   Word16        throttle;
               }));
#endif


    /* len of output signal */
    len = add(frame_length, overlap);

    nbLostCmpt_loc = nbLostFramesInRow; move16();
    nextInc = 1;  move16();
    beforeNextInc = 1;  move16();
    SWITCH (frame_dms)
    {
        case 25: 
            nbLostCmpt_loc = shr(add(nbLostFramesInRow, 3), 2);
            nextInc = (nbLostFramesInRow & 0x0003) == 1;  move16();
            beforeNextInc = (nbLostFramesInRow & 0x0003) == 0;  move16();
            BREAK;
        case 50:
            nbLostCmpt_loc = shr(add(nbLostFramesInRow, 1), 1);
            nextInc = (nbLostFramesInRow & 0x0001) == 1;  move16();
            beforeNextInc = (nbLostFramesInRow & 0x0001) == 0;  move16();
            BREAK;
    }

    IF (sub(nbLostCmpt_loc, PLC_FADEOUT_IN_MS / 10) > 0)
    {
        gain_p_fx  = 0; move16();
        *gain_c_fx = 0; move32();
        *Q_syn     = 0; move16();
        *alpha     = 0; move16();
        basop_memset(synth_fx, 0, len * sizeof(Word16));
#ifdef DYNMEM_COUNT
        Dyn_Mem_Out();
#endif
        return;
    }

    frame_length_2 = shr_pos(frame_length, 1);

    Tc = pitch_int; move16();
    if (sub(*fract, 0) > 0)
    {
        Tc = add(Tc, 1);
    }

    len_pi_lf_2 = add(Tc, frame_length_2);

    /*----------------------------------------------------------------
     * Buffer Initialization
     *
     *                exc_fx       synth_mem_fx
     * |--exc_buf_past--|--exc_buf_curr--|--syn_mem--|--x_pre--|
     *                                               |--exc2--|
     *                                               |--syn--|
     *
     *---------------------------------------------------------------*/

    /* pointer inits */
    exc_fx       = (Word16 *)scratchAlign(scratchBuffer,
                                    sizeof(Word16) * len_pi_lf_2); /* MAX_PITCH+MAX_LEN/2 + MAX_LEN+MDCT_MEM_LEN_MAX */
    synth_mem_fx = (Word16 *)scratchAlign(exc_fx, sizeof(*exc_fx) * len);           /* M */
    x_pre_fx     = (Word16 *)scratchAlign(synth_mem_fx, sizeof(*synth_mem_fx) * lpc_order); /* MAX_PITCH+MAX_LEN/2+M+1 */
    exc2_fx      = (Word16 *)scratchAlign(synth_mem_fx, sizeof(*synth_mem_fx) * lpc_order); /* MAX_LEN+MDCT_MEM_LEN_MAX+TDC_L_FIR_HP-1 */
    synth_tmp_fx = (Word16 *)scratchAlign(synth_mem_fx, sizeof(*synth_mem_fx) * lpc_order); /* MAX_LEN+MDCT_MEM_LEN_MAX */
    /* Buffers 'overlap' since they are not used at the same time */

    /*---------------------------------------------------------------*
     * LPC Residual                                                  *
     *---------------------------------------------------------------*/
    IF (sub(nbLostFramesInRow, 1) == 0)
    {
 
    /* copy buffer to pre-emphasis buffer */
        cnt = add(len_pi_lf_2, lpc_order + 1);
        basop_memmove(&x_pre_fx[0], &pcmbufHist_fx[max_len_pcm_plc - cnt], cnt * sizeof(Word16));

        /* apply pre-emphasis to the signal; x_pre = x_pre_flt * 2^(q_fx_old_exp-15-Q_exc+1) */
        Q_exc = TDC_preemph(&(x_pre_fx[1]), preemphFac_fx, sub(cnt, 1));

        /* copy memory for LPC synth */
        basop_memmove(&synth_mem_fx[0], &x_pre_fx[len_pi_lf_2 + 1], lpc_order * sizeof(Word16));

        /* LPC Residual; exc = exc_fx * 2^(q_fx_old_exp-15-Q_exc) */
        TDC_LPC_residu_fx(A_fx, &(x_pre_fx[lpc_order + 1]), &(exc_fx[-len_pi_lf_2]), len_pi_lf_2, lpc_order);
    }

    /*---------------------------------------------------------------*
     * Calculate gains                                               *
     *---------------------------------------------------------------*/

    IF (sub(nbLostFramesInRow, 1) == 0)
    {
        IF (sub(pitch_int, Tc) == 0)
        {
            gain_p_fx =
                round_fx_sat(L_shl_sat(TDC_calcGainp(&(x_pre_fx[lpc_order + Tc + 1]), &(x_pre_fx[lpc_order + 1]), frame_length_2), 15));
        }
        ELSE
        {
            tmp32   = TDC_calcGainp(&(x_pre_fx[lpc_order + Tc + 1]), &(x_pre_fx[lpc_order + 2]), frame_length_2);
            tmp32_2 = TDC_calcGainp(&(x_pre_fx[lpc_order + Tc + 1]), &(x_pre_fx[lpc_order + 1]), frame_length_2);

            IF (L_sub(tmp32, tmp32_2) > 0)
            {
                Tc         = pitch_int; move16();
                gain_p_fx = round_fx_sat(L_shl_sat(tmp32, 15));
                *fract     = 0; move16();
            }
            ELSE
            {
                gain_p_fx = round_fx_sat(L_shl_sat(tmp32_2, 15));
            }
        }

        if (gain_p_fx < 0)
        {
            gain_p_fx = 0; move16();
        }

        IF (sub(pitch_int, Tc) == 0)
        {
            /* gain_c = gain_c_32_fx * 2^(q_fx_old_exp-31) */
            TDC_calcGainc(exc_fx, Q_exc, Tc, frame_length_2, frame_dms, gain_p_fx, &gain_c_32_fx);
        }
        ELSE
        {
            TDC_calcGainc(exc_fx, Q_exc, pitch_int, frame_length_2, frame_dms, gain_p_fx, &tmp32);
            TDC_calcGainc(exc_fx, Q_exc, Tc       , frame_length_2, frame_dms, gain_p_fx, &gain_c_32_fx);

            gain_c_32_fx = L_min(gain_c_32_fx, tmp32); move32();
        }
    }
    ELSE
    {
        gain_c_32_fx = *gain_c_fx; move32();
        gain_p_fx = *alpha;
    }

    /*---------------------------------------------------------------*
     * Damping factor                                                *
     *---------------------------------------------------------------*/

    alphaPrev_fx = 0x7FFF; move16();
    IF (sub(nbLostFramesInRow,1) > 0)
    {
        alphaPrev_fx = *alpha; move16();
    }
  
    IF (nextInc != 0)
    {
        IF (sub(nbLostCmpt_loc, 1) == 0)
        {
            /* Threshold 31470 is 0.98^2 in Q15 format */
            IF (sub(gain_p_fx, 31470) > 0)
            {
                *alpha = 0x7D71; /*0.98f*/
                move16();
            }
            /* Threshold 28037 is 0.925^2 in Q15 format */
            ELSE IF (sub(gain_p_fx, 28037) < 0)
            {
                *alpha = 0x7666; /*0.925f*/
                move16();
            }
            ELSE
            {
                exp_scale = 0;
                *alpha    = Sqrt16(gain_p_fx, &exp_scale); move16();
                *alpha    = shl(*alpha, exp_scale);
            }
        }
        ELSE
        {
            SWITCH (nbLostCmpt_loc)
            {
            case 2:
                c1 = 0x50A4; /*0.630f*/
                move16();
                c2 = 0x2CCD; /*0.350f*/
                move16();
                BREAK;
            default:
                c1 = 0x5375; /*0.652f*/
                move16();
                c2 = 0x29FC; /*0.328f*/
                move16();
                BREAK;
            }

            *alpha = mult_r(stabFac_fx, c2);
            *alpha = add(*alpha, c1);

            *alpha = mult(gain_p_fx, *alpha);


            IF (sub(nbLostCmpt_loc, 2) == 0)
            {
                if (sub(*alpha, 0x75A2 /*0.919f*/) < 0)
                {
                    *alpha = 0x75A2; move16();
                }
            }
            ELSE IF (sub(nbLostCmpt_loc, 5) > 0)
            {
                gain_p_fx = *alpha; move16();
            }
        }
    }

    IF (sub(nbLostCmpt_loc,3) > 0)
    {
      SWITCH (frame_dms)
      {
      case  25: *alpha = mult(*alpha, PLC34_ATTEN_FAC_025_FX); BREAK;
      case  50: *alpha = mult(*alpha, PLC34_ATTEN_FAC_050_FX); BREAK;
      case 100: *alpha = mult(*alpha, PLC34_ATTEN_FAC_100_FX); BREAK;
      }
    }
    if (sub(nbLostCmpt_loc, 5) > 0)
    {
        gain_p_fx = *alpha; move16();
    }



    /*---------------------------------------------------------------*
     * Construct the harmonic part                                   *
     *  Last pitch cycle of the previous frame is repeatedly copied. *
     *---------------------------------------------------------------*/

    pt_exc  = harmonicBuf_fx; move16();
    pt1_exc = exc_fx - Tc;    move16();
    s       = s_min(len, Tc); move16();
    test();
    IF (sub(nbLostFramesInRow, 1) == 0)
    {
        *harmonicBuf_Q = Q_exc; move16();
        IF (sub(stabFac_fx, 32767 /*1.f Q15*/) >= 0)
        {
            basop_memmove(pt_exc, pt1_exc, Tc * sizeof(Word16));
        }
        ELSE
        {
            /* These values are necessary for the last five filtered samples */
            basop_memmove(exc_fx, &exc_fx[-Tc], (TDC_L_FIR_HP-1)/2 * sizeof(Word16));
            TDC_high_harm = TDC_high_32_harm;
            if (sub(fs_idx, 1) <= 0)
            {
                TDC_high_harm = TDC_high_16_harm;
            }
            FOR (i = 0; i < Tc; i++)
            {
                pt_exc[i] = TDC_Dot_product(&pt1_exc[i-(TDC_L_FIR_HP-1)/2], TDC_high_harm, TDC_L_FIR_HP);
            }
        }
    }
    ELSE
    {
        Q_exc = *harmonicBuf_Q; move16();
    }

    
    /*---------------------------------------------------------------*
     * Construct the random part of excitation                       *
     *---------------------------------------------------------------*/

    TDC_random_fx(seed_fx, add(len, sub(TDC_L_FIR_HP, 1)), exc2_fx);

    /* high pass noise */
    hp_filt_fx = TDC_high_32;
    if (sub(fs_idx, 1) <= 0)
    {
        hp_filt_fx = TDC_high_16;
    }

    IF (sub(nbLostFramesInRow, 1) == 0)
    {
        TDC_highPassFiltering_fx(len, exc2_fx, TDC_L_FIR_HP, hp_filt_fx);
    }
    ELSE
    {
        throttle = div_s(nbLostCmpt_loc, add(nbLostCmpt_loc, PLC3_HPBLENDTHROTTLE));
        hpBlendFac = mult(sub(0x7FFF, *alpha), throttle);
        c1 = sub(0x7FFF, hpBlendFac);
        FOR (i = 0; i < len; i++)
        {
            /* Return value of dot product is Q1 */
            tmp32      = Mpy_32_16(TDC_L_Dot_product(&exc2_fx[i], hp_filt_fx, TDC_L_FIR_HP), c1 /*Q15*/);
            exc2_fx[i] = round_fx(L_mac0(tmp32, hpBlendFac, exc2_fx[i+TDC_L_FIR_HP/2])); move16();
        }
    }

    /* normalize energy */
    TDC_normalize_energy_fx(&gain_inov_fx, &gain_inov_fx_exp, exc2_fx, frame_length);
    tmp32 = Mpy_32_16(
        L_sub(590558016l /*1.1 Q29*/, Mpy_32_16(L_shr_pos(L_deposit_h(gain_p_fx), 2), 24576 /*0.75*/)) /*Q29*/,
        gain_inov_fx /*Q15,gain_inov_e*/); /*Q29,gain_inov_e*/
    s                = norm_l(tmp32);
    tmp32            = L_shl_pos(tmp32, s);
    tmp32            = L_min(tmp32, 0x7FFEFFFF);
    gain_inov_fx_exp = add(sub(gain_inov_fx_exp, s), 31 - 29); /*->Q31*/
    gain_inov_fx     = round_fx(tmp32);                        /*Q15,gain_inov_e*/

    /* gains */
    gain_h_fx = alphaPrev_fx; move16();

    /* update steps */
    ilen = BASOP_Util_Divide1616_Scale((Word16)1, frame_length, &ilen_exp);
    step_fx = round_fx(L_shl(L_mult(sub(gain_h_fx, *alpha), ilen), ilen_exp));

    s     = norm_l(gain_c_32_fx);
    tmp32 = L_shl_pos(gain_c_32_fx, s);

    gain_c_16_fx     = extract_h(tmp32);
    gain_c_16_fx_exp = sub(15, s);

    gainc_tmp = L_mult (gain_c_16_fx, div_s(*alpha, alphaPrev_fx));
    tmp32     = L_sub (tmp32, gainc_tmp);
    step_n_fx = round_fx(L_shl(Mpy_32_16(tmp32, ilen), ilen_exp));

    /*---------------------------------------------------------------*
     * Construct the total excitation                                *
     *---------------------------------------------------------------*/

    harmonicBufPtr = harmonicBuf_fx + ((nbLostFramesInRow - 1) * frame_length) % Tc;

    s1  = add(Q_exc, add(gain_inov_fx_exp, gain_c_16_fx_exp));
    cnt = add(frame_length, TDC_L_FIR_HP / 2);

    g_fx = mult_r(gain_c_16_fx, gain_inov_fx);

    FOR (i = 0; i < len; i++)
    {
        /* harmonic */
        if (harmonicBufPtr - harmonicBuf_fx >= Tc) {
            harmonicBufPtr = harmonicBuf_fx;
        }
        exc_fx[i] = *harmonicBufPtr++;
        tmp32 = L_mult(exc_fx[i], gain_h_fx);
        /* random */
        tmp32_2 = L_shl_sat(L_mult(exc2_fx[i], g_fx), s1);
        /* total */
        exc_fx[i] = round_fx_sat(L_add_sat(tmp32, tmp32_2)); move16();
        /* update */
        gain_h_fx    = s_max(sub(gain_h_fx, step_fx), 0);
        gain_c_16_fx = s_max(sub(gain_c_16_fx, step_n_fx), 0);
        g_fx = mult_r(gain_c_16_fx, gain_inov_fx);
    }

    /* update gain */
    *gain_c_fx = L_shl(gainc_tmp, sub(gain_c_16_fx_exp, 15)); move32();


    /*----------------------------------------------------------*
     * Compute the synthesis speech                             *
     *----------------------------------------------------------*/

    /* introduce some headroom to avoid Overflows, 2 bit seem to be sufficient */
    *Q_syn = sub(Q_exc, 2);
    *Q_syn = s_max(*Q_syn, -3);

    exp_scale = sub(*Q_syn, Q_exc - 1);

    IF (sub(nbLostFramesInRow, 1) != 0)
    {
        synth_mem_fx = synthHist_fx;
    }
    Copy_Scale_sig(synth_mem_fx, &synth_tmp_fx[-lpc_order], lpc_order, exp_scale);
    TDC_LPC_synthesis_fx(sub(Q_exc, *Q_syn), A_fx, exc_fx, synth_tmp_fx, len, lpc_order);

    FOR (i=0; i<lpc_order; i++)
    {
        synthHist_fx[i] = shr_sat(synth_tmp_fx[frame_length-lpc_order+i], exp_scale);
    }
  
    /*----------------------------------------------------------*
     * Deemphasis                                               *
     *----------------------------------------------------------*/

    mem_deemph = shl(pcmbufHist_fx[max_len_pcm_plc - 1], *Q_syn);
    TDC_deemph_fx(synth_tmp_fx, synth_fx, preemphFac_fx, len, mem_deemph);

    /*----------------------------------------------------------*
     * Fade to zero                                             *
     *----------------------------------------------------------*/

    IF (beforeNextInc != 0)
    {
        IF (sub(nbLostCmpt_loc, PLC_FADEOUT_IN_MS / 10) == 0)
        {
            gain_h_fx = (Word16)0x7FFF; move16();
            step_fx   = round_fx(L_shl(L_mult(gain_h_fx, ilen), ilen_exp));
            FOR (i = 0; i < frame_length; i++)
            {
                assert(gain_h_fx >= 0);
                synth_fx[i] = mult(synth_fx[i], gain_h_fx);
                gain_h_fx   = sub(gain_h_fx, step_fx);
            }
            basop_memset(&synth_fx[frame_length], 0, overlap * sizeof(Word16));
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


/*****************************************************************************/

static Word32 syn_kern_2(Word32 L_tmp, const Word16 a[], const Word16 y[])
{
    L_tmp = L_msu_sat(L_tmp, y[-1], a[1]);
    L_tmp = L_msu_sat(L_tmp, y[-2], a[2]);
    return L_tmp;
}

static Word32 syn_kern_4(Word32 L_tmp, const Word16 a[], const Word16 y[])
{
    L_tmp = syn_kern_2(L_tmp, a, y);
    return syn_kern_2(L_tmp, a + 2, y - 2);
}

static Word32 syn_kern_8(Word32 L_tmp, const Word16 a[], const Word16 y[])
{
    L_tmp = syn_kern_4(L_tmp, a, y);
    return syn_kern_4(L_tmp, a + 4, y - 4);
}

static Word32 syn_kern_16(Word32 L_tmp, const Word16 a[], const Word16 y[])
{
    L_tmp = syn_kern_8(L_tmp, a, y);
    return syn_kern_8(L_tmp, a + 8, y - 8);
}

/*
 * TDC_Dot_product
 *
 * Parameters:
 *   x     i: x vector       Q0
 *   y     i: y vector       Q0
 *   lg    i: vector length  Q0
 *
 * Function:
 *   dot product
 *
 * Returns:
 *   dot product             Q0
 */
static Word16 TDC_Dot_product(const Word16 x[], const Word16 y[], const Word16 lg)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  sum;
    );

    sum = L_mult(x[0], y[0]);
    FOR (i = 1; i < lg-1; i++)
    {
        sum = L_mac_sat(sum, x[i], y[i]);
    }
    
    sum = mac_r_sat(sum, x[lg-1], y[lg-1]);

    Dyn_Mem_Deluxe_Out();

    return sum;
}

/*
 * TDC_L_Dot_product
 *
 * Parameters:
 *   x     i: x vector       Q0
 *   y     i: y vector       Q0
 *   lg    i: vector length  Q0
 *
 * Function:
 *   dot product
 *
 * Returns:
 *   dot product              Q1
 */
static Word32 TDC_L_Dot_product(const Word16 x[], const Word16 y[], const Word16 lg)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  sum;
    );

    sum = L_mac0(1L, x[0], y[0]);
    FOR (i = 1; i < lg; i++)
    {
        sum = L_mac0(sum, x[i], y[i]);
    }

    Dyn_Mem_Deluxe_Out();

    return sum;
}

/*
 * TDC_highPassFiltering_fx
 *
 * Parameters:
 *   L_buffer     i: buffer length
 *   exc2         i: unvoiced excitation before the high pass filtering  Qx/Qx+1
 *   l_fir_fer    i: length highpass filter
 *   hp_filt      i: highpass filter coefficients                        Q15
 *
 * Function:
 *   Highpass filter
 *
 * Returns:
 *   void
 */
static void TDC_highPassFiltering_fx(const Word16 L_buffer, Word16 exc2[], const Word16 l_fir_fer,
                                     const Word16 *hp_filt)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    FOR (i = 0; i < L_buffer; i++)
    {
        exc2[i] = round_fx(L_sub(TDC_L_Dot_product(&exc2[i], hp_filt, l_fir_fer), 1)); move16();
    }

    Dyn_Mem_Deluxe_Out();
}

/*
 * TDC_calcGainc
 *
 * Parameters:
 *   exc        i: pointer to excitation buffer
 *   Q_exc      i: Q format of excitation buffer
 *   old_fpitch i: pitch_int
 *   lg         i: length
 *   lp_gainp   i: gain p
 *   lp_gainc   o: pointer to gain (15Q16)
 *
 * Function:
 *   Gain calculation
 *
 * Returns:
 *   void
 */
static void TDC_calcGainc(Word16 *exc, Word16 Q_exc, Word16 old_fpitch, Word16 lg, Word16 frame_dms, Word16 lp_gainp, Word32 *lp_gainc)
{
    Dyn_Mem_Deluxe_In(
        Word16  tmp16, tmp_e, tmp2_e;
        Word32  L_tmp, L_tmp_max;
        Counter i;
    );

    L_tmp = L_deposit_l(0);

    FOR (i = 0; i < lg; i++)
    {
        /* gain_c += ( exc[-i-1] - *gain_p * exc[-i-1-pitch_int] ) * ( exc[-i-1] - *gain_p * exc[-i-1-pitch_int] ); */
        tmp16 = sub_sat(exc[i - lg] /*Q1*/, mult_r(lp_gainp /*Q15*/, exc[i - lg - old_fpitch] /*Q1*/) /*Q1*/);
        L_tmp = L_mac0_sat(L_tmp, tmp16, tmp16); /*Q3*/
    }

    IF (sub(frame_dms, 100) < 0)
    {
        L_tmp_max = L_deposit_l(0);
        FOR (i = 0; i < lg; i++)
        {
            L_tmp_max = L_mac0_sat(L_tmp_max, exc[i - lg], exc[i - lg]); /*Q3*/
        }
        L_tmp = L_min(L_tmp, L_tmp_max);
    }

    tmp_e = norm_l(L_tmp);
    L_tmp = L_shl(L_tmp, tmp_e);
    tmp_e = sub(sub(31, shl_pos(Q_exc, 1)), tmp_e); /*L_tmp is Q31, now*/
    tmp16 = BASOP_Util_Divide3216_Scale(L_tmp /*Q31,norm,tmp_e*/, lg /*Q15,15*/, &tmp2_e) /*Q15,tmp2_e+tmp_e-15*/;
    tmp_e = sub(add(tmp2_e, tmp_e), 15);

    IF (tmp16 != 0)
    {
        tmp16     = Sqrt16(tmp16, &tmp_e); /*Q15,norm,tmp_e*/
        tmp_e     = L_min(tmp_e, 15);
        *lp_gainc = L_shl_pos(L_deposit_l(tmp16), add(tmp_e, 1)); /*15Q16*/
        move32();
    }
    ELSE
    {
        *lp_gainc = 0;
    }

    Dyn_Mem_Deluxe_Out();
}

/*
 * TDC_calcGainp
 *
 * Parameters:
 *   x      i: input signal
 *   y      i: shifted input signal
 *   lg     i: vector length
 *
 * Function:
 *   Gain calculation
 *
 * Returns:
 *   gain (15Q16)
 */
static Word32 TDC_calcGainp(Word16 x[], Word16 y[], Word16 lg)
{
    Dyn_Mem_Deluxe_In(
        Word32  tcorr, tener, Lgain, L_tmp1, L_tmp2;
        Word16  m_corr, m_ener, negative, Q_corr, Q_ener;
        Counter i;
    );

    negative = 0; move16();

    L_tmp1 = L_deposit_l(0);
    L_tmp2 = L_deposit_l(0);
    FOR (i = 0; i < lg; i += 2)
    {
        L_tmp1 = L_mac0(L_tmp1, x[i], y[i]);
        L_tmp2 = L_mac0(L_tmp2, x[i + 1], y[i + 1]);
    }
    tcorr  = L_add(L_shr_pos(L_tmp1, 1), L_shr_pos(L_tmp2, 1));
    Q_corr = norm_l(tcorr);
    tcorr  = L_shl(tcorr, Q_corr);
    Q_corr = sub(2, Q_corr);

    L_tmp1 = L_deposit_l(0);
    L_tmp2 = L_deposit_l(0);
    FOR (i = 0; i < lg; i += 2)
    {
        L_tmp1 = L_mac0(L_tmp1, y[i], y[i]);
        L_tmp2 = L_mac0(L_tmp2, y[i + 1], y[i + 1]);
    }
    tener  = L_add(L_shr_pos(L_tmp1, 1), L_shr_pos(L_tmp2, 1));
    Q_ener = norm_l(tener);
    tener  = L_shl(tener, Q_ener);
    Q_ener = sub(2, Q_ener);

    tener = L_max(tener, 1);

    if (tcorr <= 0)
    {
        negative = 1; move16();
    }
    tcorr = L_abs_sat(tcorr);

    m_corr = extract_h(tcorr);

    m_ener = extract_h(tener);

    IF (sub(m_corr, m_ener) > 0)
    {
        m_corr = shr_pos(m_corr, 1);
        Q_corr = add(Q_corr, 1);
    }
    if (m_ener == 0)
    {
        move16();
        m_corr = 0x7FFF;
    }
    if (m_ener != 0)
    {
        m_corr = div_s(m_corr, m_ener);
    }

    Q_corr = sub(Q_corr, Q_ener);

    Lgain = L_shl(L_deposit_l(m_corr), add(Q_corr, 1));

    if (negative != 0)
    {
        Lgain = L_negate(Lgain);
    }

    Dyn_Mem_Deluxe_Out();

    return Lgain;
}

/*
 * TDC_LPC_synthesis_fx
 *
 * Parameters:
 *   sh          i  : scaling to apply for a[0]                 Q0
 *   a[]         i  : LP filter coefficients                    Qx
 *   x[]         i  : input signal                              Qx
 *   y[]         o  : output signal                             Qx-s
 *   lg          i  : size of filtering                         Q0
 *   m           i  : order of LP filter                        Q0
 *
 * Function:
 *    Apply LP filtering to obtain synthesis signal.
 *    Memory size is always m.
 *
 * Returns:
 *    void
 */
static void TDC_LPC_synthesis_fx(const Word16 sh, const Word16 a[], const Word16 x[], Word16 y[], const Word16 lg,
                                 const Word16 m)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  a0;
        Word16  q;
        Word32(*syn_kern)(Word32 L_tmp, const Word16 a[], const Word16 y[]
    ););

    ASSERT(m == 16 || m == 8);

    if (sub(m, 16) == 0)
    {
        syn_kern = syn_kern_16;
    }
    if (sub(m, 8) == 0)
    {
        syn_kern = syn_kern_8;
    }

    q        = add(norm_s(a[0]), 1);
    a0       = shr(a[0], sh);

    FOR (i = 0; i < lg; i++)
    {
        y[i] = round_fx_sat(L_shl_sat(syn_kern(L_mult(a0, x[i]), a, &y[i]), q)); move16();
    }

    Dyn_Mem_Deluxe_Out();
}

/* TDC_LPC_residu_fx
 *
 * Parameters:
 *    a           I: LP filter coefficients (Qx)
 *    x           I: input signal
 *    y           O: output signal
 *    lg          I: size of filtering
 *    m           I: lpc order
 *
 * Function:
 *    Apply inverse filtering to obtain LP residual signal.
 *
 * Returns:
 *    void
 */
static void TDC_LPC_residu_fx(const Word16 *a, Word16 *x, Word16 *y, Word16 lg, Word16 m)
{
    Dyn_Mem_Deluxe_In(
        Word16  a_exp;
        Word32  s;
        Counter i;
    );

    ASSERT(m == 16 || m == 8);

    a_exp = add(norm_s(a[0]), 1);
    a_exp = add(a_exp, 1);

    IF (sub(m, 16) == 0)
    {
        FOR (i = 0; i < lg; i++)
        {
            s = L_mult(x[i], a[0]);
            s = L_mac_sat(s, x[i - 1], a[1]);
            s = L_mac_sat(s, x[i - 2], a[2]);
            s = L_mac_sat(s, x[i - 3], a[3]);
            s = L_mac_sat(s, x[i - 4], a[4]);
            s = L_mac_sat(s, x[i - 5], a[5]);
            s = L_mac_sat(s, x[i - 6], a[6]);
            s = L_mac_sat(s, x[i - 7], a[7]);
            s = L_mac_sat(s, x[i - 8], a[8]);
            s = L_mac_sat(s, x[i - 9], a[9]);
            s = L_mac_sat(s, x[i - 10], a[10]);
            s = L_mac_sat(s, x[i - 11], a[11]);
            s = L_mac_sat(s, x[i - 12], a[12]);
            s = L_mac_sat(s, x[i - 13], a[13]);
            s = L_mac_sat(s, x[i - 14], a[14]);
            s = L_mac_sat(s, x[i - 15], a[15]);
            s = L_mac_sat(s, x[i - 16], a[16]);

            s    = L_shl_sat(s, a_exp);
            y[i] = round_fx_sat(s); move16();
        }
    }
    IF (sub(m, 8) == 0)
    {
        FOR (i = 0; i < lg; i++)
        {
            s = L_mult(x[i], a[0]);
            s = L_mac_sat(s, x[i - 1], a[1]);
            s = L_mac_sat(s, x[i - 2], a[2]);
            s = L_mac_sat(s, x[i - 3], a[3]);
            s = L_mac_sat(s, x[i - 4], a[4]);
            s = L_mac_sat(s, x[i - 5], a[5]);
            s = L_mac_sat(s, x[i - 6], a[6]);
            s = L_mac_sat(s, x[i - 7], a[7]);
            s = L_mac_sat(s, x[i - 8], a[8]);

            s    = L_shl_sat(s, a_exp);
            y[i] = round_fx_sat(s); move16();
        }
    }

    Dyn_Mem_Deluxe_Out();
}

/* TDC_random_fx
 *
 * Parameters:
 *    seed        i/o: seed for random number
 *    lg          i  : vector length
 *    y           o  : output values
 *
 * Function:
 *    Uniform distributed random generator.
 *
 * Returns:
 *    random number
 */
static void TDC_random_fx(Word16 *seed, Word16 lg, Word16 *y)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    FOR (i = 0; i < lg; i++)
    {
        *seed = extract_l(L_mac0(16831L, *seed, 12821));
        *y++  = *seed; move16();
    }

    Dyn_Mem_Deluxe_Out();
}

/*
 * TDC_preemph
 *
 * Parameters:
 *    x              i/o: signal             Qx
 *    fac            i:   preemphasis factor Q15
 *    lg             i:   vector length
 *
 * Function:
 *    Filtering through 1 - fac z^-1
 *
 * Returns:
 *    Q-factor
 */
static Word16 TDC_preemph(Word16 *x, const Word16 fac, const Word16 lg)
{
    Dyn_Mem_Deluxe_In(
        Word16  fac_sh, Q_max_value, Q_out;
        Word32  max_val;
        Counter i;
    );

    fac_sh      = shr(fac, 3);
    Q_max_value = 4096; move16();
    Q_out       = 12;   move16();

    max_val = 0; move32();
    FOR (i = sub(lg, 1); i >= 0; i--)
    {
        max_val = L_max(L_abs(L_msu(L_mult(x[i], Q_max_value), x[i - 1], fac_sh)), max_val);
    }

    IF (extract_h(max_val) != 0)
    {
        Q_out = s_min(s_max(sub(norm_s(extract_h(max_val)), 3), 0), 12);
    }

    FOR (i = sub(lg, 1); i >= 0; i--)
    {
        x[i] = round_fx(L_shl(L_msu(L_mult(x[i], Q_max_value), x[i - 1], fac_sh), Q_out)); move16();
    }

    Dyn_Mem_Deluxe_Out();

    return sub(Q_out, 2);
}

/*
 * TDC_deemph_fx
 *
 * Parameters:
 *    x              i: input signal        Qx
 *    y              o: output signal       Qx
 *    fac            i: deemphasis factor   Q15
 *    lg             i: size of filtering   Q0
 *    mem            i: memory (x[-1])
 *
 * Function:
 *    Filtering through 1/(1-fac z^-1)
 *
 * Returns:
 *    void
 */
static void TDC_deemph_fx(const Word16 *x, Word16 *y, const Word16 fac, const Word16 lg, const Word16 mem)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );

    y[0] = round_fx_sat(L_mac_sat(L_deposit_h(x[0]), mem, fac)); move16();
    FOR (i = 1; i < lg; i++)
    {
        y[i] = round_fx_sat(L_mac_sat(L_deposit_h(x[i]), y[i - 1], fac)); move16();
    }

    Dyn_Mem_Deluxe_Out();
}

/*
 * TDC_normalize_energy_fx
 *
 * Parameters:
 *   gain          o: gain
 *   gain_exp      o: exponent of gain
 *   x             i: input signal
 *   lg            i: length of input signal
 *
 * Function:
 *    Normalizes the energy.
 *
 * Returns:
 *    void
 */
static void TDC_normalize_energy_fx(Word16 *gain, Word16 *gain_exp, const Word16 *x, const Word16 lg)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  c;
        Word16  e;
        Word16  e1;
        Word16  e2;
        Word32  tmp;
        Word16  tmp16;
    );

    tmp = 0; move32();
    FOR (i = 0; i < lg; i++)
    {
        tmp16 = mult_r(x[i], 2048);
        tmp   = L_mac(tmp, tmp16, tmp16);
    }

    e     = norm_l(tmp);
    tmp   = L_shl_pos(tmp, e);
    e1    = sub(sub(30, e), -8); move16();
    tmp16 = BASOP_Util_Divide3216_Scale(tmp, lg, &e2);

    e = 0; move16();
    if (tmp16 != 0)
    {
        e = sub(add(e1, e2), 15);
    }

    c = 0x0148; /* 0.01 */
    move16();
    IF (e > 0)
    {
        c = shr(c, s_min(e, 15));
    }
    ELSE
    {
        tmp16 = shr(tmp16, s_min(negate(e), 15));
        e     = 0; move16();
    }

    e2 = 2; move16();
    if (s_and(e, 1))
    {
        e2 = 1; move16();
    }

    tmp16 = add(shr_pos(tmp16, e2), shr_pos(c, e2));
    e     = add(e, e2);

    tmp16 = Sqrt16(tmp16, &e);

    *gain     = BASOP_Util_Divide1616_Scale((Word16)0x7FFF, tmp16, &e1); move16();
    *gain_exp = sub(e1, e);                                              move16();

    Dyn_Mem_Deluxe_Out();
}

