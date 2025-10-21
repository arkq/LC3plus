/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef CR9_C_ADD_1p25MS
static Word32 limitShaping (Word32* d3_fx) {
    Word32 score = L_shl_pos(L_sub(d3_fx[0], d3_fx[1]), 1);
    score = L_add(score, L_sub(d3_fx[0], d3_fx[2]));
    score = L_shr_pos(score, 1);
    
    score = L_max(L_min(score, 536870912), 335544320); /* 5*2^26 and 8*2^26 */
    
    score = L_sub(536870912, score); /* 8 - score */
        
    score = L_shr_pos(Mpy_32_16(score, 7646), 0); /* 7646 = 0.7/3 * 2^15 */

    score = L_add(score, 20132660); /* +0.3*2^26 */
    
    score = L_shl_sat(score, 5); /* Bring score into Q31 Format for Mpy32_16() */
    /* Alex: Changed from L_shl_pos to L_shl_sat */
    
    return score;
}
#endif


#ifdef FIX_SNS_BASOP_MEAN64_CALC
/*improved precision in  mean64 accumulation will increase min SNR by 15 dB */
/* output is max energy band location  [0 ... 63]  */
static Word16  sns_compute_mean64_ip(Word32 *en_fx_m, Word16 * en_fx_exp, Word32 *L_mean64_fx, Word16* mean64_fx_exp, Word32* L_en_upd_fx, Word16* en_upd_fx_exp )
{
    Word32 i,L_tmp;
    Word16 tmp, max_exp;

    /* maximize all exponents and mantissas  , except for zero mantissas  */
    max_exp = -128; move16();
    FOR (i = 0; i < 64; i++) 
    {  
         tmp = norm_l(en_fx_m[i]);
        
         L_en_upd_fx[i] = L_shl_pos(en_fx_m[i], tmp);     /* factor 1/64  __NOT__  applied here   */
         en_upd_fx_exp[i] = sub(en_fx_exp[i], tmp);

         if (L_en_upd_fx[i] == 0) 
         {
             en_upd_fx_exp[i] = 0; move16();/* all zero mantissa -->  correct   max_exp */
         }
 
         max_exp = s_max(max_exp, en_upd_fx_exp[i]);
    
    }

    L_tmp = L_deposit_l(0); 

    /* sum up at max_exp + 6 for 1/64  */ 
    Word16  glob_shift = add(max_exp, 6);  /* add div by 64 , as margin */
  
    glob_shift = s_max(glob_shift, -31);   

    for (i = 0; i < 64; i++)
    {   
         tmp = s_min(31, sub(glob_shift, en_upd_fx_exp[i])); /* right shift including 1/64  can become larger than 31  */
         L_tmp = L_add(L_tmp, L_shr(L_en_upd_fx[i], tmp));     
    }

    //double mean64_m = round(mean64 * pow(2.0, 31.0-(double)max_exp)); /* dbg */
    //double mean64_new = (double)L_tmp * pow(2.0, (double)max_exp - 31.0);

    *L_mean64_fx = L_tmp;
    *mean64_fx_exp = max_exp;

    return -1;

}
#endif 


void processSnsComputeScf_fx(Word32 *d2_fx, Word16 d2_fx_exp, Word16 fs_idx, Word16 n_bands, Word16 *scf,
                             Word16 scf_smoothing_enabled, Word16 attdec_damping_factor, Word8 *scratchBuffer, Word16 sns_damping
#ifdef CR9_C_ADD_1p25MS
                             , LC3PLUS_FrameDuration frame_dms, Word16 norm_corr, Word16 *LT_normcorr
#endif
                             )
{
    Dyn_Mem_Deluxe_In(
        Word16  i, s, s2, nf, tmp;
        Word32  L_mean, L_tmp;
        Word32 *d3_fx;
        Word16 *d3_fx_exp;
        Word16 *d4_fx;
        Word16 *scf_smooth;
    );

    UNUSED(attdec_damping_factor);

    d3_fx     = scratchAlign(scratchBuffer, 0);                         /* Size = 4 * MAX_BANDS_NUMBER = 256 bytes */
    d3_fx_exp = scratchAlign(d3_fx, sizeof(*d3_fx) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    d4_fx = scratchAlign(d3_fx_exp, sizeof(*d3_fx_exp) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128bytes */
    scf_smooth = scratchAlign(d4_fx, sizeof(*d4_fx) * MAX_BANDS_NUMBER);    /* Size = 2 * 16 */
    
    const Word16 *preemp;
    const Word16 *preemp_e;
#ifdef CR9_C_ADD_1p25MS
    Word16 *pre_em_buf;
    Word16 *pre_em_buf_e;

    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
        pre_em_buf  = scratchAlign(scf_smooth, sizeof(*pre_em_buf) * MAX_BANDS_NUMBER);
        pre_em_buf_e  = scratchAlign(pre_em_buf, sizeof(*pre_em_buf_e) * MAX_BANDS_NUMBER);
    }
    
    Word32 limiterGain = 0x7FFFFFFF; /* 1.0f */
    Word16 ncorr_fac = 0;

#ifdef FIX_SNS_BASOP_MEAN64_CALC
    Word32 L_d5_fx[64], L_mean_ip = 1;
    Word16 d5_fx_exp[64], mean_ip_exp = 0;
#endif 
#endif

/* Smoothing and Pre-emphasis */
    IF (sub(n_bands, 32) < 0)
    {
        L_tmp = sub(32, n_bands);
        FOR (i = sub(n_bands, 1); i >= L_tmp; i--)
        {
            d2_fx[(i + L_tmp) * 2 + 1] = d2_fx[i]; move32();
            d2_fx[(i + L_tmp) * 2 + 0] = d2_fx[i]; move32();
        }
        FOR (i = sub(L_tmp, 1); i >= 0; i--)
        {
            d2_fx[i * 4 + 3] = d2_fx[i]; move32();
            d2_fx[i * 4 + 2] = d2_fx[i]; move32();
            d2_fx[i * 4 + 1] = d2_fx[i]; move32();
            d2_fx[i * 4 + 0] = d2_fx[i]; move32();
        }
        n_bands = 64; move16();
    }
    ELSE
    IF (sub(n_bands, 64) < 0)
    {
        L_tmp = sub(64, n_bands);
        FOR (i = sub(n_bands, 1); i >= L_tmp; i--)
        {
            d2_fx[i + L_tmp] = d2_fx[i]; move32();
        }
        FOR (i = sub(L_tmp, 1); i >= 0; i--)
        {
            d2_fx[i * 2 + 1] = d2_fx[i]; move32();
            d2_fx[i * 2 + 0] = d2_fx[i]; move32();
        }
        n_bands = 64; move16();
    }
    
#ifdef CR9_C_ADD_1p25MS
    /* calc long termn normcorr */
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS && lpc_pre_adapt_emphasis[fs_idx] != NULL) {
#ifdef FIX_BASOP_LT_NORMCORR_AR1  
        #ifdef ER_DEBUG  
                Word16 LT_normcorr_mem = *LT_normcorr;
        #endif 
                L_tmp = L_mult0(1, norm_corr);            /* 1/8 */
                L_tmp = L_mac0(L_tmp, 7, *LT_normcorr);   /* 7/8 */
                *LT_normcorr = round_fx(L_shl_pos(L_tmp, 16 - 3) );
                assert(*LT_normcorr >= 0);
        
                /* ncorr_fac = max(*LT_normcorr - 0.8, 0) * 5 */        
               /* improved  precision of n_corr_fac  esp. close to the 0.8 border    */
                L_tmp = L_mult0(*LT_normcorr, 5 * (1L << (15 - 3)));      /* work in the 5 * upscaled domain from start */
                L_tmp = L_sub(L_tmp, 4 * (1L << ((15 - 3) + 15)));        /* sub by   0.8   ,  limit now exactly 0.8      */
                L_tmp = L_max(L_tmp, 0L);
                ncorr_fac = round_fx(L_shl(L_tmp, 4));               /* extract in Q15 */ 
              
        #ifdef ER_DEBUG  
                LT_normcorr_mem = add(mult(LT_normcorr_mem, 28672), mult(norm_corr, 4096));
                Word16 ncorr_fac_orig = i_mult(s_max(sub(LT_normcorr_mem, 26214), 0), 5); 
                UNUSED(ncorr_fac_orig);
        #endif 
#else        
        /* LT_normcorr * (1-0.125) + normcorr * 0.125 */
        *LT_normcorr = add(mult(*LT_normcorr,28672), mult(norm_corr,4096));
        /* ncorr_fac = max(LT_normcorr - 0.8, 0) * 5 */
        ncorr_fac = i_mult(s_max(sub(*LT_normcorr, 26214), 0), 5);
#endif        

        if (ncorr_fac == 0)
        {
            preemp = lpc_pre_emphasis[fs_idx];
            preemp_e = lpc_pre_emphasis_e[fs_idx];
        } else {
            FOR (i = 0; i < n_bands; i++) {
                assert (lpc_pre_adapt_emphasis_e[fs_idx][i]>= lpc_pre_emphasis_e[fs_idx][i]);
                assert((((Word16)lpc_pre_adapt_emphasis_e[fs_idx][i]) + 1) <= 255);  /*check  Word8 storage */
                assert((((Word16)lpc_pre_adapt_emphasis_e[fs_idx][i]) + 1) >= -256); /*check Word8 storage */
#ifdef FIX_BASOP_PREEMPH_CALC
            /*  FLT      sns_preemph_new[i]    = (sns_preemph[i]        +  fac   * sns_preemph_adapt[i]); */
            /*  BASOP  MantNew[i]*2^expNew[i]  = (Mant16A[i]*2^expA[i]  +  facQ15 * Mant16B^expB[i]); */
                    /*determine downscaling factor for the base part */
                    Word16 tmp_shift        = sub(lpc_pre_adapt_emphasis_e[fs_idx][i], lpc_pre_emphasis_e[fs_idx][i]);
                    Word16 shift_factor     = shr_pos(-32768,  tmp_shift);   /* mantissa downshift  impl.  as multiplication on the negative side */

                    /*Add scaled mantissas in the same Q  but 1 bit rightshifted(for 1 bit WC addition margin) */
                    L_tmp = L_mult0(lpc_pre_adapt_emphasis[fs_idx][i], ncorr_fac); /*  Qx*Q15   --> Qx,   1 additional bit margin created with Mult0  */
                    L_tmp = L_msu0(L_tmp, lpc_pre_emphasis[fs_idx][i], shift_factor );  /* shift_dactor is negative,  Msu0  will lead to addition */
                    
                    tmp = norm_l(L_tmp);                                 
                    pre_em_buf[i] = round_fx_sat(L_shl_pos(L_tmp, tmp)); /* maximize mantissa and apply round */ /* sat needed due to warning due to RND additio of 0.5   */

                    tmp = sub(tmp, add(tmp_shift, 1)); ;      
                    pre_em_buf_e[i] = sub(lpc_pre_emphasis_e[fs_idx][i], tmp); /* adjust final exponent due to mantissa upshift and addition margin */ 

                    /* sum  ~10 ops,  but increased precision, and a maintained maximized mantissa  */

#else 
                pre_em_buf_e[i] = lpc_pre_adapt_emphasis_e[fs_idx][i] + 1;   
                Word16 tmp_s = lpc_pre_adapt_emphasis_e[fs_idx][i] - lpc_pre_emphasis_e[fs_idx][i];
                Word16 tmp_v = mult(lpc_pre_adapt_emphasis[fs_idx][i], ncorr_fac);
                pre_em_buf[i] = add(shr(lpc_pre_emphasis[fs_idx][i], (tmp_s+1)), shr(tmp_v,1));
#endif
            }
            preemp = pre_em_buf;
            preemp_e = pre_em_buf_e;
        }
    } else {
        preemp = lpc_pre_emphasis[fs_idx];
        preemp_e = lpc_pre_emphasis_e[fs_idx];
    }
#else
    preemp = lpc_pre_emphasis[fs_idx];
    preemp_e = lpc_pre_emphasis_e[fs_idx];
#endif
    
    L_tmp        = L_add(Mpy_32_16(d2_fx[0], 24576), L_shr_pos(d2_fx[1], 2));
    d3_fx[0]     = Mpy_32_16(L_tmp, preemp[0]); move32();
    d3_fx_exp[0] = add(d2_fx_exp, preemp_e[0]); move16();
    FOR (i = 1; i < n_bands - 1; i++)
    {
        L_tmp        = L_add(L_shr_pos(d2_fx[i], 1), L_add(L_shr_pos(d2_fx[i - 1], 2), L_shr_pos(d2_fx[i + 1], 2)));
        d3_fx[i]     = Mpy_32_16(L_tmp, preemp[i]); move32();
        d3_fx_exp[i] = add(d2_fx_exp, preemp_e[i]); move16();
    }
    L_tmp                  = L_add(Mpy_32_16(d2_fx[n_bands - 1], 24576), L_shr_pos(d2_fx[n_bands - 2], 2));
    d3_fx[n_bands - 1]     = Mpy_32_16(L_tmp, preemp[n_bands - 1]); move32();
    d3_fx_exp[n_bands - 1] = add(d2_fx_exp, preemp_e[n_bands - 1]); move16();

    /* Mean */
    s  = d3_fx_exp[MAX_BANDS_NUMBER - 1];
    s2 = add(s, 6);

    L_mean = L_shr(d3_fx[0], sub(s2, d3_fx_exp[0]));
    FOR (i = 1; i < MAX_BANDS_NUMBER; i++)
    {
        L_mean = L_add(L_mean, L_shr(d3_fx[i], sub(s2, d3_fx_exp[i])));
    }

#ifdef FIX_SNS_BASOP_MEAN64_CALC
    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) 
    {
        /* improved precision in  mean64 accumulation will increase minimum SNR vs FLT by 15 dB */
        Word16 max_loc = sns_compute_mean64_ip(&d3_fx[0], &d3_fx_exp[0], &L_mean_ip, &mean_ip_exp, &L_d5_fx[0], &d5_fx_exp[0]);  /*ToDo::inplace update  d3_fx, d3_fx_exp */
        /* high precision accumluation , divides all entries in d_fx by 1/64  */
        UNUSED(max_loc);

#ifdef ER_DEBUG
        //float mean64new_fx_fl = (((double)L_mean_ip))*pow(2.0, (double)mean_ip_exp - 31);
        Word16 d3_max_log2 = BASOP_Util_Log2_16(d3_fx[max_loc], d3_fx_exp[max_loc]);
        double  d3_maxlog2_dbl = (double)d3_max_log2 / 512.0;
        UNUSED(d3_maxlog2_dbl);
#endif 

#ifdef ER_DEBUG
        if (d3_fx[max_loc] != 0)   /*  max_loc is -1 when not debugging */
#endif 
        {   /* maximize the q in the  coming log2 function call  */

            basop_memcpy(d3_fx, L_d5_fx, 64 * sizeof(Word32));          /* use  upscaled mantissa values */

            basop_memcpy(d3_fx_exp, d5_fx_exp, 64 * sizeof(Word16));    /* use  corresponding  exponents values */

#ifdef ER_DEBUG
         /* s,s2  are no longer valid in relation to the previous d3_fx */
            s = d3_fx_exp[max_loc];
            s2 = add(s, 6);  /* div by 64 */

            /* debug:: redo  L_mean calc using new s, s2 and upscaled mantissa  parameters */
            tmp = s_min(31, sub(s2, d3_fx_exp[0]));
            L_mean = L_shr_pos(d3_fx[0], tmp);
            FOR(i = 1; i < MAX_BANDS_NUMBER; i++)
            {
                tmp = s_min(31, sub(s2, d3_fx_exp[i]));
                L_mean = L_add(L_mean, L_shr_pos(d3_fx[i], tmp));
            }
#endif 
        }
    }
#endif 

#ifdef ER_DEBUG
    if ((dbgflag("r") || dbgflag("ru")) && scratch->max_scratch_calculation_only == 0)
    {
        double sns64_fx_fl_all[64];
        double sns64_fx_fl_all_log2[64];  /* find range of smoothed input signal */  

        double sns64_fx_fl_opt_all[64];
        double sns64_fx_fl_opt_all_log2[64];  /* find range of smoothed input signal */

        Word32 d3_opt_fx[64];
        Word16 d3_opt_fx_exp[64];
        float fl_mean64, fl_mean64_m, fl_mean64_e;
        float L_mean_fx_fl;
        float fl_mean64_recalc;

        UNUSED(sns64_fx_fl_opt_all);
        UNUSED(sns64_fx_fl_opt_all_log2);
        UNUSED(sns64_fx_fl_all_log2);
        UNUSED(sns64_fx_fl_all);
       

        UNUSED(fl_mean64_m);
        UNUSED(fl_mean64_e);

        dbgread(&fl_mean64, sizeof(fl_mean64), 1, "..//..//..//fl_x_mean64.float.1.dat");

#ifdef FIX_SNS_BASOP_MEAN64_CALC
        float mean64new_fx_fl = (((double)L_mean_ip))*pow(2.0, (double)mean_ip_exp - 31);

        snr_diff(&mean64new_fx_fl, &fl_mean64, 1, 0, "SNRc4b_x_mean64ip_FXvsFL");

#endif 

        fl_mean64_e = s;
        fl_mean64_m = fl_mean64 * pow(2.0, 31.0 -(double)fl_mean64_e );


        fl_mean64_recalc = 0.0;
        for (int i = 0; i < 64; i++)
        {
            sns64_fx_fl_all[i] = (double)d3_fx[i] * pow(2.0, (double)d3_fx_exp[i] - 31.0);
            sns64_fx_fl_all_log2[i] = log2(sns64_fx_fl_all[i]);
            fl_mean64_recalc += sns64_fx_fl_all[i];

            int tmp_up = norm_l(d3_fx[i]);
            d3_opt_fx[i] = L_shl_pos(d3_fx[i], tmp_up);
            d3_opt_fx_exp[i] = sub(d3_fx_exp[i], tmp_up);

            sns64_fx_fl_opt_all[i] = (double)d3_opt_fx[i] * pow(2.0, (double)d3_opt_fx_exp[i] - 31.0);
            sns64_fx_fl_opt_all_log2[i] = log2(sns64_fx_fl_all[i]);
        }



        fl_mean64_recalc /= 64.0;  /* downscaling by 1/64 already applied */

        L_mean_fx_fl = (((double)L_mean))*pow(2.0, (double)s - 31.0);
        snr_diff(&L_mean_fx_fl, &fl_mean64, 1, 0, "SNRc4b_x_mean64_FXvsFL");

        snr_diff(&fl_mean64_recalc, &fl_mean64, 1, 0, "SNRc4b_x_mean64recalc_FXvsFL"); /* 15 dB better !! */


        if (dbgflag("ru"))
        {
#ifdef FIX_SNS_BASOP_MEAN64_CALC
            L_mean = round(fl_mean64*pow(2.0, 31 - mean_ip_exp));
            assert(fabs(fl_mean64*pow(2.0, 31 - mean_ip_exp)) <= (double)INT32_MAX);
#else 
            L_mean = round(fl_mean64*pow(2.0, 31 - s));
            assert(fabs(fl_mean64*pow(2.0, 31 - s)) <= (double)INT32_MAX);
#endif 

        }
    }
#endif 

#ifdef FIX_SNS_BASOP_NF_APPL 
    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
#ifdef ER_DEBUG
        Word16 nf_post = BASOP_Util_Log2_16(L_mean_ip, mean_ip_exp);
        nf_post = sub(s_max(nf_post, -25965), 6803);       /*  log2(-25965)- log2(6803) -->  -50.7  -13.28  =>  -64.00 in Q9  */
        /*if nf > -50.7  nf is applied by subtracting   2^-13.28   9996  */
#endif 

    /* calulate the  nf level in W32 before going to the less exact log2() W16 domain */
        Word16 tmp_nf = norm_l(L_mean_ip);
        Word32 L_nf_fx = L_shl_pos(L_mean_ip, tmp_nf);                   /* shift up  1/64 to maximize mantissa */

        L_nf_fx = Mpy_32_16(L_nf_fx, 26844);                   /* 26844=round(2^13/10000)  */
        Word16 nf_fx_exp = sub(mean_ip_exp, add(13, tmp_nf));  /*  division by 10000.0  as  L_mant = L_mant*(2^13/10000), nf_exp=mean_exp-13 */

        nf = BASOP_Util_Log2_16(L_nf_fx, nf_fx_exp);    /* log2_16 function ouput is in Q9 , range is values  [-64.xxx to 63.xxxx]  */

#ifdef ER_DEBUG
        /*  all zero input will get nf of -32768 */
        if (L_mean_ip == 0)
        {
            assert(nf == -32768);
        }
#endif 

    }
    ELSE
    {  /* 2.5ms ...  10ms  (nf  is at  approximated at  1/9998 below mean) */
        nf = BASOP_Util_Log2_16(L_mean, s);
        nf = sub(s_max(nf, -25965), 6803);
    }
#else 
    nf = BASOP_Util_Log2_16(L_mean, s);
    nf = sub(s_max(nf, -25965), 6803);
#endif 


/* Log-domain */
    FOR (i = 0; i < MAX_BANDS_NUMBER; i++)
    {
        d4_fx[i] = s_max(nf, BASOP_Util_Log2_16(d3_fx[i], d3_fx_exp[i])); move16();
    }

    /* Downsampling */
    L_tmp    = L_mult(d4_fx[0], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[1], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[2], 8192);
    L_tmp    = L_mac(L_tmp, d4_fx[3], 5461);
    d3_fx[0] = L_mac(L_tmp, d4_fx[4], 2731); move32();
    FOR (i = 1; i < M - 1; i++)
    {
        L_tmp    = L_mult(d4_fx[i * 4 - 1], 2731);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 0], 5461);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 1], 8192);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 2], 8192);
        L_tmp    = L_mac(L_tmp, d4_fx[i * 4 + 3], 5461);
        d3_fx[i] = L_mac(L_tmp, d4_fx[i * 4 + 4], 2731); move32();
    }
    L_tmp        = L_mult(d4_fx[59], 2731);
    L_tmp        = L_mac(L_tmp, d4_fx[60], 5461);
    L_tmp        = L_mac(L_tmp, d4_fx[61], 8192);
    L_tmp        = L_mac(L_tmp, d4_fx[62], 8192);
    d3_fx[M - 1] = L_mac(L_tmp, d4_fx[63], 8192); move32();
    
#if defined (CR9_C_ADD_1p25MS)
    /* limit shaping */
    IF (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
        limiterGain = limitShaping(d3_fx);
    }
#endif

/* Remove mean and scaling */
    L_mean = L_shr_pos(d3_fx[0], 4);
    FOR (i = 1; i < M; i++)
    {
        L_mean = L_add(L_mean, L_shr_pos(d3_fx[i], 4));
    }

#ifdef  FIX_SNS_BASOP_MEAN16_APPLY
    Word16 limGainDamping = 0;
    Word16 up_shift = 0, dn_shift = 0;

    IF(frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        up_shift = 31;  move16();

        FOR(i = 0; i < M; i++)
        {
            L_tmp = L_sub(d3_fx[i], L_mean);
            up_shift = s_min(norm_l(L_tmp), up_shift);
            d3_fx[i + M] = L_tmp; move32();  /* store mean removed variable */
        }
        dn_shift = s_max(0, sub(up_shift, 1));   /* conditional downshift from q10 to Q11, to maintain precision  */

        limGainDamping = round_fx(Mpy_32_16(limiterGain, sns_damping));
    }
#endif 

    FOR (i = 0; i < M; i++)
    {
#if defined (CR9_C_ADD_1p25MS)
        if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
        {
#ifdef  FIX_SNS_BASOP_MEAN16_APPLY
            /* maximize precision */
            L_tmp = Mpy_32_16(L_shl_pos(d3_fx[i + M], up_shift), limGainDamping );
            scf[i] = round_fx(L_shr_pos(L_tmp, dn_shift));
#else 
            scf[i] = mult_r(extract_h(Mpy_32_16(limiterGain, sns_damping)), round_fx(L_shl_pos(L_sub(d3_fx[i], L_mean), 1)));
#endif
            move16();
        } else {
            scf[i] = mult_r(sns_damping, round_fx(L_shl_pos(L_sub(d3_fx[i], L_mean), 1)));
            move16();
        }
#else
        scf[i] = mult_r(sns_damping, round_fx(L_shl_pos(L_sub(d3_fx[i], L_mean), 1)));
        move16();
#endif
    }

    /* scale factor smoothing */
    #ifdef CR9_C_ADD_1p25MS_LRSNS
    test(); test();
    IF( (scf_smoothing_enabled != 0 ) ||  (sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)  )
    #else 
    IF (scf_smoothing_enabled)
    #endif 
    {
        #ifdef CR9_C_ADD_1p25MS_LRSNS
        IF(sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)
        {

#           define A0Q15 6144                /* 0.1875 */
#           define A1Q15 (32768 - 2 * A0Q15) /* 0.675 */
#           define A2Q15 A0Q15               /* 0.1875 */ 

            L_tmp = L_mult0(scf[0], A0Q15 >> 1);
            L_tmp = L_mac0(L_tmp, scf[1], A0Q15 >> 1);/*incl a mult by 0.5)*/  /* a virtual scf[-1]*.1875 created */
            L_tmp = L_mac0(L_tmp, scf[0], A1Q15);     /*scf[0] * .675 */
            L_tmp = L_mac0(L_tmp, scf[1], A2Q15);      /*scf[1] * .1875 */
            L_tmp = L_add(L_tmp, L_tmp); /* added to realize one upshift due to margin created by mac0 */

            scf_smooth[0] = round_fx(L_tmp);
            L_mean = L_deposit_l(scf_smooth[0]);

            FOR(i = 1; i < M - 1; i++)
            {  /* normal unity gain three tap FIR filter */
                L_tmp = L_mult0(scf[i - 1], A0Q15);  /*Q(11+16+1)  Q26 */
                L_tmp = L_mac0(L_tmp, scf[i], A1Q15);
                L_tmp = L_mac0(L_tmp, scf[i + 1], A2Q15);
                L_tmp = L_add(L_tmp, L_tmp); /* added to realize one upshift due to margin created by mac0 */
                scf_smooth[i] = round_fx(L_tmp);   /* back to Q11 */
                L_mean = L_add(L_mean, L_deposit_l(scf_smooth[i])); /* mean  accumulated on the low side */
            }

            /* final coeff */
            L_tmp = L_mult0(scf[M - 2], A2Q15);       /* .1875*/
            L_tmp = L_mac0(L_tmp, scf[M - 1], A1Q15);       /* .675 */
            L_tmp = L_mac0(L_tmp, scf[M - 2], A0Q15 >> 1); /* now the  final virtual scf[M]*.1875   */
            L_tmp = L_mac0(L_tmp, scf[M - 1], A0Q15 >> 1);
            L_tmp = L_add(L_tmp, L_tmp);            /* added to realize one upshift due to margin created by mac0 */
            scf_smooth[M - 1] = round_fx(L_tmp);

            L_mean = L_add(L_mean, L_deposit_l(scf_smooth[M - 1]));
            
            L_mean = L_shr(L_add(L_mean, 1L << (4 - 1)), 4);    /* excplicit rounding in the 1/M =  0.1625   application  */
        }
        ELSE IF(sub(scf_smoothing_enabled, 1) == 0)
        { /* legacy smoothing logic */
#endif 
        scf_smooth[0] = L_shr(L_mult0(L_add(L_add(scf[0], scf[1]), scf[2]), 10923), 15);
        L_mean        = scf_smooth[0]; move16();
        scf_smooth[1] = L_shr(L_add(L_add(L_add(scf[0], scf[1]), scf[2]), scf[3]), 2);
        L_mean        = L_add(L_mean, scf_smooth[1]);
        FOR (i = 2; i < M - 2; i++)
        {
            L_tmp         = L_add(L_add(L_add(L_add(scf[i - 2], scf[i - 1]), scf[i]), scf[i + 1]), scf[i + 2]);
            tmp = norm_s(L_shr(L_abs(L_tmp), 15));
            if (tmp > 0) {
                tmp = sub(16, tmp);
                L_tmp = L_shr(L_tmp, tmp);
            }
            scf_smooth[i] = L_shr(L_mult0(L_tmp, 13107), sub(16, tmp));
            L_mean        = L_add(L_mean, scf_smooth[i]);
        }
        scf_smooth[M - 2] = L_shr(L_add(L_add(L_add(scf[M - 4], scf[M - 3]), scf[M - 2]), scf[M - 1]), 2);
        L_mean            = L_add(L_mean, scf_smooth[M - 2]);
        scf_smooth[M - 1] = L_shr(L_mult0(L_add(L_add(scf[M - 3], scf[M - 2]), scf[M - 1]), 10923), 15);
        L_mean            = L_add(L_mean, scf_smooth[M - 1]);

        L_mean = L_shr(L_mean, 4);

        
        #ifdef CR9_C_ADD_1p25MS_LRSNS
        }
        #endif 
        
        #ifdef CR9_C_ADD_1p25MS_LRSNS
         IF(sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)
        {
            ASSERT(L_mean <= 32767L && L_mean >= -32768L);
            Word16  mean =  extract_l(L_mean);
            FOR(i = 0; i < M; i++)
            {
                scf[i] = sub(scf_smooth[i], mean);  /* only mean subtraction required, no damping/scaling  */
            }
        }
        ELSE
        {
            FOR(i = 0; i < M; i++)
            {
                scf[i] = mult_r(attdec_damping_factor, sub(scf_smooth[i], L_mean));
            }   
        }
#   else 
        FOR (i = 0; i < M; i++)
        {
            scf[i] = mult_r(attdec_damping_factor, sub(scf_smooth[i], L_mean));
        }
        #endif
    }

    Dyn_Mem_Deluxe_Out();
}

