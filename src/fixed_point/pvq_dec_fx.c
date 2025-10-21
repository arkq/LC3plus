/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

Word16 pvq_dec_deidx_fx(                          /* out BER detected 1 , ok==0 */
                        Word16 *      y,          /* o:   decoded vector (non-scaled int)  */
                        const Word16  k_val,      /* i:   number of allocated pulses       */
                        const Word16  dim,        /* i:   Length of vector                 */
                        const Word16  LS_ind,     /* i; lS index              1 bit        */
                        const UWord32 UL_MPVQ_ind /* i; MPVQ  index                        */
)
{
    Dyn_Mem_Deluxe_In(
        Word16      BER_flag;
        UWord32     h_mem[1 + KMAX_FX + 1];
        PvqEntry_fx entry;
    );

    BER_flag = 0; move16();

    /* get_size will likely be called before this function,     as the range decoder needs the size to fetch the index
     */
    entry = get_size_mpvq_calc_offset_fx(dim, k_val, h_mem); /* TBD should be made into tables for N=16,10,6  */

    entry.lead_sign_ind = LS_ind;         move16();
    entry.index         = L_deposit_l(0); /* only  in case dim == 1 */
    IF (sub(dim, 1) != 0)
    {
        entry.index = UL_MPVQ_ind;

        /* safety check in case of bit errors */
        IF (L_sub(entry.index, entry.size) >= 0)
        {
            BER_flag    = 1;            move16();
            entry.index = 0; move16(); /* return something deterministic/valid, and LOW complex  */
        }
    }
    mpvq_deindex_fx(&entry, h_mem, y); /* actual deindexing  */

    Dyn_Mem_Deluxe_Out();
    return BER_flag;
}



#ifdef ENABLE_HR_MODE
void pvq_dec_scale_vec_fx(const Word32 *inQ29, Word16 adjGainQ13, Word32 *outQ)
#else
void pvq_dec_scale_vec_fx(const Word16 *inQ14, Word16 adjGainQ13, Word16 *outQ14)
#endif
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );
    FOR (i = 0; i < M; i++)
    {
#ifdef ENABLE_HR_MODE
        outQ[i] = L_shr(L_add(outQ[i], Mpy_32_16(inQ29[i], adjGainQ13)), 1);
        move16();
#else
        outQ14[i] = add(outQ14[i], mult_r(adjGainQ13, inQ14[i])); move16();
#endif
    }
    Dyn_Mem_Deluxe_Out();
}

#ifdef CR9_C_ADD_1p25MS_LRSNS
void lrsns_pvq_dec_scale_W16vec_fx(
    const Word16 *inQ14,
    Word16 adjGainQ12,
    Word16 *inQ11outQ11
)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    );
    FOR(i = 0; i < M; i++)
    {     /* Q12*Q14 +1  -16  => Q11  */
        inQ11outQ11[i] = add(inQ11outQ11[i], mult_r(adjGainQ12, inQ14[i]));   move16();
        /*   Q11   =           Q11   +   adjGainQ12 *inQ14 +1-16  ,  note  12+14+1-16=11 */
    }
    Dyn_Mem_Deluxe_Out();
}

#ifdef ENABLE_HR_MODE
void lrsns_pvq_dec_scale_W32vec_fx(
    const Word32 *inQ30,
    Word16 adjGainQ12,
    Word32 *inQ27outQ26,
    Word16 *outQ11
)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    Word32 L_tmp;
    );
    FOR(i = 0; i < M; i++)
    {
        L_tmp = L_shr_pos(L_add(inQ27outQ26[i], Mpy_32_16(inQ30[i], adjGainQ12)), 1);
        /*   Q26   =      (  Q27  +        inQ30*adjGainQ12 * +1-16 ) >>1 ,  note    30+12+1-16 =43-16 = 27   */
        inQ27outQ26[i] = L_tmp;      move32();
        L_tmp = L_add(L_tmp, 1 << 14);                    /* manual round  + 0.5*2^15     */
        L_tmp = L_min(16383L << 16, L_tmp);               /* pre_saturate, to not exceed 32767 , for extract_l() below */
        L_tmp = L_max(-(16384L << 16), L_tmp);            /* pre_saturate, to not exceed -32768. for extract_l() below */
        outQ11[i] = extract_h(L_shl_pos(L_tmp, 1));  move16(); /*26-16 +1= 26-15 = 11 -->  Q26 to Q11*/

        /* the effect of this Word32 to Word16 rounding does not result in exactly the same W16Q11 vector  as DISABLE_HR  */
    }
    Dyn_Mem_Deluxe_Out();
}
#endif  /* ENABLE_HR_MODE */

void pvq_fess_dec_en1_normQ30andQ14_fx(
    const Word16 *y /*Q0*/,
    Word16 y_up_bits,
    Word32 L_norm_factor,
    Word16 norm_factorQ,
    Word16 len,
    Word32* L_y_norm,
    Word16 *y_norm)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    Word32 L_tmp;
    Word16 shift_tot;
    );

        BASOP_sub_sub_start("pvq_fess_dec_en1_normQ30andQ14_fx");

    /* lrsns_norm_factorQ_L[N_SCF_SEARCH_SHAPES_ST2_LR] = { 31, 31, 19 + 16 }; */  /* 0_split, 1_full, 2_fixenv */
    /* maxamps = [8, 5, 12 ];  8*2^11 = 2^14 =16384 , 5*2^12 = 20480, 12*2^11= 24576   */
    /* lrsns_y_up_bits[N_SCF_SEARCH_SHAPES_ST2_LR] = { 11, 12, 11 };*/

    shift_tot = sub((30 + 16), add(norm_factorQ, y_up_bits));
    /* 0: lf_split    (30+16) - (31+11) = 4  */
    /* 1: full        (30+16) - (31+12) = 3  */
    /* 2: fixenv      (30+16) - (35+11) = 0   */
    ASSERT(shift_tot >= 0);
    FOR(i = 0; i < len; i++)
    {
        L_tmp = Mpy_32_16_0_0(L_norm_factor, shl_pos(y[i], y_up_bits));       /* Qfactor*Q_y_up    => Qfactor*Qyup + 0  - 16 . e.g. Q31*Q11-16 = Q26  */
        L_y_norm[i] = L_shl_pos(L_tmp, shift_tot);              move32();     /* shift up  the target to Q30 */

        y_norm[i] = round_fx(L_y_norm[i]);  move16(); /* 30-16=>  Q14,      used in enc side gainQ loop, and for  DISABLE_HR synthesis  */
            
    }
    BASOP_sub_sub_end(); /*  pvq_fess_dec_en1_normQ30andQ14_fx*/
    Dyn_Mem_Deluxe_Out();
}

void FESSdeenum_fx(Word16 dim_in_fx,    /* i :  dimension of vec_out   typically (M-1 == 15)  */
    Word16 n_env_fx,     /* i :  number of envelopes    */
    Word16 n_shift_fx,   /* i :  number shifts        */
    Word16 n_signs_fx,   /* i :  number signs          */
    Word16 env_ind_fx,   /*  i:indx */
    Word16 shift_ind_fx, /*  i:indx */
    Word16 sign_ind_fx,  /*  i:indx */
    Word16* vec_out_fx /* o :  FESS  integer  pulse train , with signs applied  */)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    Word16  sign_val_fx;
    );
    assert(n_env_fx >= 1 && n_env_fx <= 4);
    assert(env_ind_fx >= 0 && env_ind_fx < n_env_fx);
    assert(shift_ind_fx >= 0 && shift_ind_fx < n_shift_fx);

    UNUSED(n_env_fx);
    UNUSED(n_shift_fx);

    BASOP_sub_sub_start("FESSdeenum_fx");

    basop_memset(vec_out_fx, 0, sizeof(*vec_out_fx)*dim_in_fx);

    FOR(i = (shift_ind_fx + n_signs_fx - 1); i >= shift_ind_fx; i--)
    {
        /* low numbered coeff  signs are in the msb's */
        /* high  numbered coeff  signs are in the lsb's */
        ASSERT(i < dim_in_fx);

        sign_val_fx = sub(1, shl_pos(s_and(sign_ind_fx, 0x01), 1)); /* 1 - 2*signbit_value */
        sign_ind_fx = shr_pos_pos(sign_ind_fx, 1);

        vec_out_fx[i] = extract_l(L_mult0(sign_val_fx, lrsns_fix_env_fx[env_ind_fx][i]));
    }

    BASOP_sub_sub_end(); /*  */
    Dyn_Mem_Deluxe_Out();
}
#endif 

void pvq_dec_en1_normQ14_fx(/*  Have to be used EXACTLY the same way in both  both encoder and decoder */
#ifdef ENABLE_HR_MODE
                            Word32 *      xq, /* o:   en1 normalized decoded vector (Q14)   */
#else
                            Word16 *      xq, /* o:   en1 normalized decoded vector (Q14)   */
#endif
                            const Word16 *y,  /* i:   decoded vector (non-scaled int)  */
                            const Word16  k_val_max,
                            /* i:   max possible K   in Q0 kO or kA   */ /* OPT:  not BE , use dynamic max  pulse
                                                                            amplitude */
                            const Word16 dim                             /* i:   Length of vector                 */
)
{
#ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  L_tmp;
        Word16  shift_num, shift_tot;
        Word32  isqrtQ31_local; 
        Word16  tmp, exp, exp_shift;
        Word32  L_yy;
    );
#else
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  L_tmp;
        Word16  shift_num, shift_tot;
        Word16  isqrtQ16_local, tmp, exp, exp_shift;
        Word32  L_yy;
    );
#endif

/* energy normalization starts here */
    L_yy = L_mult0(y[0], y[0]);
    FOR (i = 1; i < dim; i++)
    {
        L_yy = L_mac0(L_yy, y[i], y[i]); /* stay in Q0 */ /* OPT: reuse some energies from PVQ linear search */
    }
    /* 16 bit */
    IF (L_sub(L_yy, SQRT_EN_MAX_FX) < 0)
    {
        ASSERT(L_yy > 4);                               /* Q16 isqrt table lookup not valid below  5  */
#ifdef ENABLE_HR_MODE
        isqrtQ31_local = isqrt_Q31tab[L_yy]; move16(); /* 1 cycle */
#else
        isqrtQ16_local = isqrt_Q16tab[L_yy]; move16(); /* 1 cycle */
#endif
    }
    ELSE
    {
        /* about 8-9 cycles */
        exp            = 15; move16(); /* set ISqrt16() exp_in to get delta exp out near 0  when Lyy is in Q0  */
        tmp            = ISqrt16(extract_l(L_yy),
                      &exp); /* exp  out is now a delta shift with a later tmp Q15 multiplication in mind  */
#ifdef ENABLE_HR_MODE
        exp_shift      = add(exp, 16);   /*  up to Q16 */
        isqrtQ31_local = L_shl(L_deposit_l(tmp), exp_shift); /* new mantissa in a fixed  Q16  */
#else
        exp_shift      = add(exp, 16 - 15);   /*  up to Q16 */
        isqrtQ16_local = shl(tmp, exp_shift); /* new mantissa in a fixed  Q16  */
#endif
    }

    shift_num = norm_s(k_val_max);      /* simply account for the preknown fixed max possible pulseamp in y */
    shift_tot = sub(14 - 1, shift_num); /* upshift  to get  to Q14 */
    FOR (i = 0; i < dim; i++) /*  upshifted y[i]  used    */
    {
#ifdef ENABLE_HR_MODE
        L_tmp = Mpy_32_16(isqrtQ31_local, shl_pos(y[i], shift_num)); /* Q(16+0+shift_num +1    =   shift_num+1  */
        xq[i] = L_shl(L_tmp, shift_tot + 1); move32();     /* Q30 ,      */
#else
        L_tmp = L_mult(isqrtQ16_local, shl_pos(y[i], shift_num)); /* Q(16+0+shift_num +1    =   shift_num+1  */
        xq[i] = round_fx(L_shl(L_tmp, shift_tot)); move16();     /* Q14 ,      */
#endif
    }

    Dyn_Mem_Deluxe_Out();
}

