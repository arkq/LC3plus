/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

typedef struct
{
    Word16 inv_bin;
    Word16 numbytes;
    Word16 c_bp;
    Word16 c_bp_side;
    Word16 bytes;
    Word16 b_left;
    Word16 b_right;
    Word16 enc;
    Word16 sim_dec;
    Word16 bfi;
    Word16 be_bp_left;
    Word16 be_bp_right;
} Pc_State_fx;

typedef struct
{
    UWord32 ac_low_fx;
    UWord32 ac_range_fx;
    Word16  ac_cache_fx;
    Word16  ac_carry_fx;
    Word16  ac_carry_count_fx;
} Encoder_State_fx;

typedef struct
{
    UWord32 ac_low_fx;
    UWord32 ac_range_fx;
    UWord32 ac_help_fx;
    Word16  BER_detect;
    Pc_State_fx pc;
} Decoder_State_fx;

static void ac_dec_init_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side,
                           Decoder_State_fx *st_fx
);

#ifdef LL_INCL_HPVC
typedef struct
{
    UWord8 *ptr;
    Word16  bp_side;
    Word16  mask_side;
} EncBitBuffState_bw;

typedef struct
{
    UWord8 *ptr;
    Word16  bp;
    UWord32 ac_low_fx;
    UWord32 ac_range_fx;
    Word16  ac_cache_fx;
    Word16  ac_carry_fx;
    Word16  ac_carry_count_fx;
} EncBitBuffState_fw;

static void ac_enc_init_fx(Encoder_State_fx *st_fx);
static void ac_enc_shift_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx);
static void ac_encode_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx,
                         UWord32 cum_freq, UWord32 sym_freq);
static Word16 ac_enc_finish_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx);
static Word16 ac_dec_update_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side,
                               Word16 cur_bin, Decoder_State_fx *st_fx,
                               UWord32 cum_freq, UWord32 sym_freq);

static inline void write_bit_backward_st(EncBitBuffState_bw *s, Word16 bit)
{
    write_bit_backward(s->ptr, &s->bp_side, &s->mask_side, bit);
}

static inline void write_uint_backward_st(EncBitBuffState_bw *s, Word16 idx, Word16 nbits)
{
    write_indice_backward(s->ptr, &s->bp_side, &s->mask_side, idx, nbits);
}

static inline void ac_enc_init_fx_st(EncBitBuffState_fw *s, UWord8 *bytes)
{
    Encoder_State_fx tmp;
    s->ptr = bytes;
    s->bp = 0;
    ac_enc_init_fx(&tmp);
    s->ac_low_fx         = tmp.ac_low_fx;
    s->ac_range_fx       = tmp.ac_range_fx;
    s->ac_cache_fx       = tmp.ac_cache_fx;
    s->ac_carry_fx       = tmp.ac_carry_fx;
    s->ac_carry_count_fx = tmp.ac_carry_count_fx;
}

static inline void ac_encode_fx_st(EncBitBuffState_fw *s, UWord32 cum, UWord32 sym)
{
    Encoder_State_fx tmp;
    tmp.ac_low_fx         = s->ac_low_fx;
    tmp.ac_range_fx       = s->ac_range_fx;
    tmp.ac_cache_fx       = s->ac_cache_fx;
    tmp.ac_carry_fx       = s->ac_carry_fx;
    tmp.ac_carry_count_fx = s->ac_carry_count_fx;
    ac_encode_fx(s->ptr, &s->bp, &tmp, cum, sym);
    s->ac_low_fx         = tmp.ac_low_fx;
    s->ac_range_fx       = tmp.ac_range_fx;
    s->ac_cache_fx       = tmp.ac_cache_fx;
    s->ac_carry_fx       = tmp.ac_carry_fx;
    s->ac_carry_count_fx = tmp.ac_carry_count_fx;
}

static inline Word16 ac_enc_finish_fx_st(EncBitBuffState_fw *s)
{
    Encoder_State_fx tmp;
    Word16 bits;
    tmp.ac_low_fx         = s->ac_low_fx;
    tmp.ac_range_fx       = s->ac_range_fx;
    tmp.ac_cache_fx       = s->ac_cache_fx;
    tmp.ac_carry_fx       = s->ac_carry_fx;
    tmp.ac_carry_count_fx = s->ac_carry_count_fx;
    bits = ac_enc_finish_fx(s->ptr, &s->bp, &tmp);
    s->ac_low_fx         = tmp.ac_low_fx;
    s->ac_range_fx       = tmp.ac_range_fx;
    s->ac_cache_fx       = tmp.ac_cache_fx;
    s->ac_carry_fx       = tmp.ac_carry_fx;
    s->ac_carry_count_fx = tmp.ac_carry_count_fx;
    return bits;
}

/* HPVC-specific helpers (forward decls; defined later in this file) */
static Word16 ac_uni_get_step_fx(Word16 Ntot, Word16 cdf_bits);
static void ac_encode_uni_fx_st(EncBitBuffState_fw *s, Word16 val, Word16 Ntot,
                                Word16 bit_range_min, Word16 bit_range_max);
static void ac_encode_W32N_uni_fx_st(EncBitBuffState_fw *s, EncBitBuffState_bw *s_bw,
                                     Word32 L_cw, Word32 L_Ntot);

static Word16 ac_decode_generic_fx(Decoder_State_fx *st_fx, const UWord16 *cumfreq,
                                   Word16 num_sym);
static Word16 ac_decode_uni_fx(Decoder_State_fx *st_fx, Word16 Ntot,
                               Word16 cdf_bits_min, Word16 cdf_bits_max,
                               UWord32 *UL_cumfreq_upd_ptr, UWord32 *UL_symfreq_upd_ptr);
static Word32 ac_decode_W32N_uni_fx(Decoder_State_fx *st_fx, UWord8 *ptr, Word16 *bp,
                                    Word16 *bp_side, Word16 *mask_side, Word32 L_Ntot);
static Word16 read_uint_hpvc(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 numbits);
static Word16 post_HPVC_update_TCX_context(Word32 *L_x_tail);
#endif /* LL_INCL_HPVC */



static void pc_init_fx(Word16 n_pc, Word16 numbytes, Word16 be_bp_left, Word16 be_bp_right, Word16 L_spec,
                                        Word16 enc, Word16 sim_dec, Word16 bfi, Pc_State_fx *pc /* i/o: Pc State */
);
static Word16 check_pc_bytes(Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin,
                                           Word16 from_left, Pc_State_fx *pc /* i/o: Pc State */
);

static void ac_enc_init_fx(Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static void ac_enc_shift_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static void write_indice_forward(UWord8 *ptr, Word16 bp, Word16 indice, Word16 numbits);

static void ac_encode_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx, /* i/o: Encoder state */
                         UWord32 cum_freq, /* i  : Cumulative frequency up to symbol   */
                         UWord32 sym_freq  /* i  : Symbol probability                  */
);

static Word16 ac_enc_finish_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx /* i/o: Encoder state       */
);

static Word16 ac_decode_fx(                         /* o  : Decoded cumulative frequency    */
                           Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                           Word16            pki
#ifdef CR13_B_FIX_PC_BINS
                           ,
                           Word16 *bp,
                           Word16 *bp_side,
                           Word16 *mask_side,
                           Word16 cur_bin,
                           Word16 from_left
#endif
                           );
static Word16 ac_decode_tns_order(                         /* o  : Decoded cumulative frequency    */
                                  Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                  Word16            enable_lpc_weighting);
static Word16 ac_decode_tns_coef(                         /* o  : Decoded cumulative frequency    */
                                 Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                 Word16            pki);
static Word16 ac_dec_update_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin,
                               Decoder_State_fx *st_fx,    /* i/o: Decoder State           */
                               UWord32           cum_freq, /* i  : Cumulative frequency    */
                               UWord32           sym_freq  /* i  : Symbol frequency        */
);

/*************************************************************************/

#  ifdef ENABLE_HR_MODE

Word16 processAriEncoder_fx(UWord8 *bytes, Word16 bp_side_in, Word16 mask_side_in,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            Word32 nbbits,
#else
                            Word16 nbbits,
#endif
                            Word32 xq[],
                            Word16 *tns_order, Word16 tns_numfilters, Word16 *tns_idx,
                            Word16 lastnz,
                            Word16 *codingdata, UWord8 *resBits,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            Word32 numResBits,
#else
                            Word16 numResBits,
#endif
                            Word16 lsbMode,
                            Word16 enable_lpc_weighting,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            Word16 ll_side[],
                            UWord8 det_curve[],
                            Word32 x_for_resbits[],
                            Word32 max_resBits_len,
                            Word16 tns_lsb_num_remove,
                            Word16 ll_cbr,
                            Word32 d_fx_orig[],
                            UWord8 tns_lsb_remove[],
                            Word16 ll_adap_flag,
#ifdef LL_INCL_HPVC
                            HpvcEncCfg* hpvcEncCfgPtr,
#endif
                            Word16 b_relative,
                            const Word16* bands_offset,
                            Word16 bands_number,
                            Word16 L_spec,
#endif
                            lc3_scratch_t scratch)
{
    Word16 resbit, i1, i2;

#ifdef CR14_A_ADD_LOSSLESS_MODE
    Dyn_Mem_Deluxe_In(Encoder_State_fx st; Word16 bp, bp_side, mask_side, extra_bits;
                      Word32 a1, b1, a1_i, b1_i, a1_msb, b1_msb; Word16 lev1; Word32 nbits_side; Word16 tmp;
                      Word32 fill_bits; UWord8 * ptr; Word32 numResBitsEnc; Word16 * lsb, nlsbs; Counter i, k, lev; Word32 n;);
#else
    Dyn_Mem_Deluxe_In(Encoder_State_fx st; Word16 bp, bp_side, mask_side, extra_bits;
                      Word32 a1, b1, a1_i, b1_i, a1_msb, b1_msb; Word16 lev1; Word16 nbits_side; Word16 tmp;
                      Word16 fill_bits; UWord8 * ptr; Word16 numResBitsEnc; Word16 * lsb, nlsbs; Counter i, n, k, lev;);
#endif

#ifdef LL_INCL_HPVC
    EncBitBuffState_bw bitBuffState_bw;
    EncBitBuffState_fw bitBuffState_fw;
#endif

    lsb = (Word16*) lc3_scratch_push( scratch, 2 * lastnz * sizeof( *lsb ) );

    if ( scratch->max_scratch_calculation_only )
    {
        lsb = (Word16*) lc3_scratch_pop( scratch, lsb );
        Dyn_Mem_Deluxe_Out();
        return 0;
    }

    /* Init */
    a1_i = 0;
    move16();
    b1_i = 1;
    move16();
    bp = 0;
    move16();
    numResBitsEnc = 0;
    move16();
    nlsbs = 0;
    move16();
    ptr     = bytes;
    bp_side = bp_side_in;
    move16();
    mask_side = mask_side_in;
    move16();

    /*Start Encoding*/
    ac_enc_init_fx(&st);

    /* TNS data */
    FOR (n = 0; n < tns_numfilters; n++)
    {
        IF (tns_order[n] > 0)
        {
            ac_encode_fx(ptr, &bp, &st, ac_tns_order_cumfreq[enable_lpc_weighting][tns_order[n] - 1],
                         ac_tns_order_freq[enable_lpc_weighting][tns_order[n] - 1]);
            FOR (k = 0; k < tns_order[n]; k++)
            {
                ac_encode_fx(ptr, &bp, &st, ac_tns_coef_cumfreq[k][tns_idx[MAXLAG * n + k]],
                             ac_tns_coef_freq[k][tns_idx[MAXLAG * n + k]]);
            }
        }
    }

#ifdef LL_INCL_HPVC
#ifdef LL_HPVC_GLOBAL_FRAC
    IF (hpvcEncCfgPtr != NULL && hpvcEncCfgPtr->active_flag != 0 && ll_adap_flag != 0)
    {
        ASSERT(hpvcEncCfgPtr->mode >= -1 && hpvcEncCfgPtr->mode <= 1);
        {
            Word16 mode_sig = add(hpvcEncCfgPtr->mode, 1);  /* -1,0,1 -> 0,1,2 */
            ac_encode_fx(ptr, &bp, &st, hpvc_GlobalTab_cumfreq[mode_sig], hpvc_GlobalTab_freq[mode_sig]);
        }
    }
#endif
#endif

    IF (lsbMode == 0)
    {

#ifdef LL_INCL_HPVC_ARICODEC
        IF (hpvcEncCfgPtr->mode < 0)
        {
#endif
        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, 0,
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                write_bit_backward(ptr, &bp_side, &mask_side, L_and(xq[a1_i], 1));
                write_bit_backward(ptr, &bp_side, &mask_side, L_and(xq[b1_i], 1));
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }
            ELSE
            {
                a1 = L_abs(xq[a1_i]);
                b1 = L_abs(xq[b1_i]);
                FOR (lev = 0; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(L_shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(L_shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
#ifdef LL_INCL_HPVC_ARICODEC
        }
        ELSE
        {
            /* HPVC + TCX joint loop (Np = 2 or 8..128). */
            Word16 lastnz_local, k_incr, k_hpvc, f, k_tmp;
            Word16 Np, Kp, NsSafe, Ns, NsHdrSafe, NsHdr, treeCount, NpIdx;
            Word16 splitRule;
            BASOP_sub_sub_start("aricod_TCX&HPVC");
            UNUSED(splitRule);

            bitBuffState_bw.ptr       = ptr;
            bitBuffState_bw.bp_side   = bp_side;
            bitBuffState_bw.mask_side = mask_side;
            bitBuffState_fw.ptr       = ptr;
            bitBuffState_fw.bp        = bp;
            bitBuffState_fw.ac_low_fx         = st.ac_low_fx;
            bitBuffState_fw.ac_range_fx       = st.ac_range_fx;
            bitBuffState_fw.ac_cache_fx       = st.ac_cache_fx;
            bitBuffState_fw.ac_carry_fx       = st.ac_carry_fx;
            bitBuffState_fw.ac_carry_count_fx = st.ac_carry_count_fx;

            ASSERT(hpvcEncCfgPtr->mode == 0 || hpvcEncCfgPtr->mode == 1);
            lastnz_local = lastnz;
            move16();
            treeCount = 0;
            k_incr = 2;
            FOR (k = 0; k < lastnz_local; k += k_incr)
            {
                Np = 2;
                f = -1;
                k_hpvc = sub(k, hpvcEncCfgPtr->startCoef);

                ASSERT(((lastnz_local - hpvcEncCfgPtr->startCoef) % LL_HPVC_N_SIGNAL) == 0);
                if (k_hpvc >= 0)
                {
                    f = shr_pos(k_hpvc, N_SIGNAL_LOG);
                }

                k_tmp = -1;
                move16();
                IF (f >= 0)
                {
                    Word16 k_trunc = shl_pos(shr_pos(k_hpvc, N_SIGNAL_LOG), N_SIGNAL_LOG);
                    k_tmp = sub(k_hpvc, k_trunc);
                }

                IF (k_tmp == 0)
                {
                    Np = hpvcEncCfgPtr->Tx_dec[f];
                    ASSERT(Np == 2 || Np == hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->Np);
                    NpIdx = s_max(0, sub(15 - N_SIGNAL_LOG, norm_s(Np)));
                    ASSERT(Np == 2 || Np == (1 << (NpIdx + N_SIGNAL_LOG - 1)));
                    ac_encode_fx_st(&bitBuffState_fw, NpTabCDF[NpIdx], NpTabPDF[NpIdx]);
                }

                IF (sub(Np, 2) != 0)
                {
                    ASSERT(Np > 0);
                    Np = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->Np;
                    ASSERT(Np == hpvcEncCfgPtr->Tx_dec[f]);
                    k_incr = Np;
                    move16();

                    Kp = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->Kp;
                    NsSafe = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->NsSafe;
                    Ns = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->Ns;
                    NsHdr = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->NsHdr;
                    NsHdrSafe = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->NsHdrSafe;

                    ac_encode_fx_st(&bitBuffState_fw,
                                    UL_deposit_l(hpvc_KpTab_cumfreq[Kp]),
                                    UL_deposit_l(hpvc_KpTab_freq[Kp]));

                    IF (Kp > 0)
                    {
                        write_bit_backward_st(&bitBuffState_bw, hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->LS);

                        IF (NsHdrSafe > 0)
                        {
                            ASSERT(NsSafe >= 1);
                            Word16 tmp_sr = sub(15 - 4, norm_s(Np));
                            Word16 splitRuleMax = splitRuleMaxPerNp[tmp_sr];
                            Word16 splitRuleNsMin = splitRuleNsMinPerNp[tmp_sr];

                            if (NsSafe == LL_HPVC_SPLITRULE_GLOBAL_NSMIN)
                            {
                                splitRuleMax = s_min(2, splitRuleMax);
                            }
                            ASSERT(splitRuleMax >= -1 && splitRuleMax <= LL_HPVC_SPLITRULE_GLOBAL_MAX);

                            IF (splitRuleMax > 0 && NsSafe >= splitRuleNsMin)
                            {
                                splitRule = hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->splitRule;
                                ASSERT(splitRule >= 0 && splitRule <= LL_HPVC_SPLITRULE_GLOBAL_MAX);

                                IF (splitRuleMax == 3)
                                {
                                    write_bit_backward_st(&bitBuffState_bw, s_and(splitRule, 0x0001));
                                    write_bit_backward_st(&bitBuffState_bw, s_and(shr_pos(splitRule, 1), 0x0001));
                                }
                                IF (splitRuleMax == 1)
                                {
                                    write_bit_backward_st(&bitBuffState_bw, splitRule);
                                }
                                IF (splitRuleMax == 2)
                                {
                                    write_bit_backward_st(&bitBuffState_bw, s_min(splitRule, 1));
                                    IF (splitRule != 0)
                                    {
                                        write_bit_backward_st(&bitBuffState_bw, sub(splitRule, 1));
                                    }
                                }
                            }
                            Ns = splitRuleAdapt_Ns_NsHdr(hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->splitRule, Np, Kp, NsSafe, &NsHdr);

                            ac_encode_W32N_uni_fx_st(&bitBuffState_fw, &bitBuffState_bw,
                                                     hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->hdrIdx,
                                                     hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->hdrSz);

                            IF (NsHdr > 1)
                            {
                                ASSERT(NsHdr == LL_HPVC_NSHDR_MAX);
                                FOR (i = 0; i < NsHdr; i++)
                                {
                                    IF (hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->splitHdrLeafSz[i] > 0)
                                    {
                                        ac_encode_W32N_uni_fx_st(&bitBuffState_fw, &bitBuffState_bw,
                                                                 hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->splitHdrLeafIdx[i],
                                                                 hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->splitHdrLeafSz[i]);
                                    }
                                }
                            }
                        }

                        /* single leaf (Ns==1) or Tx several flat main leaves */
                        FOR (i = 0; i < Ns; i++)
                        {
                            IF (hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->flatLeafSz[i] > 0)
                            {
                                ac_encode_W32N_uni_fx_st(&bitBuffState_fw, &bitBuffState_bw,
                                                         hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->flatLeafIdx[i],
                                                         hpvcEncCfgPtr->HpvcTreeEnumCfgPtr->flatLeafSz[i]);
                            }
                        }
                    }

                    hpvcEncCfgPtr->HpvcTreeEnumCfgPtr++;
                    treeCount = add(treeCount, 1);

                    a1_i += Np;
                    b1_i += Np;
                    codingdata += (3 * shr_pos(Np, 1));
#ifdef LL_INCL_HPVC_OPT_AR_ENCODE_FIX
/* copied section from TCX loop  */

#ifdef LOSSLESS_192kHz
                    /* Side bits (in sync with the decoder) */
                    nbits_side = L_sub(nbbits, L_add(L_shl_pos(bitBuffState_bw.bp_side, 3), (Word32)sub(norm_s(bitBuffState_bw.mask_side), 6)));
#else
                    nbits_side = sub(nbbits, add(shl_pos(bitBuffState_bw.bp_side, 3), sub(norm_s(bitBuffState_bw.mask_side), 6)));
#endif
                    /* Residual bits (in sync with the decoder) */
                    extra_bits = sub(norm_ul(bitBuffState_fw.ac_range_fx), 6);
                    if (bitBuffState_fw.ac_cache_fx >= 0)
                    {
                        extra_bits = add(extra_bits, 8);
                    }
                    if (bitBuffState_fw.ac_carry_count_fx > 0)
                    {
                        extra_bits = add(extra_bits, shl_pos(bitBuffState_fw.ac_carry_count_fx, 3));
                    }

#ifdef LOSSLESS_192kHz
                    n = L_sub(nbbits, L_add(shl_pos(bitBuffState_fw.bp, 3), L_add(extra_bits, nbits_side)));
#else
                    n = sub(nbbits, add(shl_pos(bitBuffState_fw.bp, 3), add(extra_bits, nbits_side)));
#endif
#endif /* LL_INCL_HPVC_OPT_AR_ENCODE_FIX */
                }
                ELSE
                {
                    /* legacy TCX 2-tuple */
                    k_incr = 2;
                    move16();
                    IF (codingdata[1] < 0)
                    {
                        ac_encode_fx_st(&bitBuffState_fw, 0,
                                        ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
                    }
                    ELSE IF (codingdata[1] == 0)
                    {
                        ac_encode_fx_st(&bitBuffState_fw,
                                        ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                                        ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                        IF (xq[a1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[a1_i], 31));
                        }
                        IF (xq[b1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[b1_i], 31));
                        }
                    }
                    ELSE IF (sub(codingdata[1], 1) == 0)
                    {
                        ac_encode_fx_st(&bitBuffState_fw,
                                        ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                                        ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                        ac_encode_fx_st(&bitBuffState_fw,
                                        ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                                        ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                        write_bit_backward_st(&bitBuffState_bw, L_and(xq[a1_i], 1));
                        write_bit_backward_st(&bitBuffState_bw, L_and(xq[b1_i], 1));
                        IF (xq[a1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[a1_i], 31));
                        }
                        IF (xq[b1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[b1_i], 31));
                        }
                    }
                    ELSE
                    {
                        a1 = L_abs(xq[a1_i]);
                        b1 = L_abs(xq[b1_i]);
                        FOR (lev = 0; lev < codingdata[1]; lev++)
                        {
                            lev1 = s_min(lev, 3);
                            ac_encode_fx_st(&bitBuffState_fw,
                                            ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                            ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                            write_bit_backward_st(&bitBuffState_bw, s_and(L_shr_pos(a1, lev), 1));
                            write_bit_backward_st(&bitBuffState_bw, s_and(L_shr_pos(b1, lev), 1));
                        }
                        lev1 = s_min(codingdata[1], 3);
                        ac_encode_fx_st(&bitBuffState_fw,
                                        ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                                        ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                        IF (xq[a1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[a1_i], 31));
                        }
                        IF (xq[b1_i] != 0)
                        {
                            write_bit_backward_st(&bitBuffState_bw, L_lshr(xq[b1_i], 31));
                        }
                    }
                    a1_i += 2;
                    b1_i += 2;
                    codingdata += 3;
                }
            } /* FOR k */
            BASOP_sub_sub_end();

            bp        = bitBuffState_fw.bp;
            bp_side   = bitBuffState_bw.bp_side;
            mask_side = bitBuffState_bw.mask_side;
            st.ac_low_fx         = bitBuffState_fw.ac_low_fx;
            st.ac_range_fx       = bitBuffState_fw.ac_range_fx;
            st.ac_cache_fx       = bitBuffState_fw.ac_cache_fx;
            st.ac_carry_fx       = bitBuffState_fw.ac_carry_fx;
            st.ac_carry_count_fx = bitBuffState_fw.ac_carry_count_fx;
        }
#endif /* LL_INCL_HPVC_ARICODEC */
    }
    ELSE
    {
        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, 0,
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                a1_msb       = s_and(codingdata[2], 0x3);
                tmp          = L_and(xq[a1_i], 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                b1_msb       = shr_pos(codingdata[2], 2);
                tmp          = L_and(xq[b1_i], 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }
            ELSE
            {
                a1           = L_abs(xq[a1_i]);
                b1           = L_abs(xq[b1_i]);
                a1_msb       = L_shr_pos(a1, 1);
                tmp          = L_and(a1, 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                b1_msb       = L_shr_pos(b1, 1);
                tmp          = s_and(b1, 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC]);
                FOR (lev = 1; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(L_shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(L_shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[a1_i], 31));
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, L_lshr(xq[b1_i], 31));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
    }

    /* Side bits (in sync with the decoder) */
#ifdef CR14_A_ADD_LOSSLESS_MODE
    nbits_side = L_sub( nbbits, L_add( L_shl_pos(bp_side, 3), (Word32)sub( norm_s( mask_side ), 6 ) ) );
#else
    nbits_side = sub(nbbits, add(shl_pos(bp_side, 3), sub(norm_s(mask_side), 6)));
#endif

    /* Residual bits (in sync with the decoder) */
    extra_bits = sub(norm_ul(st.ac_range_fx), 6);
    if (st.ac_cache_fx >= 0)
    {
        extra_bits = add(extra_bits, 8);
    }
    if (st.ac_carry_count_fx > 0)
    {
        extra_bits = add(extra_bits, shl_pos(st.ac_carry_count_fx, 3));
    }

#ifdef CR14_A_ADD_LOSSLESS_MODE
    n = L_max(L_sub( nbbits, L_add( shl_pos(bp, 3), L_add( extra_bits, nbits_side))), 0);
#else
    n = s_max(sub(nbbits, add(shl_pos(bp, 3), add(extra_bits, nbits_side))), 0);
    move16();
#endif

    IF (lsbMode == 0)
    {
        numResBitsEnc = s_min(numResBits, n);

#ifdef CR14_A_ADD_LOSSLESS_MODE
        UNUSED(resBits);

        IF (ll_adap_flag)
        {
            numResBits = 0;
            {
                UWord8* eff_det_curve_buf = (UWord8*) lc3_scratch_push( scratch, sizeof( *eff_det_curve_buf ) * L_spec );
                IF (ll_cbr)
                {
                    compute_resbits_priority( x_for_resbits, det_curve, bands_offset, bands_number, L_spec, n, b_relative, eff_det_curve_buf, scratch );
                }
                ELSE
                {
                    basop_memcpy(eff_det_curve_buf, det_curve, sizeof(*det_curve) * L_spec);
                }
                residual_encoder_lossless( x_for_resbits, resBits, det_curve, eff_det_curve_buf, L_spec, max_resBits_len, &numResBits, scratch);
                eff_det_curve_buf = (UWord8*) lc3_scratch_pop( scratch, eff_det_curve_buf );
            }

            IF (tns_lsb_num_remove && ll_cbr)
            {
                basop_memset(tns_lsb_remove, tns_lsb_num_remove, L_spec * sizeof(tns_lsb_remove[0]));
                residual_encoder_lossless(d_fx_orig, resBits, tns_lsb_remove, tns_lsb_remove, L_spec, max_resBits_len, &numResBits, scratch);
            }
        }

        numResBitsEnc = L_min( numResBits, n );

        IF (ll_adap_flag)
        {
            FOR( i = 0; i < numResBitsEnc; i++ )
            {
                write_bit_backward( ptr, &bp_side, &mask_side, resBits[i] );
            }
        } ELSE {
            FOR( i = 0; i < numResBitsEnc; i++ )
            {
                resbit = 0;
                move16();
                i1 = shr( i, RESBITS_PACK_SHIFT );
                i2 = s_and( i, RESBITS_PACK_MASK );
                if ( s_and( resBits[i1], shl( 1, i2 ) ) )
                {
                    resbit = 1;
                }

                write_bit_backward( ptr, &bp_side, &mask_side, resbit );
            }
        }
#else /* CR14_A_ADD_LOSSLESS_MODE */
        FOR (i = 0; i < numResBitsEnc; i++)
        {
            resbit = 0; move16();
            i1 = shr(i, RESBITS_PACK_SHIFT);
            i2 = s_and(i, RESBITS_PACK_MASK);
            if (s_and(resBits[i1], shl(1, i2)))
            {
                resbit = 1;
            }
            write_bit_backward(ptr, &bp_side, &mask_side, resbit);
        }
#endif
    }
    ELSE
    {
        nlsbs = s_min(nlsbs, n);
        FOR (k = 0; k < nlsbs; k++)
        {
            write_bit_backward(ptr, &bp_side, &mask_side, lsb[k]);
        }
    }

    /* End arithmetic coder, overflow management */
    extra_bits = ac_enc_finish_fx(ptr, &bp, &st);

    /* Fill bits (for debugging, the exact number of fill bits cannot be computed in the decoder)*/
    fill_bits = nbbits - (bp * 8 + extra_bits + nbits_side + nlsbs + numResBitsEnc);

#ifdef CR14_A_ADD_LOSSLESS_MODE
    ll_side[0] = bp;
    ll_side[1] = bp_side;
#endif

    Dyn_Mem_Deluxe_Out();

    lsb = (Word16*) lc3_scratch_pop( scratch, lsb );
    return fill_bits;
}

#  else /* ENABLE_HR_MODE */

Word16 processAriEncoder_fx(UWord8 *bytes, Word16 bp_side_in, Word16 mask_side_in, Word16 nbbits, Word16 xq[],
                            Word16 *tns_order, Word16 tns_numfilters, Word16 *tns_idx, Word16 lastnz,
                            Word16 *codingdata, UWord8 *resBits, Word16 numResBits, Word16 lsbMode,
                            Word16 enable_lpc_weighting, lc3_scratch_t scratch)
{
#ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(Encoder_State_fx st; Word16 bp, bp_side, mask_side, extra_bits;
                      Word16 a1, b1, a1_i, b1_i, a1_msb, b1_msb; Word16 lev1; Word16 nbits_side; Word16 tmp;
                      Word16 fill_bits; UWord8 * ptr; Word16 numResBitsEnc; Word16 * lsb, nlsbs; Counter i, n, k, lev;
                      );
#else
    Dyn_Mem_Deluxe_In(Encoder_State_fx st; Word16 bp, bp_side, mask_side, extra_bits;
                      Word16 a1, b1, a1_i, b1_i, a1_msb, b1_msb; Word16 lev1; Word16 nbits_side; Word16 tmp;
                      Word16 fill_bits; UWord8 * ptr; Word16 numResBitsEnc; Word16 * lsb, nlsbs; Counter i, n, k, lev;
                      );
#endif


    lsb = (Word16*) lc3_scratch_push( scratch, 2 * lastnz * sizeof( *lsb ) );

    if ( scratch->max_scratch_calculation_only )
    {
        lsb = (Word16*) lc3_scratch_pop( scratch, lsb );
        Dyn_Mem_Deluxe_Out();
        return 0;
    }

    /* Init */
    a1_i = 0;
    move16();
    b1_i = 1;
    move16();
    bp = 0;
    move16();
    numResBitsEnc = 0;
    move16();
    nlsbs = 0;
    move16();
    ptr     = bytes;
    bp_side = bp_side_in;
    move16();
    mask_side = mask_side_in;
    move16();

    /*Start Encoding*/
    ac_enc_init_fx(&st);

    /* TNS data */
    FOR (n = 0; n < tns_numfilters; n++)
    {
        IF (tns_order[n] > 0)
        {
            ac_encode_fx(ptr, &bp, &st, ac_tns_order_cumfreq[enable_lpc_weighting][tns_order[n] - 1],
                         ac_tns_order_freq[enable_lpc_weighting][tns_order[n] - 1]);
            FOR (k = 0; k < tns_order[n]; k++)
            {
                ac_encode_fx(ptr, &bp, &st, ac_tns_coef_cumfreq[k][tns_idx[MAXLAG * n + k]],
                             ac_tns_coef_freq[k][tns_idx[MAXLAG * n + k]]);
            }
        }
    }

    IF (lsbMode == 0)
    {

        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][0],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                write_bit_backward(ptr, &bp_side, &mask_side, s_and(xq[a1_i], 1));
                write_bit_backward(ptr, &bp_side, &mask_side, s_and(xq[b1_i], 1));
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE
            {
                a1 = abs_s(xq[a1_i]);
                b1 = abs_s(xq[b1_i]);
                FOR (lev = 0; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
    }
    ELSE
    {
        /*Main Loop through the 2-tuples*/
        FOR (k = 0; k < lastnz; k += 2)
        {
            IF (codingdata[1] < 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][0],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][0]);
            }
            ELSE IF (codingdata[1] == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][codingdata[2]]);
                IF (xq[a1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (xq[b1_i] != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE IF (sub(codingdata[1], 1) == 0)
            {
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0]]][VAL_ESC]);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[1]]][codingdata[2]]);
                a1_msb       = s_and(codingdata[2], 0x3);
                tmp          = s_and(xq[a1_i], 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                b1_msb       = shr_pos(codingdata[2], 2);
                tmp          = s_and(xq[b1_i], 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }
            ELSE
            {
                a1           = abs_s(xq[a1_i]);
                b1           = abs_s(xq[b1_i]);
                a1_msb       = shr_pos(a1, 1);
                tmp          = s_and(a1, 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (a1_msb == 0 && tmp > 0)
                {
                    if (xq[a1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[a1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                b1_msb       = shr_pos(b1, 1);
                tmp          = s_and(b1, 1);
                lsb[nlsbs++] = tmp;
                move16();
                test();
                IF (b1_msb == 0 && tmp > 0)
                {
                    if (xq[b1_i] > 0)
                    {
                        lsb[nlsbs++] = 0;
                        move16();
                    }
                    if (xq[b1_i] < 0)
                    {
                        lsb[nlsbs++] = 1;
                        move16();
                    }
                }
                ac_encode_fx(ptr, &bp, &st, ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[0]]][VAL_ESC]);
                FOR (lev = 1; lev < codingdata[1]; lev++)
                {
                    lev1 = s_min(lev, 3);
                    ac_encode_fx(ptr, &bp, &st,
                                 ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC],
                                 ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][VAL_ESC]);
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(a1, lev), 1));
                    write_bit_backward(ptr, &bp_side, &mask_side, s_and(shr_pos(b1, lev), 1));
                }
                lev1 = s_min(codingdata[1], 3);
                ac_encode_fx(ptr, &bp, &st,
                             ari_spec_cumfreq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]],
                             ari_spec_freq[ari_spec_lookup[codingdata[0] + Tab_esc_nb[lev1]]][codingdata[2]]);
                IF (a1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[a1_i], 15));
                }
                IF (b1_msb != 0)
                {
                    write_bit_backward(ptr, &bp_side, &mask_side, lshr(xq[b1_i], 15));
                }
            }

            a1_i += 2;
            b1_i += 2;
            codingdata += 3;

        } /*end of the 2-tuples loop*/
    }

    /* Side bits (in sync with the decoder) */
    nbits_side = sub(nbbits, add(shl_pos(bp_side, 3), sub(norm_s(mask_side), 6)));

    /* Residual bits (in sync with the decoder) */
    extra_bits = sub(norm_ul(st.ac_range_fx), 6);
    if (st.ac_cache_fx >= 0)
    {
        extra_bits = add(extra_bits, 8);
    }
    if (st.ac_carry_count_fx > 0)
    {
        extra_bits = add(extra_bits, shl_pos(st.ac_carry_count_fx, 3));
    }

    n = s_max(sub(nbbits, add(shl_pos(bp, 3), add(extra_bits, nbits_side))), 0);
    move16();

    IF (lsbMode == 0)
    {
        numResBitsEnc = s_min(numResBits, n);
        FOR (i = 0; i < numResBitsEnc; i++)
        {
            FOR (i = 0; i < numResBitsEnc; i++)
            {
                write_bit_backward(ptr, &bp_side, &mask_side, resBits[i]);
            }
        }
    }
    ELSE
    {
        nlsbs = s_min(nlsbs, n);
        FOR (k = 0; k < nlsbs; k++)
        {
            write_bit_backward(ptr, &bp_side, &mask_side, lsb[k]);
        }
    }

    /* End arithmetic coder, overflow management */
    extra_bits = ac_enc_finish_fx(ptr, &bp, &st);

    /* Fill bits (for debugging, the exact number of fill bits cannot be computed in the decoder)*/
    fill_bits = nbbits - (bp * 8 + extra_bits + nbits_side + nlsbs + numResBitsEnc);

    Dyn_Mem_Deluxe_Out();

    lsb = (Word16*) lc3_scratch_pop( scratch, lsb );
    return fill_bits;
}

#endif /* ENABLE_HR_MODE */

void processAriDecoder_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side,
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32 nbbits,
#else
    Word16 nbbits,
#endif
    Word16 L_spec,
    Word16 fs_idx, Word16 enable_lpc_weighting, Word16 tns_numfilters, Word16 lsbMode,
    Word16 lastnz, Word16 *bfi, Word16 *tns_order, Word16 fac_ns_idx, Word16 gg_idx,
    LC3PLUS_FrameDuration frame_dms,
                          Word16 n_pc, Word16 be_bp_left, Word16 be_bp_right, Word16 mode, Word16 *spec_inv_idx,
                          Word16 *b_left,
#ifdef CR14_A_ADD_LOSSLESS_MODE
                          Word32* resBits,
#else
                          Word16 *resBits,
#endif
#  ifdef ENABLE_HR_MODE
                          Word32 *x,
#  else
                          Word16 *x,
#  endif
                          Word16 *nf_seed, UWord8 *resQdata, Word16 *tns_idx, Word16 *zero_frame, lc3_scratch_t scratch
#  ifdef ENABLE_HR_MODE
                          , Word16 hrmode
#  endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                          , UWord8 *det_curve, Word16 ll_flag
                          , Word16 lossless
#    endif
#ifdef LL_INCL_HPVC
                          , HpvcDecCfg* hpvcDecCfgPtr
#endif
)
{
    Decoder_State_fx st;
#ifdef LL_INCL_HPVC
    /* hpvcDecCfgPtr used in HPVC branch below */
#endif
#  ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32           resbit, i1, i2;
#else
    Word16           resbit, i1, i2;
#endif
    Word32           a, b;
#  else
    Word16           a, b;
#  endif
    Word16           t, a1, b1, a1_i, b1_i, bp;
    Word16           esc_nb;
    Word16           rateFlag;
    Word16           r;
    Word16           nt_half;
    Word16           c;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32           nbits_side;
#else
    Word16           nbits_side;
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Word32           extra_bits, nbits_ari;
#else
    Word16           extra_bits, nbits_ari;
#endif
    UWord8 *         ptr;
    Word32           tmp32;
    Word16           lsb_ind_c;
    Word16 *         lsb_ind;
    Word16           tmp;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Counter          k, lev;
    Counter          i;
    Word32           n;
#else
    Counter          n, k, lev;
    Counter          i;
#endif
    Word16           max_lev = 14;

#ifdef DYNMEM_COUNT
struct _dynmem
{
Decoder_State_fx st;
        Pc_State_fx      pc;
        Word16           resbit, i1, i2;
#    ifdef ENABLE_HR_MODE
        Word32           a, b;
#    else
        Word16           a, b;
#    endif
        Word16           t, a1, b1, a1_i, b1_i, bp;
        Word16           esc_nb;
        Word16           rateFlag;
        Word16           r;
        Word16           nt_half;
        Word16           c;
        Word16           nbits_side, extra_bits, nbits_ari;
        UWord8 *         ptr;
        Word32           tmp32;
        Word16           lsb_ind_c;
        Word16 *         lsb_ind;
        Word16           tmp;
        Counter          i, n, k, lev;
    };
    Dyn_Mem_In("processAriDecoder_fx", sizeof(struct _dynmem));
#endif

#  ifdef ENABLE_HR_MODE
    if (hrmode == 1)
    {
        max_lev = max_lev + 8;
    }
#  endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    /* match encoder maxlevs=21 + 8-bit norm */
    if ( ll_flag == 1 )
    {
        max_lev = 31;
    }
#endif

    lsb_ind = (Word16*) lc3_scratch_push( scratch, L_spec * sizeof( *lsb_ind ) );

    if ( scratch->max_scratch_calculation_only )
    {
        lsb_ind = (Word16*) lc3_scratch_pop( scratch, lsb_ind);
#ifdef DYNMEM_COUNT
        Dyn_Mem_Out();
#endif
        return;
    }

    /* Rate flag */
    rateFlag = 0;
    move16();
#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF ( fs_idx != 5 || (fs_idx==5 && ll_flag) )
    {
        IF ((L_sub(nbbits, add(160, i_mult(fs_idx, 160))) > 0) || ll_flag )
#else
    IF (fs_idx != 5)
    {
        IF (sub(nbbits, add(160, i_mult(fs_idx, 160))) > 0)
#endif
        {
            rateFlag = 2 << NBITS_CONTEXT;
            move16();
        }
    }

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF ( lossless && !ll_flag && fs_idx >= 5 )
    {
        IF ((L_sub(nbbits, add(160, i_mult(fs_idx, 160))) > 0) || ll_flag )
        {
            rateFlag = 2 << NBITS_CONTEXT;
            move16();
        }
    }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    pc_init_fx(n_pc, L_shr_pos(nbbits, 3), be_bp_left, be_bp_right, L_spec, mode==1, mode==2, *bfi, &st.pc);
#else
    pc_init_fx(n_pc, shr_pos(nbbits, 3), be_bp_left, be_bp_right, L_spec, mode==1, mode==2, *bfi, &st.pc);
#endif

    /* Init */
    nt_half = shr_pos(L_spec, 1);
    c       = 0;
    move16();
    t = 0;
    move16();
    a1_i = 0;
    move16();
    b1_i = 1;
    move16();
    bp = 0;
    move16();
    if (mode != 1)
    {
        bp = add(bp, st.pc.bytes);
        move16();
    }
    *spec_inv_idx = L_spec;
    move16();
    *b_left = -1;
    move16();
    lsb_ind_c = 0;
    move16();

ptr = bytes;

/* Start Decoding */
ac_dec_init_fx(ptr, &bp, bp_side, mask_side, &st);

    /* Decode TNS data */
    tmp = MAXLAG;
IF (sub(frame_dms, LC3PLUS_FRAME_DURATION_2p5MS) == 0)
{
tmp = shr_pos(tmp, 1);
}
IF (sub(frame_dms, LC3PLUS_FRAME_DURATION_5MS) == 0)
{
tmp = shr_pos(tmp, 1);
}

    FOR (n = 0; n < tns_numfilters; n++)
    {
        IF (tns_order[n] > 0)
        {
            tns_order[n] = ac_decode_tns_order(&st, enable_lpc_weighting);
            move16();
            tns_order[n] = add(tns_order[n], 1);
            move16();
            IF (tns_order[n] > tmp)
            {
                GOTO ber_detect;
            }
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, 0, &st,
                                 ac_tns_order_cumfreq[enable_lpc_weighting][tns_order[n] - 1],
                                 ac_tns_order_freq[enable_lpc_weighting][tns_order[n] - 1]) != 0)
            {
                GOTO ber_detect;
            }
            FOR (k = 0; k < tns_order[n]; k++)
            {
                IF (sub(*bp_side, bp) < 0)
                {
                    GOTO ber_detect;
                }
                tns_idx[MAXLAG * n + k] = ac_decode_tns_coef(&st, k);
                move16();
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, 0, &st,
                                     ac_tns_coef_cumfreq[k][tns_idx[MAXLAG * n + k]],
                                     ac_tns_coef_freq[k][tns_idx[MAXLAG * n + k]]) != 0)
                {
                    GOTO ber_detect;
                }
            }
        }
    }
    IF (st.BER_detect > 0)
    {
        GOTO ber_detect;
    }

#ifdef LL_INCL_HPVC
#ifdef LL_HPVC_GLOBAL_FRAC
    IF (hpvcDecCfgPtr != NULL)
    {
        hpvcDecCfgPtr->mode = -1;
        move16();
    }
    IF (hpvcDecCfgPtr != NULL && hpvcDecCfgPtr->active_flag != 0 && lossless != 0 && ll_flag != 0)
    {
        Word16 mixed_mode_pre = ac_decode_generic_fx(&st, hpvc_GlobalTab_cumfreq, 3);
        UWord32 UL_cumFreq = UL_deposit_l(hpvc_GlobalTab_cumfreq[mixed_mode_pre]);
        UWord32 UL_symFreq = UL_deposit_l(hpvc_GlobalTab_freq[mixed_mode_pre]);
        IF (st.BER_detect != 0)
        {
            GOTO ber_detect;
        }
        IF (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, 0, &st, UL_cumFreq, UL_symFreq) != 0)
        {
            GOTO ber_detect;
        }
        hpvcDecCfgPtr->mode = sub(mixed_mode_pre, 1);  /* 0,1,2 -> -1,0,1 */

        IF (hpvcDecCfgPtr->mode >= 0)
        {
#ifdef LL_HPVC_ALIGN_STARTCOEFF_TO_LASTNZ
            hpvcDecCfgPtr->startCoefNom = hpvcDecCfgPtr->startCoefListNom[hpvcDecCfgPtr->mode];
            hpvc_adjust_startcoefs(hpvcDecCfgPtr->startCoefListNom, lastnz, LL_HPVC_N_SIGNAL, N_SIGNAL_LOG, hpvcDecCfgPtr->startCoefList);
            hpvcDecCfgPtr->startCoef = hpvcDecCfgPtr->startCoefList[hpvcDecCfgPtr->mode];
#else
            hpvcDecCfgPtr->startCoefList[0] = hpvcDecCfgPtr->startCoefListNom[0];
            hpvcDecCfgPtr->startCoefList[1] = hpvcDecCfgPtr->startCoefListNom[1];
            hpvcDecCfgPtr->startCoefNom = hpvcDecCfgPtr->startCoefListNom[hpvcDecCfgPtr->mode];
            hpvcDecCfgPtr->startCoef = hpvcDecCfgPtr->startCoefList[hpvcDecCfgPtr->mode];
#endif
        }
    }
#endif
#endif

IF (lsbMode == 0)
{

#ifdef LL_INCL_HPVC_ARICODEC_DEC
    IF (hpvcDecCfgPtr == NULL || hpvcDecCfgPtr->mode < 0)
    {
#endif
/*Main Loop through the 2-tuples*/
FOR (k = 0; k < lastnz; k += 2)
{

/* Get context */
t = add(c, rateFlag);
if (sub(k, nt_half) > 0)
{
t = add(t, 1 << NBITS_CONTEXT);
}

            r = ac_decode_fx(&st, ari_spec_lookup[t]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
            );
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t]][r],
                                 ari_spec_freq[ari_spec_lookup[t]][r]) != 0)
            {
                GOTO ber_detect;
            }

            IF (r == 0)
            {
                x[a1_i] = 0;
                move16();
                x[b1_i] = 0;
                move16();
                c = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(r, VAL_ESC) < 0)
            {
                a = s_and(r, 0x3);
                b = shr_pos(r, 2);
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a, b), 1));
                IF (a > 0)
                {
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        a = negate(a);
                    }
                }
                x[a1_i] = a;
                move16();
                IF (b > 0)
                {
if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
                }
                x[b1_i] = b;
                move16();
            }
            ELSE
            {
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
GOTO ber_detect;
}
a = read_bit(ptr, bp_side, mask_side);
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
GOTO ber_detect;
}
                b = read_bit(ptr, bp_side, mask_side);
                r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[1]]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
                );
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                     ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[1]]][r],
                                     ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[1]]][r]) != 0)
                {
                    GOTO ber_detect;
                }
                IF (sub(r, VAL_ESC) < 0)
                {
                    a1 = s_and(r, 0x3);
                    b1 = shr_pos(r, 2);
                    a  = add(shl_pos(a1, 1), a);
                    b  = add(shl_pos(b1, 1), b);
                    IF (a > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a;
                    move16();
                    IF (b > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i] = b;
                    move16();
                    c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1, b1), 1), 1));
                }
                ELSE
                {
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
a = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), a);
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    b = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), b);
                    FOR (lev = 2; lev < max_lev; lev++)
                    {
                        esc_nb = s_min(lev, 3);
                        r      = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[esc_nb]]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
                        );
                        if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                             ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r],
                                             ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r]) != 0)
                        {
                            GOTO ber_detect;
                        }
                        IF (sub(r, VAL_ESC) < 0)
                        {
                            BREAK;
                        }
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }

#  ifdef ENABLE_HR_MODE
                        a = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  else
                        a = add(shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  endif

                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
                        {
                            GOTO ber_detect;
                        }
#  ifdef ENABLE_HR_MODE
                        b = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  else
                        b = add(shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  endif
                    }
                    /* check for bitflip */
                    IF (sub(lev, max_lev) == 0)
                    {
                        GOTO ber_detect;
                    }

                    b1 = shr_pos(r, 2);
                    a1 = s_and(r, 0x3);

#  ifdef ENABLE_HR_MODE
                    a  = L_add(L_shl(a1, lev), a);
                    b  = L_add(L_shl(b1, lev), b);
#  else
                    a = add(shl(a1, lev), a);
                    b = add(shl(b1, lev), b);
#  endif
                    IF (a > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
#  ifdef ENABLE_HR_MODE
                            a = L_negate(a);
#  else
                            a = negate(a);
#  endif
                        }
                    }
                    x[a1_i] = a;
                    move16();
                    IF (b > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
#  ifdef ENABLE_HR_MODE
                            b = L_negate(b);
#  else
                            b = negate(b);
#  endif
                        }
                    }
                    x[b1_i] = b;
                    move16();
                    c = add(shl_pos(s_and(c, 0xf), 4), add(esc_nb, 12));
                }
            }

            test();
            test();
            IF ((sub(sub(bp, *bp_side), 3) > 0 && sub(st.pc.c_bp, st.pc.c_bp_side) == 0) || st.BER_detect > 0)
{
GOTO ber_detect;
}

a1_i += 2;
b1_i += 2;
}
#ifdef LL_INCL_HPVC_ARICODEC_DEC
    } /* IF mode<0 */
    ELSE
    {
        /* HPVC + TCX joint decoder loop. */
        Word16 lastnz_local, k_incr, k_hpvc, f, k_tmp, Ntmp;
        Word16 Np, Kp, NsSafe, Ns, NsHdrSafe, NsHdr, treeCount, NpIdx, splitRule;
        UWord32 UL_cumFreq, UL_symFreq;
        Word16 hdrVals[LL_HPVC_NS_MAX];
        Word16 splitHdrVals[LL_HPVC_NSHDR_MAX];
        HpvcTreeEnumCfg *hpvcTreeEnumCfg;
        Word16 ii_, accPos;
        Word16 LS_ind_sub;
        Word16 v_tmp[LL_HPVC_NP_MAX];

        ASSERT(hpvcDecCfgPtr->mode == 0 || hpvcDecCfgPtr->mode == 1);
        lastnz_local = lastnz;
        move16();
        treeCount = 0;
        k_incr = 2;
        move16();
        FOR (k = 0; k < lastnz_local; k += k_incr)
        {
            Np = 2;
            f = -1;
            k_hpvc = sub(k, hpvcDecCfgPtr->startCoef);

            if (k_hpvc >= 0)
            {
                f = shr_pos(k_hpvc, N_SIGNAL_LOG);
            }
            k_tmp = -1;
            move16();
            IF (f >= 0)
            {
                Word16 k_trunc = shl_pos(shr_pos(k_hpvc, N_SIGNAL_LOG), N_SIGNAL_LOG);
                k_tmp = sub(k_hpvc, k_trunc);
            }

            k_incr = 2;
            move16();

            IF (k_tmp == 0)
            {
                /* multiplex HPVC tree e.g. Np = 64,128 */
                NpIdx = ac_decode_generic_fx(&st, NpTabCDF, LL_HPVC_NB_BLOCKTYPE);
                IF (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, NpTabCDF[NpIdx], NpTabPDF[NpIdx]) != 0)
                {
                    GOTO ber_detect;
                }
                Np = shl_pos(1, add(NpIdx, N_SIGNAL_LOG - 1));
                if (NpIdx == 0)
                {
                    Np = 2;
                    move16();
                }
                ASSERT(Np == 2 || Np == (LL_HPVC_N_SIGNAL << (NpIdx - 1)));
                k_incr = Np;
                move16();
            }

            IF (sub(Np, 2) != 0)
            {
                ASSERT(Np > 0);
                hpvcTreeEnumCfg = hpvcDecCfgPtr->HpvcTreeEnumCfgPtr;
                hpvcTreeEnumCfg->Np = Np;
                k_incr = Np;
                move16();

                Kp = ac_decode_generic_fx(&st, hpvc_KpTab_cumfreq, LL_HPVC_KP_MAX + 1);
                UL_cumFreq = UL_deposit_l(hpvc_KpTab_cumfreq[Kp]);
                UL_symFreq = UL_deposit_l(hpvc_KpTab_freq[Kp]);
                IF (st.BER_detect != 0)
                {
                    GOTO ber_detect;
                }
                IF (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, UL_cumFreq, UL_symFreq) != 0)
                {
                    GOTO ber_detect;
                }

                hpvcTreeEnumCfg->Kp = Kp;
                move16();

                NsSafe = MPVQ_HPVC_SplitSetup(Np, Kp, &NsHdrSafe);
                hpvcTreeEnumCfg->NsSafe = NsSafe;
                move16();
                hpvcTreeEnumCfg->splitRule = -1;
                move16();

                Ns = NsSafe;
                hpvcTreeEnumCfg->Ns = Ns;
                move16();

                NsHdr = NsHdrSafe;
                hpvcTreeEnumCfg->NsHdr = NsHdr;
                move16();
                hpvcTreeEnumCfg->NsHdrSafe = NsHdrSafe;
                move16();

                {
                    Word32 *L_x_top = &x[a1_i];
                    basop_memset(L_x_top, 0, Np * sizeof(Word32));

                    IF (Kp > 0)
                    {
                        hpvcTreeEnumCfg->LS = read_bit(ptr, bp_side, mask_side);

                        IF (NsHdrSafe > 0)
                        {
                            ASSERT(hpvcTreeEnumCfg->NsHdr > 0);
                            hpvcTreeEnumCfg->splitRule = -1;
                            move16();
                            splitRule = -1;
                            move16();

                            Word16 tmp_sr = sub(15 - 4, norm_s(Np));
                            Word16 splitRuleMax = splitRuleMaxPerNp[tmp_sr];
                            Word16 splitRuleNsMin = splitRuleNsMinPerNp[tmp_sr];

                            if (NsSafe == LL_HPVC_SPLITRULE_GLOBAL_NSMIN)
                            {
                                splitRuleMax = s_min(2, splitRuleMax);
                            }
                            ASSERT(splitRuleMax >= -1 && splitRuleMax <= LL_HPVC_SPLITRULE_GLOBAL_MAX);

                            IF (splitRuleMax > 0 && NsSafe >= splitRuleNsMin)
                            {
                                if (splitRuleMax == 3)
                                {
                                    Word16 tmp_b;
                                    tmp_b = read_bit(ptr, bp_side, mask_side);
                                    splitRule = read_bit(ptr, bp_side, mask_side);
                                    splitRule = add(shl_pos(splitRule, 1), tmp_b);
                                }
                                if (splitRuleMax == 1)
                                {
                                    splitRule = read_bit(ptr, bp_side, mask_side);
                                }
                                if (splitRuleMax == 2)
                                {
                                    Word16 tmp_b;
                                    splitRule = read_bit(ptr, bp_side, mask_side);
                                    IF (splitRule != 0)
                                    {
                                        tmp_b = read_bit(ptr, bp_side, mask_side);
                                        splitRule = add(splitRule, tmp_b);
                                    }
                                }
                                hpvcTreeEnumCfg->splitRule = splitRule;
                                move16();
                                ASSERT(splitRule >= 0 && splitRule <= LL_HPVC_SPLITRULE_GLOBAL_MAX);
                            }

                            /* Ns, NsHdr may be reduced by splitRules (1,2,3) */
                            Ns = splitRuleAdapt_Ns_NsHdr(hpvcTreeEnumCfg->splitRule, Np, Kp, NsSafe, &NsHdr);
                            hpvcTreeEnumCfg->Ns = Ns;
                            hpvcTreeEnumCfg->NsHdr = NsHdr;

                            Ntmp = LL_HPVC_NSHDR_MAX;
                            Word16 NsHdrMinusOne;
                            NsHdrMinusOne = sub(hpvcTreeEnumCfg->NsHdr, 1);
                            if (NsHdrMinusOne == 0)
                            {
                                Ntmp = Ns;
                                move16();
                            }

                            hpvcTreeEnumCfg->hdrSz = hpvc_leaf_sz_Nmpvq(Ntmp, Kp);
                            hpvcTreeEnumCfg->hdrIdx = ac_decode_W32N_uni_fx(&st, ptr, &bp, bp_side, mask_side, hpvcTreeEnumCfg->hdrSz);

                            IF (st.BER_detect != 0)
                            {
                                GOTO ber_detect;
                            }
                            Word16 *hdrValsPtr = splitHdrVals;
                            basop_memset(splitHdrVals, 0, LL_HPVC_NSHDR_MAX * sizeof(Word16));
                            basop_memset(hdrVals, 0, Ns * sizeof(Word16));
                            ASSERT(Ns <= LL_HPVC_NS_MAX);
                            if (NsHdrMinusOne == 0)
                            {
                                hdrValsPtr = hdrVals;
                            }
                            st.BER_detect = hpvc_leaf_dec_deidx_fx(hdrValsPtr, Kp, Ntmp, hpvcTreeEnumCfg->LS, (UWord32)hpvcTreeEnumCfg->hdrIdx);
                            IF (st.BER_detect != 0)
                            {
                                GOTO ber_detect;
                            }

                            IF (NsHdrMinusOne > 0)
                            {
                                ASSERT((NsHdr == LL_HPVC_NSHDR_MAX) && (Ns > 1) && (Ns <= LL_HPVC_NS_MAX));
                                Word16 nHdrLeaves[LL_HPVC_NSHDR_MAX];
                                hpvc_get_quasi_uniform_sizes(Ns, NsHdr, nHdrLeaves);

                                Word16 accHdrPos = 0;
                                FOR (i = 0; i < NsHdr; i++)
                                {
                                    hpvcTreeEnumCfg->splitHdrNsDbg[i] = nHdrLeaves[i];
                                    hpvcTreeEnumCfg->splitHdrLeafSz[i] = hpvc_leaf_sz_Nmpvq(nHdrLeaves[i], abs_s(splitHdrVals[i]));
                                    IF (hpvcTreeEnumCfg->splitHdrLeafSz[i] != 0)
                                    {
                                        hpvcTreeEnumCfg->splitHdrLeafIdx[i] = ac_decode_W32N_uni_fx(&st, ptr, &bp, bp_side, mask_side, hpvcTreeEnumCfg->splitHdrLeafSz[i]);
                                        IF (st.BER_detect != 0)
                                        {
                                            GOTO ber_detect;
                                        }
                                        Word16 LS_sub = lshr(splitHdrVals[i], 15);
                                        st.BER_detect |= hpvc_leaf_dec_deidx_fx(&hdrVals[accHdrPos], abs_s(splitHdrVals[i]), nHdrLeaves[i], LS_sub, (UWord32)hpvcTreeEnumCfg->splitHdrLeafIdx[i]);
                                        IF (st.BER_detect != 0)
                                        {
                                            GOTO ber_detect;
                                        }
                                    }
                                    accHdrPos = add(accHdrPos, nHdrLeaves[i]);
                                }
                            }
                        } /* NsHdrSafe > 0 */

                        /* single leaf (Ns==1) or Rx several flat main leaves */
                        Word16 nFlatLeaves[16];
                        IF (NsHdr == 0)
                        {
                            nFlatLeaves[0] = Np;
                            move16();
                            hdrVals[0] = Kp;
                            move16();
                            if (hpvcTreeEnumCfg->LS != 0)
                            {
                                hdrVals[0] = negate(Kp);
                            }
                            ASSERT(Ns == 1);
                        }
                        ELSE
                        {
                            IF (splitRule <= 0 || sub(splitRule, 3) == 0)
                            {
                                hpvc_get_quasi_uniform_sizes(Np, Ns, nFlatLeaves);
                            }
                            ELSE
                            {
                                hpvc_get_log_sizes(splitRule, Np, Ns, nFlatLeaves);
                            }
                        }

                        accPos = 0;
                        move16();
                        for (i = 0; i < Ns; i++)
                        {
                            Word32 *L_x = &x[add(a1_i, accPos)];
                            IF (hdrVals[i] != 0)
                            {
                                hpvcTreeEnumCfg->flatLeafSz[i] = hpvc_leaf_sz_Nmpvq(nFlatLeaves[i], abs_s(hdrVals[i]));
                                hpvcTreeEnumCfg->flatLeafIdx[i] = ac_decode_W32N_uni_fx(&st, ptr, &bp, bp_side, mask_side, hpvcTreeEnumCfg->flatLeafSz[i]);
                                IF (st.BER_detect != 0)
                                {
                                    GOTO ber_detect;
                                }
                                ASSERT(abs(hdrVals[i]) <= LL_HPVC_KP_MAX);
                                LS_ind_sub = lshr(hdrVals[i], 15);
                                hpvc_leaf_dec_deidx_fx(v_tmp, abs_s(hdrVals[i]), nFlatLeaves[i], LS_ind_sub, (UWord32)hpvcTreeEnumCfg->flatLeafIdx[i]);
                                FOR (ii_ = 0; ii_ < nFlatLeaves[i]; ii_++)
                                {
                                    L_x[ii_] = L_deposit_l(v_tmp[ii_]);
                                }
                            }
                            accPos = add(accPos, nFlatLeaves[i]);
                        }
                    } /* Kp > 0 */

                    treeCount = add(treeCount, 1);

                    a1_i += Np;
                    b1_i += Np;

                    /* update TCX context based on last 4 coeffs of X */
                    c = post_HPVC_update_TCX_context(&x[sub(a1_i, 2 * 2)]);
                }
            }
            ELSE /* Np == 2 — legacy TCX 2-tuple decode (inline) */
            {
                t = add(c, rateFlag);
                if (sub(k, nt_half) > 0)
                {
                    t = add(t, 1 << NBITS_CONTEXT);
                }
                r = ac_decode_fx(&st, ari_spec_lookup[t]
#ifdef CR13_B_FIX_PC_BINS
                                 , &bp, bp_side, mask_side, k, 1
#endif
                );
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                     ari_spec_cumfreq[ari_spec_lookup[t]][r],
                                     ari_spec_freq[ari_spec_lookup[t]][r]) != 0)
                {
                    GOTO ber_detect;
                }
                IF (r == 0)
                {
                    x[a1_i] = 0; move16();
                    x[b1_i] = 0; move16();
                    c = add(shl_pos(s_and(c, 0xf), 4), 1);
                }
                ELSE IF (sub(r, VAL_ESC) < 0)
                {
                    a = s_and(r, 0x3);
                    b = shr_pos(r, 2);
                    c = add(shl_pos(s_and(c, 0xf), 4), add(add(a, b), 1));
                    IF (a > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) { GOTO ber_detect; };
                        if (read_bit(ptr, bp_side, mask_side) != 0) a = negate(a);
                    }
                    x[a1_i] = a; move16();
                    IF (b > 0)
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0) { GOTO ber_detect; };
                        if (read_bit(ptr, bp_side, mask_side) != 0) b = negate(b);
                    }
                    x[b1_i] = b; move16();
                }
                ELSE
                {
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;};
                    a = read_bit(ptr, bp_side, mask_side);
                    if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) { GOTO ber_detect;};
                    b = read_bit(ptr, bp_side, mask_side);
                    r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[1]]
#ifdef CR13_B_FIX_PC_BINS
                                     , &bp, bp_side, mask_side, k, 1
#endif
                    );
                    if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                         ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[1]]][r],
                                         ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[1]]][r]) != 0)
                    {
                        GOTO ber_detect;
                    }
                    IF (sub(r, VAL_ESC) < 0)
                    {
                        a1 = s_and(r, 0x3);
                        b1 = shr_pos(r, 2);
                        a = add(shl_pos(a1, 1), a);
                        b = add(shl_pos(b1, 1), b);
                        IF (a > 0)
                        {
                            if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;}
                            if (read_bit(ptr, bp_side, mask_side) != 0) a = negate(a);
                        }
                        x[a1_i] = a; move16();
                        IF (b > 0)
                        {
                            if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0) {GOTO ber_detect;}
                            if (read_bit(ptr, bp_side, mask_side) != 0) b = negate(b);
                        }
                        x[b1_i] = b; move16();
                        c = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1, b1), 1), 1));
                    }
                    ELSE
                    {
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;}
                        a = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), a);
                        if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;}
                        b = add(shl_pos(read_bit(ptr, bp_side, mask_side), 1), b);
                        FOR (lev = 2; lev < max_lev; lev++)
                        {
                            esc_nb = s_min(lev, 3);
                            r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[esc_nb]]
#ifdef CR13_B_FIX_PC_BINS
                                             , &bp, bp_side, mask_side, k, 1
#endif
                            );
                            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                                 ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r],
                                                 ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r]) != 0)
                            {
                                GOTO ber_detect;
                            }
                            IF (sub(r, VAL_ESC) < 0) { BREAK; };
                            if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) { GOTO ber_detect;};
#  ifdef ENABLE_HR_MODE
                            a = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  else
                            a = add(shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  endif
                            if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;};
#  ifdef ENABLE_HR_MODE
                            b = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  else
                            b = add(shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  endif
                        }
                        IF (sub(lev, max_lev) == 0) {GOTO ber_detect;}
                        b1 = shr_pos(r, 2);
                        a1 = s_and(r, 0x3);
#  ifdef ENABLE_HR_MODE
                        a = L_add(L_shl(a1, lev), a);
                        b = L_add(L_shl(b1, lev), b);
#  else
                        a = add(shl(a1, lev), a);
                        b = add(shl(b1, lev), b);
#  endif
                        IF (a > 0)
                        {
                            if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0) {GOTO ber_detect;};
                            if (read_bit(ptr, bp_side, mask_side) != 0)
                            {
#  ifdef ENABLE_HR_MODE
                                a = L_negate(a);
#  else
                                a = negate(a);
#  endif
                            }
                        }
                        x[a1_i] = a; move16();
                        IF (b > 0)
                        {
                            if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0) {GOTO ber_detect;};
                            if (read_bit(ptr, bp_side, mask_side) != 0)
                            {
#  ifdef ENABLE_HR_MODE
                                b = L_negate(b);
#  else
                                b = negate(b);
#  endif
                            }
                        }
                        x[b1_i] = b; move16();
                        c = add(shl_pos(s_and(c, 0xf), 4), add(esc_nb, 12));
                    }
                }
                test();
                test();
                IF ((sub(sub(bp, *bp_side), 3) > 0 && sub(st.pc.c_bp, st.pc.c_bp_side) == 0) || st.BER_detect > 0)
                {
                    GOTO ber_detect;
                }
                a1_i += 2;
                b1_i += 2;
            } /* Np == 2 */
        } /* FOR k joint loop */
    } /* IF mode>=0 (HPVC) */
#endif /* LL_INCL_HPVC_ARICODEC_DEC */
}
ELSE
{
/*Main Loop through the 2-tuples*/
FOR (k = 0; k < lastnz; k += 2)
{

/* Get context */
t = add(c, rateFlag);
if (sub(k, nt_half) > 0)
{
t = add(t, 1 << NBITS_CONTEXT);
}

            r = ac_decode_fx(&st, ari_spec_lookup[t]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
            );
            if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st, ari_spec_cumfreq[ari_spec_lookup[t]][r],
                                 ari_spec_freq[ari_spec_lookup[t]][r]) != 0)
            {
                GOTO ber_detect;
            }

            IF (r == 0)
            {
                x[a1_i] = 0;
                move16();
                x[b1_i] = 0;
                move16();
                c = add(shl_pos(s_and(c, 0xf), 4), 1);
            }
            ELSE IF (sub(r, VAL_ESC) < 0)
            {
                a = s_and(r, 0x3);
                b = shr_pos(r, 2);
                c = add(shl_pos(s_and(c, 0xf), 4), add(add(a, b), 1));
                IF (a > 0)
                {
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        a = negate(a);
                    }
                }
                x[a1_i] = a;
                move16();
                IF (b > 0)
                {
if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
                }
                x[b1_i] = b;
                move16();
            }
            ELSE
            {
                r = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[1]]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
                );
                if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                     ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[1]]][r],
                                     ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[1]]][r]) != 0)
                {
                    GOTO ber_detect;
                }
                IF (sub(r, VAL_ESC) < 0)
                {
                    a1 = s_and(r, 0x3);
                    b1 = shr_pos(r, 2);
                    a  = shl_pos(a1, 1);
                    b  = shl_pos(b1, 1);
                    IF (a > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            a = negate(a);
                        }
                    }
                    x[a1_i] = a;
                    move16();
                    IF (b > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
                            b = negate(b);
                        }
                    }
                    x[b1_i] = b;
                    move16();
                    c                    = add(shl_pos(s_and(c, 0xf), 4), add(shl_pos(add(a1, b1), 1), 1));
                    lsb_ind[lsb_ind_c++] = k;
                    move16();
                }
                ELSE
                {
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
a = shl_pos(read_bit(ptr, bp_side, mask_side), 1);
if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
{
  GOTO ber_detect;
}
                    b = shl_pos(read_bit(ptr, bp_side, mask_side), 1);
                    FOR (lev = 2; lev < max_lev; lev++)
                    {
                        esc_nb = s_min(lev, 3);
                        r      = ac_decode_fx(&st, ari_spec_lookup[t + Tab_esc_nb[esc_nb]]
#ifdef CR13_B_FIX_PC_BINS
                             ,
                             &bp,
                             bp_side,
                             mask_side,
                             k,
                             1
#endif
                        );
                        if (ac_dec_update_fx(ptr, &bp, bp_side, mask_side, k, &st,
                                             ari_spec_cumfreq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r],
                                             ari_spec_freq[ari_spec_lookup[t + Tab_esc_nb[esc_nb]]][r]) != 0)
                        {
                            GOTO ber_detect;
                        }
                        IF (sub(r, VAL_ESC) < 0)
                        {
                            BREAK;
                        }
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
#  ifdef ENABLE_HR_MODE
                        a = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  else
                        a = add(shl(read_bit(ptr, bp_side, mask_side), lev), a);
#  endif
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
#  ifdef ENABLE_HR_MODE
                        b = L_add(L_shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  else
                        b = add(shl(read_bit(ptr, bp_side, mask_side), lev), b);
#  endif
                    }
                    /* check for bitflip */
                    IF (sub(lev, max_lev) == 0)
                    {
                        GOTO ber_detect;
                    }

                    b1 = shr_pos(r, 2);
                    a1 = s_and(r, 0x3);
#  ifdef ENABLE_HR_MODE
                    a  = L_add(L_shl(a1, lev), a);
                    b  = L_add(L_shl(b1, lev), b);
#  else
                    a = add(shl(a1, lev), a);
                    b = add(shl(b1, lev), b);
#  endif
                    IF (a > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, a1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
#  ifdef ENABLE_HR_MODE
                            a = L_negate(a);
#  else
                            a = negate(a);
#  endif
                        }
                    }
                    x[a1_i] = a;
                    move16();
                    IF (b > 0)
                    {
  if (check_pc_bytes(&bp, bp_side, mask_side, b1_i, 0, &st.pc) != 0)
  {
      GOTO ber_detect;
  }
                        if (read_bit(ptr, bp_side, mask_side) != 0)
                        {
#  ifdef ENABLE_HR_MODE
                            b = L_negate(b);
#  else
                            b = negate(b);
#  endif
                        }
                    }
                    x[b1_i] = b;
                    move16();
                    c                    = add(shl_pos(s_and(c, 0xf), 4), add(esc_nb, 12));
                    lsb_ind[lsb_ind_c++] = k;
                    move16();
                }
            }

            test();
            test();
            IF ((sub(sub(bp, *bp_side), 3) > 0 && sub(st.pc.c_bp, st.pc.c_bp_side) == 0) || st.BER_detect > 0)
{
GOTO ber_detect;
}

a1_i += 2;
b1_i += 2;
}
}

IF (L_spec > k)
{
basop_memset(&x[k], 0, (L_spec - k) * sizeof(*x));
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF (ll_flag)
    {
        process_lsb_add(x, det_curve,  L_spec,  x);
    }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
    nbits_side = L_sub(nbbits, L_add(L_shl_pos((Word32)*bp_side, 3), L_sub(norm_s(*mask_side), 6)));
#else
    nbits_side = sub(nbbits, add(shl_pos(*bp_side, 3), sub(norm_s(*mask_side), 6)));
#endif
    extra_bits = sub(norm_ul(st.ac_range_fx), 6);
#ifdef CR14_A_ADD_LOSSLESS_MODE
    nbits_ari  = L_shl_pos(L_sub(bp, 3), 3);
#else
    nbits_ari  = shl_pos(sub(bp, 3), 3);
#endif
IF (mode != 1)
{
IF (st.pc.c_bp == 0)
{
#ifdef CR14_A_ADD_LOSSLESS_MODE
nbits_ari = L_shl_pos(L_sub(L_sub(bp, st.pc.bytes), 3), 3);
#else
nbits_ari = shl_pos(sub(sub(bp, st.pc.bytes), 3), 3);
#endif
}
ELSE
{
#ifdef CR14_A_ADD_LOSSLESS_MODE
nbits_ari = L_shl_pos(L_add(bp, L_sub(L_sub(st.pc.b_left, st.pc.bytes), 3)), 3);
#else
nbits_ari = shl_pos(add(bp, sub(sub(st.pc.b_left, st.pc.bytes), 3)), 3);
#endif
}

        IF (st.pc.c_bp_side != 0)
        {
            nbits_side = sub(add(sub(nbbits, shl_pos(st.pc.b_left, 3)), shl_pos(sub(st.pc.bytes, *bp_side), 3)),
                             sub(norm_s(*mask_side), 6));
        }
    }

#ifdef CR14_A_ADD_LOSSLESS_MODE
    n = L_sub(nbbits, L_add(nbits_ari, L_add(extra_bits, nbits_side)));
#else
    n = sub(nbbits, add(nbits_ari, add(extra_bits, nbits_side)));
#endif
    move16();

IF (n < 0)
{
GOTO ber_detect;
}

    IF (lsbMode == 0)
    {
        *resBits = n;
        move16();
        i=0;

#  ifdef CR9_C_ADD_1p25MS
        Counter l, lMax = 1;
        if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
            lMax = 3;
        }
        FOR (l = 0; l < lMax; l++)
        {
#  endif
#ifdef ENABLE_HR_MODE
        FOR (k = 0; k < L_spec; k++)
        {
            IF (x[k] != 0)
            {
                IF (n == 0)
                {
                    BREAK;
                }
                if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                {
                    GOTO ber_detect_res;
                }
                i1 = shr(i, RESBITS_PACK_SHIFT);
                i2 = s_and(i, RESBITS_PACK_MASK);
                resbit = read_bit(ptr, bp_side, mask_side);
                if (resbit)
                {
                    resQdata[i1] = (UWord8) s_or(resQdata[i1], shl(1, i2));
                }
                i = add(i, 1);
                move16();
#ifdef CR14_A_ADD_LOSSLESS_MODE
                n = L_sub(n, 1);
#else
                n = sub(n, 1);
#endif
            }
        }
#else
        FOR (k = 0; k < L_spec; k++)
        {
            IF (x[k] != 0)
            {
                IF (n == 0)
                {
                    BREAK;
                }

                if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                {
                    GOTO ber_detect_res;
                }

                *resQdata++ = (UWord8) read_bit(ptr, bp_side, mask_side);
                move16();
                n = sub(n, 1);
            }
        }
#endif
#  ifdef CR9_C_ADD_1p25MS
        }
#  endif
#  ifdef ENABLE_HR_MODE
        if (hrmode)
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            Word32 idx_len     = L_sub(*resBits, n); /* Number of nonzero bits */
            Word32 idx_len_lim = idx_len * EXT_RES_ITER_MAX;
#else
            Word16 idx_len     = sub(*resBits, n); /* Number of nonzero bits */
            Word16 idx_len_lim = idx_len * EXT_RES_ITER_MAX;
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
            Word32 res_bits_hrmode = 0;
            IF (ll_flag)
            {
                res_bits_hrmode = *resBits - idx_len;
            }
            ELSE
            {
                res_bits_hrmode = s_min( idx_len_lim, *resBits ) - idx_len;
            }
#else
            Word16 res_bits_hrmode = s_min(idx_len_lim, *resBits) - idx_len;
#endif
            /* idx_len bits have been read in the previous loop */

            for (k = 0; k < res_bits_hrmode; k++)
            {
                if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
                {
                    GOTO ber_detect_res;
                }

#ifdef CR14_A_ADD_LOSSLESS_MODE
                i1 = L_shr(i, RESBITS_PACK_SHIFT);
                i2 = L_and(i, RESBITS_PACK_MASK);
#else
                i1 = shr(i, RESBITS_PACK_SHIFT);
                i2 = s_and(i, RESBITS_PACK_MASK);
#endif
                resbit = read_bit(ptr, bp_side, mask_side);
                if (resbit)
                {
                    resQdata[i1] = (UWord8) s_or(resQdata[i1], shl(1, i2));
                }
#ifdef CR14_A_ADD_LOSSLESS_MODE
                i = L_add(i, 1);
                move16();
                n = L_sub(n, 1);
#else
                i = add(i, 1);
                move16();
                n = sub(n, 1);
#endif
            }
        }
#  endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
        *resBits = L_sub(*resBits, n);
#else
        *resBits = sub(*resBits, n);
#endif
    }
    ELSE
    {
        *resBits = 0;
        FOR (k = 0; k < lsb_ind_c; k++)
        {
            a = x[lsb_ind[k]];
            move16();
            IF (n == 0)
            {
                BREAK;
            }
if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
{
GOTO ber_detect_res;
}
            tmp = read_bit(ptr, bp_side, mask_side);
#ifdef CR14_A_ADD_LOSSLESS_MODE
            n   = L_sub(n, 1);
#else
            n   = sub(n, 1);
#endif
            IF (tmp > 0)
            {
#  ifdef ENABLE_HR_MODE
                if (a > 0)
                {
                    a = L_add(a, 1);
                }
                if (a < 0)
                {
                    a = L_sub(a, 1);
                }
#  else
                if (a > 0)
                {
                    a = add(a, 1);
                }
                if (a < 0)
                {
                    a = sub(a, 1);
                }
#  endif
                IF (a == 0)
                {
                    IF (n == 0)
                    {
                        BREAK;
                    }
                    a = 1;
if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
{
  GOTO ber_detect_res;
}
if (read_bit(ptr, bp_side, mask_side) != 0)
{
  a = negate(a);
}
#ifdef CR14_A_ADD_LOSSLESS_MODE
n = L_sub(n, 1);
#else
n = sub(n, 1);
#endif
}
}

            x[lsb_ind[k]] = a;
            move16();
            b = x[lsb_ind[k] + 1];
            move16();
            IF (n == 0)
            {
                BREAK;
            }
if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
{
GOTO ber_detect_res;
}
            tmp = read_bit(ptr, bp_side, mask_side);

#ifdef CR14_A_ADD_LOSSLESS_MODE
            n   = L_sub(n, 1);
#else
            n   = sub(n, 1);
#endif

            IF (tmp > 0)
            {
#  ifdef ENABLE_HR_MODE
                if (b > 0)
                {
                    b = L_add(b, 1);
                }
                if (b < 0)
                {
                    b = L_sub(b, 1);
                }
#  else
                if (b > 0)
                {
                    b = add(b, 1);
                }
                if (b < 0)
                {
                    b = sub(b, 1);
                }
#  endif
                IF (b == 0)
                {
                    IF (n == 0)
                    {
                        BREAK;
                    }
                    b = 1;
if (check_pc_bytes(&bp, bp_side, mask_side, st.pc.inv_bin, 0, &st.pc) != 0)
{
  GOTO ber_detect_res;
}
                    if (read_bit(ptr, bp_side, mask_side) != 0)
                    {
                        b = negate(b);
                    }
#ifdef CR14_A_ADD_LOSSLESS_MODE
                    n = L_sub(n, 1);
#else
                    n = sub(n, 1);
#endif
                }
            }
            x[lsb_ind[k] + 1] = b;
            move16();
        }
    }

/* Noise Filling seed */
    tmp32 = L_deposit_l(0);

#ifdef CR14_A_ADD_LOSSLESS_MODE
    IF( ll_flag == 0 )
#endif
    {
      FOR (i = 0; i < L_spec; i++)
      {
#     ifdef ENABLE_HR_MODE
          tmp32 = L_mac0(tmp32, L_and(L_abs(x[i]), 32767), i);
#     else
          tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
#     endif
      }
    }
    *nf_seed = extract_l(tmp32);
    move16();

    /* Detect zero frame */
    test();
    test();
    test();
    test();
    IF (sub(lastnz, 2) == 0 && sub(x[0], 0) == 0 && sub(x[1], 0) == 0 && sub(gg_idx, 0) == 0 && sub(fac_ns_idx, 7) == 0)
    {
        *zero_frame = 1;
        move16();
    }
    ELSE
    {
        *zero_frame = 0;
        move16();
    }

IF (mode == 1)
    {
        IF (st.pc.bytes > 0)
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
            IF (L_sub(st.pc.b_left, L_shr_pos(nbbits, 3)) > 0)
#else
            IF (sub(st.pc.b_left, shr_pos(nbbits, 3)) > 0)
#endif
            {
                *b_left = sub(*bp_side, st.pc.bytes);
            }
        }
    }
IF (mode == 2)
{
IF (st.pc.bytes > 0)
{
#ifdef CR14_A_ADD_LOSSLESS_MODE
IF (L_sub(st.pc.b_left, L_shr_pos(nbbits,3)) > 0)
#else
IF (sub(st.pc.b_left, shr_pos(nbbits,3)) > 0)
#endif
{
*b_left = *bp_side;
}
}
}

IF (sub(*bfi, 2) == 0)
{
IF (sub(*spec_inv_idx, L_spec) == 0)
{
*bfi = 0;
}
}
GOTO bail;

/* goto for bit error handling */
ber_detect:
    *bfi = 1;
    move16();
    *b_left = st.pc.b_left;
    move16();
    test();
    IF (st.pc.inv_bin > 0 && sub(st.pc.inv_bin, L_spec) <= 0)
    {
        *spec_inv_idx = st.pc.inv_bin;
        move16();
        *bfi = 2;
        move16();
        *resBits = 0;
        move16();
        *zero_frame = 0;
        move16();
        /* Noise Filling seed */
        tmp32 = L_deposit_l(0);
        FOR (i = 0; i < *spec_inv_idx; i++)
        {
            tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
        }
        *nf_seed = extract_l(tmp32);
        move16();
    }
GOTO bail;

/* goto for bit error handling in residual signal */
ber_detect_res:
    *b_left = st.pc.b_left;
    move16();
    *resBits = 0;
    move16();
    *bfi = 0;
    move16();
    *zero_frame = 0;
    move16();
    /* Noise Filling seed */
    tmp32 = L_deposit_l(0);
    FOR (i = 0; i < *spec_inv_idx; i++)
    {
        tmp32 = L_mac0(tmp32, abs_s(x[i]), i);
    }
    *nf_seed = extract_l(tmp32);
    move16();
    GOTO bail;

/* goto, because of dynmem out */
bail:
Dyn_Mem_Deluxe_Out();

    lsb_ind = (Word16*) lc3_scratch_pop( scratch, lsb_ind );
}


void processAriDecoderScaling_fx(
#  ifdef ENABLE_HR_MODE
    Word32 *datain,
#  else
    Word16 *data16,
#  endif
    Word16 dataLen, Word32 *data32, Word16 *data_e)
{
    Counter i;

#  ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(Word16 shift; Word32 tmp, x_min, x_max;);
#  else
    Dyn_Mem_Deluxe_In(Word16 shift; Word16 tmp, x_min, x_max;);
#  endif


#ifdef ENABLE_HR_MODE
    x_max = 0;
    move32();
    x_min = 0;
    move32();
#else
    x_max = 0;
    move16();
    x_min = 0;
    move16();
#endif

    FOR (i = 0; i < dataLen; i++)
    {
#ifdef ENABLE_HR_MODE
        if (datain[i] > 0)
            x_max = L_max(x_max, datain[i]);
        if (datain[i] < 0)
            x_min = L_min(x_min, datain[i]);
#else
        if (data16[i] > 0)
            x_max = s_max(x_max, data16[i]);
        if (data16[i] < 0)
            x_min = s_min(x_min, data16[i]);
#endif
    }

#ifdef ENABLE_HR_MODE
    tmp   = L_max(x_max, L_negate(x_min));
    shift = norm_l(tmp);
    if (tmp == 0)
    {
        shift = 31;
        move32();
    }
#else
    tmp = s_max(x_max, negate(x_min));
    shift = norm_s(tmp);
    if (tmp == 0)
    {
        shift = 15;
        move16();
    }
#endif

    FOR (i = 0; i < dataLen; i++)
    {
#ifdef ENABLE_HR_MODE
        data32[i] = L_shl_pos(datain[i], shift);
#else
        data32[i] = L_shl_pos(L_deposit_h(data16[i]), shift);
#endif
    }

#ifdef ENABLE_HR_MODE
    *data_e = sub(31, shift);
    move16();
#else
    *data_e = sub(15, shift);
    move16();
#endif

    Dyn_Mem_Deluxe_Out();
}


/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/

static UWord32 UL_addNs24(UWord32 UL_var1, UWord32 UL_var2, UWord16 *wrap)
{
    return UL_lshr(UL_addNs(UL_lshl(UL_var1, 8), UL_lshl(UL_var2, 8), wrap), 8);
}
#ifdef ENABLE_HR_MODE
Word16 find_last_nz_pair(const Word32 x[], Word16 length)
#else
Word16 find_last_nz_pair(const Word16 x[], Word16 length)
#endif
{
    Dyn_Mem_Deluxe_In(Word16 last_nz, lobs[4]; Counter stage, i;);

    lobs[0] = 4;
    move16();
    lobs[1] = shr_pos(length, 1); /* length/2 */
    move16();
    lobs[2] = add(lobs[1], shr_pos(length, 2));
    move16();
    lobs[3] = add(lobs[2], shr_pos(length, 3));
    move16();

    last_nz = 0;
    move16();
    i = length;
    move16();
    FOR (stage = 3; stage >= 0; --stage)
    {
        /* unmapped kernel */
        FOR (; i >= lobs[stage]; i -= 2)
        {
            if (x[i - 2] != 0)
            {
                last_nz = s_max(last_nz, i);
            }
            if (x[i - 1] != 0)
            {
                last_nz = s_max(last_nz, i);
            }
        }
        IF (last_nz > 0)
        {
            BREAK;
        }
    }

    Dyn_Mem_Deluxe_Out();
    return s_max(last_nz, 2);
}


void write_bit_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 bit)
{
    if (bit > 0)
    {
        ptr[*bp] = (UWord8)s_or((Word16)ptr[*bp], *mask);
        move16();
    }
    *mask = lshl_pos(*mask, 1);
    move16();
    if (sub(*mask, 0x100) == 0)
    {
        *mask = 1;
        move16();
    }
    if (sub(*mask, 1) == 0)
    {
        *bp = sub(*bp, 1);
        move16();
    }
}


void write_indice_backward(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 indice, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(Counter k; Word16 bit;);

    FOR (k = 0; k < numbits; k++)
    {
        bit = s_and(indice, 1);
        write_bit_backward(ptr, bp, mask, bit);
        indice = lshr(indice, 1);
    }

    Dyn_Mem_Deluxe_Out();
}


static void write_indice_forward(UWord8 *ptr, Word16 bp, Word16 indice, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(Counter k; Word16 bit, mask, tmp;);

    tmp = (Word16)ptr[bp];
    move16();
    mask = 0x80;
    move16();
    FOR (k = 0; k < numbits; k++)
    {
        bit = s_and(indice, mask);
        tmp = s_or(tmp, mask);
        if (bit == 0)
        {
            tmp = sub(tmp, mask);
        }
        mask = lshr(mask, 1);
    }
    ptr[bp] = (UWord8)tmp;
    move16();

    Dyn_Mem_Deluxe_Out();
}

static void ac_enc_init_fx(Encoder_State_fx *st_fx) /* i/o: Encoder state       */
{
    st_fx->ac_low_fx = L_deposit_l(0);
    move32();
    st_fx->ac_range_fx = 0x00ffffff;
    move32();
    st_fx->ac_cache_fx = -1;
    move16();
    st_fx->ac_carry_fx = 0;
    move16();
    st_fx->ac_carry_count_fx = 0;
    move16();
}

static void ac_enc_shift_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx) /* i/o: Encoder state */
{
    test();
    L_sub(0, 0); /* For comparision in if */
    IF (st_fx->ac_low_fx < (0x00ff0000UL) || sub(st_fx->ac_carry_fx, 1) == 0)
    {
        IF (st_fx->ac_cache_fx >= 0)
        {
            ptr[(*bp)++] = (UWord8)add(st_fx->ac_cache_fx, st_fx->ac_carry_fx);
            move16();
        }

        WHILE (st_fx->ac_carry_count_fx > 0)
        {
            ptr[(*bp)++] = (UWord8)s_and(add(st_fx->ac_carry_fx, 0xff), 255);
            move16();
            st_fx->ac_carry_count_fx = sub(st_fx->ac_carry_count_fx, 1);
            move16();
        }

        st_fx->ac_cache_fx = u_extract_l(UL_lshr_pos(st_fx->ac_low_fx, 16));
        move16();
        st_fx->ac_carry_fx = 0;
        move16();
    }
    ELSE
    {
        st_fx->ac_carry_count_fx = add(st_fx->ac_carry_count_fx, 1);
        move16();
    }
    st_fx->ac_low_fx = UL_and(UL_lshl_pos(st_fx->ac_low_fx, 8), 0x00ffffff);
    move32();
}

static void ac_encode_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx, /* i/o: Encoder state */
                                       UWord32 cum_freq, /* i  : Cumulative frequency up to symbol   */
                                       UWord32 sym_freq) /* i  : Symbol probability                  */
{
    Dyn_Mem_Deluxe_In(UWord32 r, tmp; UWord16 carry;);

    r   = UL_lshr_pos(st_fx->ac_range_fx, 10);
    tmp = UL_Mpy_32_32(r, cum_freq);

    assert(r < (1U << 24));
    assert(cum_freq < (1U << 24));
    assert(tmp < (1U << 24));
    assert(st_fx->ac_low_fx < (1U << 24));
    st_fx->ac_low_fx = UL_addNs24(st_fx->ac_low_fx, tmp, &carry);
    move32();

    if (carry != 0)
    {
        st_fx->ac_carry_fx = carry;
        move16();
    }

    st_fx->ac_range_fx = UL_Mpy_32_32(r, sym_freq);
    move32();

    assert(cum_freq < (1U << 24));
    assert(st_fx->ac_range_fx < (1U << 24));
    WHILE (st_fx->ac_range_fx < (1U << 16))
    {
        L_sub(0, 0); /* Comparison in while */
        st_fx->ac_range_fx = UL_lshl_pos(st_fx->ac_range_fx, 8);
        move32();

        assert(st_fx->ac_range_fx < (1U << 24));

        ac_enc_shift_fx(ptr, bp, st_fx);
    }

    Dyn_Mem_Deluxe_Out();
}

static Word16 ac_enc_finish_fx(UWord8 *ptr, Word16 *bp, Encoder_State_fx *st_fx) /* i/o: Encoder state */
{
    Dyn_Mem_Deluxe_In(UWord32 val, mask, high; Word16 bits; UWord16 over1, over2;);

    /*bits = 24 - log2_i(st->ac_range); */
    bits = sub(norm_ul(st_fx->ac_range_fx), 7);

    mask = UL_lshr(0x00ffffff, bits);

    val  = UL_addNs24(st_fx->ac_low_fx, mask, &over1);
    high = UL_addNs24(st_fx->ac_low_fx, st_fx->ac_range_fx, &over2);

    L_xor(0, 0);    /* For bit not */
    UL_and(1U, 1U); /* added counters */
    val = L_and(val, (~mask) & 0x00ffffff);

    L_xor(0, 0); /* For bit not */
    IF ((L_xor(over1, over2)) == 0)
    {
        L_sub(0, 0); /* For comparision in if */
        IF (UL_addNsD(val, mask) >= high)
        {
            bits = add(bits, 1);
            mask = UL_lshr_pos(mask, 1);
            val  = UL_and(UL_addNsD(st_fx->ac_low_fx, mask), (~mask) & 0x00ffffff);
            L_xor(0, 0);
            UL_and(1, 1); /* For bit not , mask */
        }

        if (val < st_fx->ac_low_fx)
        {
            st_fx->ac_carry_fx = 1;
            move16();
        }
    }

    st_fx->ac_low_fx = val;
    move32();

    FOR (; bits > 0; bits -= 8)
    {
        ac_enc_shift_fx(ptr, bp, st_fx);
    }
    bits = add(bits, 8);

    assert(st_fx->ac_carry_fx == 0);

    IF (st_fx->ac_carry_count_fx > 0)
    {
        ptr[(*bp)++] = (UWord8)st_fx->ac_cache_fx;
        move16();

        FOR (; st_fx->ac_carry_count_fx > 1; st_fx->ac_carry_count_fx--)
        {
            ptr[(*bp)++] = 0xff;
            move16();
        }
        write_indice_forward(ptr, *bp, lshr(0xff, sub(8, bits)), bits);
    }
    ELSE
    {
        write_indice_forward(ptr, *bp, st_fx->ac_cache_fx, bits);
    }

    Dyn_Mem_Deluxe_Out();
    return bits;
}


Word16 read_bit(UWord8 *ptr, Word16 *bp, Word16 *mask)
{
    Dyn_Mem_Deluxe_In(Word16 bit;);

    bit = 0;
    move16();
    if (s_and((Word16)ptr[*bp], *mask) > 0)
    {
        bit = 1;
        move16();
    }
    *mask = lshl_pos(*mask, 1);
    move16();
    if (sub(*mask, 0x100) == 0)
    {
        *mask = 1;
        move16();
    }
    if (sub(*mask, 1) == 0)
    {
        *bp = sub(*bp, 1);
        move16();
    }

    Dyn_Mem_Deluxe_Out();
    return bit;
}


static void ac_dec_init_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side,
                                         Decoder_State_fx *st_fx) /* i/o: Decoder State */
{
    Dyn_Mem_Deluxe_In(Counter i;);


    st_fx->ac_low_fx = L_deposit_l(0);
    move32();

    st_fx->ac_range_fx = 0x00ffffff;
    move32();
    FOR (i = 0; i < 3; i++)
    {
        if (check_pc_bytes(bp, bp_side, mask_side, 0, 1, &st_fx->pc) != 0)
        {
            Dyn_Mem_Deluxe_Out();
            return;
        }
        st_fx->ac_low_fx = UL_addNsD(UL_lshl_pos(st_fx->ac_low_fx, 8), UL_deposit_l((Word16)ptr[(*bp)++]));
        move32();
        assert(st_fx->ac_low_fx < (1U << 24));
    }

    st_fx->BER_detect = 0;
    move16();

    Dyn_Mem_Deluxe_Out();
}

#ifdef CR13_B_FIX_PC_BINS
static Word16 pc_check_bytes_ac_decode_fx(
    Word16 *bp,
    Word16 *bp_side,
    Word16 *mask_side,
    Word16 cur_bin,
    Word16 from_left,
    Pc_State_fx *pc
)
{
    Word16 bp_local, bp_side_local;

    if (pc->bytes > 0)
    {
        bp_local = *bp;
        bp_side_local = *bp_side;

        if (from_left)
        {
            if (*mask_side == 1)
            {
                bp_side_local = add(bp_side_local, 1);
            }
        }
        else
        {
            bp_local = sub(bp_local, 1);
        }

        if (!pc->enc && pc->b_right > -1)
        {
            if (pc->bfi == 2)
            {
                if (pc->c_bp && bp_local > pc->be_bp_left)
                {
                    pc->inv_bin = cur_bin;
                    return 1;
                }
            }
        }
    }

    return 0;
}
#endif

/* o  : Decoded cumulative frequency    */
static Word16 ac_decode_fx(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                         Word16            pki
#ifdef CR13_B_FIX_PC_BINS
                                         ,
                                         Word16 *bp,
                                         Word16 *bp_side,
                                         Word16 *mask_side,
                                         Word16 cur_bin,
                                         Word16 from_left
#endif
                                         )
{
    Dyn_Mem_Deluxe_In(UWord16 sgn; Word16 val, r;);

#ifdef CR13_B_FIX_PC_BINS
    IF (pc_check_bytes_ac_decode_fx(bp, bp_side, mask_side, cur_bin, from_left, &st_fx->pc))
    {
        st_fx->BER_detect = 1;
        return 0;
    }
#endif

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10);
    move32();
    val = 0;
    move16();

    r = add(val, 8);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][r]), &sgn);
    IF (sgn == 0)
    {
        val = r;
        move16();
        IF (sub(val, 15) == 0)
        {
            UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ari_spec_cumfreq[pki][16]), &sgn);
            if (sgn == 0)
            {
                val = 16;
                move16();
            }
            UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
            if (sgn == 0)
            {
                st_fx->BER_detect = 1;
                move16();
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

/* o  : Decoded cumulative frequency    */
static Word16 ac_decode_tns_order(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                                Word16            enable_lpc_weighting)
{
    Dyn_Mem_Deluxe_In(UWord16 sgn; Word16 val, r;);

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10);
    move32();
    val = 0;
    move16();

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_order_cumfreq[enable_lpc_weighting][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
    if (sgn == 0)
    {
        st_fx->BER_detect = 1;
        move16();
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

/* o  : Decoded cumulative frequency    */
static Word16 ac_decode_tns_coef(Decoder_State_fx *st_fx, /* i/o: Decoder State                   */
                                               Word16            pki)
{
    Dyn_Mem_Deluxe_In(UWord16 sgn; Word16 val, r;);

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10);
    move32();
    val = 0;
    move16();

    r = add(val, 8);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 4);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 2);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
    }

    r = add(val, 1);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][r]), &sgn);
    if (sgn == 0)
    {
        val = r;
        move16();
        IF (sub(val, 15) == 0)
        {
            UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, ac_tns_coef_cumfreq[pki][16]), &sgn);
            if (sgn == 0)
            {
                val = 16;
                move16();
            }
            UL_subNs(st_fx->ac_low_fx, UL_lshl(st_fx->ac_help_fx, 10), &sgn);
            if (sgn == 0)
            {
                st_fx->BER_detect = 1;
                move16();
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return val;
}

static Word16 ac_dec_update_fx(UWord8 *ptr, Word16 *bp, Word16 *bp_side, Word16 *mask_side,
                                             Word16 cur_bin, Decoder_State_fx *st_fx, /* i/o: Decoder State */
                                             UWord32 cum_freq,                        /* i  : Cumulative frequency    */
                                             UWord32 sym_freq                         /* i  : Symbol frequency        */
)
{
    UWord32 UL_tmp;


    assert(st_fx->ac_help_fx < (1U << 24));
    assert(cum_freq < (1U << 24));

    UL_tmp = UL_Mpy_32_32(cum_freq, st_fx->ac_help_fx);
    assert(UL_tmp < (1U << 24));

    st_fx->ac_low_fx = UL_subNsD(st_fx->ac_low_fx, UL_tmp);
    move32(); /*0+0*/
    assert(st_fx->ac_low_fx < (1U << 24));

    st_fx->ac_range_fx = UL_Mpy_32_32(st_fx->ac_help_fx, sym_freq);
    move32();

    assert(st_fx->ac_range_fx < (1U << 24));
    /* updated to 16 from 24 */
    WHILE (st_fx->ac_range_fx < (1U << 16))
    {
        L_sub(0, 0); /* For comparision in while*/

        st_fx->ac_low_fx =
            UL_and(st_fx->ac_low_fx, 0x0000ffFF); /*  make sure upshift doe not lead to more than 24 bits */
        assert(st_fx->ac_low_fx < 1U << 16);

        if (check_pc_bytes(bp, bp_side, mask_side, cur_bin, 1, &st_fx->pc) != 0)
            return 1;

        /*shift in 8 bits */
        st_fx->ac_low_fx = UL_addNsD(UL_lshl_pos(st_fx->ac_low_fx, 8), UL_deposit_l((Word16)ptr[(*bp)++]));
        move32();

        assert(st_fx->ac_low_fx < (1U << 24));
        st_fx->ac_range_fx = UL_lshl_pos(st_fx->ac_range_fx, 8);
        move32();
        assert(st_fx->ac_range_fx < (1U << 24));
    }
    return 0;
}

static void pc_init_fx(Word16 n_pc, Word16 numbytes, Word16 be_bp_left, Word16 be_bp_right, Word16 L_spec,
    Word16 enc, Word16 sim_dec, Word16 bfi, Pc_State_fx *pc /* i/o: Pc State */
)
{
pc->inv_bin     = add(L_spec, 1);       move16();
pc->numbytes    = numbytes;             move16();
pc->c_bp        = 0;                    move16();
pc->c_bp_side   = 0;                    move16();
pc->bytes       = shr(add(n_pc, 1),1);  move16();
pc->b_left      = add(numbytes,1);      move16();
pc->b_right     = -1;                   move16();
pc->enc         = enc;                  move16();
pc->sim_dec     = sim_dec;              move16();
pc->bfi         = bfi;                  move16();
pc->be_bp_left  = shr(be_bp_left, 3);   move16();
pc->be_bp_right = shr(be_bp_right, 3);  move16();
    assert(pc->be_bp_right < pc->bytes || pc->bytes == 0);
}

static Word16 check_pc_bytes(Word16 *bp, Word16 *bp_side, Word16 *mask_side, Word16 cur_bin,
                                           Word16 from_left, Pc_State_fx *pc /* i/o: Pc State */)
{
    Dyn_Mem_Deluxe_In(Word16 bp_local, bp_side_local, offset;);

    IF (pc->bytes > 0)
    {
        test();
        IF (from_left == 0 && sub(*mask_side, 1) != 0)
        {
            Dyn_Mem_Deluxe_Out();
            return 0;
        }
        test();
        IF (pc->c_bp_side > 0 && *bp_side < 0)
        {
            assert(*mask_side == 1);
            assert(pc->b_right != -1);
            *bp_side = pc->b_right;
            Dyn_Mem_Deluxe_Out();
            return 0;
        }
        bp_local      = *bp;
        bp_side_local = *bp_side;

        IF (from_left != 0)
        {
            if (sub(*mask_side, 1) == 0)
            {
                bp_side_local = add(bp_side_local, 1);
            }
        }
        ELSE
        {
            bp_local = sub(bp_local, 1);
        }

        IF (pc->b_right < 0)
        {
            offset = -1;
            move16();
            if (pc->enc == 0)
            {
                offset = add(offset, pc->bytes);
            }

            IF (add(bp_side_local, sub(offset, bp_local)) == pc->bytes)
            {
                pc->b_left  = add(bp_local, 1);
                pc->b_right = sub(bp_side_local, 1);
                IF (pc->enc != 0)
                {
                    assert(pc->b_right - pc->b_left + 1 == pc->bytes);
                    Dyn_Mem_Deluxe_Out();
                    return 1;
                }
            }
        }

        test();
        IF (pc->enc == 0 && pc->b_right >= 0)
        {
            test();
            IF (from_left != 0 && sub(*bp, pc->b_left) == 0)
            {
                IF (pc->sim_dec == 1)
                {
                    pc->b_left = *bp;
                    Dyn_Mem_Deluxe_Out();
                    return 1;
                }
                *bp = 0;  move16();
                pc->c_bp = 1;  move16();
            }
            test();
            IF (from_left == 0 && sub(bp_side_local, pc->b_right) == 0)
            {
                *bp_side = sub(pc->bytes, 1);
                move16();
                pc->c_bp_side = 1;
                move16();
            }
            IF (sub(pc->bfi, 2) == 0)
            {
                test();
                test();
                IF ((pc->c_bp != 0 && sub(*bp, pc->be_bp_left) >= 0) ||
                    (pc->c_bp_side != 0 && sub(*bp_side, pc->be_bp_right) <= 0))
                {
                    pc->inv_bin = cur_bin;
                    move16();
                    Dyn_Mem_Deluxe_Out();
                    return 1;
                }
                ELSE IF ((pc->c_bp != 0 && *bp >= 0) || (pc->c_bp_side != 0 && sub(*bp_side, sub(pc->bytes, 1)) <= 0))
                {
                    pc->inv_bin = s_min(pc->inv_bin, cur_bin);
                    Dyn_Mem_Deluxe_Out();
                    return 0;
                }
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
    return 0;
}


#ifdef LL_INCL_HPVC

static Word16 read_uint_hpvc(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(Word16 indice, bit; Counter i;);

    ASSERT(numbits > 0);
    indice = read_bit(ptr, bp, mask);
    FOR (i = 1; i < numbits; i++)
    {
        bit = read_bit(ptr, bp, mask);
        indice = add(indice, lshl_pos(bit, i));
    }
    Dyn_Mem_Deluxe_Out();
    return indice;
}

static Word16 ac_uni_get_step_fx(Word16 Ntot, Word16 cdf_bits)
{
    Word16 uniform_stepQ6;
    Word16 invNtotNegQx;
    Word16 Ntrail_z, Nodd;
    Word16 inv_top_bits;
    Word16 tot_r_shift;

    ASSERT(Ntot > 0 && Ntot <= 256);

    Ntrail_z = sub(14, norm_s(s_and(Ntot, negate(Ntot))));
    Nodd = shr_pos(Ntot, Ntrail_z);
    ASSERT((Nodd) & 0x0001);

    inv_top_bits = sub(14, norm_s(negate(Nodd)));
    invNtotNegQx = hpvc_odd_invNtotNegQx[shr_pos(Nodd, 1)];

    tot_r_shift = sub(add((15 - 6), add(inv_top_bits, Ntrail_z)), cdf_bits);
    uniform_stepQ6 = shr(invNtotNegQx, tot_r_shift);
    uniform_stepQ6 = negate(uniform_stepQ6);

    ASSERT(uniform_stepQ6 > 0 && uniform_stepQ6 <= (1 << (cdf_bits + 6 - 1)));
    return uniform_stepQ6;
}

static void ac_encode_uni_fx_st(EncBitBuffState_fw *s, Word16 val, Word16 Ntot,
                                Word16 bit_range_min, Word16 bit_range_max)
{
    UWord32 UL_cumfreq;
    UWord32 UL_symfreq;
    UWord32 UL_cdf_pow2_range;
    Word16  uniform_stepQ6;
    Word16  cdf_bits;
    Word16  cdf_bits1;
    Word16  high_val;
    UWord32 UL_cumfreq_hr;

    Dyn_Mem_Deluxe_In(UWord32 UL_r, UL_tmp; UWord16 carry;);
    BASOP_sub_sub_start("ac_encode_uni_fx_st");

    cdf_bits1 = sub(15, norm_s(Ntot));
    cdf_bits = add(cdf_bits1, 4);
    cdf_bits = s_min(cdf_bits, bit_range_max);
    cdf_bits = s_max(cdf_bits, bit_range_min);

    ASSERT(Ntot >= 2);
    UL_cdf_pow2_range = UL_lshl(1UL, cdf_bits);
    ASSERT((Ntot << 2) < (Word16)UL_cdf_pow2_range);

    uniform_stepQ6 = ac_uni_get_step_fx(Ntot, cdf_bits);

    high_val = add(val, 1);
    UL_cumfreq_hr = UL_lshr((UWord32)L_mult0(high_val, uniform_stepQ6), 6);
    UL_cumfreq_hr = (UWord32)L_min((Word32)UL_cdf_pow2_range, (Word32)UL_cumfreq_hr);
    if (sub(high_val, Ntot) == 0)
    {
        UL_cumfreq_hr = UL_cdf_pow2_range;
        move32();
    }
    UL_cumfreq = UL_lshr((UWord32)L_mult0(val, uniform_stepQ6), 6);
    UL_symfreq = UL_subNsD(UL_cumfreq_hr, UL_cumfreq);

    UL_r = UL_lshr_pos(s->ac_range_fx, cdf_bits);
    UL_tmp = UL_Mpy_32_32(UL_r, UL_cumfreq);

    /* update low, potentially setting the carry */
    {
        Encoder_State_fx tmp;
        tmp.ac_low_fx         = s->ac_low_fx;
        tmp.ac_range_fx       = s->ac_range_fx;
        tmp.ac_cache_fx       = s->ac_cache_fx;
        tmp.ac_carry_fx       = s->ac_carry_fx;
        tmp.ac_carry_count_fx = s->ac_carry_count_fx;

        tmp.ac_low_fx = UL_addNs24(tmp.ac_low_fx, UL_tmp, &carry);
        if (carry != 0)
        {
            tmp.ac_carry_fx = carry;
        }
        tmp.ac_range_fx = UL_Mpy_32_32(UL_r, UL_symfreq);

        WHILE (tmp.ac_range_fx < (1U << 16))
        {
            tmp.ac_range_fx = UL_lshl_pos(tmp.ac_range_fx, 8);
            ac_enc_shift_fx(s->ptr, &s->bp, &tmp);
        }
        s->ac_low_fx         = tmp.ac_low_fx;
        s->ac_range_fx       = tmp.ac_range_fx;
        s->ac_cache_fx       = tmp.ac_cache_fx;
        s->ac_carry_fx       = tmp.ac_carry_fx;
        s->ac_carry_count_fx = tmp.ac_carry_count_fx;
    }

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
}

static void ac_encode_W32N_uni_fx_st(EncBitBuffState_fw *s, EncBitBuffState_bw *s_bw,
                                     Word32 L_cw, Word32 L_Ntot)
{
    Dyn_Mem_Deluxe_In(
        Word16 bits_ceil;
        Word16 bits_int;
        Word16 Ntail, Ntail_msb, Ntail_lsb, lsb_bits, msb_bits;
        Word16 Nrc;
        Word32 L_mask, L_lsbs;);

    BASOP_sub_sub_start("ac_encode_W32N_uni_fx_st");

    bits_ceil = L_sub(31, norm_l(L_negate(L_Ntot)));
    bits_int = sub(bits_ceil, 1);
    if (L_add(L_Ntot, L_shl_pos(-1L, bits_ceil)) == 0)
    {
        bits_int = bits_ceil;
    }
#ifndef LL_INCL_HPVC_UNIW32_OPT
    UNUSED(bits_int);
#endif

#ifdef LL_INCL_HPVC_UNIW32_OPT
    IF (sub(bits_int, bits_ceil) == 0 && sub(bits_int, 28) <= 0)
    {
        /* exact pow2 - transmit cw as 1 to 28 whole bits  */
        Ntail_lsb = s_min(14, bits_ceil);
        ASSERT(Ntail_lsb > 0);
        Ntail_msb = s_max(0, sub(bits_ceil, Ntail_lsb));
        L_mask = L_sub(L_shl_pos(1L, bits_ceil), 1L);
        L_lsbs = L_and(L_cw, L_mask);
        ASSERT(L_lsbs <= (1L << (14 + 14)));
        lsb_bits = extract_l(L_lsbs);
        write_uint_backward_st(s_bw, lsb_bits, Ntail_lsb);
        msb_bits = 0;
        move16();
        IF (Ntail_msb > 0)
        {
            msb_bits = extract_l(L_shr_pos(L_lsbs, Ntail_lsb));
            write_uint_backward_st(s_bw, msb_bits, Ntail_msb);
        }
    }
    ELSE
#endif
    {
        IF (sub(bits_ceil, W32_UNI_RANGE_BITS_MAX) > 0)
        {
            Ntail = sub(bits_ceil, W32_UNI_RANGE_BITS_MAX);
            Ntail_lsb = s_min(14, Ntail);
            ASSERT(Ntail_lsb > 0);
            Ntail_msb = s_max(0, sub(Ntail, Ntail_lsb));
            L_mask = L_sub(L_shl_pos(1L, Ntail), 1L);
            L_lsbs = L_and(L_cw, L_mask);
            lsb_bits = extract_l(L_lsbs);
            write_uint_backward_st(s_bw, lsb_bits, Ntail_lsb);
            msb_bits = 0;
            move16();
            IF (Ntail_msb > 0)
            {
                msb_bits = extract_l(L_shr_pos(L_lsbs, Ntail_lsb));
                write_uint_backward_st(s_bw, msb_bits, Ntail_msb);
            }
            Nrc = extract_l(L_shr_pos(L_Ntot, Ntail));
            if (L_and(L_Ntot, L_mask) != 0L)
            {
                Nrc = add(Nrc, 1);
            }
            {
                Word16 cw = extract_l(L_shr_pos(L_cw, Ntail));
#ifdef LL_INCL_HPVC_UNIW32_OPT
                IF (s_and(Nrc, sub(Nrc, 1)) == 0)
                {
                    bits_ceil = sub(15, norm_s(negate(Nrc)));
                    ASSERT((1 << bits_ceil) == Nrc);
                    write_uint_backward_st(s_bw, cw, bits_ceil);
                }
                ELSE
                {
                    ac_encode_uni_fx_st(s, cw, Nrc, LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE);
                }
#else
                ac_encode_uni_fx_st(s, cw, Nrc, LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE);
#endif
            }
        }
        ELSE
        {
            ASSERT(L_Ntot >= 2 && L_Ntot <= 256);
            ac_encode_uni_fx_st(s, extract_l(L_cw), extract_l(L_Ntot),
                                LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE);
        }
    }
    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
}


static Word16 ac_decode_generic_fx(Decoder_State_fx *st_fx, const UWord16 *cumfreq, Word16 num_sym)
{
    Dyn_Mem_Deluxe_In(Word16 val;);
    BASOP_sub_sub_start("ac_decode_generic_fx");

    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, 10);
    if (st_fx->ac_low_fx >= UL_lshl(st_fx->ac_help_fx, 10))
    {
        st_fx->BER_detect = 1;        move16();
    }
    val = sub(num_sym, 1);
    WHILE (st_fx->ac_low_fx < UL_Mpy_32_32(st_fx->ac_help_fx, UL_deposit_l(cumfreq[val])))
    {
        val = sub(val, 1);
    }
    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return val;
}

static Word16 ac_decode_uni_fx(Decoder_State_fx *st_fx, Word16 Ntot,
                               Word16 cdf_bits_min, Word16 cdf_bits_max,
                               UWord32 *UL_cumfreq_upd_ptr, UWord32 *UL_symfreq_upd_ptr)
{
    Dyn_Mem_Deluxe_In(UWord16 sgn; Word16 val, tst_val; UWord32 UL_tmp, UL_cumfreq, UL_cdf_pow2_range;);

    BASOP_sub_sub_start("ac_decode_uni_fx");
    UNUSED(tst_val);

    Word16 cdf_bits;
    Word16 cdf_bits1 = sub(15, norm_s(Ntot));
    cdf_bits = add(cdf_bits1, 4);
    cdf_bits = s_max(cdf_bits_min, cdf_bits);
    cdf_bits = s_min(cdf_bits_max, cdf_bits);

    UL_cdf_pow2_range = UL_lshl(1UL, cdf_bits);
    st_fx->ac_help_fx = UL_lshr_pos(st_fx->ac_range_fx, (Word16)cdf_bits);
    val = Ntot; move16();

    Word16 uniform_stepQ6 = ac_uni_get_step_fx(Ntot, cdf_bits);

    UL_tmp = UL_lshl(st_fx->ac_help_fx, cdf_bits);
    ASSERT(st_fx->ac_low_fx <= UL_tmp);

    Word16 sym0;
    UWord32 UL_num = UL_Mpy_32_32(st_fx->ac_low_fx, UL_deposit_l(Ntot));
    Word16 num_shift = sub(norm_ul(UL_num), 1);
    UWord32 UL_num_up = UL_lshl(UL_num, num_shift);
    Word16 den_shift = norm_ul(UL_tmp);
    UWord32 UL_tmp_up = UL_lshl(UL_tmp, den_shift);

    UL_subNsD(0, 0);
    if (UL_lshr(UL_tmp_up, 32 - 8) < 255)
    {
        UL_tmp_up = UL_addNsD(UL_tmp_up, 0x00800000U);
    }
    Word16 low_res_den = (Word16)extract_l(L_lshr((Word32)UL_tmp_up, 24));

    #define INV_TOP_BITS 7
    Word32 L_sym = Mpy_32_16_0_0((Word32)UL_num_up, hpvc_invDenNegQx[sub(low_res_den, 128)]);
    L_sym = L_negate(L_sym);

    Word16 total_shift = sub(num_shift, den_shift);
    total_shift = add(total_shift, INV_TOP_BITS + (32 - 24) + 15);
    total_shift = s_min(31, total_shift);
    ASSERT(total_shift >= 0);

    sym0 = L_shr_pos(L_sym, total_shift);

    Word16 sym1 = add(sym0, 1);
    sym1 = s_max(2, sym1);
    sym1 = s_min(sym1, sub(Ntot, 1));

    val = sub(sym1, 2); move16();
    tst_val = add(val, 1);
    UL_cumfreq = UL_lshr((UWord32)L_mult0(tst_val, uniform_stepQ6), 6);
    UL_tmp = UL_Mpy_32_32(st_fx->ac_help_fx, UL_cumfreq);
    UL_subNs(st_fx->ac_low_fx, UL_tmp, &sgn);
    if (sgn == 0) { val = tst_val; move16(); }

    tst_val = add(tst_val, 1);
    UL_cumfreq = UL_lshr((UWord32)L_mult0(tst_val, uniform_stepQ6), 6);
    UL_subNs(st_fx->ac_low_fx, UL_Mpy_32_32(st_fx->ac_help_fx, UL_cumfreq), &sgn);
    if (sgn == 0) { val = tst_val; move16(); }

    Word16 high_val = add(val, 1);
    UL_cumfreq = UL_lshr((UWord32)L_mult0(high_val, uniform_stepQ6), 6);
    UL_cumfreq = (UWord32)L_min((Word32)UL_cdf_pow2_range, (Word32)UL_cumfreq);
    if (sub(high_val, Ntot) == 0)
    {
        UL_cumfreq = UL_cdf_pow2_range;
        move32();
    }
    *UL_cumfreq_upd_ptr = UL_lshr((UWord32)L_mult0(val, uniform_stepQ6), 6);
    *UL_symfreq_upd_ptr = UL_subNsD(UL_cumfreq, *UL_cumfreq_upd_ptr);

    if (sub(val, Ntot) >= 0)
    {
        st_fx->BER_detect = 1;
        move16();
    }
    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return val;
    #undef INV_TOP_BITS
}

static Word32 ac_decode_W32N_uni_fx(Decoder_State_fx *st_fx, UWord8 *ptr, Word16 *bp,
                                    Word16 *bp_side, Word16 *mask_side, Word32 L_Ntot)
{
    Dyn_Mem_Deluxe_In(
        Word32 L_cw;
        Word16 cw;
        Word16 bits_ceil;
        Word16 bits_int;
        Word16 Ntail, Ntail_lsb, Ntail_msb;
        Word16 lsb_bits, msb_bits;
        Word16 Nrc;
        Word32 L_mask, L_lsbs;
        UWord32 UL_cumfreq_upd;
        UWord32 UL_symfreq_upd;);

    BASOP_sub_sub_start("ac_decode_W32N_uni_fx");

    bits_ceil = sub(31, norm_l(L_negate(L_Ntot)));
    bits_int = sub(bits_ceil, 1);
    if (L_add(L_Ntot, L_shl_pos(-1L, bits_ceil)) == 0)
    {
        bits_int = bits_ceil;
    }
#ifndef LL_INCL_HPVC_UNIW32_OPT
    UNUSED(bits_int);
#endif

#ifdef LL_INCL_HPVC_UNIW32_OPT
    IF (sub(bits_int, bits_ceil) == 0 && sub(bits_int, 28) <= 0)
    {
        /* exact pow2 - read cw as 1 to 28 whole bits  */
        Ntail_lsb = s_min(14, bits_ceil);
        Ntail_msb = s_max(0, sub(bits_ceil, Ntail_lsb));
        lsb_bits = read_uint_hpvc(ptr, bp_side, mask_side, Ntail_lsb);
        msb_bits = 0;
        move16();
        L_lsbs = L_deposit_l(lsb_bits);
        IF (Ntail_msb > 0)
        {
            msb_bits = read_uint_hpvc(ptr, bp_side, mask_side, Ntail_msb);
            L_lsbs = L_add(L_shl_pos(L_deposit_l(msb_bits), Ntail_lsb), L_lsbs);
        }
        ASSERT(L_lsbs <= (1L << (14 + 14)));
        L_cw = L_add(L_lsbs, 0);
    }
    ELSE
#endif
    {
        IF (sub(bits_ceil, W32_UNI_RANGE_BITS_MAX) > 0)
        {
            Ntail = sub(bits_ceil, W32_UNI_RANGE_BITS_MAX);
            Ntail_lsb = s_min(14, Ntail);
            Ntail_msb = s_max(0, sub(Ntail, Ntail_lsb));
            L_mask = L_sub(L_shl_pos(1L, Ntail), 1L);

            lsb_bits = read_uint_hpvc(ptr, bp_side, mask_side, Ntail_lsb);
            msb_bits = 0;
            move16();
            L_lsbs = L_deposit_l(lsb_bits);
            IF (Ntail_msb > 0)
            {
                msb_bits = read_uint_hpvc(ptr, bp_side, mask_side, Ntail_msb);
                L_lsbs = L_add(L_shl_pos(L_deposit_l(msb_bits), Ntail_lsb), L_lsbs);
            }

            Nrc = extract_l(L_shr_pos(L_Ntot, Ntail));
            if (L_and(L_Ntot, L_mask) != 0L)
            {
                Nrc = add(Nrc, 1);
            }
#ifdef LL_INCL_HPVC_UNIW32_OPT
            IF (s_and(Nrc, sub(Nrc, 1)) == 0)
            {
                bits_ceil = sub(15, norm_s(negate(Nrc)));
                ASSERT((1 << bits_ceil) == Nrc);
                cw = read_uint_hpvc(ptr, bp_side, mask_side, bits_ceil);
            }
            ELSE
            {
                cw = ac_decode_uni_fx(st_fx, Nrc, LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE,
                                      &UL_cumfreq_upd, &UL_symfreq_upd);
                ac_dec_update_fx(ptr, bp, bp_side, mask_side, 0, st_fx, UL_cumfreq_upd, UL_symfreq_upd);
            }
#else
            cw = ac_decode_uni_fx(st_fx, Nrc, LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE,
                                  &UL_cumfreq_upd, &UL_symfreq_upd);
            ac_dec_update_fx(ptr, bp, bp_side, mask_side, 0, st_fx, UL_cumfreq_upd, UL_symfreq_upd);
#endif
            L_cw = L_add(L_shl_pos(L_deposit_l(cw), Ntail), L_lsbs);
        }
        ELSE
        {
            ASSERT(L_Ntot >= 2 && L_Ntot <= 256);
            cw = ac_decode_uni_fx(st_fx, extract_l(L_Ntot), LL_HPVC_MIN_RANGE, LL_HPVC_MAX_RANGE,
                                  &UL_cumfreq_upd, &UL_symfreq_upd);
            L_cw = L_deposit_l(cw);
            ac_dec_update_fx(ptr, bp, bp_side, mask_side, 0, st_fx, UL_cumfreq_upd, UL_symfreq_upd);
        }
    }

    if (L_sub(L_cw, L_Ntot) >= 0)
    {
        st_fx->BER_detect = 1;
        move16();
    }
    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return L_cw;
}

static Word16 post_HPVC_update_TCX_context(Word32 *L_x_tail)
{
    Counter i;
    Word16 c, ctx[2];
    Word32 L_a, L_b, L_abMax, L_ab0L1norm;
    Word16 dn_shift;
    Word16 lev, esc_nb;

    FOR (i = 0; i < 2; i++)
    {
        L_a = L_abs(L_x_tail[2 * i]);
        L_b = L_abs(L_x_tail[2 * i + 1]);
        L_abMax = L_max(L_a, L_b);
        IF (L_abMax == 0)
        {
            lev = 0;
            move16();
            ctx[i] = 0;
            move16();
        }
        ELSE
        {
            ASSERT(L_abMax > 0);
            dn_shift = sub(30, norm_l(L_abMax));
            lev = s_max(0, sub(dn_shift, 1));
            L_ab0L1norm = L_add(L_shr_pos(L_a, lev), L_shr_pos(L_b, lev));
            esc_nb = s_min(3, lev);

            ctx[i] = extract_l(L_ab0L1norm);
            if (esc_nb > 0)
            {
                ctx[i] = extract_l(L_shl_pos(L_ab0L1norm, 1));
            }
            if (sub(esc_nb, 1) > 0)
            {
                ctx[i] = add(11, esc_nb);
            }
        }
        ctx[i] = add(ctx[i], 1);
        ASSERT(ctx[i] > 0 && ctx[i] <= 0xf);
    }
    c = s_or(shl_pos(ctx[0], 4), ctx[1]);
    ASSERT(c >= 0 && c <= 0xff);
    return c;
}
#endif /* LL_INCL_HPVC */
