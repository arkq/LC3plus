/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void ac_shift_fl(Encoder_State_fl* st);
static void ac_encode_fl(Encoder_State_fl* st, LC3_INT sym_freq, LC3_INT cum_freq);
static void ac_finalize_fl(Encoder_State_fl* st);
static void write_uint_forward_fl(Encoder_State_fl* st, LC3_INT val, LC3_INT numbits);
static void ari_enc_init(Encoder_State_fl* st, LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side);
static LC3_INT  sign(LC3_INT x);

static void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit);

static void ac_dec_init_fl(LC3_UINT8* ptr, LC3_INT* bp, Decoder_State_fl* st_fl, LC3_INT from_left, LC3_INT mask_side, LC3_INT *bp_side);

static LC3_INT32 ac_decode_fl(Decoder_State_fl* st, const LC3_INT16* sym_freq, LC3_INT32 num_sym, LC3_UINT8* ptr, LC3_INT32* bp, LC3_INT32 from_left, LC3_INT32 mask_side, LC3_INT32 *bp_side, LC3_INT16 cur_bin);

static LC3_INT16 pc_check_bytes(LC3_INT32* bp, Decoder_State_fl* st_fl, LC3_INT32 from_left, LC3_INT32 mask_side, LC3_INT32 *bp_side, LC3_INT16 cur_bin);

static void calculate_nfseed(LC3_INT *x, LC3_INT L_spec, LC3_INT *nf_seed);
static void findNonZero(LC3_INT* in, LC3_INT len, LC3_INT* outLen);

void findNonZero(LC3_INT* in, LC3_INT len, LC3_INT* outLen)
{
    LC3_INT i = 0, j = 0;

    for (i = 0; i < len; i++) {
        if (in[i] != 0) {
            j++;
        }
    }

    *outLen = j;
}

void calculate_nfseed(LC3_INT *x, LC3_INT L_spec, LC3_INT *nf_seed)
{
    LC3_INT k;
    
    *nf_seed = 0;
    
    for (k = 0; k < L_spec; k++) {
        *nf_seed = *nf_seed + (abs(x[k]) & 32767) * k;
    }
    *nf_seed = *nf_seed & 65535;

    if (*nf_seed >= 32768) {
        *nf_seed = *nf_seed - 65536;
    }
}

static LC3_INT16 pc_check_bytes(LC3_INT32* bp, Decoder_State_fl* st_fl, LC3_INT32 from_left, LC3_INT32 mask_side, LC3_INT32 *bp_side, LC3_INT16 cur_bin)
{
    LC3_INT32 bp_local, bp_side_local, offset;
#ifdef WMOPS
    push_wmops("pc_check_bytes");
#endif

    if (st_fl->pc_bytes > 0)
    {
        if (!from_left && mask_side != 1)
        {
            return 0;
        }

        if (st_fl->pc_c_bp_side > 0 && *bp_side < 0)
        {
            assert(mask_side == 1);
            assert(st_fl->pc_b_right != -1);
            *bp_side = st_fl->pc_b_right;
            
            return 0;
        }

        bp_local = *bp;
        bp_side_local = *bp_side;
        
        if (from_left)
        {
            if (mask_side == 1)
            {
                bp_side_local = bp_side_local + 1;
            }
        } else {
            bp_local = bp_local - 1;
        }
        
        if (st_fl->pc_b_right == -1)
        {
            offset = -1;
            if (!st_fl->pc_enc)
            {
                offset = offset + st_fl->pc_bytes;
            }
            
            if ((bp_side_local + offset - bp_local) == st_fl->pc_bytes)
            {
                st_fl->pc_b_left = bp_local + 1;
                st_fl->pc_b_right = bp_side_local - 1;
                
                if (st_fl->pc_enc)
                {
                    assert(st_fl->pc_b_right - st_fl->pc_b_left + 1 == st_fl->pc_bytes);
                    return 1;
                }
            }
        }
        
        if (!st_fl->pc_enc && st_fl->pc_b_right > -1)
        {
            if (from_left && *bp == st_fl->pc_b_left)
            {
                *bp = 0;
                st_fl->pc_c_bp = 1;
            }
            
            if (!from_left && bp_side_local == st_fl->pc_b_right)
            {
                *bp_side = st_fl->pc_bytes - 1;
                st_fl->pc_c_bp_side = 1;
            }
            
            if (st_fl->pc_bfi == 2)
            {
                
                if ((st_fl->pc_c_bp && (*bp + 1) >= st_fl->pc_be_bp_left) || (st_fl->pc_c_bp_side && (*bp_side + 1) <= st_fl->pc_be_bp_right))
                {
                    st_fl->pc_inv_bin = cur_bin;
                    return 1;
                } else if ((st_fl->pc_c_bp && *bp >= 0) || (st_fl->pc_c_bp_side && *bp_side <= (st_fl->pc_bytes - 1)))
                {
                    st_fl->pc_inv_bin = MIN(st_fl->pc_inv_bin, cur_bin);
                    return 0;
                }
            }
        }   
    }

#ifdef WMOPS
    pop_wmops();
#endif   
    return 0;
}

void ac_dec_init_fl(LC3_UINT8* ptr, LC3_INT* bp, Decoder_State_fl* st_fl, LC3_INT from_left, LC3_INT mask_side, LC3_INT *bp_side)
{
    LC3_INT i;

    if (!st_fl->pc_enc)
    {
        *bp = *bp + st_fl->pc_bytes;
    }

    st_fl->ac_low_fl = 0;

    st_fl->ac_range_fl = (LC3_UINT32) 16777215;  /* 2^24 -1 */
    for (i = 0; i < 3; i++) {
        if(pc_check_bytes(bp, st_fl, from_left, mask_side, bp_side, 0) != 0)
        {
            return;
        }
        
        st_fl->ac_low_fl = (st_fl->ac_low_fl << 8) + (LC3_UINT32)ptr[*bp];
        *bp              = *bp + 1;
    }

    st_fl->BER_detect = 0;
}

/* Returns val */
LC3_INT32 ac_decode_fl(Decoder_State_fl* st, const LC3_INT16* freq, LC3_INT32 num_sym, LC3_UINT8* ptr, LC3_INT32* bp, LC3_INT32 from_left, LC3_INT32 mask_side, LC3_INT32 *bp_side, LC3_INT16 cur_bin)
{
    LC3_INT val, tmp, symfreq_loc;
#ifdef WMOPS
    push_wmops("ac_decode_fl");    
#endif

    tmp = st->ac_range_fl >> 10;

    if (st->ac_low_fl >= (LC3_UINT32)(tmp << 10)) {
        st->BER_detect = 1;
    }

    val = num_sym - 1;

    while (st->ac_low_fl < (LC3_UINT32)(tmp * freq[val])) {
        val--;
    }
    
    symfreq_loc = freq[val + 1] - freq[val];

    st->ac_low_fl   = st->ac_low_fl - tmp * freq[val];
    st->ac_range_fl = tmp * symfreq_loc;

    while (st->ac_range_fl < 65536) {
        st->ac_low_fl   = ((LC3_INT32)st->ac_low_fl) & ((LC3_INT32)(16777215));

        if(pc_check_bytes(bp, st, from_left, mask_side, bp_side, cur_bin) != 0)
        {
            st->BER_detect = 1;
            return 1;
        }
        
        st->ac_low_fl   = st->ac_low_fl << 8;
        st->ac_low_fl   = st->ac_low_fl + ptr[*bp];
        *bp             = *bp + 1;
        st->ac_range_fl = st->ac_range_fl << 8;
    }
    
#ifdef WMOPS
    pop_wmops();
#endif
    return val;
}

void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit)
{
    if (ptr[*bp_side] & *mask_side) {
        *bit = 1;
    } else {
        *bit = 0;
    }

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side * 2;
    }
}

void processAriDecoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT L_spec, LC3_INT fs_idx, LC3_INT enable_lpc_weighting,
                          LC3_INT tns_numfilters, LC3_INT lsbMode, LC3_INT lastnz, LC3_INT* bfi, LC3_INT* tns_order, LC3_INT fac_ns_idx,
                          LC3_INT gg_idx, uint8_t * resBits, LC3_INT* x, LC3_INT* nf_seed, LC3_INT* tns_idx, LC3_INT* zero_frame, LC3_INT numbytes,
                          LC3_INT* nbits_residual, LC3_INT* residualPresent, LC3_INT frame_dms,
                          LC3_INT32 n_pc, LC3_INT32 be_bp_left, LC3_INT32 be_bp_right, LC3_INT32 enc, LC3_INT32 *b_left, LC3_INT32 *spec_inv_idx,
                          LC3_INT hrmode
)
{
    Decoder_State_fl st;
    LC3_INT              a, b, t, bp;
    LC3_INT              c;
    LC3_INT              nbits_side, extra_bits;
    LC3_UINT8*           ptr;
    LC3_INT              n, k, lev;
    LC3_INT              max_lev, tmp;
    LC3_INT              bit, lev1, pki, sym, save_lev[MAX_LEN], idx_len, total_bits, nbits_ari, rateFlag;

#ifdef WMOPS
    push_wmops("processAriDecoder_fl");
#endif

    total_bits = 8 * numbytes;
    rateFlag = 0;
    
    memset(&st, 0, sizeof(st));
    
    st.pc_bytes = (n_pc + 1) >> 1;
    st.pc_b_left = numbytes + 1;
    st.pc_b_right = -1;
    st.pc_enc = enc;
    st.pc_bfi = *bfi;
    st.pc_be_bp_left = floor(be_bp_left / 8);
    st.pc_be_bp_right = floor(be_bp_right / 8) - 1;
    *spec_inv_idx = L_spec + 1;
    assert(st.pc_be_bp_right < st.pc_bytes || st.pc_bytes == 0);

    /* Rate flag */
    if (fs_idx != 5)
    {
        if (total_bits > (160 + fs_idx * 160)) {
            rateFlag = 512;
        }
    }

    /* Init */
    c  = 0;
    t  = 0;
    bp = 0;
    
    *b_left = -1;

    ptr = bytes;

    /* Start Decoding */
    ac_dec_init_fl(ptr, &bp, &st, 1, mask_side, &bp_side);
    
    /* Decode TNS data */
    tmp = MAXLAG;
    

    if (frame_dms <= 50)
    {
        tmp /= 2;
    }

    /* Decode TNS data */
    for (n = 0; n < tns_numfilters; n++) {
    
        if (tns_order[n] > 0) {
            tns_order[n] = ac_decode_fl(&st, &ari_tns_order_cf[enable_lpc_weighting][0], 8, ptr, &bp, 1, mask_side, &bp_side, 0);
            
            tns_order[n] = tns_order[n] + 1;
            
            if (tns_order[n] > tmp || st.BER_detect > 0)
            {
                goto ber_detect;
            }

            for (k = 0; k < tns_order[n]; k++) {
                if (bp_side < bp)
                {
                    *bfi = 1;
                    return;
                }
            
                tns_idx[n * 8 + k] = ac_decode_fl(&st, &ari_tns_freq_cf[k][0], 17, ptr, &bp, 1, mask_side, &bp_side, 0);
                
                if (st.BER_detect > 0)
                {
                    goto ber_detect;
                }
            }
        }
    }

    /* Spectral data */
    for (k = 0; k < lastnz; k = k + 2) {
        /* Context */
        t = c + rateFlag;

        if (k > (L_spec >> 1)) {
            t = t + 256;
        }

        /* Decode amplitude */
        x[k]     = 0;
        x[k + 1] = 0;

        if (hrmode == 1) {
            max_lev = 13 + 8;
        } else {
            max_lev = 13;
        }

        for (lev = 0; lev <= max_lev; lev++) {
            lev1 = MIN(lev, 3);
            pki  = ari_spec_lookup_fl[t + lev1 * 1024];

            sym = ac_decode_fl(&st, &ari_spec_cumfreq_fl[pki][0], 17, ptr, &bp, 1, mask_side, &bp_side, k);

            if (sym < 16) {
                break;
            }

            if (lsbMode == 0 || lev > 0) {
                if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                {
                    goto ber_detect;
                }
                read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                
                x[k] = x[k] + (bit << lev);
                if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                {
                    goto ber_detect;
                }
                read_bit_fl(ptr, &mask_side, &bp_side, &bit);

                x[k + 1] = x[k + 1] + (bit << lev);
            }
        }
        
        if ((lev - 1) == 13 && sym == 16)
        {
            goto ber_detect;
        }
        
        if (hrmode == 0) {
            lev = MIN(lev, 13);
        }

        if (lsbMode == 1) {
            save_lev[k] = lev;
        }

        a = sym & 3;
        b = sym >> 2;

        x[k]     = x[k] + (a << lev);
        x[k + 1] = x[k + 1] + (b << lev);

        /* Decode signs */
        if (x[k] > 0) {
            if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
            {
                goto ber_detect;
            }
            read_bit_fl(ptr, &mask_side, &bp_side, &bit);

            if (bit == 1) {
                x[k] = -x[k];
            }
        }

        if (x[k + 1] > 0) {
            if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
            {
                goto ber_detect;
            }
            read_bit_fl(ptr, &mask_side, &bp_side, &bit);
            
            if (bit == 1) {
                x[k + 1] = -x[k + 1];
            }
        }

        /* Context */
        lev1 = MIN(lev, 3);
        if (lev1 <= 1) {
            t = 1 + (a + b) * (lev1 + 1);
        } else {
            t = 12 + lev1;
        }

        c = (c & 15) * 16 + t;

        if (((bp - bp_side) > 3 && (st.pc_c_bp == st.pc_c_bp_side))) {

            if ((0 < *spec_inv_idx) && (*spec_inv_idx < (L_spec + 1)))
            {
                *bfi = 2;
                calculate_nfseed(x, k, nf_seed);
                return; 
            }

            *bfi = 1;
            return;
        }
        
        if (st.BER_detect > 0)
        {
            goto ber_detect;
        }
    }

    /* Residual bits */
    nbits_side      = total_bits - (8 * bp_side + 8 - (31 - clz_func(mask_side)));
    nbits_ari       = (bp - 3) * 8;
    extra_bits      = 25 - (31 - clz_func(st.ac_range_fl));

    if (enc == 0)
    {
        if (st.pc_c_bp == 0)
        {
            nbits_ari = (bp - st.pc_bytes - 3) * 8;
        } else {
            nbits_ari = (bp + st.pc_b_left - st.pc_bytes - 3) * 8;
        }
        
        if (st.pc_c_bp_side != 0)
        {
            nbits_side = total_bits - 8 * (st.pc_b_left) + 8 * (st.pc_bytes - bp_side) - (8 - LC3_LOGTWO(mask_side));
        }
    }

    
    *nbits_residual = total_bits - (nbits_side + nbits_ari + extra_bits);

    if (*nbits_residual < 0) {
        if ((0 < *spec_inv_idx) && (*spec_inv_idx < (L_spec + 1)))
        {
            *bfi = 2;
            calculate_nfseed(x, k, nf_seed);
            return;
        }
    
        *bfi = 1;
        return;
    }

    if (lsbMode == 0) {
        findNonZero(x, L_spec, &idx_len);

        if (hrmode)
        {
            idx_len *= EXT_RES_ITER_MAX;
        }
        *nbits_residual  = MIN(*nbits_residual, idx_len);
        *residualPresent = 1;

        memset(resBits, 0, MAX_RESBITS_LEN);

        for (k = 0; k < *nbits_residual; k++) {
            if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
            {
                goto ber_detect_res;
            }
            read_bit_fl(ptr, &mask_side, &bp_side, &tmp);
            
            resBits[k >> 3] |= tmp << (k & 7);
        }
    } else {
        for (k = 0; k < lastnz; k = k + 2) {
            if (save_lev[k] > 0) {
                if (*nbits_residual == 0) {
                    break;
                }

                if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                {
                    goto ber_detect_res;
                }
                read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                
                *nbits_residual = *nbits_residual - 1;

                if (bit == 1) {
                    if (x[k] > 0) {
                        x[k] = x[k] + 1;
                    } else if (x[k] < 0) {
                        x[k] = x[k] - 1;
                    } else {
                        if (*nbits_residual == 0) {
                            break;
                        }

                        if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                        {
                            goto ber_detect_res;
                        }
                        read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                        
                        *nbits_residual = *nbits_residual - 1;

                        if (bit == 0) {
                            x[k] = 1;
                        } else {
                            x[k] = -1;
                        }
                    }
                }

                if (*nbits_residual == 0) {
                    break;
                }

                if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                {
                    goto ber_detect_res;
                }
                read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                
                *nbits_residual = *nbits_residual - 1;

                if (bit == 1) {
                    if (x[k + 1] > 0) {
                        x[k + 1] = x[k + 1] + 1;
                    } else if (x[k + 1] < 0) {
                        x[k + 1] = x[k + 1] - 1;
                    } else {
                        if (*nbits_residual == 0) {
                            break;
                        }

                        if(pc_check_bytes(&bp, &st, 0, mask_side, &bp_side, k) != 0)
                        {
                            goto ber_detect_res;
                        }
                        read_bit_fl(ptr, &mask_side, &bp_side, &bit);
                        
                        *nbits_residual = *nbits_residual - 1;

                        if (bit == 0) {
                            x[k + 1] = 1;
                        } else {
                            x[k + 1] = -1;
                        }
                    }
                }
            }
        }
    }

    /* Noise-filling seed */
    calculate_nfseed(x, L_spec, nf_seed);

    /* Zero frame flag */
    if (lastnz == 2 && x[0] == 0 && x[1] == 0 && gg_idx == 0 && fac_ns_idx == 7) {
        *zero_frame = 1;
    } else {
        *zero_frame = 0;
    }
    
    if (enc)
    {
        if (st.pc_bytes > 0)
        {
            if (st.pc_b_left > numbytes)
            {
                *b_left = bp_side - st.pc_bytes;
            }
        }
    } else {
        if (st.pc_bytes > 0)
        {
            if (st.pc_b_left > numbytes)
            {
                *b_left = bp_side;
            }
        }
    }
    
    if ((*bfi == 2) && (*spec_inv_idx == (L_spec + 1)))
    {
        *bfi = 0;
    }
    
    *spec_inv_idx = *spec_inv_idx - 1;

    goto bail;

/* goto for bit error handling */
ber_detect:
    *bfi = 1;
    *b_left = st.pc_b_left;

    if (st.pc_inv_bin > 0 && (st.pc_inv_bin - L_spec) <= 0)
    {
        *spec_inv_idx = st.pc_inv_bin;
        *bfi = 2;
        *resBits = 0;
        *zero_frame = 0;
        /* Noise Filling seed */
        calculate_nfseed(x, *spec_inv_idx, nf_seed);
    }
    goto bail;

/* goto for bit error handling in residual signal */
ber_detect_res:
    *b_left = st.pc_b_left;
    *resBits = 0;
    *bfi = 0;
    *zero_frame = 0;
    /* Noise Filling seed */
    calculate_nfseed(x, *spec_inv_idx, nf_seed);
    goto bail;

/* goto, because of dynmem out */
bail:

#ifdef WMOPS
    pop_wmops();
#endif
    /* Avoid warning "label at end of compound statement" when WMOPS is inactive */
    (void)0;
}

void ac_encode_fl(Encoder_State_fl* st, LC3_INT sym_freq, LC3_INT cum_freq)
{
    LC3_INT r;

    r       = st->range >> 10;
    st->low += r * cum_freq;

    if ((st->low >> 24) == 1) {
        st->carry = 1;
    }

    st->low   &= (16777215); /* 2^24 -1 */
    st->range = r * sym_freq;

    while (st->range < 65536) {  /* 2^16 */
        st->range <<= 8;
        ac_shift_fl(st);
    }
}

void ac_shift_fl(Encoder_State_fl* st)
{
    if (st->low < 16711680 || st->carry == 1) {
        if (st->cache >= 0) {
            st->ptr[st->bp] = st->cache + st->carry;
            st->bp          = st->bp + 1;
        }

        while (st->carry_count > 0) {
            st->ptr[st->bp] = (st->carry + 255) & 255;
            st->bp          = st->bp + 1;
            st->carry_count = st->carry_count - 1;
        }

        st->cache = st->low >> 16;
        st->carry = 0;
    } else {
        st->carry_count = st->carry_count + 1;
    }

    st->low = st->low << 8;
    st->low = (st->low) & (16777215); /* 2^24 - 1 */
}


void ac_finalize_fl(Encoder_State_fl* st)
{
    LC3_INT bits = 0, mask = 0, val = 0, over1 = 0, high = 0, over2 = 0, c = 0, b = 0;

    bits = 24 - (31 - clz_func(st->range));
    mask  = 16777215 >> bits;
    val   = st->low + mask;
    over1 = val >> 24;

    val   = (val) & 16777215;
    high  = st->low + st->range;
    over2 = high >> 24;
    high  = high & 16777215;
    val   = val & (16777215 - mask);

    if (over1 == over2) {
        if (val + mask >= high) {
            bits = bits + 1;
            mask = mask >> 1;
            val  = ((st->low + mask) & (16777215)) & (16777215 - mask);
        }

        if (val < st->low) {
            st->carry = 1;
        }
    }

    st->low = val;

    b = bits;

    if (bits > 8) {
        for (; b >= 1; b = b - 8) {
            ac_shift_fl(st);
        }
    } else {
        ac_shift_fl(st);
    }

    bits = b;
    if (bits < 0) {
        bits += 8;
    }

    if (st->carry_count > 0) {
        st->ptr[st->bp] = st->cache;
        st->bp          = st->bp + 1;

        for (c = st->carry_count; c >= 2; c--) {
            st->ptr[st->bp] = 255;
            st->bp          = st->bp + 1;
        }

        write_uint_forward_fl(st, 255 << (bits - 8), bits);
    } else {
        write_uint_forward_fl(st, st->cache, bits);
    }
}

void write_uint_forward_fl(Encoder_State_fl* st, LC3_INT val, LC3_INT numbits)
{
    LC3_INT k, bit, mask = 128;

    for (k = 0; k < numbits; k++) {
        bit = val & mask;

        if (bit == 0) {
            st->ptr[st->bp] = st->ptr[st->bp] & (255 - mask);
        } else {
            st->ptr[st->bp] = st->ptr[st->bp] | mask;
        }

        mask = mask >> 1;
    }
}

void ari_enc_init(Encoder_State_fl* st, LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side)
{
    st->ptr         = bytes;
    st->bp_side     = bp_side;
    st->mask_side   = mask_side;
    st->bp          = 0;
    st->low         = 0;
    st->range       = 16777215;
    st->cache       = -1;
    st->carry       = 0;
    st->carry_count = 0;
}

LC3_INT sign(LC3_INT x)
{
    if (x > 0)
        return 1;

    if (x < 0)
        return -1;

    return 0;
}

void processAriEncoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT* x, LC3_INT* tns_order, LC3_INT tns_numfilters,
                          LC3_INT* tns_idx, LC3_INT lastnz, LC3_INT* codingdata, uint8_t* res_bits, LC3_INT resBitsLen, LC3_INT lsbMode,
                          LC3_INT nbbits, LC3_INT enable_lpc_weighting)
{
    LC3_INT              total_bits, cumfreq, symfreq, k, i, j, lev, lev1;
    LC3_INT              bit1, bit2, lsb1, lsb2, a, b, bit, pki, nbits_side;
    LC3_INT              nbits_residual_enc, nbits_ari, lsbs[MAX_LEN], lsbsLen = 0;
	LC3_INT              abs_x_k, abs_x_kp1;
    Encoder_State_fl st;

#ifdef WMOPS
    push_wmops("processAriEncoder_fl");
#endif

    ari_enc_init(&st, bytes, &bp_side, &mask_side);

    total_bits = nbbits;

    /* TNS data  */
    for (i = 0; i < tns_numfilters; i++) {
        if (tns_order[i] > 0) {
            symfreq = tns_freq_cf[enable_lpc_weighting][tns_order[i]] - tns_freq_cf[enable_lpc_weighting][tns_order[i] - 1];
            cumfreq = tns_freq_cf[enable_lpc_weighting][tns_order[i] - 1];
            ac_encode_fl(&st, symfreq, cumfreq);

            for (j = 0; j < tns_order[i]; j++) {
                symfreq = tns_cf[j][tns_idx[i * 8 + j] + 1] - tns_cf[j][tns_idx[i * 8 + j]];
                cumfreq = tns_cf[j][tns_idx[i * 8 + j]];
                ac_encode_fl(&st, symfreq, cumfreq);
            }
        }
    }

    /* Spectral data */
    for (k = 0; k < lastnz; k = k + 2) {
        abs_x_k = abs(x[k]);
        abs_x_kp1 = abs(x[k + 1]);
        for (lev = 0; lev < codingdata[1]; lev++) {
            lev1 = MIN(lev, 3);
            pki  = ari_spec_lookup_fl[codingdata[0] + lev1 * 1024];
            symfreq = ari_spec_cumfreq_fl[pki][17] - ari_spec_cumfreq_fl[pki][16];
            cumfreq = ari_spec_cumfreq_fl[pki][16];

            ac_encode_fl(&st, symfreq, cumfreq);
            bit1 = (abs_x_k >> lev) & 1;
            bit2 = (abs_x_kp1 >> lev) & 1;


            if (lsbMode == 1 && lev == 0) {
                lsb1 = bit1;
                lsb2 = bit2;
            } else {
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit1);
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit2);
            }
        }

        lev1 = MIN(MAX(codingdata[1], 0), 3);
        pki  = ari_spec_lookup_fl[codingdata[0] + lev1 * 1024];

        symfreq = ari_spec_cumfreq_fl[pki][codingdata[2] + 1] - ari_spec_cumfreq_fl[pki][codingdata[2]];
        cumfreq = ari_spec_cumfreq_fl[pki][codingdata[2]];
        ac_encode_fl(&st, symfreq, cumfreq);

        a = abs_x_k;
        b = abs_x_kp1;

        if (lsbMode == 1 && codingdata[1] > 0) {
            a             = a >> 1;
            lsbs[lsbsLen] = lsb1;
            lsbsLen++;

            if (a == 0 && x[k] != 0) {
                bit           = MAX(0, -sign(x[k]));
                lsbs[lsbsLen] = bit;
                lsbsLen++;
            }

            b             = b >> 1;
            lsbs[lsbsLen] = lsb2;
            lsbsLen++;

            if (b == 0 && x[k + 1] != 0) {
                bit           = MAX(0, -sign(x[k + 1]));
                lsbs[lsbsLen] = bit;
                lsbsLen++;
            }
        }

        if (a != 0) {
            bit = MAX(0, -sign(x[k]));
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit);
        }

        if (b != 0) {
            bit = MAX(0, -sign(x[k + 1]));
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, bit);
        }

        codingdata += 3;
    }

    /* Residual bits */
        nbits_side = total_bits - (8 * (*(st.bp_side) + 1) +  8 - (31 - clz_func(*(st.mask_side))));
        nbits_ari  =               8 * (st.bp      + 1) + 25 - (31 - clz_func(st.range    )) ;

    if (st.cache >= 0) {
        nbits_ari = nbits_ari + 8;
    }

    if (st.carry_count > 0) {
        nbits_ari = nbits_ari + st.carry_count * 8;
    }

    nbits_residual_enc = MAX(total_bits - (nbits_side + nbits_ari), 0);
    /* the max operation avoids in very rare cases, that
    * nbits_residual_enc becomes negative; having overwritten
    * the last bit(s) of the side information is in this case
    * assumed to be not critical, since no spectral data bits
    * were written */

    if (lsbMode == 0) {
        nbits_residual_enc = MIN(nbits_residual_enc, resBitsLen);
        for (k = 0; k < nbits_residual_enc; k++) {
            if (res_bits[k >> 3] & (1 << (k & 7)))
            {
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, 1);
            }
            else
            {
                write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, 0);
            }
        }
    } else {
        nbits_residual_enc = MIN(nbits_residual_enc, lsbsLen);

        for (k = 0; k < nbits_residual_enc; k++) {
            write_bit_backward_fl(st.ptr, st.bp_side, st.mask_side, lsbs[k]);
        }
    }

    ac_finalize_fl(&st);
#ifdef WMOPS
    pop_wmops();
#endif
}

