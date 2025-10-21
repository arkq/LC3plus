/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef CR9_C_ADD_1p25MS_LRSNS

static LC3_INT16 get_lead_sign( LC3_UINT32* ind_in );
static void mind2vec_one( LC3_UINT8 k_val_in, LC3_INT16 leading_sign, LC3_INT32* vec_out );

LC3_INT16 get_lead_sign( LC3_UINT32* ind_in )
{
    LC3_INT16 leading_sign;
    leading_sign = +1;
    if ( ( ( *ind_in ) & 0x1 ) != 0 )
    {
        leading_sign = -1;
    }
    ( *ind_in ) = ( *ind_in >> 1 );

    return leading_sign;
}

void mind2vec_one( LC3_UINT8 k_val_in,     /* i:  nb unit pulses  */
                   LC3_INT16 leading_sign, /* i: leading sign  -1, 1 */
                   LC3_INT32* vec_out /* o:  updated pulse train   */ )
{
        LC3_INT32 amp;
    amp = k_val_in;
    if ( leading_sign < 0 )
    {
        amp = -k_val_in;
    }
    *vec_out = amp;
}

static LC3_UINT8 setval_update_sign( LC3_INT16 k_delta, LC3_UINT8 k_max_local_in, LC3_INT16* leading_sign, LC3_UINT32* ind_in, LC3_INT32* vec_out );

static void mind2vec_tab( LC3_UINT8 dim_in, LC3_UINT8 k_max_local, LC3_INT16 leading_sign, LC3_UINT32 ind, LC3_INT32* vec_out );

#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
static void MPVQdeenum( LC3_UINT8 dim_in, LC3_UINT8 k_val_in, LC3_INT32 LS_ind, LC3_INT32 MPVQ_ind, LC3_INT32* vec_out );
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS

LC3_UINT8 setval_update_sign( LC3_INT16 k_delta,        /* i  */
                              LC3_UINT8 k_max_local_in, /* i  */
                              LC3_INT16* leading_sign,  /* i/o */
                              LC3_UINT32* ind_in,       /* i/o */
                              LC3_INT32* vec_out /* i/o */ )
{
    LC3_UINT8 k_max_local_out;
    k_max_local_out = k_max_local_in;
    if ( k_delta != 0 )
    {
        mind2vec_one( k_delta, *leading_sign, vec_out );
        *leading_sign = get_lead_sign( ind_in );
        k_max_local_out -= k_delta;
    }

    return k_max_local_out;
}

#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS

void mind2vec_tab( LC3_UINT8 dim_in,       /* i:  dimension        */
                   LC3_UINT8 k_max_local,  /* i:  nb unit pulses   */
                   LC3_INT16 leading_sign, /* i:  leading sign     */
                   LC3_UINT32 ind,         /* i:  MPVQ-index       */
                   LC3_INT32* vec_out /* o:  pulse train      */ )
{
     /* pvq_enc_A = MPVQ_offsets */
        LC3_UINT8 pos, k_acc;
        LC3_UINT32 UL_tmp_offset;
        LC3_INT32 UL_diff;
        LC3_INT16 wrap_flag, k_delta;
        const LC3_UINT32* h_row_ptr;
    /* init */
    h_row_ptr = &( pvq_enc_A[( dim_in - 1 )][0] );
    k_acc = k_max_local;

    /* loop over positions */
    for ( pos = 0; pos < dim_in; pos++ )
    {
        if ( ind != 0 )
        {
            k_acc = k_max_local;
            UL_tmp_offset = h_row_ptr[k_acc];

            wrap_flag = ( ind < UL_tmp_offset );
            UL_diff = (LC3_INT32) ( (LC3_INT32) ind - (LC3_INT32) UL_tmp_offset );

            while ( wrap_flag != 0 )
            {
                k_acc--;
                wrap_flag = ( ind < h_row_ptr[k_acc] );
                UL_diff = (LC3_INT32) ( (LC3_INT32) ind - (LC3_INT32) h_row_ptr[k_acc] );
            }
            ind = UL_diff;
            k_delta = k_max_local - k_acc;
        }
        else
        {
            mind2vec_one( k_max_local, leading_sign, &vec_out[pos] );
            break;
        }
        k_max_local = setval_update_sign( k_delta, k_max_local, &leading_sign, &ind, &vec_out[pos] );
        h_row_ptr -= 11; /* reduce dimension in MPVQ_offsets table */
    }
}

#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS

void MPVQdeenum( LC3_UINT8 dim_in,   /* i :  dimension of vec_out     */
                 LC3_UINT8 k_val_in, /* i :  number of unit pulses    */
                 LC3_INT32 LS_ind,   /* i :  leading sign index       */
                 LC3_INT32 MPVQ_ind, /* i :  MPVQ shape index         */
                 LC3_INT32* vec_out /* o :  PVQ integer pulse train  */ )
{
        LC3_INT32 leading_sign;


    leading_sign = 1;
    if ( LS_ind != 0 )
    {
        leading_sign = -1;
    }

    mind2vec_tab( dim_in, k_val_in, leading_sign, MPVQ_ind, vec_out );
}

/* local funcs*/

static void FESSdeenum(LC3_UINT8 dim_in,    /* i :  dimension of vec_out     */
    LC3_UINT8 n_env,     /* i :  number envelopes    */
    LC3_UINT8 n_shift,   /* i :  number shifts        */
    LC3_UINT8 n_signs,   /* i :  number signs          */
    LC3_INT32 env_ind,   /*  i:indx */
    LC3_INT32 shift_ind, /*  i:indx */
    LC3_INT32 sign_ind,  /*  i:indx */
    LC3_INT32* vec_out /* o :  FESS  integer  pulse train  */);

LC3_INT32 snsQuantScfEncLRSt1ABC(LC3_FLOAT* env, LC3_INT32* L_index, LC3_FLOAT *min_mse_saveBCA_ptr,
    LC3_INT32* ind_saveB_ptr, LC3_FLOAT* st1_vectors,
    LC3_INT32 pitch_rx, LC3_INT32 ltpf_rx);
 
#endif /*  CR9_C_ADD_1p25MS_LRSNS*/

static void pvq_dec(LC3_INT k, LC3_INT m, LC3_INT LS_ind, LC3_INT MPVQ_ind, LC3_INT* pulses);
static LC3_INT  find_last_indice_le(LC3_INT compare, const LC3_UINT32* array, LC3_INT len);
static void idct_II(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT len);

void idct_II(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT len)
{
    LC3_INT i;
    LC3_FLOAT norm1, sum;

    norm1 = 0.353553390593274; /* sqrt(2 / 16) */
    
    for (i = 0; i < len; i++) {
        sum = mac_loop(in, idct_lookup[i], len);
        out[i] = norm1 * sum;
    }
}

#ifdef  FIX_FLOAT_ENC_PVQ_PULSE_LOOP 
static LC3_INT32 pvq_pulse_search_lc(LC3_FLOAT* xabs, LC3_FLOAT* ener, LC3_FLOAT* corr, LC3_INT32* y, LC3_INT32 start, LC3_INT32 end)
{
    Dyn_Mem_Deluxe_In(
    LC3_INT32 i;
    LC3_INT32 nBest;
    LC3_FLOAT max_norm_ratio, cand_norm_ratio;
    LC3_FLOAT currCorr;
    LC3_INT32 currEnInt;
    const LC3_FLOAT* isqrt_QxTab;
    );

    const LC3_INT16 isqrt_Q16tab[1 + 64] = {     /* table generated using ISqrt16 function + shift to Q16 */
                                                32767, 32767, 32767, 32767, 32766, 29308, 26754, 24770, 23169, 21844,
                                                20723, 19759, 18918, 18176, 17515, 16921, 16383, 15894, 15446, 15034,
                                                14654, 14300, 13972, 13665, 13377, 13107, 12852, 12612, 12385, 12169,
                                                11965, 11770, 11584, 11407, 11238, 11077, 10922, 10773, 10631, 10493,
                                                10361, 10234, 10112, 9993, 9879, 9769, 9662, 9559, 9459, 9362,
                                                9268, 9176, 9088, 9001, 8918, 8836, 8757, 8680, 8605, 8532,
                                                8460, 8390, 8323, 8256, 8191
    };
     LC3_FLOAT isqrt_Q16tabFlt[1 + 64];

    const LC3_INT16 isqrt_Q15tab[1 + 6] = {
        /* onepulse_search_tab based search , Q15 table generated using isqrt_Q16tab table for idx=[8(4*2),16(4*4), 20(4*5),24(4*6)  ]    */
        /* this slighly suboptimal inv_sqrt(x) table is used to enable optional exact use of the ROM saving BASOP ISqrt16() function */

                                                   32767/*0*/, 32766/*1*/, 23169 /*2*/, 18918 /*3*/ , 16384/*4*/ , 14654 /*5*/,  13377/*6*/
                                                   /* i.e. value in isqrt_Q15tab[n=0..6]  ==  isqrt_Q16tab[4*n]  */
    };
   LC3_FLOAT isqrt_Q15tabFlt[1 + 6];

    for ( i = 0; i <= 64; i++ ) 
    {
        isqrt_Q16tabFlt[i] = isqrt_Q16tab[i] / 65536.0; 
    }
    for (i = 0; i <= 6; i++)
    {
        isqrt_Q15tabFlt[i] = isqrt_Q15tab[i] / 32768.0;
    }


    isqrt_QxTab = isqrt_Q16tabFlt;    /* Q16 table valid for energies 4...64  */
  
    assert( *ener >= 0.0 );
    *ener += 1.0;  /* Added once for the entire loop */
    if (*ener  <= 3.0 ) /* +1 for the inloop value of *ener  was preadded before */
    {
        isqrt_QxTab = isqrt_Q15tabFlt;      /* energies:  1...6   */
    }

    nBest = -1;
    max_norm_ratio = -1.0; 
    /* Iterative max search using tabulated  inv sqrt  */
    for (  i = start; i < end; i++)
    {
        currCorr = *corr + xabs[i];
        currEnInt  = ((LC3_INT32) *ener)  + (2 * y[i]);

        cand_norm_ratio =  currCorr *  isqrt_QxTab[currEnInt];

        if ( cand_norm_ratio >=  max_norm_ratio ) 
        {
            nBest = i;
        }
        max_norm_ratio = LC3_FMAX(max_norm_ratio, cand_norm_ratio); /* always update */
    }

    *corr += xabs[nBest];
    *ener += (2 * y[nBest]);

    assert(nBest >= 0);
    assert(nBest <=M);
    assert(*ener <= 64.0);
    y[nBest] += 1; /* Add the selected unit pulse */

    Dyn_Mem_Deluxe_Out();
    return nBest;
}
#endif 

static LC3_INT pvq_pulse_search(LC3_FLOAT *xabs, LC3_FLOAT *ener, LC3_FLOAT *corr, LC3_INT *y, LC3_INT start, LC3_INT end)
{
    LC3_INT i;
    LC3_INT nBest;
    LC3_FLOAT bestCorrSq, bestEn;
    LC3_FLOAT corrSq, currCorr, currEn;

    nBest = 0;
    bestCorrSq = 0.0;
    bestEn = 0.0;

    *ener += 1; /* Added once for the entire loop */

    i = start;

    currCorr = *corr + xabs[i];
    currEn = *ener + (2 * y[i]);

    corrSq = currCorr * currCorr;

    bestEn = currEn;
    bestCorrSq = corrSq;
    nBest = i;

    /* Iterative max search as recommended in the spec */
    for (; i < end; i++)
    {
        currCorr = *corr + xabs[i];
        currEn = *ener + (2 * y[i]);

        corrSq = currCorr * currCorr;

        if ((corrSq * bestEn) > (bestCorrSq * currEn))
        {
            bestEn = currEn;
            bestCorrSq = corrSq;
            nBest = i;
        }
    }

    *corr += xabs[nBest];
    *ener += (2 * y[nBest]);

    y[nBest] += 1; /* Add the selected unit pulse */

    return nBest;
}

static void pvq_enc_vec_normalize(LC3_FLOAT *vec, LC3_INT N)
{
    LC3_FLOAT mag = 0.0, norm_fac;
    LC3_INT i;

    for (i = 0; i < N; i++)
    {
        mag += (vec[i] * vec[i]);
    }

    norm_fac = 1.0 / LC3_SQRT(mag);

    for (i = 0; i < N; i++)
    {
        vec[i] = vec[i] * norm_fac;
    }

    return;
}

static void pvq_enc_search(LC3_FLOAT* x_in, LC3_INT y[4][M])
{
    LC3_INT i, N, K, pulse_total, N_setA;
    LC3_FLOAT abs_sum, projfac;
    LC3_FLOAT xabs[16];
    LC3_FLOAT yy, xy;

    abs_sum = 0.0;

    /* Step 1 : Projection to pyramid N=16, K=6 */
    N = 16;
    K = 6;
    pulse_total = 0;
    N_setA = 10;

    yy = xy = 0.0f;

    for (i = 0; i < N; i++)
    {
        xabs[i] = LC3_FABS(x_in[i]);
        abs_sum += xabs[i];
    }

    projfac = (K - 1) / abs_sum;

    for (i = 0; i < N; i++)
    {
        y[3][i] = floor(xabs[i] * projfac);

        pulse_total += y[3][i];

        yy += (y[3][i] * y[3][i]);
        xy += (xabs[i] * y[3][i]);
    }

    /* Step 2: Adding unit pulses up to K = 6 */
    for (; pulse_total < K; pulse_total++)
    {
        pvq_pulse_search(xabs, &yy, &xy, y[3], 0, N);
    }

    /* Step 3: Adding unit pulses up to K = 8 */
    memcpy(y[2], y[3], sizeof(LC3_INT)*N);
    K = 8;

    for (; pulse_total < K; pulse_total++)
    {
        pvq_pulse_search(xabs, &yy, &xy, y[2], 0, N);
    }

    memcpy(y[1], y[2], sizeof(LC3_INT)*N_setA);

    /* Step 4: Remove unit pulses not belonging to set A */
    for (i = N_setA; i < N; i++)
    {
        y[1][i] = 0;
    }

    /* Step 5: Update yy and xy terms to reflect y1 */
    yy = 0;
    xy = 0;
    pulse_total = 0;

    for (i = 0; i < N_setA; i++)
    {
        yy += (y[1][i] * y[1][i]);
        xy += (xabs[i] * y[1][i]);

        pulse_total += y[1][i];
    }

    /* Step 6: Add unit pulses until K = 10 over N = 10 */
    K = 10;
    for (; pulse_total < K; pulse_total++)
    {
        pvq_pulse_search(xabs, &yy, &xy, y[1], 0, N_setA);
    }

    memcpy(y[0], y[1], sizeof(LC3_INT)*N);

    /* Step 7: Add unit pulses until K = 1 over N = 6 in set B*/
    pvq_pulse_search(xabs, &yy, &xy, y[0], N_setA, N);

    /* Step 8: Add signs to each of the 4 vectors from x */
    for (i = 0; i < N; i++)
    {
        if (x_in[i] < 0)
        {
            y[0][i] = -y[0][i];
            y[1][i] = -y[1][i];
            y[2][i] = -y[2][i];
            y[3][i] = -y[3][i];
        }
    }

    return;
}

static LC3_FLOAT calc_mse(LC3_FLOAT *t2rot, LC3_FLOAT *y, LC3_FLOAT gain, LC3_INT N)
{
    LC3_FLOAT mse;
    LC3_INT i;

    mse = 0.0;

    for (i = 0; i < N; i++)
    {
        LC3_FLOAT err = (t2rot[i] - gain * y[i]);
        mse += (err * err);
    }

    return mse;
}

static void sns_quant_adj_gain_shape_search(LC3_FLOAT *t2rot, LC3_INT y[4][M] ,
    LC3_INT *gain_idx, LC3_INT *shape_idx, LC3_FLOAT *y_norm, LC3_FLOAT *scq_gain)
{
    LC3_INT gidx, sidx;
    LC3_FLOAT min_mse, mse;
    LC3_INT N;
    LC3_FLOAT yCur[4][16];
    LC3_INT i;

    const LC3_INT gain_levels[4] = { 2, 4, 4, 8 };
    const LC3_FLOAT *sns_vq_gains[4] = { sns_vq_reg_adj_gains_fl , sns_vq_reg_lf_adj_gains_fl ,
        sns_vq_near_adj_gains_fl , sns_vq_far_adj_gains_fl };

    min_mse = -1.0;
    N = 16;


    *gain_idx = *shape_idx = 0;

    for (sidx = 0; sidx < 4; sidx++)
    {
        for (i = 0; i < N; i++)
        {
            yCur[sidx][i] = (LC3_FLOAT)y[sidx][i];
        }

        /* Step 9: Normalize the vectors */
        pvq_enc_vec_normalize(yCur[sidx], N);

        for (gidx = 0; gidx < gain_levels[sidx]; gidx++)
        {
            mse = calc_mse(t2rot, yCur[sidx], sns_vq_gains[sidx][gidx], N);

            if ((mse < min_mse)  || (min_mse < 0))
            {
                *gain_idx = gidx;
                *shape_idx = sidx;
                min_mse = mse;
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        y_norm[i] = yCur[*shape_idx][i];
    }

    *scq_gain = sns_vq_gains[*shape_idx][*gain_idx];

    return;
}

static void enc_push_sign(LC3_FLOAT val, LC3_UINT32 *next_sign_ind, LC3_INT *index)
{
    if (((*next_sign_ind & 0x80000000U) == 0) && (val != 0)) {
        *index = 2 * (*index) + *next_sign_ind;
    }
    if (val < 0) {
        *next_sign_ind = 1;
    }
    if (val > 0) {
        *next_sign_ind = 0;
    }

    return;
}

static void MPVQ_enum(LC3_INT dim, LC3_INT *sns_vec, LC3_INT *index_val, LC3_INT *lead_sign_ind)
{
    LC3_UINT32 next_sign_ind;
    LC3_INT k_val_acc;
    LC3_INT pos;
    LC3_INT index, n;
    LC3_INT const *row_ptr;

    /* MPVQ-index composition loop */
    LC3_INT tmp_h_row;
    LC3_INT tmp_val;

    next_sign_ind = 0x80000000U;
    k_val_acc = 0;
    pos = dim;
    index = 0;
    n = 0;

    row_ptr = (LC3_INT const *)&(pvq_enc_A[n]);
    tmp_h_row = row_ptr[0];

    for (pos--; pos >= 0; pos--)
    {
        tmp_val = sns_vec[pos];
        enc_push_sign(tmp_val, &next_sign_ind, &index);

        index += tmp_h_row;
        k_val_acc += abs(tmp_val);
        if (pos != 0) {
            n += 1; /* switch row in offset table MPVQ_offsets(n, k) */
        }
        row_ptr = (LC3_INT const *)&(pvq_enc_A[n]);

        tmp_h_row = row_ptr[k_val_acc];
    }

    *index_val = index;
    *lead_sign_ind = next_sign_ind;

    return;
}

static LC3_INT MSEsearch (LC3_FLOAT *scf, const LC3_FLOAT sns_CB[8][32])
{
    LC3_FLOAT distance, mse;
    LC3_INT i, n, ind;

    ind = 0;

    distance = (LC3_FLOAT) LC3_CONST_POW_2_100;
    for (i = 0; i < 32; i++) {
        mse = 0;
        for (n = 0; n < 8; n++) {
            mse += (scf[n] - sns_CB[n][i]) * (scf[n] - sns_CB[n][i]);
        }

        if (mse < distance) {
            distance = mse;
            ind      = i;
        }
    }
    return ind;
}

void process_snsQuantizesScf_Enc(LC3_FLOAT* env, LC3_INT* index, LC3_FLOAT* envq, Dct2 dct2structSNS)
{
    LC3_FLOAT stage2_en1_norm_sub[M];
    LC3_INT i, j;
    LC3_FLOAT st1_vector[M];
    LC3_FLOAT pvq_target_pre[M];
    LC3_FLOAT pvq_target[M];
    LC3_FLOAT stage2_en1_norm_pre_sub[M];
    LC3_INT gain, shape;
    LC3_FLOAT scfq_gain;
    LC3_INT y[4][M];

    /* Stage 1 split VQ */
    index[0] = MSEsearch(&env[0], sns_LFCB);  /* ind_LF */
    index[1] = MSEsearch(&env[8], sns_HFCB);  /* ind_HF */

    j = 8;
    for (i = 0; i < 8; i++, j++) {
        st1_vector[i] = sns_LFCB[i][index[0]];
        st1_vector[j] = sns_HFCB[i][index[1]];
    }

    /* STAGE 2 */
    for (i = 0; i < 16; i++) {
        pvq_target_pre[i] = env[i] - st1_vector[i];
    }

    dct2_apply(&dct2structSNS, pvq_target_pre, pvq_target);
    pvq_enc_search(pvq_target, y);
    sns_quant_adj_gain_shape_search(pvq_target, y, &gain, &shape, stage2_en1_norm_pre_sub, &scfq_gain);

    /* Inverse transform */
    idct_II(stage2_en1_norm_pre_sub, stage2_en1_norm_sub, M);

    index[2] = shape;
    index[3] = gain;

    if (shape < 2) {
        MPVQ_enum(10, y[shape], &index[5], &index[4]);
    }
    else {
        MPVQ_enum(M, y[shape], &index[5], &index[4]);
    }

    if (shape == 0) {
        LC3_INT ls_ind, ind;
        MPVQ_enum(6, &y[shape][10], &ind, &ls_ind);
        index[6] = ind * 2 + ls_ind;
    }
    else if (shape == 2) {
        index[6] = -1;
    }
    else {
        index[6] = -2;
    }

    for (i = 0; i < M; i++) {
        envq[i] = st1_vector[i] + (stage2_en1_norm_sub[i] * scfq_gain);
    }
}

LC3_INT find_last_indice_le(LC3_INT compare, const LC3_UINT32* array, LC3_INT len)
{
    LC3_INT idx = 0, i = 0;

    for (i = 0; i < len; i++) {
        if ((LC3_UINT32)compare >= array[i]) {
            idx++;
        }
    }

    if (idx > 0) {
        idx--;
    }

    return idx;
}

void pvq_dec(LC3_INT k, LC3_INT m, LC3_INT LS_ind, LC3_INT MPVQ_ind, LC3_INT* pulses)
{
    LC3_INT leading_sign, idx, k_delta = 0, pos;

    leading_sign = 1 - 2 * LS_ind;

    /* Decoding loop */

    for (pos = 0; pos < m; pos++) {
        if (MPVQ_ind != 0) {
            /* Find last indice */
            idx      = find_last_indice_le(MPVQ_ind, &pvq_enc_A[m - pos - 1][0], k + 1);
            MPVQ_ind = MPVQ_ind - pvq_enc_A[m - pos - 1][idx];
            k_delta  = k - idx;
        } else {
            pulses[pos] = leading_sign * k;
            break;
        }

        if (k_delta != 0) {
            pulses[pos] = leading_sign * k_delta;
            if ((MPVQ_ind % 2) != 0) {
                leading_sign = -1;
            } else {
                leading_sign = 1;
            }

            MPVQ_ind = floor(MPVQ_ind / 2);
            k        = k - k_delta;
        }
    }
}

void process_snsQuantizesScf_Dec(LC3_INT* scf_idx, LC3_FLOAT* scf_q)
{
    LC3_INT   i, submode;
    LC3_INT   pulses2[6] = {0}, pulses[M] = {0};
    LC3_FLOAT st2_vector[M], st2_vector_idct[M], sum = 0;

    /* Decode first stage */

    for (i = 0; i < 8; i++) {
        scf_q[i]     = sns_LFCB[i][scf_idx[0]];
        scf_q[i + 8] = sns_HFCB[i][scf_idx[1]];
    }

    /* STAGE 2 */
    /* Decode submode */

    submode = scf_idx[2];

    /* Decode pulses */

    if (submode < 2) {
        pvq_dec(10, 10, scf_idx[4], scf_idx[5], pulses);

        if (submode == 0) {
            pvq_dec(1, 6, (scf_idx[6] % 2), floor(scf_idx[6] / 2), pulses2);

            move_int(&pulses[10], pulses2, 6);

        } else {
            pulses[15] = 0;
        }
    } else if (submode == 2) {
        pvq_dec(8, 16, scf_idx[4], scf_idx[5], pulses);
    } else {
        pvq_dec(6, 16, scf_idx[4], scf_idx[5], pulses);
    }

    /* Normalization */

    for (i = 0; i < M; i++) {
        sum += pulses[i] * pulses[i];
    }

    sum = 1.0 / LC3_SQRT(sum);

    for (i = 0; i < M; i++) {
    st2_vector[i] = pulses[i] * sum;
    }

    /* Inverse transform */
    idct_II(st2_vector, st2_vector_idct, M);

    /* Gain */
    /* Add stage 1 and stage 2 */
    for (i = 0; i < M; i++) {
        scf_q[i] += st2_vector_idct[i] * sns_dec_gains[submode][scf_idx[3]];
    }
}

#ifdef  CR9_C_ADD_1p25MS_LRSNS

/* 29/30 bits  optimized search functions for PVQ and FESS  */
/* stage 2 submode shape  0:  "splitLF"  (N=5,K=6)(N=8,K=2)   or (N=5,K=8)(N=8,K=0) , 4 gains */
/* stage 2 submode shape  1:  "full" FB (N=15,K=5), 4xfixed  , 8 gains */
/* stage 2 submode shape  2-5: "fixed env "  (N=13-15,K=10-12), 4xfixed  , 8 gains */
 
static void pvq_fess_enc_search(LC3_FLOAT* x_in    /* i:   0...M-1 */
    , LC3_INT32 y[SNSLR_MAX_PVQ_SEARCH_CAND][M],   /* o:  [3]*[0...M-1] */
    LC3_INT32  *fixShapeNb                         /* o:   [-1, 0...3]  */
)
{
    LC3_INT32 i, j, curr_cand;
    LC3_INT32 NcandLF, Ngrp, Kgrp, Ktop_LF, Ntop_LF;
    LC3_INT32 NcandFB, Kfull;
    LC3_INT32 NcandFix, Kfix, Nsigns_fix;

    LC3_INT32 N, K, pulse_total, top_pulses;
    LC3_FLOAT abs_sum, projfac;
    LC3_FLOAT yy, xy;    /* in-loop energy  and corr  */

    LC3_FLOAT  yy_n5k8, xy_n5k8;
    LC3_FLOAT  xabs[M];
    LC3_INT32  y_tmp_n5k8[M];
    LC3_INT32  y_tmp[M];
    LC3_INT16  n, fix_ind, shift_ind;
    LC3_INT16 best_env_ind;
    LC3_INT16 best_shift_ind;
    LC3_INT16 best_ind;
    LC3_INT16 pulse_total_n5k8;
    LC3_FLOAT f_normcorr_fixenv[SNSLR_N_FIXENV][SNSLR_N_FIXENV_SHIFTS];
    LC3_FLOAT *p;
    const LC3_INT32* env_ptr;

#define SC 1 /*Starting coeff position of PVQ target  , (we have no DC) */
    {
        /*  no DC coeff , only coeff 1..15 is quantized in  the  range 0..15  */
        NcandLF = 1;

        Ngrp = 5;      /* PVQ(5, 6) = 10.94b , + P(8,2)=7b  + P(2,0)=0b */
        Kgrp = 6;
        Ntop_LF = 8;   /*   allow to encode less than M-sc -Ngrp  */
        Ktop_LF = 2;   /*   fill a few unit pulses in the HF region   */
                       /*   P(5,8)+P(10,0)  will also fit in  13 bits */
        /*  Fullband safety net */
        NcandFB = 1;
        Kfull = 5; ;      /* PVQ(15,5)= 16.65   bits  ,   */

        /*fixed env */
        NcandFix = 4;         /* 4th only has 10 signs */
        Kfix = 12;            /* number of signs  */
        Nsigns_fix = 12;       /* area to cover and  number of signs */
    }

    for (i = SC; i < (M); i++)
    {
        xabs[i] = LC3_FABS(x_in[i]);
    }
    xabs[0] = 0;

    /*   splitLF LFfocus subshape idx 0   */
    if (Kgrp > 0 && NcandLF > 0)
    {
        curr_cand = 0;
        /*   Projection to LF pyramid N=Ngrp, K=Kgrp  */
        abs_sum = 0.0;
        N = Ngrp;
        K = Kgrp;
        pulse_total = 0;

        for (i = SC; i < (SC + N); i++)
        {
            abs_sum += xabs[i];
        }
        projfac = (K - 1) / abs_sum;

        if (abs_sum == 0.0)
        {
            projfac = 0.0;
        }

        yy = xy = 0.0f;
        for (i = SC; i < (SC + N); i++)
        {
            y_tmp[i] = LC3_FLOOR(xabs[i] * projfac);    /* local projection within group over N coeffs used here  */

            pulse_total += y_tmp[i];

            yy += (y_tmp[i] * y_tmp[i]);
            xy += (xabs[i] * y_tmp[i]);
        }

        /* LF  Adding unit pulses up to K in LF */
        for (; pulse_total < K; pulse_total++)
        {
            pvq_pulse_search(xabs, &yy, &xy, y_tmp, SC, SC + N);
        }

        {
            /* allow (N=5,K=(6+2)=8) in low band only,  higher part is then zeroed  */
            /* can be multiplexed put in 13 bits */
            yy_n5k8 = yy;
            xy_n5k8 = xy;
            pulse_total_n5k8 = pulse_total;

            memset(y_tmp_n5k8, 0, sizeof(*y_tmp_n5k8)*M);
            memcpy(y_tmp_n5k8, y_tmp, sizeof(*y_tmp)*(SC + N));  /* cpy LF initial search result for n5k6 to n5k8 output  */

            for (; pulse_total_n5k8 < (K + 2); pulse_total_n5k8++)
            {
                pvq_pulse_search(xabs, &yy_n5k8, &xy_n5k8, y_tmp_n5k8, SC, SC + N);
            }
        }


        /* split  HF    :   longer M vector  */
        memset(y[curr_cand], 0, sizeof(*y_tmp)*M);         /*zero full candidate vector  */
        memcpy(&(y[curr_cand]), y_tmp, sizeof(*y_tmp)*(SC + N));  /* cpy LF-part   to output  */


        if (Ktop_LF && (M - Ngrp - SC >= Ntop_LF))
        {

            /* xy = 0.0; yy = 0.0;   */
            memset(&(y_tmp[Ngrp]), 0, (M - SC - Ngrp) * sizeof(*y_tmp));
            top_pulses = 0;

            for (; top_pulses < Ktop_LF; top_pulses++)
            {
                pvq_pulse_search(xabs, &yy, &xy, y_tmp, SC + N, SC + N + Ntop_LF);
            }

            for (j = SC + N; j < (SC + N + Ntop_LF); j++)
            {
                y[curr_cand][j] = y_tmp[j];
            }
        }

        {
            /* shape 0 : decide on best split alternative  "0a" n5k6 + n8k2 + n2k0   or  "0b"  n5k8 + n10k0,
              shape only assessed , i.e. assume same  gain quantizer */
            if ((xy_n5k8*xy_n5k8)*yy > (xy*xy)*yy_n5k8)
            {
                /* use  n5k8, with zeroed high band */
                /* i.e. better to put two additional pulses in the n5 low band than in the higher band */
                memcpy(&(y[curr_cand]), y_tmp_n5k8, sizeof(*y_tmp_n5k8)*(M));   /* cpy LFonly   to output  */
            }
        }

        /*  assign signs to the cand  vector  from x_in */
        for (i = SC; i < (M); i++)
        {
            if (x_in[i] < 0)
            {
                y[curr_cand][i] *= -1;
            }
        }
        y[curr_cand][0] = 0;  /* DC */

    }

    /* shape 1 :  Full band shape analysis      */
    if (Kfull > 0 && NcandFB > 0)
    {
        curr_cand = 1;
        N = M - SC;
        K = Kfull;

        /* project */
        pulse_total = 0;
        yy = xy = 0.0f;
        abs_sum = 0.0;
        for (i = SC; i < SC + N; i++)
        {
            abs_sum += xabs[i];/* over N */
        }

        projfac = (K - 1) / abs_sum;
        if (abs_sum == 0.0)
        {
            projfac = 0.0;
        }


        memset(y_tmp, 0, sizeof(*y_tmp)*M);
        for (i = SC; i < SC + N; i++)
        {
            y_tmp[i] = LC3_FLOOR(xabs[i] * projfac);    /*  projection over N coeffs  */
            pulse_total += y_tmp[i];
            yy += (y_tmp[i] * y_tmp[i]);
            xy += (xabs[i] * y_tmp[i]);
        }

        for (; pulse_total < K; pulse_total++)
        {
            pvq_pulse_search(xabs, &yy, &xy, y_tmp, SC, SC + N);
        }
        memcpy(&(y[curr_cand]), y_tmp, sizeof(*y_tmp)*M);  /* cpy to output  */

       /* Step: Add signs to each of the full  vector from input  */
        for (i = SC; i < M; i++)
        {
            if (x_in[i] < 0)
            {
                y[curr_cand][i] *= -1; ;
            }
        }
        y[curr_cand][0] = 0;
    } /* Fullband  */


    { /* fixed envelopes shifted signs shapes 2 ... 5 */
        best_env_ind = -1;
        best_shift_ind = -1;
 
        /* search fixed env shapes first */
         /* search shape of all fixed flat-ish envelopes using an optimal  nb_env*nbshifts shape analysis */
        if (Kfix > 0 && NcandFix > 0)
        {
            curr_cand = 2;
            N = Nsigns_fix;
            K = Kfix;

            /* find the minimum shape error across all possible 4  fixed envelopes and all 4 shifts */
            /* maximise normalized cross correlation  target*y  *(1/sqrt(y^2))  to minimize shape error */

            for (fix_ind = 0; fix_ind < SNSLR_N_FIXENV; fix_ind++)
            {
                for (shift_ind = 0; shift_ind < SNSLR_N_FIXENV_SHIFTS; shift_ind++)
                {
                    f_normcorr_fixenv[fix_ind][shift_ind] = 0.0;
                    for (n = 0; n < signs_fix[fix_ind]; n++)
                    {
                        f_normcorr_fixenv[fix_ind][shift_ind] += (xabs[SC + n + shift_ind] * env_ptrs[fix_ind][n + shift_ind]);
                    }
                    /*   energy normalize, for the specific shift and env,
                        otherwise it is not possible to compare among the shifted envelopes properly */

                    f_normcorr_fixenv[fix_ind][shift_ind] *= shift_en_norm_factors[fix_ind][shift_ind];
                }
            }

            /* now actually pick the env[0,1,2,3] and env shift[0,1,2,3] option which maximizes the normcorr
               (minimizes shape mse as all fixed envelopes have the same gain quantizer) */
            p = &(f_normcorr_fixenv[0][0]);
            best_ind = 0;
            for (n = 1; n < (SNSLR_N_FIXENV * SNSLR_N_FIXENV_SHIFTS); n++)
            {
                if (p[n] > p[best_ind])
                {
                    best_ind = n;
                }
            }
            best_env_ind = best_ind / SNSLR_N_FIXENV_SHIFTS;
            best_shift_ind = best_ind - best_env_ind * SNSLR_N_FIXENV;
            assert(best_env_ind >= 0 && best_env_ind < SNSLR_N_FIXENV);
            assert(best_shift_ind >= 0 && best_shift_ind < SNSLR_N_FIXENV_SHIFTS);
        }
        /* o: best_shift_ind; o: best_env_ind */

        *fixShapeNb = best_env_ind; /* best normalized correlation out of the 4x4 = 16 envelope options */
        /*best shift ind, later re-established/found via output vector y[2] */

        /* Fixed envelope "flat-ish" signband coding , including  sign coding  of shifted  block  */
        /* submode idx 2,3,4     , "1"/env    1(s0)+ 2(shift)+ 11bits(s1..s11)   +   always a 3 bits gain  */
        /* submode idx  5        , "1"/env    1(s0)+ 2(shift)+  9bits(s1..s9)    +   always a 3 bits gain  */
        /*
         2,  init_bell_12signs , [ 8,8,8, 7,7 ... ]
         3,  decaying envelope 12 signs, [ 12,12,11,11, ... ]
         4,  start_bell_12signs ,[ 7,7,8,8,8, 7,... ]
         5,  early_bell_10signs,[ 6,6, 7,7,8,8,8,7,... ]
       */

        {
            /* construct the selected fix shape with its corresponding shift  */
            memset(y[curr_cand], 0, sizeof(LC3_INT32)*M);    /* zero  output  */
            env_ptr = env_ptrs[best_env_ind];

            Nsigns_fix = signs_fix[best_env_ind];

            /*  now assign unit amplitude and signs to each of the  fixed env  vector from x_in  in the band */
            for (i = best_shift_ind; i < (Nsigns_fix + best_shift_ind); i++)
            {
                y[curr_cand][SC + i] = env_ptr[i];
            }
            for (i = (SC + best_shift_ind); i < (SC + Nsigns_fix + best_shift_ind); i++)
            {
                if (x_in[i] < 0)
                {
                    y[curr_cand][i] *= -1;
                }
            }
            y[curr_cand][0] = 0;
        }
    }

    return;
}

static void lrsns_quant_shape_gain(
    LC3_FLOAT* t2rot, LC3_INT32 y[SNSLR_MAX_PVQ_SEARCH_CAND][M],
    LC3_INT32* gain_idx, LC3_INT32 *shape_idx, LC3_FLOAT* y_norm, LC3_FLOAT* scq_gain, LC3_INT32 n_cand )
{
        LC3_INT32 gidx, sidx;
    LC3_FLOAT min_mse, mse;
 
    LC3_INT32 i;
    LC3_INT16  start_shape;
    LC3_INT16  end_shape;
    LC3_INT32 gain_levels[SNSLR_MAX_PVQ_SEARCH_CAND];   
    LC3_FLOAT yCur[SNSLR_MAX_PVQ_SEARCH_CAND][M];
    const LC3_FLOAT *lrsns_vq_gains[SNSLR_MAX_PVQ_SEARCH_CAND];

    gain_levels[0] = 4;     /* splitLF */
    gain_levels[1] = 8;     /* full     */
    gain_levels[2] = 8;     /* Fix-env{ 0,1,2,3,4}  */

    lrsns_vq_gains[0] = &(lrsns_gains_Q11[0][0]);
    lrsns_vq_gains[1] = &(lrsns_gains_Q11[1][0]);
    lrsns_vq_gains[2] = &(lrsns_gains_Q11[2][0]);
   
    min_mse = LC3_CONST_FLOATMAX;
   

    *gain_idx = 0;
    *shape_idx = 0;
    start_shape = 0;
    end_shape = n_cand;

    for (sidx = start_shape; sidx < end_shape; sidx++)
    {
        /*  Normalize the vectors */
        for (i = 0; i < M; i++)
        {
            yCur[sidx][i] = (LC3_FLOAT)y[sidx][i];
        }
        pvq_enc_vec_normalize(yCur[sidx], M);

        for (gidx = 0; gidx < gain_levels[sidx]; gidx++)
        {
            mse = calc_mse(t2rot, yCur[sidx], lrsns_vq_gains[sidx][gidx], M);

            if (mse < min_mse)
            {
                *gain_idx = gidx;
                *shape_idx = sidx;
                min_mse = mse;
            }
        }

    }

    for (i = 0; i < M; i++)
    {
        y_norm[i] = yCur[*shape_idx][i];
    }

    *scq_gain = lrsns_vq_gains[*shape_idx][*gain_idx];

    return;
}

LC3_INT32 MSEsearchGeneric(LC3_FLOAT *scf, const LC3_FLOAT *sns_CB, LC3_INT32 v_len, LC3_INT32 cb_len, LC3_FLOAT* min_mse)
{

        LC3_FLOAT f_tmp, mse;
    LC3_INT32 i, n, ind;

    ind = -1;

    *min_mse = (LC3_FLOAT)LC3_CONST_POW_2_100;
    for (i = 0; i < cb_len; i++)
    {
        mse = 0;
        for (n = 0; n < v_len; n++)
        {
            f_tmp = (scf[n] - sns_CB[i * v_len + n]);
            mse += (f_tmp * f_tmp);
        }

        if (mse < *min_mse)
        {
            *min_mse = mse;
            ind = i;
        }
    }
    assert(ind >= 0 && ind < cb_len);

    return ind;
}

void snslr_st1B_vector_dec(LC3_INT16 idx, const LC3_INT16* LFCB, const LC3_INT16 *HFCB, const LC3_INT16* seg_cnt_cum, const LC3_INT16* idx12b_cb, LC3_FLOAT *st1B_vector)
{
    /* decompose the received 0 ... 169 index , into the correct intger  and float st1B vector  */
    LC3_INT32 i, seg; /*counters*/
    const LC3_INT16 *lf_cb, *hf_cb;
    LC3_INT16 idx_12b, lf_sign, hf_sign;
    LC3_INT16 idx_LF, idx_HF;
    LC3_INT16 buf[M];
    LC3_INT16 st1B_W16Q11[M];
    LC3_FLOAT f_Q11_scale = (1.0 / 2048.0);

    assert(idx >= 0 && idx < 170);
    seg = 0;
    while (seg_cnt_cum[seg + 1] <= idx) {
        seg++;
    }
    assert(seg >= 0 && seg < 4);

    idx_12b = idx12b_cb[idx];  /* from sequential value to a coded 12b index*/

    lf_sign = 1; /*   assume a 0 bit -> "+"  */
    if ((0x0800 & idx_12b) != 0) {
        lf_sign = -1; /*   assume a 1 bit -> " -"  */
    }

    hf_sign = 1; /*   assume a 0 bit -> "+"  */
    if ((0x0400 & idx_12b) != 0) {
        hf_sign = -1; /*   assume a 1 bit -> "-"  */
    }
    idx_LF = (0x03e0 & idx_12b) >> 5;
    idx_HF = (0x001f & idx_12b);

    /* extseg0  f,f  */
    lf_cb = &(LFCB[idx_LF*(M / 2)]);
    hf_cb = &(HFCB[idx_HF*(M / 2)]);
    for (i = 0; i < (M / 2); i++)
    {
        st1B_W16Q11[i] = lf_sign * lf_cb[i];     /* imult()  or negate */
        st1B_W16Q11[M / 2 + i] = hf_sign * hf_cb[i];
    }
    memcpy(buf, st1B_W16Q11, sizeof(*buf)*M); /* buffer cpy needed for reversal sections */

    if ((seg & 0x0002) != 0)
    {   /* r,*  */     /* flip LF */
        for (i = 0; i < (M / 2); i++)
        {
            st1B_W16Q11[i] = buf[(M / 2 - 1) - i];
        }
    }

    if ((seg & 0x0001) != 0)
    {  /* *,r */ /* flip HF */
        for (i = 0; i < (M / 2); i++)
        {
            st1B_W16Q11[(M / 2) + i] = buf[(M - 1) - i];
        }
    }

    /* Cfloat: convert  the 16bit integer  Q11 values from LFCB, and HFCB into  floats */
    for (i = 0; i < M; i++)
    {
        st1B_vector[i] = ((LC3_FLOAT)st1B_W16Q11[i]) * f_Q11_scale;
    }
}

LC3_INT32 MSEsearchCbBIdxMap(const LC3_FLOAT *scf, const LC3_INT16 *LFCB, const LC3_INT16 *HFCB, const LC3_INT16 *seg_cnt_cum, const LC3_INT16* idx12b_cb, LC3_INT32 v_len, LC3_INT32 cb_len, LC3_FLOAT* min_mse)
{
    LC3_INT32  seg, i, j; /*counters */
    LC3_INT32  L_mse;
    LC3_INT16   scfLF_Q11[M], tmp_buf[M];
    LC3_INT16*  scfHF_Q11;
    const LC3_INT16 *lf_cb, *hf_cb;
    LC3_INT16  err, best_ind, idx_12b, signbitLF, signbitHF, idx_LF, idx_HF;
    LC3_INT32 L_mse_best = INT_MAX; /*a huge positive number */
    
    UNUSED(cb_len); /*cb_len only used for assert/verification */

    assert(v_len == M);

    scfHF_Q11 = (&scfLF_Q11[v_len / 2]); /* ptr init */

    /*set up   fwd,fwd search */
    for (i = 0; i < M; i++)
    {
        tmp_buf[i] = (LC3_INT16)LC3_ROUND(scf[i] * 2048.0);  /*   scf  target  moved to  BASOP signed Word16Q11 integer domain */
    }

    best_ind = -1;   /*for debug*/
    for (seg = 0; seg < 4; seg++)
    {
        memcpy(scfLF_Q11, tmp_buf, M * sizeof(*scfLF_Q11));  /*  M * move16() */
        /*seg==0: fwd, fwd */
        /*seg==1: fwd, rev */
        /*seg==2: fwd, fwd */
        /*seg==3: rev, rev */

        if ((seg & 0x0002) != 0)
        {   /* {r,*}  */     /* flip LF */
            for (i = 0; i < (M / 2); i++)
            {
                scfLF_Q11[i] = tmp_buf[(M / 2 - 1) - i];
            }
        }
        if ((seg & 0x0001) != 0)
        {  /* {*,r} */    /* flip HF */
            for (i = 0; i < (M / 2); i++)
            {
                scfLF_Q11[M / 2 + i] = tmp_buf[(M - 1) - i];
            }
        }

        for (i = seg_cnt_cum[seg]; i < seg_cnt_cum[seg + 1]; i++)
        {
            idx_12b = idx12b_cb[i];           /* indirect adressing lookup of  12b index pointing to LF and HF + individual sign swaps  */
            signbitLF = (0x0800 & idx_12b);   /* b11 logical          */
            signbitHF = (0x0400 & idx_12b);   /* b10 logical          */
            idx_LF = (0x03e0 & idx_12b) >> 5; /* b9...b5              */
            idx_HF = (0x001f & idx_12b);      /* b4..b0 lowest 5 bits */

            lf_cb = &(LFCB[idx_LF * M / 2]);    /* ptr init */
            hf_cb = &(HFCB[idx_HF * M / 2]);    /* ptr init */

            L_mse = 0; /* move32() */ /* we accumulate energies on the positive  side   */
            for (j = 0; j < (M / 2); j++)
            {
                err = (scfLF_Q11[j] - lf_cb[j]);  /* a "+"sign,  LF err in Q11 */
                if (signbitLF != 0)
                {  /* negate LF*/
                    err = (scfLF_Q11[j] + lf_cb[j]); /* "-""-"  --> "+"  LF err in Q11 ,single BASOP */
                }
                L_mse += (err * err);  /* simulate L_mac0 , accumulate towards max positive side */

                err = (scfHF_Q11[j] - hf_cb[j]); /* a "+"sign, HF err in Q11 */
                if (signbitHF != 0)
                {   /* negate HF */
                    err = (scfHF_Q11[j] + hf_cb[j]);  /* -- --> "+",  LF err in Q11 ,single BASOP */
                }
                L_mse += (err * err);   /*simulate L_mac0 */ /* now total error */
            }

            L_mse_best = MIN(L_mse, L_mse_best); /* 1 cycle BASOP preupdate best  */
            if ((L_mse - L_mse_best) == 0)
            {
                best_ind = i;      /* update winner, single BASOP */
            }
        }
    }/* segment seg */

    *min_mse = (LC3_FLOAT)((double)L_mse_best) / (2048.0*2048.0);  /* Word32 mse -> Cfloat output calculation  */

    assert(best_ind >= 0 && best_ind < cb_len);

    return (LC3_INT32)best_ind;
}

void snslr_st1C_vector_dec(LC3_INT16 idx, const LC3_INT8* CBW8, LC3_INT16 scaleQ4, LC3_INT16 inv_scaleQ15, LC3_INT32 v_len, LC3_INT32 cb_len, LC3_FLOAT *st1C_vector)
{
    /* decompose the received 0... 169 index , into the correct (integer and)  float st1C vector  */
    /* even in C-float the st1C coeffs  are put into a S16Q11 final integers domain  */
    /* Enables BE compatibility between {BASOP, float, double}  arithmetic implmentations */


    LC3_INT32 i; /*counter*/
    const LC3_INT8 *cb;
    LC3_INT32 L_tmp;
    LC3_INT16  s_tmp, st1C_Q11[M];

    UNUSED(v_len);
    UNUSED(cb_len);
    UNUSED(inv_scaleQ15);  /* req for debugging only */
    assert(idx >= 0 && idx < cb_len);
    assert(v_len == M);

    cb = &(CBW8[idx*M]);                 /* pointer init */
    for (i = 0; i < (M); i++)
    {
        /* BASOP: L_tmp = L_mult0((int16_t)cb[i], scaleQ4 ); */ /*S8Q7 * S15Q4 */ /*sign+7bit, sign+4 bits --> sign+11bit  .lt  sign+23 bits*/
        L_tmp = ((int16_t)cb[i] /*S8Q7*/ * scaleQ4 /* S16Q4 */);  /* S32Q11 */

        /* Cfloat: convert to  Word16 Q11  integer  as  Word16 Q11 SNS domain bits, and then into floats */    
        assert(L_tmp >= -32768L && L_tmp <= 32767L); /* INT16 domain check*/
        s_tmp = (int16_t)(L_tmp);
        st1C_Q11[i] = s_tmp;
    }

    /* Cfloat: convert  the Word16 Q11 SNS domain bits, and then into floats */
    for (i = 0; i < M; i++)
    {
        st1C_vector[i] = ((LC3_FLOAT)st1C_Q11[i])*(1.0 / 2048.0);
    }
}

LC3_INT32 MSEsearchGenericScaledW8(LC3_FLOAT *scf, const LC3_INT8 *sns_CBW8, LC3_INT16 scaleQ4, LC3_INT16 inv_scaleQ15, LC3_INT32 v_len, LC3_INT32 cb_len, LC3_FLOAT* min_mse)
{
    /*   scf float input values are  typically in the range  +12.0  to -12.0.
        ROM table stored in WORD8 [+127,-128], format  corresponding to   ]+1.0  .. -1.0 ]
        inv_scaleQ15, [downscaling  value in Q15]  applied before search
        scaleQ12   upscaling value quantized in Q12,  used in the mse calulation and in the common float and BASOP  synthesis routines

        L_mse evaluated here in a positive integer Word32 domain to match   BASOP

            Fuzzing/saturation considerations
            max M==16 values (4 bits) yields  16*(256*256)=>2^(4+8+8)=2^20  .lt  2^31,
            IntegerWord32  is a safe search domain

           a WC input analysis  would be when half  target  entries are  "-256" and half +255

     */
        LC3_INT32 i, n, best_ind; /*counters*/
    LC3_INT32 L_mse, L_mse_best;
    LC3_INT16  targetW16[M], err;
    LC3_FLOAT  f_tmp, f_scale;
    const LC3_INT8  * cbW8;

    f_scale = ((LC3_FLOAT)inv_scaleQ15) / 32768.0;
    for (i = 0; i < v_len; i++) {
        f_tmp = (scf[i] * f_scale);
        assert(f_tmp >= -4.0 && f_tmp < 4.0);    /* check for about 2  bit  integer margin in the target */
        targetW16[i] = (LC3_INT16)LC3_FLOOR(f_tmp*128.0);    /*cast to INT16   W8Q7   */
    };
    /* in BASOP  1/32768 and  *128 is handled by one single shift) */

    L_mse_best = INT_MAX;  /* largest possible positive number in INT32 */
    best_ind = -1;

    for (i = 0; i < cb_len; i++) {
        L_mse = (0L);   /* move32() */
        /* we accumulate energies on the positive side Word8 vectors do not have any issue with  saturation  (16*256*256) ==2^(4+8+8)  .lt 2^31  MAX_INT */

        cbW8 = &sns_CBW8[i*v_len];  /*ptr init */
        for (n = 0; n < v_len; n++)
        {
            err = (targetW16[n] - ((LC3_INT16)cbW8[n]));     /* cast from Word8 to Word16 is not for free,
                                                                 actual cost in a Word16 architecture is a  L_and(x,0x00ffff)  or  a shr(x,8) */

            L_mse += (err*err);    /* L_mse = L_mac0(L_mse, err, err); */
        }

        if ((L_mse - L_mse_best) <= 0)
        {  /* a value closer to 0  */
            best_ind = i;     /* single BASOP for best idx update */
        }
        L_mse_best = MIN(L_mse, L_mse_best); /* always update best MSE using L_min() in the idx loop, reduces WC WMOPS  */
    }
    assert(best_ind >= 0 && best_ind < cb_len);


    *min_mse = (LC3_FLOAT)((double)(L_mse_best)) / (128.0*128.0);  /* Word32 L_mse in Q7   ->  Cfloat output calculation  */

    f_tmp = (((double)scaleQ4) / 16.0);
    *min_mse *= (f_tmp*f_tmp);              /* make gain scaling a part of  Word32 L_mse  ->  Cfloat output calculation  */

    return best_ind;
}

/*   LRSNS stage 1 functionality */
LC3_INT32 snsQuantScfEncLRSt1ABC(LC3_FLOAT* env, LC3_INT32* L_index, LC3_FLOAT *min_mse_saveBCA_ptr,
    LC3_INT32* ind_saveB_ptr, LC3_FLOAT* st1_vectors, LC3_INT32 pitch_rx, LC3_INT32 ltpf_rx)
{

    LC3_INT32 i;
    LC3_FLOAT *st1_vector, *st1_vectorA, *st1_vectorB, *st1_vectorC, target[M];
    LC3_FLOAT  st1_vectorBC[M];
    LC3_FLOAT  min_mse_saveA, min_mse_saveB, min_mse_saveC, min_mse_saveBC;
    LC3_INT32 ind_saveA, ind_saveC;
    const LC3_FLOAT *cb;
    LC3_INT16  stage1_mode; /*0=A, 1=B, 2=C. -1==fail*/
    LC3_FLOAT min_mse_saveB_fxlike;
    LC3_INT32 ind_saveB_fxlike = -1;
    LC3_FLOAT st1_vectorB_idx[M];
    LC3_FLOAT st1C_lim;
    LC3_FLOAT f_mse_tmp;
    LC3_INT32 ind_saveC_ScaledW8 = -1;
    LC3_FLOAT min_mse_saveC_ScaledW8;

#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
    UNUSED(ltpf_rx);
#endif 

    stage1_mode = -1;         /* output mode */

    st1_vectorA = &(st1_vectors[0]);
    st1_vectorB = &(st1_vectors[1 * M]);
    st1_vectorC = &(st1_vectors[2 * M]);
    st1_vector = &(st1_vectors[3 * M]);

    /* snslr stage1  B(170) and C(170), A(2)  evaluation */
    /*        b0-b8     b9
       segm ,  idx9b , stop bit,  comment use
       -----+--------+---------
        A   | 510,511| n/a,   2 entries,  9 bit total
      ------+--------+--------
      ------+--------+--------+-------
        B   | 0--169 |   1  ,  170 entries,  10 bit total
      ------+--------+--------
        C   | 170-339|   1  , 170 entries,   10 bit total
      ------+--------+--------+------------
      ------+--------+--------+------------
        B*  | 340-509|   1     --> aux=1, 170, 3b+17b for stage2 'LR_full/LR_fix', 30 bit total
      ------+--------+--------+-------
        B*  | 0--169 |   0  ,  --> aux=0, 170,  2b+17b for stage2 'LR_SplitLF' ,  29 bit total
      ------+--------+--------+-------
        B*  | 170-339|   0  ,  --> aux=1, 170,  2b+17b for stage2 'LR_SplitLF',29 bit total
      ------+--------+--------+-------
        B*  | 340-509|   0     --> aux=0, 170, 3b+17b for stage2 'LR_full/LR_fix', 30 bit total
      ------+--------+--------+-------
     */

    {   /* stage 1 section A(2),  a very small 2xM entry  cb  */
        cb = lrsns_st1A_topTab_1bitNoDC;

        memcpy(target, env, sizeof(LC3_FLOAT)*M);

        ind_saveA = MSEsearchGeneric(target, cb, M, 2, &min_mse_saveA);
        memcpy(st1_vectorA, &(cb[ind_saveA*M]), sizeof(LC3_FLOAT)*M);
    }

    min_mse_saveB = LC3_CONST_FLOATMAX; /*safety init*/

    {  /* stage1 section B(170) MSE analysis */
        memcpy(target, env, sizeof(LC3_FLOAT)*M);
        *ind_saveB_ptr = -1;

        ind_saveB_fxlike = MSEsearchCbBIdxMap(target, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11,
            lrsns_st1B_merged170orderSortedSegmCum, lrsns_st1B_merged170orderSort12bitIdx, M, 170, &min_mse_saveB_fxlike); /*st1B LF,HF idx lookup search 170 ,170 Word16s ,  0.34kB ROM */
        snslr_st1B_vector_dec(ind_saveB_fxlike, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, lrsns_st1B_merged170orderSortedSegmCum, lrsns_st1B_merged170orderSort12bitIdx,
            st1_vectorB_idx);

        *ind_saveB_ptr = ind_saveB_fxlike;
        min_mse_saveB = min_mse_saveB_fxlike;
        memcpy(st1_vectorB, st1_vectorB_idx, sizeof(LC3_FLOAT)*M); /*use low ROM version */

        {
            /* remove DC that can remain in the LF,HF index  stored stage cb  B structure  */
            /* a very slight offline decrease in perf 0.001 dB in AvSD when searching with DC above,
               but it allows much better stage1 ROM-reuse performance
            */
            LC3_FLOAT dc = snslr_remove_st1_DC_fQ11(st1_vectorB, M); /* inplace removal of dc in st1_vectorB */
            LC3_FLOAT f_tmp;
            assert(min_mse_saveB >= 0.0);

            f_tmp = min_mse_saveB - M * (dc*dc);
            min_mse_saveB = MAX(f_tmp, 0.0); /* truncation in DC removal can cause negative MSE, limit to 0 */
        }
    } /* end of stage 1 section B , search */



    st1C_lim = 3.97; /*  corresponding  to an  SD limit of  = 1.5  */
    f_mse_tmp = MIN(min_mse_saveA*0.875, min_mse_saveB);
    
    if (f_mse_tmp < st1C_lim)  /*  skip C search if SD is already low enough < 1.5dB)  to save average WMOPS */

    {
        min_mse_saveC = 2 * 5.579999923706055;  /*disable C selection in consecutive  logic */
        ind_saveC = 0;    /* a valid index that will not be used */
    }
    else
    {   /* search stage1 C */
        /* search another set (pitch dependent section C) of 170 mean residual vectors  */
        /* pitch_info rx selects a  mean   and then employ a trained  residual 1x170W8 based harmonic outlier table */

            /* float means are based on  Word16 S16Q11 values,  so that BASOP and float  always may become BE in synthesis*/
        int16_t CBCmeanp_ind = pitch_rx; /* 0 or 1 */

#ifndef  LRSNS_CBC_NO_LTPF_DEPENDENCY 
        if (pitch_rx != 0 && ltpf_rx != 0)
        {
            CBCmeanp_ind = CBCmeanp_ind + 1;  /*LTPF_rx is also active */
        }
#else 
    /* CbC only dependent on LTP transmission on/off */
#endif 

        for (i = 0; i < M; i++)
        {
            target[i] = env[i] - lrsns_st1CTrainedMapMeans[CBCmeanp_ind][i];
        }

        ind_saveC_ScaledW8 = MSEsearchGenericScaledW8(target, lrsns_st1C_Both_Word8,
            lrsns_st1C_Both_scaleQ4_7p4bits_fx[1], lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1],
            M, 170, &(min_mse_saveC_ScaledW8));

        snslr_st1C_vector_dec(ind_saveC_ScaledW8, lrsns_st1C_Both_Word8, lrsns_st1C_Both_scaleQ4_7p4bits_fx[1], lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1],
            M, 170, st1_vectorC);
        /* integer based decoder side synthesis of scaled W8 table for best possible interop */

        ind_saveC = ind_saveC_ScaledW8;
        min_mse_saveC = min_mse_saveC_ScaledW8;

        for (i = 0; i < M; i++) {
            st1_vectorC[i] += lrsns_st1CTrainedMapMeans[CBCmeanp_ind][i]; /*  Q11 means*/
        }
    }

    /* BC stage1    comparison   */
    /* initially assume B as stage 1 winner */
    min_mse_saveBC = min_mse_saveB;
    L_index[0] = *ind_saveB_ptr;  /* 0...169   */

    memcpy(st1_vector, st1_vectorB, sizeof(LC3_FLOAT)*M);   /* stage 1 segmentB result without DC  copied as base for st2 */


    if (min_mse_saveC < min_mse_saveBC)
    {  /* C better than B */
        min_mse_saveBC = min_mse_saveC;
        L_index[0] = 170 + ind_saveC;  /* [2x 170] (9+1)b             [170-339] ,  ">=170"  is a signal to multiplexor  */
        memcpy(st1_vector, st1_vectorC, sizeof(LC3_FLOAT)*M);
    }

    memcpy(st1_vectorBC, st1_vector, sizeof(LC3_FLOAT)*M);
    L_index[1] = -10;
    assert(min_mse_saveBC >= 0.0);

    /* (9(A)<->10(BC) bit weighted comparison   */
    *min_mse_saveBCA_ptr = min_mse_saveBC;
    if (min_mse_saveA*0.875 < min_mse_saveBC)  /* a minor favouring of the 9b vector results  sqrt(0.875)  =>  approx  0.6dB level domain  */
    {
        L_index[0] = 510 + ind_saveA;  /*   only [510, 511]  possible */
        L_index[1] = -9;

        cb = lrsns_st1A_topTab_1bitNoDC;
        memcpy(st1_vectorA, &(cb[ind_saveA*M]), sizeof(LC3_FLOAT)*M);
        memcpy(st1_vector, st1_vectorA, sizeof(LC3_FLOAT)*M);

        *min_mse_saveBCA_ptr = min_mse_saveA * 0.875;;
    }

    /* index0_saveBCA = index[0];*/    /*  0 ... 511 */
    /* index1_saveBCA = index[1];*/    /* -9 or -10 */

    stage1_mode = 2; /* C , assumed */
    if (L_index[0] >= 510)
    {
        stage1_mode = 0; /* A */
    }
    if (L_index[0] <= 169) {
        stage1_mode = 1; /* B */
    }

    return stage1_mode; /* return best  stage1  mode */
}

LC3_INT16 snsQuantScfEncLR(LC3_FLOAT* env, LC3_INT32* L_index, LC3_FLOAT* envq, Dct2 dct2structSNS, LC3_INT32 pitch_rx, LC3_INT32 ltpf_rx)
{
    LC3_INT32 i;
    LC3_FLOAT min_mse_saveBCA;
    LC3_INT32 ind_saveB;
    LC3_INT16 st1_mode;

    LC3_FLOAT stage2_en1_norm_sub[M];
    LC3_FLOAT st1_vectors[(SNSLR_MAX_PVQ_SEARCH_CAND+1)*M], *st1_vectorA, *st1_vectorB, *st1_vectorC, *st1_vector;
    LC3_FLOAT pvq_target_pre[M];
    LC3_FLOAT pvq_target[M];
#ifdef  LRSNS_WMC_FIX
    LC3_INT32 y[SNSLR_MAX_PVQ_SEARCH_CAND][M];  /* o:  [3]*[0...M-1] */
#else
    LC3_INT32 y_tmp[3*M];
#endif 
    LC3_FLOAT stage2_en1_norm_pre_sub[M];
    LC3_FLOAT envq_st1B_st2[M];
    LC3_FLOAT mse_st1B_st2;
    LC3_FLOAT mse_st1B_st2_dct_domain;
    LC3_INT32 gain_idx, shape;
    LC3_FLOAT scfq_gain;
    LC3_INT32 fix_shape_idx;
    LC3_INT16 envelope_bits; /* function output */  
    LC3_INT32 fix_shift_ind, fix_shift_bits, fix_end_sign, LS_tmp_ind;
    LC3_INT32 shape_local;

    UNUSED(fix_shift_bits); /* used for assert */

#ifndef  LRSNS_WMC_FIX
    LC3_INT32(*y)[M];   /* C-construct to allow matrix adressing into a scratch area  */
#endif 
    envelope_bits = -1;   /* output : 9,10, 29/30 */
    gain_idx = 1;         /* stage 2 gain idx range 0 ... 7 , or 0 ... 3 */
    shape = -1;           /* stage 2 shape range    0 ... 5  */

    st1_vectorA = &(st1_vectors[0]);
    st1_vectorB = &(st1_vectors[1 * M]);
    st1_vectorC = &(st1_vectors[2 * M]);
    st1_vector = &(st1_vectors[3 * M]);


#ifndef  LRSNS_WMC_FIX
    y = (LC3_INT32(*)[M]) &(y_tmp[0]); /*  y is an NxM Matrix Ptr */
#endif 
    { /*  stage1 A,B,C   */

        ind_saveB = -1;
        min_mse_saveBCA = M * 32.0*32.0;
        st1_mode = snsQuantScfEncLRSt1ABC(
            env, L_index, &min_mse_saveBCA, &ind_saveB, st1_vectors, pitch_rx, ltpf_rx);

        if (ind_saveB < 0 || ind_saveB > 169)
        {
            assert(ind_saveB >= 0 && ind_saveB <= 169); /* idxB always needed in case stage2 is  activated  */
        }

        if (st1_mode == 0)
        {
            envelope_bits = 9;
            memcpy(st1_vector, st1_vectorA, sizeof(LC3_FLOAT)*M);
            assert(L_index[0] >= 510 && L_index[0] <= 511);
        }

        if (st1_mode == 1)
        {
            envelope_bits = 10;
            memcpy(st1_vector, st1_vectorB, sizeof(LC3_FLOAT)*M);
            assert(L_index[0] >= 0 && L_index[0] <= 169);
        }

        if (st1_mode == 2)
        {
            envelope_bits = 10;
            memcpy(st1_vector, st1_vectorC, sizeof(LC3_FLOAT)*M);
            assert(L_index[0] >= 170 && L_index[0] < 2 * 170);
        }

        L_index[1] = -910; /* aux */
        L_index[2] = -envelope_bits; /* signal shape  */

    }
    /* only run stage 2 when necessary  */

    {  
        LC3_FLOAT  mse_lim_smooth;
        mse_lim_smooth = (5.41);   /*    1.75 SD */

        mse_st1B_st2 = 2.0* min_mse_saveBCA + 1.0;   /*  indicate that st1B+st2 is not used by setting a higher MSE than st1BCA  */

        if (min_mse_saveBCA > mse_lim_smooth)
        {
            /* run  and evaluate STAGE 2,  using vector B as  stage 1 */
            for (i = 0; i < M; i++)
            {
                pvq_target_pre[i] = env[i] - st1_vectorB[i];
            }

            dct2_apply(&dct2structSNS, pvq_target_pre, pvq_target);

            pvq_target[0] = 0.0;  /* DC always zero  */

            fix_shape_idx = -1;
            pvq_fess_enc_search(pvq_target, y, &fix_shape_idx);  /* best shape search  splitLF, full, best_fess_env      */

            assert(y[0][0] == 0 && y[1][0] == 0 && y[2][0] == 0);
            lrsns_quant_shape_gain(pvq_target, y, &gain_idx, &shape, stage2_en1_norm_pre_sub, &scfq_gain, 2 + 1);

            if (shape == 2)
            {
                shape = 2 + fix_shape_idx;
            }

            /*   check if MSE after stage 2 is better already here in dct domain, avoid unnecessary IDCT-II calls  */
            envq_st1B_st2[0] = 0;
            for (i = 1; i < M; i++)
            {
                envq_st1B_st2[i] = (stage2_en1_norm_pre_sub[i] * scfq_gain);
            }

            mse_st1B_st2_dct_domain = calc_mse(pvq_target, envq_st1B_st2, 1.0f, M);

            if (min_mse_saveBCA < mse_st1B_st2_dct_domain)
            {
                /* no need for an IDCT as stage2 was worse than only stage1  */
                mse_st1B_st2 = mse_st1B_st2_dct_domain;
            }
            else
            {
                /* Inverse transform  */
                idct_II(stage2_en1_norm_pre_sub, stage2_en1_norm_sub, M);

                for (i = 0; i < M; i++)
                {
                    envq_st1B_st2[i] = st1_vectorB[i] + (stage2_en1_norm_sub[i] * scfq_gain);
                }
                mse_st1B_st2 = calc_mse(env, envq_st1B_st2, 1.0f, M);
            }
        }
    } /*end of stage2 search */

    /* post-evaluate if one of (st1B, st1C, st1A) was actually better than st1B+stage2 */
    if (mse_st1B_st2 < min_mse_saveBCA)
    {  /*   use stage1B + st2  at  29b/30b bits total */
        L_index[0] = ind_saveB;
        L_index[1] = 2930;         /* later stage2    aux  value  LS_splitLF or LS_full or s0,   put here  as a 0 or 1 */
        L_index[2] = shape;        /* 0=splitLF, 1=full,  ( 2=fixEnv0, 3=fixEnv1, 4: fixEnv2, 5: fixEnv3 )  */
        L_index[3] = gain_idx;     /*  gain_idx  with a shape dependent number of  levels  (4  or 8  levels  ) */

        memcpy(envq, envq_st1B_st2, sizeof(LC3_FLOAT)*M);
        memcpy(st1_vector, st1_vectorB, sizeof(LC3_FLOAT)*M); /* save final  st1B result,  st1 in combination with stage 2, for verification   */
        envelope_bits = 29;  /* 'LR_splitLF' */
        if (shape > 0)
        {
            envelope_bits += 1; /*30 'LR_full/LR_fixenv' */
        }

        {
            /* DBG check values  */
            assert(shape >= 0);
            assert(envelope_bits >= 29);
            assert(L_index[0] <= 170);  /*only B allowed */
            assert(L_index[1] >= 0);

            assert(gain_idx >= 0); /*index*/
            assert(scfq_gain > 0.0); /* value */
        }
    }
    else
    {   /* stick to stage1(best of BCA)   at  9 or 10  bits */
        assert(L_index[1] < 0 && L_index[0] >= 0 && L_index[0] < 512);
        envelope_bits = ((L_index[0] >= 510) ? 9 : 10);
        shape = -envelope_bits; /* signal an invalid shape number to enc-entropy */

        memcpy(envq, st1_vector, M * sizeof(LC3_FLOAT)); /* output */
        memset(stage2_en1_norm_sub, 0, M * sizeof(*stage2_en1_norm_sub));
        scfq_gain = 0.0f;
        gain_idx = -1;   /* L_index sentinel */
        L_index[2] = shape;
    }

    /******************************************************************/
    /*  signal to enc_entropy for LRSNS semi-fractional multiplexing  */
    /******************************************************************/
    /* integer multiplexing   29/30 bit modes into intermediate  unmuxed integer indeces  0...7  */
    /* a bit of fractional multiplexing  for these  indeces  0...7  is done later,  in function enc_entropy()  */
    if (shape >= 0)
    {  /* stage 2 multiplexing manipulations */
    

        if (shape == 0)
        {   /*  splitLF   */
            LC3_INT32  n5k = 0;
            for (i = 0; i < 5; i++)
            {
                n5k += abs(y[shape][1 + i]);
            }

            if (n5k == 6)
            {
                MPVQ_enum(5, &(y[shape][1]), &(L_index[4]), &LS_tmp_ind); /*  P(N=5,K=6) (10)=10 bit L_index    */
                L_index[1] = LS_tmp_ind;     /* set the aux bit for the  splitLF path, plant the first LS as aux */
                assert((L_index[4] >= 0) && (L_index[4] < (SNSLR_NPVQ_L5K6 >> 1)));

                MPVQ_enum(8, &(y[shape][1 + 5]), &(L_index[5]), &LS_tmp_ind);
                L_index[5] = (L_index[5] << 1) + LS_tmp_ind; /* A full PVQ 7 bit index for the  P(N=8,K=2) B config*/
                assert((L_index[5] >= 0) && (L_index[5] < (1 << 7)));
            }
            else
            {
                MPVQ_enum(5, &(y[shape][1]), &(L_index[4]), &LS_tmp_ind);  /*  PVQ(N=5,K=8) (12.x   in total, i.e.  LS+ 11.x )    */
                L_index[1] = LS_tmp_ind;      /* aux bit for the  splitLF path  ,  plant the first LS as aux  */
                assert((L_index[4] >= 0) && (L_index[4] < (SNSLR_NPVQ_L5K8 >> 1)));
                L_index[5] = -8; /* signal LF PVQ(5,8)  and  zeroed HF(10,0) */
            }
        }
        if (shape == 1)
        {   /*  full (15,5) ,  LS kept separated  */
            /* indicate a stage2 path in the 9 bit stage1 index   */
            MPVQ_enum(15, &(y[shape][1]), &(L_index[4]), &LS_tmp_ind);   /*  mPVQ 16.66 bits in index[4], and LS 1 bit in index[1] */
            L_index[1] = LS_tmp_ind;   /*aux bit location,  0 or 1 , we plant the LS there  */
            assert((L_index[4] >= 0) && (L_index[4] < (SNSLR_NPVQ_L15K5 >> 1)));
        }
        if (shape >= 2 && shape <= 5)
        {
            /* fixEnv0, fixEnv1, fixEnv2, fixEnv3 */
            /* send the fixed env subshape mode to enc_entropy  */
            L_index[4] = (shape - 2); /* env shape, 0-->"1" , 1--> "env1"  */   /* L_index[2] has original shape 0...5 */

            shape_local = 2;     /*a single y shape vector for all fixed env */

            fix_shift_ind = 0;
            while (y[shape_local][1 + fix_shift_ind] == 0)
            {
                fix_shift_ind += 1;
            }
            fix_shift_bits = 2;

            assert(fix_shift_ind < (1 << fix_shift_bits));

            L_index[1] = (y[shape_local][1 + fix_shift_ind /* + 0*/] < 0);     /* aux_bit : 0 (or 1)    , will indicate the  s0 sign in the FESS fix shape */


            fix_end_sign = 12;
            if (shape == 5) 
            {
                assert(L_index[4] == 3);
                fix_end_sign = (fix_end_sign - 2);  /* shape 4 has 2 bits shift and a total of 10 signs =2^10*2^2 = 2^12 = 4096 */
            }

            L_index[5] = fix_shift_ind; /* the two shift bits will be pushed up to b11,b12 , for 11 signs s1-s11  */

            for (int sign_ind = 1; sign_ind < fix_end_sign; sign_ind++)  /* push the remaining  sequential signs s1-s11(or s1-s9),  into a single idx */
            {  /* s1 is in the MSB, and s11 is in the lsb*/
                L_index[5] = L_index[5] << 1;
                if (y[shape_local][1 + fix_shift_ind + sign_ind] < 0) /*"1" means negative, "0" means positive */
                {
                    L_index[5] += 1;
                }
            }
            assert(L_index[5] >= 0 && L_index[5] < (1 << (fix_shift_bits + (fix_end_sign - 1))));
        }
    } /* end of stage2 premultiplexing  */

    {
        assert(envelope_bits == 9 || envelope_bits == 10 || envelope_bits == 29 || envelope_bits == 30);
    }

    return envelope_bits;
}

void FESSdeenum(LC3_UINT8 dim_in,    /* i :  dimension of vec_out     */
    LC3_UINT8 n_env,     /* i :  number of envelopes    */
    LC3_UINT8 n_shift,   /* i :  number shifts        */
    LC3_UINT8 n_signs,   /* i :  number signs          */
    LC3_INT32 env_ind,   /*  i:indx */
    LC3_INT32 shift_ind, /*  i:indx */
    LC3_INT32 sign_ind,  /*  i:indx */
    LC3_INT32* vec_out /* o :  FESS  integer  pulse train  */)
{

        LC3_INT32 i;
    LC3_INT32  sign_val;

    assert(n_env >= 1 && n_env <= 4);
    assert(env_ind >= 0 && env_ind < n_env);
    assert(shift_ind >= 0 && shift_ind < n_shift);

    UNUSED(n_env);
    UNUSED(n_shift);
    memset(vec_out, 0, sizeof(*vec_out)*dim_in);

    for (i = (shift_ind + n_signs - 1); i >= shift_ind; i--)
    {
        /* low numbered coeff  signs are in the msb's */
        /* high  numbered coeff  signs are in the lsb's */
        assert(i < dim_in);
        sign_val = 1 - 2 * (sign_ind & 0x01);
        sign_ind = (sign_ind >> 1);   
        vec_out[i] = sign_val * env_ptrs[env_ind][i];    /* vec_out[i] = sign_val * amps[env_ind*(M - 1) + i]; */
    }
}

void snsQuantScfDecLR(LC3_INT32* sns_vq_idx, LC3_FLOAT* scf_q, LC3_INT32 pitch_rx, LC3_INT32 ltpf_rx)
{
    LC3_INT32   i;
    LC3_INT32 mPVQ_ind;         /* can be 16-17 bits */
    LC3_INT16   shape_idx, gain_idx, cb_idx, aux_idx, LS_ind;
    LC3_INT16   env_ind, shift_ind, sign_ind, n_signs;

    LC3_INT32 Y_shape_j[M];
    LC3_FLOAT Xq_shape_j[M], Xq_shape_j_idct[M], sum;
    const LC3_FLOAT *cb;
    const LC3_FLOAT *gainTab;
    LC3_FLOAT st1_scf_q[M];  
    LC3_INT16 CBCmeanp_ind = pitch_rx; /* 0 or 1 */
    const LC3_FLOAT *mean_cb;
    LC3_INT16 sign_mask = 0x07ff;


#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
    UNUSED(ltpf_rx);
#endif 

    sum = 0;
    memset(Y_shape_j, 0, sizeof(LC3_INT32) * M);

    gainTab = &(lrsns_gains_Q11[0][0]);  /* gcc warning  init   */
    gain_idx = 0;                     /* gcc  warning init   */

    cb_idx = sns_vq_idx[0];
    aux_idx = sns_vq_idx[1];
    shape_idx = sns_vq_idx[2];  /* analysis order shape idx  -9,-10,    0,1,  2,3,4,5 */
    gain_idx = sns_vq_idx[3];   /* stage 2 gain */

   /* Stage1 cand   */
    if (shape_idx == -9)
    {
        /* minminal   2*16 SNS codebook, no DC */
        cb = lrsns_st1A_topTab_1bitNoDC;
        memcpy(scf_q, &(cb[cb_idx * M]), sizeof(LC3_FLOAT) * M);
    }
    else if (shape_idx == -10)
    {    /* 0..339 */ /* stage 1 only, transmitted  in 9+1= 10 bits */
        if (cb_idx < 170)
        {  /*Stage 1B   */
            snslr_st1B_vector_dec(cb_idx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, lrsns_st1B_merged170orderSortedSegmCum, lrsns_st1B_merged170orderSort12bitIdx, scf_q);
            snslr_remove_st1_DC_fQ11(scf_q, M);
        }
        else
        {
            cb_idx -= 170;
            /*  Stage 1C  harm outlier  CB  with a  pitch dependent mean   */
            /*  Q11 values , so that BASOP and float  always becomes BE in synthesis*/
            assert(cb_idx >= 0 && cb_idx < 170);
            snslr_st1C_vector_dec(cb_idx, lrsns_st1C_Both_Word8, lrsns_st1C_Both_scaleQ4_7p4bits_fx[1], lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1],
                M, 170, scf_q
            );


            /* cbC add harmonic mean , based on pitch_info availability */
            pitch_rx = sns_vq_idx[3];    /* LTP active flag directly from dec_entropy */
#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
            ltpf_rx = 0;                               /* CB_C has no dependency on LTPF active flag */
#else 
            ltpf_rx = sns_vq_idx[4];   /* LTPF active flag LTP active flag directly from dec_entropy */
#endif 
            CBCmeanp_ind = pitch_rx; /* 0 or 1 */
#ifndef  LRSNS_CBC_NO_LTPF_DEPENDENCY
            if (pitch_rx != 0 && ltpf_rx != 0)
            {
                CBCmeanp_ind = CBCmeanp_ind + 1;  /* high corr ltpf_rx is also active */
            }
#endif

            mean_cb = lrsns_st1CTrainedMapMeans[CBCmeanp_ind]; /* point to pitch dependent mean */
            for (i = 0; i < M; i++)
            {
                scf_q[i] += mean_cb[i];
            }
            /* remove_DC()  call is not required for section  C  */
            /* a very small DC can still exist though,  due to Q7+Q4 quantization of values */
        }
    }
    else
    {     /* 0...169 */   /* st1B*   used with a stage 2  shape submode  */
        assert(shape_idx >= 0);
        snslr_st1B_vector_dec(cb_idx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, lrsns_st1B_merged170orderSortedSegmCum, lrsns_st1B_merged170orderSort12bitIdx, scf_q);
        snslr_remove_st1_DC_fQ11(scf_q, M);      /* DC needs removal for st1 part B */
    }
    memcpy(st1_scf_q, scf_q, sizeof(LC3_FLOAT) * M); /* keep track of stage 1 contribution */

    if (shape_idx >= 0) /* stage 2 shapes 0,1,  2,3,4,5   ( negative idx used for stage1 only   */
    {
        /* Stage 2 SNS VQ decoding */
        /* Decode shape_j */
        Y_shape_j[0] = 0;   /* no DCT-II DC-coeff decoded */

        switch (shape_idx)
        {
        case 0:  /* splitLF 29 bits total  */
            LS_ind = aux_idx;
            mPVQ_ind = sns_vq_idx[4];  /* mPVQ(5,6) or mPVQ(5,8) */

            if (sns_vq_idx[5] >= 0)
            {
                MPVQdeenum(5, 6, LS_ind, mPVQ_ind, &Y_shape_j[1]);
                LS_ind = sns_vq_idx[5] & 0x1;
                mPVQ_ind = sns_vq_idx[5] >> 1;
                MPVQdeenum(8, 2, LS_ind, mPVQ_ind, &Y_shape_j[1 + 5]);
            }
            else
            {
                MPVQdeenum(5, 8, LS_ind, mPVQ_ind, &Y_shape_j[1]);
            }
            gainTab = &(lrsns_gains_Q11[0][0]);/* 4 levels in 2 bits  */

            break;
        case 1: /* full  30 bits total */
            LS_ind = aux_idx;
            mPVQ_ind = sns_vq_idx[4];
            MPVQdeenum(15, 5, LS_ind, mPVQ_ind, &Y_shape_j[1]);
            gainTab = &(lrsns_gains_Q11[1][0]); /* 8 levels in 3 bits  */

            break;
        case 2: /* fix env 0  ,  init_bell      12 signs  */
        case 3: /* fix env 1  ,  decay 24-->13  12 signs */
        case 4: /* fix env 2  ,  start bell      12 signs */
        case 5: /* fix env 3  ,  early bell     10 signs  */
            LS_ind = aux_idx; /* s0 */
            env_ind = sns_vq_idx[4];
            assert(env_ind == (shape_idx - 2));

            n_signs = 12; /* including s0 */
            if (env_ind == 3) {
                n_signs -= 2;
            }
            sign_mask = (sign_mask >> (12 - n_signs));

            shift_ind = sns_vq_idx[5] >> (n_signs - 1);
            sign_ind = sns_vq_idx[5] & sign_mask;

            /* put s0 , right next to s1 , to make the sign decoding loop easier */
            sign_ind = (sign_ind)+(LS_ind << (n_signs - 1));   /* s0 put as MSB at 12th position  2^11 , lsb at 2^0 */

            /*FixEnvShiftSigns deenumeration */
            FESSdeenum(15, 4, 4, n_signs, env_ind, shift_ind, sign_ind, &Y_shape_j[1]);    /*30b    ,   4xenv,4xshifts, 10 or12 signs, over 15 positions,*/
            gainTab = &(lrsns_gains_Q11[2][0]);; /* 8 levels in 3 bits  */
            /* fix_envshift_nb = env_ind * 4 + shift_ind; */          /* index for fast normalization lookup */
            break;
        default:

            break;
        }

        /* Unit energy normalization of the received shape */
        for (i = 0; i < M; i++)
        {
            sum += (Y_shape_j[i] * Y_shape_j[i]);
        }

        sum = 1.0 / LC3_SQRT(sum);   /* all shapes will have tabled inv_sqrt() divisions as factors in BASOP */

        for (i = 0; i < M; i++)
        {
            Xq_shape_j[i] = Y_shape_j[i] * sum;
        }

        /* Reconstruction of the quantized SNS scale factors */
        idct_II(Xq_shape_j, Xq_shape_j_idct, M);
        for (i = 0; i < M; i++) {
            scf_q[i] += Xq_shape_j_idct[i] * gainTab[gain_idx];
        }
    }
    else
    {  /* -9, -10 */
      /* LRSNS stage 1 variations  only */
        memcpy(scf_q, st1_scf_q, sizeof(LC3_FLOAT) * M);
    }
}
 
/* LRSNS integer precision based function needed in both encoder and decoder */
LC3_FLOAT snslr_remove_st1_DC_fQ11(LC3_FLOAT *scfq, LC3_INT32  len)
{
    LC3_INT32  i; /*Counter*/
    LC3_INT32  L_dcQ11;
    LC3_FLOAT  f_dcQ11 = 0.0;
    L_dcQ11 = 0L;
    
    for (i = 0; i < len; i++) {
        L_dcQ11 += (LC3_INT32)(scfq[i] * 2048.0f);     /* BE simulation of DC Q11 summation of truncated values in BASOP, preferably  BE in synthesis in FLP/BASOP decoder  */
    }

    assert(len == M);
    {
        L_dcQ11 = L_dcQ11 >> 4;         /* make the average in integer domain ,  no rounding applied before shift,   on purpose */
        f_dcQ11 = ((LC3_FLOAT)L_dcQ11) *(1.0f / 2048.0f);  /*  now a Q11 value to match the overall generic Q11 BASOP scaling of stage1 variables  */
    }

    for (i = 0; i < len; i++)
    {
        scfq[i] -= f_dcQ11;  /* result update */
    }
    return f_dcQ11;  /* output used for encoder side  mse update*/
}

#  endif
