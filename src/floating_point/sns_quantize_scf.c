/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void pvq_dec(LC3_INT k, LC3_INT m, LC3_INT LS_ind, LC3_INT MPVQ_ind, LC3_INT* pulses);
static LC3_INT  find_last_indice_le(LC3_INT compare, const LC3_INT* array, LC3_INT len);
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

static LC3_INT pvq_pulse_search(LC3_FLOAT *xabs, LC3_FLOAT *ener, LC3_FLOAT *corr, LC3_INT *y, LC3_INT start, LC3_INT end)
{
    LC3_INT i;
    LC3_INT nBest;
    LC3_FLOAT bestCorrSq, bestEn;
    LC3_FLOAT corrSq, currCorr, currEn;

    nBest = 0;
    bestCorrSq = 0.0;
    bestEn = 0.0;

    *ener += 1; // Added once for the entire loop

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

static inline LC3_FLOAT calc_mse(LC3_FLOAT *t2rot, LC3_FLOAT *y, LC3_FLOAT gain, LC3_INT N)
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

LC3_INT find_last_indice_le(LC3_INT compare, const LC3_INT* array, LC3_INT len)
{
    LC3_INT idx = 0, i = 0;

    for (i = 0; i < len; i++) {
        if (compare >= array[i]) {
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
