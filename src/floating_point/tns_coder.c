/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static void xcorr(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT lag, LC3_INT inLen);
static void levdown(LC3_FLOAT* anxt, LC3_FLOAT* out_a, LC3_INT* len);
static void poly2rc(LC3_FLOAT* a, LC3_FLOAT* out, LC3_INT len);

void xcorr(LC3_FLOAT* in, LC3_FLOAT* out, LC3_INT lag, LC3_INT inLen)
{
	LC3_INT32   m;


    for (m = 0; m <= lag; m++) {
        /* Calculate correlation */
        out[m] = mac_loop(in, &in[m], (inLen - m));
    }
}

void levinsonDurbin(LC3_FLOAT* r, LC3_FLOAT* out_lev, LC3_FLOAT* rc_unq, LC3_FLOAT* error, LC3_INT len)
{
    LC3_INT   t, i, j;
    LC3_FLOAT g, v, sum, buf_tmp[10];

    g          = r[1] / r[0];
    out_lev[0] = g;

    v         = (1.0 - g * g) * r[0];
    rc_unq[0] = -g;

    for (t = 1; t < len; t++) {

        sum = 0;
        for (i = 1; i <= t; i++) {
            sum += out_lev[i - 1] * r[i];
        }

        g = (r[t + 1] - sum) / v;

        j = 1;
        for (i = t - 1; i >= 0; i--) {
            buf_tmp[j] = out_lev[j - 1] - g * out_lev[i];
            j++;
        }

        move_float(&out_lev[1], &buf_tmp[1], len);

        out_lev[0] = g;

        v         = v * (1 - g * g);
        rc_unq[t] = -g;
    }

    /* Reorder out_lev */
    out_lev[0] = 1;
    j          = 1;
    for (i = len - 1; i >= 0; i--) {
        buf_tmp[j] = -out_lev[i];
        j++;
    }

    move_float(&out_lev[1], &buf_tmp[1], (len - 1));

    out_lev[len] = rc_unq[len - 1];

    *error = v;
}

void levdown(LC3_FLOAT* anxt, LC3_FLOAT* out_a, LC3_INT* len)
{
    LC3_INT32   i, j;
    LC3_FLOAT tmp_buf[8], knxt;
    LC3_FLOAT norm;

	memset(tmp_buf, 0, 8 * sizeof(LC3_FLOAT));
    /* Initial length = 9 */

    /* Drop the leading 1 */

    *len = *len - 1; /* Lenght = 8 */

    /* Last coefficient */
    knxt = anxt[*len]; /* At [7] */

    *len = *len - 1; /* Lenght = 7 */

    j = 0;
    for (i = *len - 1; i >= 0; i--) {
        tmp_buf[j] = knxt * anxt[i + 1];
        j++;
    }
    
    norm = 1.0 / (1.0 - (LC3_FABS(knxt)) * (LC3_FABS(knxt)));

    out_a[0] = 1;
    for (i = 0; i < *len; i++) {
        out_a[i + 1] = (anxt[i + 1] - tmp_buf[i]) * norm;
    }

    *len = *len + 1; /* Length = 8 */
}

void poly2rc(LC3_FLOAT* a, LC3_FLOAT* out, LC3_INT len)
{
    LC3_INT   k, i, len_old;
    LC3_FLOAT buf[9];

    len_old = len;

    zero_float(out, len - 1);

    /* Length = 9 */

    /* Normalize */
    for (i = 0; i < len; i++) {
        a[i] = a[i] / a[0];
    }

    out[len - 1] = a[len - 1];

    /* Process */
    for (k = len - 2; k >= 0; k--) {
        levdown(a, buf, &len);
        out[k] = buf[len - 1]; /* Store last value */

        move_float(a, buf, len);
    }

    /* Shift output array by one to the left to lose leading 1 */
    for (i = 0; i < len_old - 1; i++) {
        out[i] = out[i + 1];
    }
}


void processTnsCoder_fl(LC3_FLOAT* x, LC3_INT bw_cutoff_idx, LC3_INT bw_fcbin, LC3_INT fs, LC3_INT N, LC3_INT frame_dms, LC3_INT nBits,
                        LC3_INT* order_out, LC3_INT* rc_idx, LC3_INT* tns_numfilters, LC3_INT* bits_out
                        , LC3_INT16 near_nyquist_flag
)
{
    LC3_INT i, stopfreq[2], startfreq[2], f, numfilters, maxOrder, bits, sub,
        subdiv_startfreq, subdiv_stopfreq, j, rc_idx_tmp[MAXLAG], order_tmp, tmp, tns;
    LC3_FLOAT minPGfac, minPredictionGain, maxPG, xcorr_out[MAXLAG + 1], sum,
          subdiv_len, nSubdivisions, r[MAXLAG + 1], rc_unq[MAXLAG + 1], error_lev, predGain,
          alpha, rc[MAXLAG], st[MAXLAG + 1] = {0}, s, tmpSave, tmp_fl;
    const LC3_INT* order;
    LC3_FLOAT inv_sum, x_val;
    LC3_FLOAT alpha_loc;
    LC3_INT32 iIndex;

    /* Init */

    if (fs >= 32000 && frame_dms >= 50) {
        numfilters = 2;
    } else {
        numfilters = 1;
    }

    /* 40 * frame_dms / 10 = 4 * frame_dms */
    if (N > 4 * frame_dms)
    {
        N = 4 * frame_dms;
        fs = 40000;
    }

    if (numfilters == 1) {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        stopfreq[0]  = N;
    } else {
        startfreq[0] = floor(600 * N * 2 / fs) + 1;
        startfreq[1] = N / 2 + 1;
        stopfreq[0]  = N / 2;
        stopfreq[1]  = N;
    }
    
    switch (frame_dms)
    {
        case 25:
            maxOrder      = 4;
            nSubdivisions = 2.0;
            break;
        case 50:
            maxOrder      = 4;
            nSubdivisions = 2.0;
            break;
        case 75:
            maxOrder      = 8;
            nSubdivisions = 3;
            break;
        case 100:
            maxOrder      = 8;
            nSubdivisions = 3.0;
            break;
    }

    minPredictionGain = 1.5;

    if (nBits >= 4.8 * frame_dms) {
        order = order1_tns;
    } else {
        order = order2_tns;
    }
    
    /* Processing */
    if (bw_cutoff_idx >= 3 && numfilters == 2) {
        numfilters   = 2;
        startfreq[1] = bw_fcbin / 2 + 1;
        stopfreq[0]  = bw_fcbin / 2;
        stopfreq[1]  = bw_fcbin;
    } else {
        numfilters  = 1;    
        stopfreq[0] = bw_fcbin;
    }

    bits = 0;

    for (f = 0; f < numfilters; f++) {
        subdiv_len = ((LC3_FLOAT)stopfreq[f] + 1.0 - (LC3_FLOAT)startfreq[f]) / nSubdivisions;

        zero_float(r, MAXLAG+1);

        for (sub = 1; sub <= nSubdivisions; sub++) {
            subdiv_startfreq = floor(subdiv_len * (sub - 1)) + startfreq[f] - 1;
            subdiv_stopfreq  = floor(subdiv_len * sub) + startfreq[f] - 1;
            
            if (fs == 32000 && frame_dms == 75)
            {
                if (subdiv_startfreq == 83)
                {
                    subdiv_startfreq = 82;
                }
                
                if (subdiv_stopfreq == 83)
                {
                    subdiv_stopfreq = 82;
                }
                
                if (subdiv_startfreq == 160)
                {
                    subdiv_startfreq = 159;
                }
                
                if (subdiv_stopfreq == 160)
                {
                    subdiv_stopfreq = 159;
                }
            }

            sum = 0;
            for (i = subdiv_startfreq; i < subdiv_stopfreq; i++) {
                sum += x[i] * x[i];
            }

            if (sum < LC3_EPS) {
                zero_float(r, MAXLAG+1);
                r[0] = 1;
                break;
            }

            xcorr(&x[subdiv_startfreq], xcorr_out, maxOrder, subdiv_stopfreq - subdiv_startfreq);

            inv_sum = 1.0 / sum;
            for (i = 0; i <= maxOrder; i++) {
                r[i] = r[i] + xcorr_out[i] * inv_sum;
            }
        }

        for (i = 0; i <= maxOrder; i++) {
            r[i] = r[i] * lagw_tns[i];
        }

        levinsonDurbin(r, xcorr_out, rc_unq, &error_lev, maxOrder);

        predGain = r[0] / error_lev;

        if (predGain > minPredictionGain && near_nyquist_flag == 0) {
            tns = 1;
        } else {
            tns = 0;
        }

        bits++;

        if (tns == 1) {
            minPGfac = 0.85;
            maxPG    = 2;
            if (nBits >= 4.8 * frame_dms) {
                maxPG = minPredictionGain;
            }

            /* LPC weighting */
            if (predGain < maxPG) {
                alpha = (maxPG - predGain) * (minPGfac - 1.0) / (maxPG - minPredictionGain) + 1.0;

                alpha_loc = 1;
                for (i = 0; i <= maxOrder; i++) {
                    xcorr_out[i] = xcorr_out[i] * alpha_loc;
                    alpha_loc *= alpha;
                }

                poly2rc(xcorr_out, rc_unq, maxOrder + 1);
            }

            /* PARCOR Quantization */
            for (i = 0; i < maxOrder; i++)
            {
                iIndex = 1;
                x_val  = rc_unq[i];

                while ((iIndex < 17) && (x_val > quants_thr_tns[iIndex - 1]))
                {
                    iIndex = (iIndex + 1);
                }
                rc_idx_tmp[i] = (iIndex - 2);
            }
            
            /* Filter Order */
            order_tmp = 0;
            for (i = 0; i < maxOrder; i++) {
                rc[i] = quants_pts_tns[rc_idx_tmp[i]];

                if (rc[i] != 0) {
                    order_tmp = i + 1;
                }
            }

            order_out[f] = order_tmp;

            // Disable TNS if order is 0:
            if (order_out[f] == 0) {
                tns = 0;

                // Jump to else statement
                goto tns_disabled;
            }
            tmp = order[order_out[f] - 1];

            /* Huffman Coding of PARCOR coefficients */
            for (i = 0; i <= order_out[f] - 1; i++) {
                tmp += huff_bits_tns[i][rc_idx_tmp[i]];
            }

            bits = bits + ceil((LC3_FLOAT)tmp / 2048.0);

            j = 0;
            for (i = f * 8; i <= f * 8 + order_out[f] - 1; i++) {
                rc_idx[i] = rc_idx_tmp[j];
                j++;
            }
        } else {
tns_disabled:
            order_out[f] = 0;
        }

        /* Filtering */
        if (tns == 1) {
            for (i = startfreq[f]; i <= stopfreq[f]; i++) {
                s       = x[i - 1];
                tmpSave = s;

                for (j = 0; j < order_out[f] - 1; j++) {
                    tmp_fl = rc[j] * s + st[j];
                    s += rc[j] * st[j];

                    st[j]   = tmpSave;
                    tmpSave = tmp_fl;
                }

                s += rc[order_out[f] - 1] * st[order_out[f] - 1];

                st[order_out[f] - 1] = tmpSave;
                x[i - 1]             = s;
            }
        }
    }

    *tns_numfilters = numfilters;
    *bits_out       = bits;
}
