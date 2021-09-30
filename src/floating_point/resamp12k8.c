/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void process_resamp12k8_fl(LC3_FLOAT x[], LC3_INT x_len, LC3_FLOAT mem_in[], LC3_INT mem_in_len, LC3_FLOAT mem_50[], LC3_FLOAT mem_out[],
                           LC3_INT mem_out_len, LC3_FLOAT y[], LC3_INT* y_len, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT fs)
{
    

    LC3_INT   len_12k8 = 0, N12k8 = 0, i = 0, k = 0;
    LC3_FLOAT mac = 0, buf_out[120 + MAX_LEN] = {0}, bufdown[128] = {0}, buf[120 + MAX_LEN] = {0};

    const LC3_FLOAT *filter;
    const LC3_FLOAT *filt_input, *filt_coeff;


    switch (frame_dms)
    {
        case 25:
            len_12k8 = LEN_12K8 / 4;
            break;
        case 50:
            len_12k8 = LEN_12K8 / 2;
            break;
        case 100:
            len_12k8 = LEN_12K8;
            break;
    }

    *y_len = len_12k8;
    N12k8  = x_len * 12800 / fs;

    /* Init Input Buffer */
    memmove(buf, mem_in, mem_in_len * sizeof(LC3_FLOAT));
    memmove(&buf[mem_in_len], x, x_len * sizeof(LC3_FLOAT));
    memmove(mem_in, &buf[x_len], mem_in_len * sizeof(LC3_FLOAT));



    filter = lp_filter[fs_idx];

    /* Upsampling & Low-pass Filtering & Downsampling */

    LC3_INT32 index_int  = 1;
    LC3_INT32 index_frac = 0;
    LC3_INT32 resamp_upfac    = resamp_params[fs_idx][0];
    LC3_INT32 resamp_delay    = resamp_params[fs_idx][1];
    LC3_INT32 resamp_off_int  = resamp_params[fs_idx][2];
    LC3_INT32 resamp_off_frac = resamp_params[fs_idx][3];

    k = 0;
    for (i = 0; i < N12k8; i++) {

        filt_input = &buf[index_int];
        filt_coeff = &filter[index_frac * resamp_delay * 2];

        mac = mac_loop(filt_input, filt_coeff, (2 * resamp_delay));
        
        bufdown[k++] = mac;
        
        index_int  = index_int + resamp_off_int;
        index_frac = index_frac + resamp_off_frac;
        
        if ((resamp_upfac - index_frac) <= 0)
        {
            index_int  = index_int + 1;
            index_frac = index_frac - resamp_upfac;
        }
    }


    /* 50Hz High-Pass */
    LC3_FLOAT u_11 = mem_50[0];
    LC3_FLOAT u_21 = mem_50[1];
    LC3_FLOAT u_1, u_2;

    for (i = 0; i < len_12k8; i++) {
        LC3_FLOAT y1  = (highpass50_filt_b[0] * bufdown[i] + u_11);
        u_1        = (highpass50_filt_b[1] * bufdown[i] + u_21) - highpass50_filt_a[1] * y1;
        u_2        = highpass50_filt_b[2] * bufdown[i] - highpass50_filt_a[2] * y1;
        u_11       = u_1;
        u_21       = u_2;
        bufdown[i] = (LC3_FLOAT)y1;
    }

    mem_50[0] = (LC3_FLOAT)u_11;
    mem_50[1] = (LC3_FLOAT)u_21;

    /* Output Buffer */
    memmove(buf_out, mem_out, mem_out_len * sizeof(LC3_FLOAT));

    memmove(&buf_out[mem_out_len], bufdown, len_12k8 * sizeof(LC3_FLOAT));

    memmove(y, buf_out, (*y_len + 1) * sizeof(LC3_FLOAT));

    memmove(mem_out, &buf_out[N12k8], mem_out_len * sizeof(LC3_FLOAT));
}
