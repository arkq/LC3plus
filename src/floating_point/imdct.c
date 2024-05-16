/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

/* Function expects already flipped window */
void ProcessingIMDCT_fl(LC3_FLOAT* y, LC3_INT yLen, const LC3_FLOAT* win, LC3_INT winLen, LC3_INT last_zeros, LC3_FLOAT* mem, LC3_FLOAT* x, Dct4* dct)
{
    LC3_FLOAT x_tda[MAX_LEN], x_ov[2 * MAX_LEN];
    LC3_INT   i, j;

    /* Flip imdct window up to down */
    i = winLen - 1;
    j = 0;

    dct4_apply(dct, y, x_tda);

    move_float(x_ov, &x_tda[yLen / 2], yLen / 2);

    j = yLen / 2;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[yLen - 1 - i];
        j++;
    }

    j = yLen;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[yLen / 2 - 1 - i];
        j++;
    }

    j = yLen + yLen / 2;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[i];
        j++;
    }

    for (i = 0; i < winLen; i++) {
        x_ov[i] = x_ov[i] * win[winLen - 1 - i];
    }

    /* Buffer update */
    j = 0;
    for (i = last_zeros; i < yLen; i++) {
        x_ov[i] = x_ov[i] + mem[j];
        j++;
    }

    move_float(&x[0], &x_ov[last_zeros], yLen);

    move_float(&mem[0], &x_ov[yLen + last_zeros], (winLen - (yLen + last_zeros)));
}

void ProcessingITDA_WIN_OLA_fl(LC3_FLOAT* x_tda, LC3_INT32 yLen, const LC3_FLOAT* win, LC3_INT32 winLen, LC3_INT32 last_zeros, LC3_FLOAT* mem, LC3_FLOAT* x)
{
    LC3_FLOAT x_ov[2 * MAX_LEN];
    LC3_INT32 i, j;

    move_float(x_ov, &x_tda[yLen / 2], yLen / 2);

    j = yLen / 2;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[yLen - 1 - i];
        j++;
    }

    j = yLen;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[yLen / 2 - 1 - i];
        j++;
    }

    j = yLen + yLen / 2;
    for (i = 0; i < yLen / 2; i++) {
        x_ov[j] = -x_tda[i];
        j++;
    }

    for (i = 0; i < winLen; i++) {
        x_ov[i] = x_ov[i] * win[winLen - 1 - i];
    }

    /* Buffer update */
    j = 0;

    for (i = last_zeros; i < yLen; i++) {
        x[j] = x_ov[i] + mem[j];
        j++;
    }

    move_float(&x[j], &x_ov[last_zeros+j], yLen-j);

    move_float(&mem[0], &x_ov[yLen + last_zeros], (winLen - (yLen + last_zeros)));
}
