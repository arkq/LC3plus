/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void dct2_init(Dct2* dct, int length)
{
    assert(length <= MAX_LEN);
    dct->length = length;
    fft_init(&dct->fft, length);
}

void dct2_free(Dct2* dct)
{
    if (dct) {
        fft_free(&dct->fft);
        memset(dct, 0, sizeof(*dct));
    }
}

void dct2_apply(Dct2* dct, const LC3_FLOAT* input, LC3_FLOAT* output)
{
    Complex   tmp1[MAX_LEN];
    Complex   tmp2[MAX_LEN];
    int       i;
    assert(input != output);

    for (i = 0; i < 8; i++) {
        tmp1[i]           = cmplx(input[i * 2], 0);
        tmp1[16 - i - 1]  = cmplx(input[i * 2 + 1], 0);
    }

    fft_apply(&dct->fft, tmp1, tmp2);

    for (i = 0; i < 16; i++) {
        output[i] = cmul(tmp2[i], dct2_16[i]).r;
    }
    output[0] /= (LC3_FLOAT)1.414213562373095; /* SQRT(2) */
}


void dct4_init(Dct4* dct, int length)
{
    int i;
    assert(length <= MAX_LEN);
    dct->length = length;
    dct->twid1  = calloc(sizeof(*dct->twid1), length / 2);
    dct->twid2  = calloc(sizeof(*dct->twid2), length / 2);
    for (i = 0; i < length / 2; i++) {
        dct->twid1[i] = cexpi(-(LC3_FLOAT)M_PI_LC3PLUS * (i + (LC3_FLOAT)0.25) / length);
        dct->twid2[i] = cexpi(-(LC3_FLOAT)M_PI_LC3PLUS * i / length);
    }
    fft_init(&dct->fft, length / 2);
}

void dct4_free(Dct4* dct)
{
    if (dct) {
        free(dct->twid1);
        free(dct->twid2);
        fft_free(&dct->fft);
        memset(dct, 0, sizeof(*dct));
    }
}

void dct4_apply(Dct4* dct, const LC3_FLOAT* input, LC3_FLOAT* output)
{
    Complex     tmp2[MAX_LEN / 2];
    int         i    = 0;
    Complex*    tmp1 = (Complex*)output;
    const int   len  = dct->length;
    const LC3_FLOAT norm = (LC3_FLOAT)1.0 / LC3_SQRT((LC3_FLOAT)(len >> 1));
    assert(input != output);

    for (i = 0; i < len >> 1; i++) {
        tmp1[i] = cmul(cmplx(input[i * 2], input[len - i * 2 - 1]), dct->twid1[i]);
    }

    fft_apply(&dct->fft, tmp1, tmp2);

    for (i = 0; i < len >> 1; i++) {
        Complex t               = cmul(tmp2[i], dct->twid2[i]);
        output[i * 2]           = t.r * norm;
        output[len - i * 2 - 1] = -t.i * norm;
    }
}
