/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               


#include <assert.h>
#include <string.h> /* for mmove */
#include <stdio.h>
#include <stdlib.h>
#include "iisfft.h"
#include "cfft.h"

/*  the fixed length fft functions have been split into sevelral headers to
    have smaller files. to give the compiler more room to optimize the ffts
    can't be in separate compilation units. the header approach seemed to be
    the best compromise. to prevent them being included from anywhere else,
    they are guarded by the INCLUDED_FROM_IISFFT_C macro.
*/
#define INCLUDED_FROM_IISFFT_C
#include "fft_2_9.h"
#include "fft_15_16.h"
#include "fft_32.h"
#include "fft_60_128.h"
#include "fft_240_480.h"
#include "fft_384_768.h"
#include "fft_generic.h"


void LC3_iisfft_apply(Iisfft* handle, LC3_FLOAT* x)
{
    if (handle->sign == -1) {
        if (!fft_n(x, handle->length))
            pfaDFT(x, handle->length, handle->scratch, handle->num_factors, handle->factors, handle->scratch2,
                   handle->isPrime);
    } else {
        assert(0);
    }
}

/* returns 1 if there is no specialized function for length or 1 if a scratch needs to be allocated.
   check the fft_n function */
static LC3_INT need_scratch(LC3_INT n)
{
    return n != 2 && n != 3 && n != 4 && n != 5 && n != 7 && n != 8 && n != 9 && n != 15 && n != 16 && n != 32 &&
           n != 60 && n != 64 && n != 128 && n != 240 && n != 256 && n != 384 && n != 480 && n != 512 && n != 768 &&
           n != 1024;
}

IIS_FFT_ERROR LC3_iisfft_plan(Iisfft* handle, LC3_INT length, LC3_INT sign)
{
    memset(handle, 0, sizeof(Iisfft));
    if (length < 2)
        return IIS_FFT_LENGTH_ERROR;
    handle->length = length;
    handle->sign   = sign;
    if (need_scratch(length)) {
        /* only needed for prime numbers bigger than BORDER_FOR_SECOND_SCRATCH */
        LC3_INT i                    = 0;
        LC3_INT lengthOfPrimeScratch = BORDER_FOR_SECOND_SCRATCH;
        if (!factorize(length, &handle->num_factors, handle->factors, handle->isPrime))
            return IIS_FFT_LENGTH_ERROR;
        handle->scratch = (LC3_FLOAT*)malloc(sizeof(LC3_FLOAT) * 2 * length);
        /* create additional scratch for primeFFT() */
        for (i = 0; i < handle->num_factors; i++) {
            if (handle->isPrime[i] == 1 && handle->factors[i] > lengthOfPrimeScratch) {
                lengthOfPrimeScratch = handle->factors[i];
            }
        }
        if (lengthOfPrimeScratch > BORDER_FOR_SECOND_SCRATCH) {
            handle->scratch2 = (LC3_INT*)malloc(sizeof(LC3_INT) * lengthOfPrimeScratch);
            if (!handle->scratch2)
                return IIS_FFT_MEMORY_ERROR;
        }
        if (!handle->scratch)
            return IIS_FFT_MEMORY_ERROR;
    }

    return IIS_FFT_NO_ERROR;
}

void LC3_iisfft_free(Iisfft* handle)
{
    handle->length = 0;
    if (handle->scratch)
        free(handle->scratch);
    if (handle->scratch2)
        free(handle->scratch2);
}

