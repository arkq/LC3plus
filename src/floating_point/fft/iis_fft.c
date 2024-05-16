/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "iis_fft.h"

/**************************************************************************************************/

/* AFFT uses two fft implementations
 * cfft is used for lengths of power of two, >= 256.
 * iisfft is used for everything else. it is optimized for certain lengths. for a list of
   fast lengths, check the fft_n function.
*/

#define FFT_COMPLEX 1
#define FFT_REAL 2


static IIS_FFT_ERROR create(HANDLE_IIS_FFT* handle, LC3_INT type, LC3_INT len, IIS_FFT_DIR sign)
{
    IIS_FFT_ERROR  err = IIS_FFT_MEMORY_ERROR;

    /* for real transforms the actual performed fft is half length */
    LC3_INT trlen = (type == FFT_COMPLEX) ? len : len / 2;

    /* check argument sanity */
    if ((sign != IIS_FFT_FWD) && (sign != IIS_FFT_BWD))
    {
        return IIS_FFT_INTERNAL_ERROR;
    }


    if (!(*handle))
    {
      (*handle) = (HANDLE_IIS_FFT)calloc(1, sizeof(IIS_FFT));
    }
    if (!(*handle))
    {
        return IIS_FFT_MEMORY_ERROR;
    }

    (*handle)->len  = len;
    (*handle)->sign = sign;

    /* create sine lookup table for real ffts */
    if (type == FFT_REAL)
    {
        LC3_create_sine_table(len, (*handle)->sine_table);
        if (!(*handle)->sine_table)
        {
           goto handle_error1;
        }
    }

    /*  set default cfft_plan to 0(length). (and default iisfft_plan to zero length)  */
    (*handle)->cfft.len = 0;  /* 0 length means that cfft should not be  called */
    (*handle)->iisfft.length = 0; /*saftey setting for  iisfft length struct */

    /* use cfft for legth of power two larger than 256. for length below iisfft is faster */
    if (trlen >= 256 && CFFT_PLAN_SUPPORT(trlen)) {
        LC3_INT s = (type == FFT_REAL) ? IIS_FFT_FWD : sign;
        err       = LC3_cfft_plan(&(*handle)->cfft, trlen, s) ? IIS_FFT_NO_ERROR : IIS_FFT_INTERNAL_ERROR;
    } else {
        LC3_INT s = (type == FFT_REAL) ? IIS_FFT_FWD : sign;
        err       = LC3_iisfft_plan(&(*handle)->iisfft, trlen, s);
    }

    return IIS_FFT_NO_ERROR;

handle_error1:
    free((*handle));

    return err;
}

IIS_FFT_ERROR LC3_IIS_RFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT32 len, IIS_FFT_DIR sign)
{
    return create(handle, FFT_REAL, len, sign);
}

static IIS_FFT_ERROR destroy(HANDLE_IIS_FFT* handle)
{
    if (handle && *handle) {
        LC3_iisfft_free(&(*handle)->iisfft);
        LC3_cfft_free(&(*handle)->cfft);
        free(*handle);
        *handle = NULL;
    }
    return IIS_FFT_NO_ERROR;
}

IIS_FFT_ERROR LC3_IIS_CFFT_Create(HANDLE_IIS_FFT* handle, LC3_INT len, IIS_FFT_DIR sign)
{
    return create(handle, FFT_COMPLEX, len, sign);
}


IIS_FFT_ERROR LC3_IIS_xFFT_Destroy(HANDLE_IIS_FFT* handle) { return destroy(handle); }

IIS_FFT_ERROR LC3_IIS_CFFT_Destroy(HANDLE_IIS_FFT* handle) { return destroy(handle); }

static IIS_FFT_ERROR real_destroy(HANDLE_IIS_FFT* handle)
{
    if (handle && *handle) {
        LC3_iisfft_free(&(*handle)->iisfft);
        *handle = NULL;
    }
    return IIS_FFT_NO_ERROR;
}

IIS_FFT_ERROR LC3_IIS_RFFT_Destroy(HANDLE_IIS_FFT* handle) { return real_destroy(handle); }

IIS_FFT_ERROR LC3_IIS_FFT_Apply_CFFT(HANDLE_IIS_FFT handle, const Complex* input, Complex* output)
{
    LC3_FLOAT* dummy;
    if (!handle)
    {
        return IIS_FFT_INTERNAL_ERROR;
    }

    /* check for inplace operation */
    memmove(output, input, sizeof(*input) * handle->len);
    dummy = (LC3_FLOAT*)output;
    if (handle->cfft.len > 0) {
        LC3_cfft_apply(&handle->cfft, dummy, dummy + 1, 2);
    } else {
        LC3_iisfft_apply(&handle->iisfft, dummy);
    }

    return IIS_FFT_NO_ERROR;
}


IIS_FFT_ERROR LC3_IIS_FFT_Apply_RFFT(HANDLE_IIS_FFT handle, const LC3_FLOAT* in, LC3_FLOAT* out)
{
   if (!handle) {
      return IIS_FFT_INTERNAL_ERROR;
   }

   memmove(out, in, sizeof(LC3_FLOAT) * handle->len);

   if (handle->sign == IIS_FFT_BWD) {
      LC3_rfft_pre(handle->sine_table, out, handle->len);
   }

   if (handle->cfft.len > 0) {
      LC3_cfft_apply(&handle->cfft, out, out + 1, 2);
   }
   else {
      LC3_iisfft_apply(&handle->iisfft, out);
   }

   if (handle->sign == IIS_FFT_FWD) {
      LC3_rfft_post(handle->sine_table, out, handle->len);
   }

   return IIS_FFT_NO_ERROR;
}
