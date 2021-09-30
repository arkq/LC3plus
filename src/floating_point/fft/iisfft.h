/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef IISFFT_H
#define IISFFT_H


#ifndef M_PIl
#define M_PIl 3.1415926535897932384626433832795029L /* pi */
#endif

/* compiler specific macros

   the restrict keyword only gives a improvelent if more than one pointers are
   passed to a function. also note that the MSVC __restrict behaves differently
   from c99, the restrict property is not transferred to aliases.

   alloca is a bit problematic because behavior is not defined in case of stack
   overflow. most probably the program will crash. it might be possible to catch
   those errors but it depends on compiler support. msvc has a safer _malloca
   but gcc has nothing similar. */
#if defined _MSC_VER || defined __INTEL_COMPILER
#include <malloc.h>
#define ALLOCA(size) _alloca(size)
#define restrict __restrict
#define inline __inline
#elif defined __GNUC__ || defined __clang__
#define ALLOCA(size) __builtin_alloca(size)
#define restrict __restrict__
#define inline __inline
#elif defined __TI_COMPILER_VERSION__
#include <assert.h>
#define ALLOCA(size) (assert(0 && "ALLOCA is not present for your compiler"), NULL)
#warn "no stack allocation for you compiler"
#else
#error "no stack allocation for your compiler"
#endif


#define IISFFT_MAXSTACKLENGTH 1024
#define IISFFT_MAXFACTORS 10

typedef struct {
    LC3_INT*   scratch2;
    LC3_INT    length;
    LC3_INT    sign;
    LC3_INT    num_factors;
    LC3_INT    factors[IISFFT_MAXFACTORS];
    LC3_INT    isPrime[IISFFT_MAXFACTORS];
} Iisfft;

typedef enum {
    IIS_FFT_NO_ERROR = 0,
    IIS_FFT_INTERNAL_ERROR, /**< a mystical error appeard */
    IIS_FFT_LENGTH_ERROR,   /**< the requested fft length is not supported */
    IIS_FFT_MEMORY_ERROR    /**< memory allocation failed */
} IIS_FFT_ERROR;

typedef enum {
    IIS_FFT_FWD = -1, /**< forward transform */
    IIS_FFT_BWD = 1   /**< inverse / backward transform */
} IIS_FFT_DIR;

/* plan, apply and free forward / backward fft */
IIS_FFT_ERROR LC3_iisfft_plan(Iisfft* handle, LC3_INT length, LC3_INT sign);
void          LC3_iisfft_apply(Iisfft* handle, LC3_FLOAT* x);
void          LC3_iisfft_free(Iisfft* handle);

/* fft related helper functions */
void       LC3_create_sine_table(LC3_INT32 len, LC3_FLOAT *sine_table);

void       LC3_rfft_pre(const LC3_FLOAT* restrict sine_table, LC3_FLOAT* restrict buf, LC3_INT len);
void       LC3_rfft_post(const LC3_FLOAT* restrict sine_table, LC3_FLOAT* restrict buf, LC3_INT len);
void       LC3_fftf_interleave(const LC3_FLOAT* restrict re, const LC3_FLOAT* restrict im, LC3_FLOAT* restrict out,
                               LC3_INT len);
void LC3_fftf_deinterleave(const LC3_FLOAT* restrict in, LC3_FLOAT* restrict re, LC3_FLOAT* restrict im, LC3_INT len);

#endif /* IISFFT_H */
