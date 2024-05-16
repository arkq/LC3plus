/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef UTIL_H
#define UTIL_H

#include "clib.h"
#include "math.h"

#ifdef _MSC_VER
/* strcasecmp is not available on visual studio */
static LC3_INT strcasecmp(const char* a, const char* b) {
  return _stricmp(a,b);
}
#endif

/* restrict is not available on visual studio */
#ifdef _MSC_VER
#define restrict __restrict
/* inline is not recognized in visual studio 13 required by matlab r2015a in win10 */
#define inline __inline
#endif

/* number of elements in array */
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

/* min max with no side effects */
static inline LC3_INT imin(LC3_INT a, LC3_INT b) { return a < b ? a : b; }
static inline LC3_INT imax(LC3_INT a, LC3_INT b) { return a > b ? a : b; }

/* restrict x to range [min, max] */
static inline LC3_INT iclamp(LC3_INT min, LC3_INT x, LC3_INT max) {
  return x < min ? min : x > max ? max : x;
}
static inline double fcmamp(double min, double x, double max) {
  return x < min ? min : x > max ? max : x;
}
static inline LC3_FLOAT fclampf(LC3_FLOAT min, LC3_FLOAT x, LC3_FLOAT max) {
  return x < min ? min : x > max ? max : x;
}

/* x² */
static inline LC3_FLOAT sqrf(LC3_FLOAT x) { return x * x; }

/* convenience wrappers around memmove */
static inline void move_float(LC3_FLOAT *dst, const LC3_FLOAT *src, LC3_INT len) {
#ifdef WMOPS
    LC3_INT i;
    for (i = 0; i < len; i++)
    {
        dst[i] = src[i];
    }
#else
  memmove(dst, src, len * sizeof(LC3_FLOAT));
#endif
}
static inline void move_int(LC3_INT *dst, const LC3_INT *src, LC3_INT len) {
#ifdef WMOPS
    LC3_INT i;
    for (i = 0; i < len; i++)
    {
        dst[i] = src[i];
    }
#else
  memmove(dst, src, len * sizeof(LC3_INT));
#endif
}

/* convenience wrappers around memset */
static inline void zero_float(LC3_FLOAT *x, LC3_INT len) {
#ifdef WMOPS
    LC3_INT i;
    for (i = 0; i < len; i++)
    {
        x[i] = 0;
    }
#else
  memset(x, 0, len * sizeof(LC3_FLOAT));
#endif
}
static inline void zero_int(LC3_INT *x, LC3_INT len) {
#ifdef WMOPS
    LC3_INT i;
    for (i = 0; i < len; i++)
    {
        x[i] = 0;
    }
#else
  memset(x, 0, len * sizeof(LC3_INT));
#endif
}

/* multiply float vectors element by element, in-place */
static inline void mult_vec(LC3_FLOAT *a, const LC3_FLOAT *b,
                            LC3_INT len) {
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    a[i] *= b[i];
  }
}

/* multiply float vector with constant, in-place */
static inline void mult_const(LC3_FLOAT *a, LC3_FLOAT b, LC3_INT len) {
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    a[i] *= b;
  }
}

/* sum of vector */
static inline LC3_FLOAT sum_vec(const LC3_FLOAT *x, LC3_INT len) {
  LC3_FLOAT sum = 0;
  LC3_INT i = 0;
  for (i = 0; i < len; i++) {
    sum += x[i];
  }
  return sum;
}

/* complex constructor */
static inline Complex cmplx(LC3_FLOAT r, LC3_FLOAT i) { return (Complex){r, i}; }

/* complex a + b */
static inline Complex cadd(Complex a, Complex b) {
  return cmplx(a.r + b.r, a.i + b.i);
}

/* complex a * b */
static inline Complex cmul(Complex a, Complex b) {
  return cmplx(a.r * b.r - a.i * b.i, a.i * b.r + a.r * b.i);
}

/* mac operator */
static inline LC3_FLOAT mac_loop(const LC3_FLOAT *array1, const LC3_FLOAT *array2, LC3_INT len)
{
    LC3_INT i;
    LC3_FLOAT sum = 0.0;

    for (i = 0; i < len; i++)
    {
        sum += (*array1++) * (*array2++);
    }

    return sum;
}

/* complex eᶦˣ */
static inline Complex cexpi(LC3_FLOAT x) { return cmplx(LC3_COS(x), LC3_SIN(x)); }

/* complex -x */
static inline Complex cneg(Complex x) { return cmplx(-x.r, -x.i); }

/* convert string to number. return true on success */
static inline bool str_to_int(const char *str, LC3_INT *value) {
  char *end = NULL;
  long v = str ? strtol(str, &end, 0) : 0;
  *value = (LC3_INT)v;
  return str && *end == 0 && v >= INT_MIN && v <= INT_MAX;
}

/* returns true if str ends with str ends with suffix. ignoring case. str may be
 * NULL */
static inline bool str_ends_with(const char *str, const char *suffix) {
  char *tmp = str ? strrchr(str, suffix[0]) : NULL;
  return tmp && !strcasecmp(tmp, suffix);
}

/* complex a - b */
static inline Complex csub(Complex a, Complex b) {
  return cmplx(a.r - b.r, a.i - b.i);
}

static inline void move_cmplx(Complex *dst, const Complex  *src, LC3_INT32 len) {
  if (len > 0) {
    memmove(dst, src, len * sizeof(Complex));
    assert(src[len - 1].r == dst[len - 1].r && src[len - 1].i == dst[len - 1].i); /*check that Cmplx is stored contiguously*/
    assert(src[0].r == dst[0].r && src[0].i == dst[0].i); /*check that Cmplx is stored contiguously*/
  }
}

static inline void zero_cmplx(Complex *x, LC3_INT32 len) {
   if(len > 0) {
     memset(x, 0, len * sizeof(Complex));
     assert(x[0].r == 0 && x[0].i == 0 &&  x[len-1].r == 0 && x[len-1].i == 0);
   }  
}

static inline Complex realtoc(LC3_FLOAT r) { return cmplx(r, 0); }

/* set float vector to constant */
static inline void set_vec(const LC3_FLOAT c, LC3_FLOAT *x, LC3_INT32 len) {
    LC3_INT32 i = 0;
    for (i = 0; i < len; i++) {
        x[i] = c;
    }
}

/* set float vector to constant */
static inline void set_vec_int(const LC3_INT32 c, LC3_INT32 *x, LC3_INT32 len) {
    LC3_INT32 i = 0;
    for (i = 0; i < len; i++) {
        x[i] = c;
    }
}

static inline LC3_INT32 clz_func(LC3_INT32 inp)
{
#if defined(__clang__) || defined(__GNUC__)
    if (inp == 0)
    {
        return 0;
    }
    return __builtin_clz(inp);

#elif defined(_WIN32) || defined(_WIN64)
    LC3_INT32 leading_zero = 0;
   
    if (_BitScanReverse(&leading_zero, inp))
    {
        return 31 - leading_zero;
    }
    else
        return 0;
    
#else
    LC3_INT32 i = 0;
    int64_t x = inp;
    
    if (inp == 0)
    {
        return 0;
    }
    
    inp = (inp < 0) ? ~inp : inp;

    while (x < (int64_t)0x80000000L)
    {
        inp <<= 1;
        i += 1;
    }

    return i;
#endif
}


#endif
