/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"



/* Union holding buffers to conserve stack memory. */

void processMdct_fx(
#ifdef ENABLE_HR_MODE
    Word32 x[],             /* i:   time input signal */
#else
    Word16 x[],             /* i:   time input signal */
#endif
    Word16 x_exp, Word16 N, /* i:   block size N */
#ifdef ENABLE_HR_MODE
    const Word32 w[],       /* i:   window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
#else
    const Word16 w[],       /* i:   window coefficients including normalization of sqrt(2/N) and scaled by 2^4 */
#endif
    Word16       wLen,      /* i:   window length */
#ifdef ENABLE_HR_MODE
    Word32       mem[],     /* i/o: last block of input samples */
#else
    Word16       mem[],     /* i/o: last block of input samples */
#endif
    Word16       memLen,    /* i:   length of last sample block */
    Word32       y[],       /* o:   spectral data */
    Word16 *     y_e,       /* o:   spectal data exponent */
    Word8 *      scratchBuffer)
{
    Counter i;
    Word16  z, s, m;
#ifdef ENABLE_HR_MODE
    Word32 *buf;
#else
    Word16 *buf;
#endif
    Word32 *workBuffer;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processMdct_fx", sizeof(struct {
                   Counter i;
                   Word16  z, s, m;
                   Word16 *buf;
                   Word32 *workBuffer;
               }));
#endif

    /* Buffers overlap since they are not used at the same time */
#ifdef ENABLE_HR_MODE
    buf        = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN */
#else
    buf        = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN */
#endif
    workBuffer = (Word32 *)scratchAlign(scratchBuffer, 0); /* Size = 4 * MAX_LEN */

    /* Init (constant per sample rate) */
    z = (N << 1) - wLen; /* number of leading zeros in window */
    m = N >> 1;          /* half block size */

#ifdef ENABLE_HR_MODE
    basop_memmove(buf, mem, memLen * sizeof(Word32));

    basop_memmove(&buf[memLen], x, (N - memLen) * sizeof(Word32));

    basop_memmove(mem, &x[N - memLen], memLen * sizeof(Word32));
#else
    basop_memmove(buf, mem, memLen * sizeof(Word16));

    basop_memmove(&buf[memLen], x, (N - memLen) * sizeof(Word16));

    basop_memmove(mem, &x[N - memLen], memLen * sizeof(Word16));
#endif

    FOR (i = 0; i < m; i++)
    {
#ifdef ENABLE_HR_MODE
        y[m + i] = Msu_32_32_0(Mpy_32_32_0(w[i], buf[i]), w[2 * m - 1 - i], buf[2 * m - 1 - i]); move32();
#else
        y[m + i] = L_msu0(L_mult0(buf[i], w[i]), buf[2 * m - 1 - i], w[2 * m - 1 - i]); move32();
#endif
    }

    FOR (i = 0; i < z; i++)
    {
#ifdef ENABLE_HR_MODE
        y[m - 1 - i] = Mpy_32_32_0(w[2 * m + i], x[2 * m - memLen + i]); move32();
#else
        y[m - 1 - i] = L_mult0(x[2 * m - memLen + i], w[2 * m + i]); move32();
#endif
    }
    
    FOR (i = i; i < m; i++)
    {
#ifdef ENABLE_HR_MODE
        y[m - 1 - i] = Mac_32_32_0(Mpy_32_32_0(w[2 * m + i], x[2 * m - memLen + i]), w[4 * m - 1 - i], x[4 * m - memLen - 1 - i]); move32();
#else
        y[m - 1 - i] = L_mac0(L_mult0(x[2 * m - memLen + i], w[2 * m + i]), x[4 * m - memLen - 1 - i],
                              w[4 * m - 1 - i]); move32();
#endif
    }

    s = s_max(0, getScaleFactor32(y, N));
    FOR (i = 0; i < N; i++)
    {
        y[i] = L_shl(y[i], s); move32();
    }

    *y_e = sub(sub(x_exp, 2), s);

    /* N=20 only for 2.5ms possible */
    /* maybe implement this a pre init of shift */
    if (sub(N, 20) <= 0)
    {
        *y_e = add(*y_e, 2);
    }
    else if (sub(N, 120) <= 0)
    {
        *y_e = add(*y_e, 1);
    }

    dct_IV(y, y_e, N, workBuffer);

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

