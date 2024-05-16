/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#ifndef __BASOP_MPY_H
#define __BASOP_MPY_H

#include "functions.h"

#ifdef ENABLE_HR_MODE
Word32 Mpy_32_16_0(Word32 x, Word16 y);
#endif

Word32 Mac_32_16_0(Word32 z, Word32 x, Word16 y);
Word32 Mpy_32_32_0(Word32 x, Word32 y);
Word32 Mac_32_32_0(Word32 z, Word32 x, Word32 y);
Word32 Msu_32_16_0(Word32 z, Word32 x, Word16 y);
Word32 Msu_32_32_0(Word32 z, Word32 x, Word32 y);

/**
 * \brief     32*16 Bit fractional Multiplication using 40 bit OPS
 *          Performs a multiplication of a 32-bit variable x by
 *          a 16-bit variable y, returning a 32-bit value.
 *
 * \param[i] x
 * \param[i] y
 *
 * \return x*y
 */
Word32 Mpy_32_16(Word32 x, Word16 y);

/**
 * \brief     32*32 Bit fractional Multiplication using 40 bit OPS
 *
 *          Performs a multiplication of a 32-bit variable x by
 *          a 32-bit variable y, returning a 32-bit value.
 *
 * \param[i] x
 * \param[i] y
 *
 * \return x*y
 */
Word32 Mpy_32_32(Word32 x, Word32 y);

#  ifdef ENABLE_HR_MODE
#  define cplxMpy32_32_32_2(re, im, a, b, c, d)                                                                          \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_sub(Mpy_32_32_0(a, c), Mpy_32_32_0(b, d));                                                              \
        im = L_add(Mpy_32_32_0(a, d), Mpy_32_32_0(b, c));                                                              \
    } while (0)
#  endif
#  define cplxMpy32_32_16_2(re, im, a, b, c, d)                                                                          \
    do                                                                                                                 \
    {                                                                                                                  \
        re = L_sub(L_shr_pos(Mpy_32_16(a, c), 1), L_shr_pos(Mpy_32_16(b, d), 1));                                      \
        im = L_add(L_shr_pos(Mpy_32_16(a, d), 1), L_shr_pos(Mpy_32_16(b, c), 1));                                      \
    } while (0)

void cplxMpy_32_16(Word32 *c_Re, Word32 *c_Im, const Word32 a_Re, const Word32 a_Im, const Word16 b_Re,
                   const Word16 b_Im);
void cplxMpy_32_32(Word32 *c_Re, Word32 *c_Im, const Word32 a_Re, const Word32 a_Im, const Word32 b_Re,
                   const Word32 b_Im);

#endif /* __BASOP_SETTINGS_H */
