/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"

#include "functions.h"

#ifdef ENABLE_HR_MODE
Word32 Mpy_32_32_0(Word32 x, Word32 y)
{
  Word32 z;

  z = L_shr(L_add(Mpy_32_32(x, y),1), 1);
  
  return (z);
}

Word32 Mpy_32_16_0(Word32 x, Word16 y)
{
  Word32 z;

  z = L_shr(L_add(Mpy_32_16(x, y),1), 1);
  
  return (z);
}
#endif

Word32 Mpy_32_16(Word32 x, Word16 y)
{
    Word32  mh;
    UWord16 ml;

    Mpy_32_16_ss(x, y, &mh, &ml);

    return (mh);
}


Word32 Mpy_32_32(Word32 x, Word32 y)
{
    Word32  mh;
    UWord32 ml;

    Mpy_32_32_ss(x, y, &mh, &ml);

    return (mh);
}


void cplxMpy_32_16(Word32 *c_Re, Word32 *c_Im, const Word32 a_Re, const Word32 a_Im, const Word16 b_Re,
                   const Word16 b_Im)
{
    *c_Re = L_sub(Mpy_32_16(a_Re, b_Re), Mpy_32_16(a_Im, b_Im)); move32();
    *c_Im = L_add(Mpy_32_16(a_Re, b_Im), Mpy_32_16(a_Im, b_Re)); move32();
}

void cplxMpy_32_32(Word32 *c_Re, Word32 *c_Im, const Word32 a_Re, const Word32 a_Im, const Word32 b_Re,
                   const Word32 b_Im)
{
    *c_Re = L_sub(Mpy_32_32(a_Re, b_Re), Mpy_32_32(a_Im, b_Im)); move32();
    *c_Im = L_add(Mpy_32_32(a_Re, b_Im), Mpy_32_32(a_Im, b_Re)); move32();
}

#ifdef ENABLE_HR_MODE
Word32 Mac_32_16_0(Word32 z, Word32 x, Word16 y)
{
  z = L_add(z, Mpy_32_16_0(x, y));

  return (z);
}

Word32 Mac_32_32_0(Word32 z, Word32 x, Word32 y)
{
  z = L_add(z, Mpy_32_32_0(x, y));

  return (z);
}

Word32 Msu_32_16_0(Word32 z, Word32 x, Word16 y)
{
  z = L_sub(z, Mpy_32_16_0(x, y));

  return (z);
}

Word32 Msu_32_32_0(Word32 z, Word32 x, Word32 y)
{
  z = L_sub(z, Mpy_32_32_0(x, y));

  return (z);
}

#endif /* #ifdef ENABLE_HR_MODE */
