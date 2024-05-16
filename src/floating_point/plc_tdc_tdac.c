/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processTdcTdac_fl(const LC3_FLOAT *synth_inp, const LC3_FLOAT *win, LC3_INT32 frame_length, LC3_INT32 la_zeroes, LC3_FLOAT *ola_mem)
{
        LC3_INT32 i, L, LD2, NZ, synth_len;
        LC3_FLOAT synth[(MAX_LEN + MDCT_MEM_LEN_MAX)], *synth1, *synth2, *ola_mem1, *ola_mem2, sz;
        const LC3_FLOAT *win1, *win2, *win3, *win4;

    assert(la_zeroes <= frame_length / 2);

    L         = frame_length;
    LD2       = L/2;
    NZ        = LD2 - la_zeroes;
    synth_len = 2*L - la_zeroes;

    move_float(synth, synth_inp, synth_len);

    /* calculate x_ov[L+la_zeroes] ... x_ov[2*L-1] */
    win1 = &win[L + LD2 - 1];
    win2 = &win[L + LD2];

    win3 = &win[LD2 - 1];
    win4 = &win[LD2];

    synth1 = &synth[L + LD2 - 1 - la_zeroes];
    synth2 = &synth[L + LD2 - la_zeroes];

    ola_mem1 = &ola_mem[LD2 - la_zeroes];
    ola_mem2 = &ola_mem[LD2 - la_zeroes - 1];

    for (i = 0; i < NZ; i++)
    {
        /* analysis windowing + 2N -> N */
        sz = *synth1 * *win1 + *synth2 * *win2;

        /* N -> 2N + synthesis windowing */
        *ola_mem1 = sz * *win3;
        *ola_mem2 = sz * *win4;

        /* pointer update */
        win1--;
        win2++;
        win3--;
        win4++;
        synth1--;
        synth2++;
        ola_mem1++;
        ola_mem2--;
    }

    for (; i < LD2; i++)
    {
        /* analysis windowing + 2N -> N */
        sz = *synth1 * *win1;

        /* N -> 2N + synthesis windowing */
        *ola_mem1 = sz * *win3;

        /* pointer update */
        win1--;
        win2++;
        win3--;
        synth1--;
        synth2++;
        ola_mem1++;
    }
}

