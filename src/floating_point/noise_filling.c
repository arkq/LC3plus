/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNoiseFilling_fl(LC3_FLOAT xq[], LC3_INT nfseed, LC3_INT fac_ns_idx, LC3_INT bw_stopband, LC3PLUS_FrameDuration frame_dms, LC3_FLOAT fac_ns_pc, LC3_INT spec_inv_idx)
{
    LC3_INT   zeroLines[MAX_LEN];
    LC3_INT   nTransWidth = 0, startOffset = 0, j, k, nzeros = 0, kZeroLines;
    LC3_FLOAT fac_ns = 0;

    switch (frame_dms)
    {
#ifdef CR9_C_ADD_1p25MS
        case LC3PLUS_FRAME_DURATION_1p25MS:
            nTransWidth = 1;
            startOffset = 6;
            break;
#endif
        case LC3PLUS_FRAME_DURATION_2p5MS:
            nTransWidth = 1;
            startOffset = 6;
            break;
        case LC3PLUS_FRAME_DURATION_5MS:
            nTransWidth = 1;
            startOffset = 12;
            break;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            nTransWidth = 2;
            startOffset = 18;
            break;
        case LC3PLUS_FRAME_DURATION_10MS:
            nTransWidth = 3;
            startOffset = 24;
            break;
        case LC3PLUS_FRAME_DURATION_UNDEFINED:
            assert(0);
    }

    fac_ns = (8.0 - fac_ns_idx) / 16.0;

    j = 0;

    for (k = startOffset - nTransWidth; k < startOffset + nTransWidth; k++)
    {
        if (xq[k] != 0)
        {
            nzeros = -2 * nTransWidth - 1;
        }
        if (xq[k] == 0)
        {
            nzeros ++;
        }
    }
    for (k = startOffset; k < bw_stopband - nTransWidth; k++)
    {
        if (xq[k + nTransWidth] != 0)
        {
            nzeros = -2 * nTransWidth - 1;
        }
        if (xq[k + nTransWidth] == 0)
        {
            nzeros ++;
        }
        if (nzeros >= 0)
        {
            zeroLines[j++] = k;
        }
    }

    for (k = bw_stopband - nTransWidth; k < bw_stopband; k++)
    {
        nzeros ++;
        if (nzeros >= 0)
        {
            zeroLines[j++] = k;
        }
    }

    kZeroLines = j;

    for (k = 0; k < kZeroLines; k++) {
        nfseed = (13849 + (nfseed + 32768) * 31821) & 65535;
        nfseed -= 32768;

        if (nfseed >= 0) {
            if (zeroLines[k] < spec_inv_idx)
            {
                xq[zeroLines[k]] = fac_ns;
            } else {
                xq[zeroLines[k]] = fac_ns_pc;
            }
        } else {
            if (zeroLines[k] < spec_inv_idx)
            {
                xq[zeroLines[k]] = -fac_ns;
            } else {
                xq[zeroLines[k]] = -fac_ns_pc;
            }
        }
    }

    return;
}
