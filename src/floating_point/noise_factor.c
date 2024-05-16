/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processNoiseFactor_fl(LC3_INT* fac_ns_idx, LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gg, LC3_INT BW_cutoff_idx, LC3_INT frame_dms,
                           LC3_INT target_bytes
)
{
    LC3_INT sumZeroLines = 0, kZeroLines = 0, startOffset = 0, nTransWidth = 0, i = 0, j = 0, k = 0, m = 0, nzeros = 0;
    LC3_FLOAT fac_ns_unq = 0, idx = 0, nsf1 = 0, nsf2 = 0;
    LC3_INT   zeroLines[MAX_LEN];

    switch (frame_dms)
    {
        case 25:
            nTransWidth = 1;
            startOffset = 6;
            break;
        case 50:
            nTransWidth = 1;
            startOffset = 12;
            break;
        case 75:
            nTransWidth = 2;
            startOffset = 18;
            break;
        case 100:
            nTransWidth = 3;
            startOffset = 24;
            break;
    }

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
    for (k = startOffset; k < BW_cutoff_idx - nTransWidth; k++)
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

    for (k = BW_cutoff_idx - nTransWidth; k < BW_cutoff_idx; k++)
    {
        nzeros ++;
        if (nzeros >= 0)
        {
            zeroLines[j++] = k;
        }
    }

    if (j == 0) {
        fac_ns_unq = 0;
    }
    else
    {
        kZeroLines = j;

        fac_ns_unq = 0;
        for (j = 0; j < kZeroLines; j++) {
            fac_ns_unq += LC3_FABS(x[zeroLines[j]]);
        }

        fac_ns_unq /= (gg) * kZeroLines;



        if (kZeroLines > 1 && target_bytes <= 20 && frame_dms == 100) {

            j = 0, k = 0, nsf1 = 0, nsf2 = 0, sumZeroLines = 0;

            for (i = 0; i < kZeroLines; i++) {
                sumZeroLines += zeroLines[i];
            }

            m = floor(sumZeroLines / kZeroLines);

            for (i = 0; i < kZeroLines; i++) {
                if (zeroLines[i] <= m) {
                    j++;
                    nsf1 += LC3_FABS(x[zeroLines[i]]);
                }
                else {
                    nsf2 += LC3_FABS(x[zeroLines[i]]);
                    k++;
                }
            }

            nsf1 /= (gg) * j;
            nsf2 /= (gg) * k; 

            fac_ns_unq = MIN(nsf1, nsf2);
        }

    }

    idx = round(8 - 16 * fac_ns_unq);
    idx = MIN(MAX(idx, 0), 7);

    *fac_ns_idx = idx;
}
