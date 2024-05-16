/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void attack_detector_fl(LC3_FLOAT* in, LC3_INT frame_size, LC3_INT fs, LC3_INT* lastAttackPosition, LC3_FLOAT* accNrg, LC3_INT* attackFlag,
                        LC3_FLOAT* attdec_filter_mem, LC3_INT attackHandlingOn, LC3_INT attdec_nblocks, LC3_INT attdec_hangover_threshold)
{
    LC3_FLOAT f_sig[160], block_nrg[4], sum, tmpEne, *ptr, tmp[162];
    LC3_INT   i, j, attackPosition;
    LC3_FLOAT mval;
    LC3_INT frame_size_16k;

    if (attackHandlingOn) {

        mval = 0; j = 0;
        frame_size_16k = attdec_nblocks * 40;
        ptr = &tmp[2];
        
        /* Decimate 48 and 32 kHz signals to 16 kHz */
        if (fs == 48000) {
            j = 0;
            for (i = 0; i < frame_size;) {
                ptr[j] = (in[i] + in[i + 1] + in[i + 2]);
                i      = i + 3;
                j++;
            }
        } else if (fs == 32000) {
            j = 0;
            for (i = 0; i < frame_size;) {
                ptr[j] = (in[i] + in[i + 1]);
                i      = i + 2;
                j++;
            }
        }

        /* Filter */
        ptr[-2] = (LC3_FLOAT)attdec_filter_mem[0];
        ptr[-1] = (LC3_FLOAT)attdec_filter_mem[1];

        attdec_filter_mem[0] = ptr[frame_size_16k - 2];
        attdec_filter_mem[1] = ptr[frame_size_16k - 1];

        for (i = 159; i >= 0; i--) {
            tmpEne = 0;

            tmpEne += ptr[i] * 0.375;
            tmpEne += ptr[i - 1] * (-0.5);
            tmpEne += ptr[i - 2] * (0.125);

            f_sig[i] = tmpEne;
        }

        for (i = 0; i < attdec_nblocks; i++) {
            sum = 0;
            for (j = 0; j < 40; j++) {
                sum += f_sig[j + i * 40] * f_sig[j + i * 40];
            }

            block_nrg[i] = sum;
        }

        *attackFlag    = 0;
        attackPosition = -1;

        for (i = 0; i < attdec_nblocks; i++) {
            tmpEne = block_nrg[i] / 8.5;

            if (tmpEne > MAX(*accNrg, mval)) {
                *attackFlag    = 1;
                attackPosition = i + 1;
            }

            *accNrg = MAX(block_nrg[i], 0.25 * (*accNrg));
        }

        if (*lastAttackPosition > attdec_hangover_threshold) {
            *attackFlag = 1;
        }

        *lastAttackPosition = attackPosition;
    }
}
