/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

LC3_FLOAT pc_peak_detector(LC3_FLOAT *q_d_prev, LC3_INT32 yLen);

void processPcClassify_fl(LC3_INT32 pitch_present, LC3_INT32 frame_dms, LC3_FLOAT *q_d_prev, LC3_FLOAT *q_old_res, LC3_INT32 yLen, LC3_INT32 spec_inv_idx, LC3_FLOAT stab_fac, LC3_INT32 *bfi)
{
        LC3_INT32 maxPitchBin, xover, i;
        LC3_FLOAT part_nrg, full_nrg;
    
    part_nrg = 0; full_nrg = 0;
    
    if (spec_inv_idx < (4 * frame_dms / 10))
    {
        if (stab_fac < 0.5)
        {
            *bfi = 1;
        } else if (pitch_present == 1)
        {
            maxPitchBin = 8;
            if (frame_dms == 50)
            {
                maxPitchBin = 4;
            }
            xover = pc_peak_detector(q_d_prev, yLen);
            if (spec_inv_idx < xover || spec_inv_idx < maxPitchBin)
            {
                *bfi = 1;
            }
        } else {
            for (i = 0; i < spec_inv_idx; i++)
            {
                part_nrg += LC3_POW(q_old_res[i], 2);
            }
            
            for (i = 0; i < yLen; i++)
            {
                full_nrg += LC3_POW(q_old_res[i], 2);
            }

            if (part_nrg < (0.3 * full_nrg))
            {
                *bfi = 1;
            }
        }
    }
}

LC3_FLOAT pc_peak_detector(LC3_FLOAT *q_d_prev, LC3_INT32 yLen)
{
        LC3_INT32 block_size, thresh1, i, peak, j, k;
        LC3_FLOAT fac, mean_block_nrg, cur_max, block_cent, maxPeak, next_max, prev_max;
    
    mean_block_nrg = 0;
    
    block_size = 3;
    thresh1 = 8;
    fac = 0.3;
    
    for (i = 0; i < yLen; i++)
    {
        mean_block_nrg += LC3_POW(q_d_prev[i], 2);
    }
    
    mean_block_nrg /= yLen;
    
    maxPeak = 0;
    peak = 0;
    
    if (LC3_FABS(q_d_prev[0]) > LC3_FABS(q_d_prev[1]))
    {
        block_cent = LC3_POW(q_d_prev[0], 2) + LC3_POW(q_d_prev[1], 2);
        
        if ((block_cent / block_size) > (thresh1 * mean_block_nrg))
        {
            cur_max = MAX(LC3_FABS(q_d_prev[0]), LC3_FABS(q_d_prev[1]));
            next_max = array_max_abs(&q_d_prev[2], 3);
            
            if (cur_max > next_max)
            {
                maxPeak = block_cent;
                peak = 1;
            }
        }
    }
    
    for (i = 0; i < block_size; i++)
    {
        if ((LC3_FABS(q_d_prev[i + 1]) >= LC3_FABS(q_d_prev[i])) && LC3_FABS(q_d_prev[i + 1]) >= LC3_FABS(q_d_prev[i + 2]))
        {
            block_cent = 0;
            for (j = i; j < i + block_size; j++)
            {
                block_cent += LC3_POW(q_d_prev[j], 2);
            }
            
            if ((block_cent / block_size) > (thresh1 * mean_block_nrg))
            {
                cur_max = array_max_abs(&q_d_prev[i], block_size);
                prev_max = 0;
                
                for (k = i - block_size; k < i; k++)
                {
                    if (k > 0)
                    {
                        prev_max = MAX(LC3_FABS(q_d_prev[k]), prev_max);
                    }
                }
                next_max = array_max_abs(&q_d_prev[i + block_size], block_size);
                
                if ((cur_max >= prev_max) && (cur_max > next_max))
                {
                    if (block_cent > (fac * maxPeak))
                    {
                        peak = i + block_size - 1;
                        if (block_cent >= maxPeak)
                        {
                            maxPeak = block_cent;
                        }
                    }
                }
            }
        }
    }
    
    for (i = block_size; i < yLen - (2 * block_size); i++)
    {
        if ((LC3_FABS(q_d_prev[i + 1]) >= LC3_FABS(q_d_prev[i])) && LC3_FABS(q_d_prev[i + 1]) >= LC3_FABS(q_d_prev[i + 2]))
        {
            block_cent = 0;
            for (j = i; j < i + block_size; j++)
            {
                block_cent += LC3_POW(q_d_prev[j], 2);
            }
            
            if ((block_cent / block_size) > (thresh1 * mean_block_nrg))
            {
                cur_max = array_max_abs(&q_d_prev[i], block_size);
                prev_max = array_max_abs(&q_d_prev[i - block_size], block_size);
                next_max = array_max_abs(&q_d_prev[i + block_size], block_size);
                
                if ((cur_max >= prev_max) && (cur_max > next_max))
                {
                    if (block_cent > (fac * maxPeak))
                    {
                        peak = i + block_size - 1;
                        if (block_cent >= maxPeak)
                        {
                            maxPeak = block_cent;
                        }
                    }
                }
            }
        }
    }
    
    return peak;
}

