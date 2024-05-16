/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"


void processNoiseFactor_fx(Word16 *fac_ns_idx, Word16 x_e, Word32 x[],
#    ifdef ENABLE_HR_MODE
                           Word32 xq[],
#    else
                           Word16 xq[],
#    endif
                           Word16 gg, Word16 gg_e, Word16 BW_cutoff_idx, Word16 frame_dms, Word16 target_bytes,
                           Word8 *scratchBuffer
#    ifdef ENABLE_HR_MODE
                           ,Word16 hrmode
#    endif
)
{
    Dyn_Mem_Deluxe_In(Counter k; Word16 nzeros, s1, s2, s3, c, idx, fac_unq, *ind;
                      Word16 noisefillwidth, noisefillstart, N; Word32 Lsum;);



    ind = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN bytes */

    noisefillwidth = 0;
    noisefillstart = 0;
    c              = 0;
    move16();
    
#ifdef ENABLE_HR_MODE
    if (hrmode)
    {
        N = BW_cutoff_bin_all_HR[BW_cutoff_idx];
        move16();    
    }
    else
#endif
    {
        N = BW_cutoff_bin_all[BW_cutoff_idx];
        move16();
    }

    SWITCH (frame_dms)
    {
    case 25:
        N              = shr_pos(N, 2);
        noisefillwidth = NOISEFILLWIDTH_2_5MS;
        noisefillstart = NOISEFILLSTART_2_5MS;
        BREAK;
    case 50:
        N              = shr_pos(N, 1);
        noisefillwidth = NOISEFILLWIDTH_5MS;
        noisefillstart = NOISEFILLSTART_5MS;
        BREAK;
    case 75:
        N              = add(shr_pos(N, 2), add(shr_pos(N, 2), shr_pos(N, 2)));
        noisefillwidth = NOISEFILLWIDTH_7_5MS;
        noisefillstart = NOISEFILLSTART_7_5MS;
        BREAK;
    case 100:
        noisefillwidth = NOISEFILLWIDTH;
        noisefillstart = NOISEFILLSTART;
        BREAK;
    }

    nzeros = -2 * noisefillwidth - 1;
    move16();

    FOR (k = noisefillstart - noisefillwidth; k < noisefillstart + noisefillwidth; k++)
    {
        if (xq[k] != 0)
        {
            nzeros = -2 * noisefillwidth - 1;
            move16();
        }
        if (xq[k] == 0)
        {
            nzeros = add(nzeros, 1);
        }
    }

    FOR (k = noisefillstart; k < N - noisefillwidth; k++)
    {
        if (xq[k + noisefillwidth] != 0)
        {
            nzeros = -2 * noisefillwidth - 1;
            move16();
        }
        if (xq[k + noisefillwidth] == 0)
        {
            nzeros = add(nzeros, 1);
        }
        if (nzeros >= 0)
        {
            ind[c++] = k;
            move16();
        }
    }

    FOR (k = N - noisefillwidth; k < N; k++)
    {
        nzeros = add(nzeros, 1);
        if (nzeros >= 0)
        {
            ind[c++] = k;
            move16();
        }
    }

    IF (c == 0)
    {
        fac_unq = 0;
        move16();
    }
    ELSE
    {

        IF (target_bytes <= 20 && frame_dms == 100)
        {
            Word32 ind_sum;
            Word16 mean_ind;

            Word16 fac_unq1, fac_unq2;

            /* calculate mean index */
            ind_sum = ind[0];
            move32();
            FOR (k = 1; k < c; k++)
            {
                ind_sum = L_add(ind_sum, ind[k]);
            }

            mean_ind = BASOP_Util_Divide3216_Scale(ind_sum, c, &s2);
            mean_ind = shl(mean_ind, s2 + 1);

            assert(0 <= mean_ind && mean_ind <= ind[c - 1]);

            /* calculate noise filling gain for low frequencies */
            s2 = 0; move16();
            if (sub(mean_ind, ind[0]) > 0)
            {
                /* calculate scale to ensure that Lsum does not overflow */
                s2 = s_max(sub(sub(14, norm_s(c)), getScaleFactor32(&x[ind[0]], sub(mean_ind, ind[0]))), 0);
            }      
            Lsum = L_shr_pos_pos(L_abs(x[ind[0]]), s2);

            FOR (k = 1; k < c && ind[k] <= mean_ind; k++)
            {
                /* scale before adding to Lsum */
                Lsum = L_add(Lsum, L_shr_pos_pos(L_abs(x[ind[k]]), s2));
            }
            fac_unq1 = BASOP_Util_Divide3216_Scale(Lsum, k, &s1);
            /* add scale applied during summing */
            s1 = add(s1, s2);
            fac_unq1 = BASOP_Util_Divide1616_Scale(fac_unq1, gg, &s2);
            s3       = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2       = norm_s(fac_unq1);
            test();
            IF (fac_unq1 != 0 && add(s3, s2) < 0)
            {
                fac_unq1 = MAX_16;
                move16();
            }
            ELSE
            {
                s3 = s_min(s_max(s3, -15), 15);
                fac_unq1 = shr_r(fac_unq1, s3);
            }

            /* calculate noise filling gain for high frequencies */
            Lsum = 0;
            move16();
            idx = sub(c, k);
            FOR (; k < c; k++)
            {
                Lsum = L_add(Lsum, L_abs(x[ind[k]]));
            }
            fac_unq2 = BASOP_Util_Divide3216_Scale(Lsum, idx, &s1);
            fac_unq2 = BASOP_Util_Divide1616_Scale(fac_unq2, gg, &s2);
            s3       = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2       = norm_s(fac_unq1);
            test();
            IF (fac_unq2 != 0 && add(s3, s2) < 0)
            {
                fac_unq2 = MAX_16;
                move16();
            }
            ELSE
            {
                s3 = s_min(s_max(s3, -15), 15);
                fac_unq2 = shr_r(fac_unq2, s3);
            }

            /* calculate noise filling gain as minimum over high and low frequencies */
            fac_unq = s_min(fac_unq1, fac_unq2);
        }
        ELSE
        {
            /* calculate scale to ensure that Lsum does not overflow */
            s2 = s_max(sub(sub(14, norm_s(c)), getScaleFactor32(&x[ind[0]], sub(N,ind[0]))), 0);
            Lsum = L_abs(x[ind[0]]);
            FOR (k = 1; k < c; k++)
            {
                /* scale before adding to Lsum */
                Lsum = L_add(Lsum, L_shr_pos_pos(L_abs(x[ind[k]]), s2));
            }
            fac_unq = BASOP_Util_Divide3216_Scale(Lsum, c, &s1);
            /* add scale applied during summing */
            s1 = add(s1, s2);
            fac_unq = BASOP_Util_Divide1616_Scale(fac_unq, gg, &s2);
            s3      = sub(15, add(x_e, add(s1, sub(s2, gg_e))));
            s2      = norm_s(fac_unq);
            test();
            IF (fac_unq != 0 && add(s3, s2) < 0)
            {
                fac_unq = MAX_16;
                move16();
            }
            ELSE
            {
                fac_unq = shr_r(fac_unq, s_max(s_min(s3, 15), -15));
            }
        }
    }

    idx = round_fx(L_sub(0x80000, L_mult(fac_unq, 16)));
    if (sub(idx, 7) > 0)
    {
        idx = 7;
        move16();
    }
    if (idx < 0)
    {
        idx = 0;
        move16();
    }
    *fac_ns_idx = idx;
    move16();

    Dyn_Mem_Deluxe_Out();
}

