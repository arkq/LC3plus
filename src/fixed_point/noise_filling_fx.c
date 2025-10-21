/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

/*************************************************************************/



void processNoiseFilling_fx(Word32 xq[], Word16 nfseed, Word16 xq_e, Word16 fac_ns_idx, Word16 BW_cutoff_idx,
                            LC3PLUS_FrameDuration frame_dms, Word16 fac_ns_pc, Word16 spec_inv_idx, Word8 *scratchBuffer
#ifdef ENABLE_HR_MODE
                            , Word16 hrmode
#endif
)
{
    Dyn_Mem_Deluxe_In(
        Counter k;
        Word16  nzeros, fac_ns, *ind, c;
        Word16  noisefillwidth, noisefillstart, N;
        Word32  L_tmp, L_tmp_neg, L_tmp_pc, L_tmp_neg_pc;
    );

    noisefillwidth = 0; move16();
    noisefillstart = 0; move16();
    ind = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * MAX_LEN bytes */

    c = 0;                                move16();
    
#ifdef ENABLE_HR_MODE
    if (hrmode == 1)
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
#ifdef CR9_C_ADD_1p25MS
    case LC3PLUS_FRAME_DURATION_1p25MS:
        N              = shr_pos(N, 3);
        noisefillwidth = NOISEFILLWIDTH_1_25MS;
        noisefillstart = NOISEFILLSTART_1_25MS;
        BREAK;
#endif
    case LC3PLUS_FRAME_DURATION_2p5MS:
        N              = shr_pos(N, 2);
        noisefillwidth = NOISEFILLWIDTH_2_5MS;
        noisefillstart = NOISEFILLSTART_2_5MS;
        BREAK;
    case LC3PLUS_FRAME_DURATION_5MS:
        N              = shr_pos(N, 1);
        noisefillwidth = NOISEFILLWIDTH_5MS;
        noisefillstart = NOISEFILLSTART_5MS;
        BREAK;
    case LC3PLUS_FRAME_DURATION_7p5MS:
        N              = add(shr_pos(N, 2), add(shr_pos(N, 2), shr_pos(N, 2)));
        noisefillwidth = NOISEFILLWIDTH_7_5MS;
        noisefillstart = NOISEFILLSTART_7_5MS;
        BREAK;
    case LC3PLUS_FRAME_DURATION_10MS:
        noisefillwidth = NOISEFILLWIDTH;
        noisefillstart = NOISEFILLSTART;
        BREAK;
    case LC3PLUS_FRAME_DURATION_UNDEFINED: assert(0);
    }

    nzeros = -2 * noisefillwidth - 1; move16();

    FOR (k = noisefillstart - noisefillwidth; k < noisefillstart + noisefillwidth; k++)
    {
        if (xq[k] != 0)
        {
            nzeros = -2 * noisefillwidth - 1; move16();
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
            nzeros = -2 * noisefillwidth - 1; move16();
        }
        if (xq[k + noisefillwidth] == 0)
        {
            nzeros = add(nzeros, 1);
        }
        if (nzeros >= 0)
        {
            ind[c++] = k; move16();
        }
    }

    FOR (k = N - noisefillwidth; k < N; k++)
    {
        nzeros = add(nzeros, 1);
        if (nzeros >= 0)
        {
            ind[c++] = k; move16();
        }
    }

    IF (c > 0)
    {
        fac_ns       = shl_pos(sub(8, fac_ns_idx), 11);
        L_tmp        = L_shr_sat(L_deposit_l(fac_ns), sub(xq_e, 16));
        L_tmp_neg    = L_negate(L_tmp);
        L_tmp_pc     = L_shr_sat(L_deposit_l(fac_ns_pc), sub(xq_e, 16));
        L_tmp_neg_pc = L_negate(L_tmp_pc);

        FOR (k = 0; k < c; k++)
        {
            nfseed = extract_l(L_mac0(13849, nfseed, 31821));
            IF (nfseed >= 0)
            {
                IF (ind[k] < spec_inv_idx)
                {
                    xq[ind[k]] = L_tmp; move32();
                }
                ELSE
                {
                    xq[ind[k]] = L_tmp_pc; move32();
                }
            }
            IF (nfseed < 0)
            {
                IF (ind[k] < spec_inv_idx)
                {
                    xq[ind[k]] = L_tmp_neg; move32();
                }
                ELSE
                {
                    xq[ind[k]] = L_tmp_neg_pc; move32();
                }
            }
        }
    }

    Dyn_Mem_Deluxe_Out();
}

