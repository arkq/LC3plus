/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "defines.h"

#include "functions.h"


void processPLCmain_fx(Word16 plcMeth, Word16 *concealMethod, Word16 *nbLostFramesInRow, Word16 bfi, Word16 prev_bfi,
                       Word16 frame_length, Word16 la_zeroes,
#ifdef ENABLE_HR_MODE
                       const Word32 w[],
#else
                       const Word16 w[],
#endif
                       Word16 x_fx[], Word16 ola_mem[],
                       Word16 *ola_mem_exp, Word16 q_old_d_fx[], Word16 *q_old_fx_exp, Word32 q_d_fx[],
                       Word16 *q_fx_exp, Word16 yLen, Word16 fs_idx, const Word16 *band_offsets, Word16 bands_number, Word16 *damping,
                       Word16 old_pitch_int, Word16 old_pitch_fr, Word16 *ns_cum_alpha, Word16 *ns_seed,
                       AplcSetup *plcAd, Word16 frame_dms, Word8 *scratchBuffer, Word16 *pc_nbLostFramesInRow
#ifdef ENABLE_HR_MODE
                       , Word16 hrmode
#endif
                       )
{
    IF (sub(bfi, 1) == 0 && plcAd )
    {
        /* FFLC increases the PFLC counter */
        *pc_nbLostFramesInRow = add(*pc_nbLostFramesInRow, 1);
    }

    processPLCclassify_fx(plcMeth, concealMethod, nbLostFramesInRow, bfi, old_pitch_int, frame_length, frame_dms,
                          fs_idx, yLen, q_old_d_fx, band_offsets, bands_number, plcAd, scratchBuffer
#ifdef ENABLE_HR_MODE
                          , hrmode
#endif
                          );

    processPLCapply_fx(*concealMethod, *nbLostFramesInRow, bfi, prev_bfi, frame_length, la_zeroes, w, x_fx, ola_mem,
                       ola_mem_exp, q_old_d_fx, q_old_fx_exp, q_d_fx, q_fx_exp, yLen, fs_idx, damping, old_pitch_int,
                       old_pitch_fr, ns_cum_alpha, ns_seed, frame_dms, plcAd, scratchBuffer
#ifdef ENABLE_HR_MODE
                       , hrmode
#endif
                       );

    IF (bfi == 0)
    {
        processPLCupdateSpec_fx(q_old_d_fx, q_old_fx_exp, q_d_fx, q_fx_exp, yLen);
    }

#ifdef ENABLE_HR_MODE
    IF (plcAd != NULL && (sub(plcAd->PhECU_frame_ms , 10) == 0) && (hrmode == 0))
#else
    IF (plcAd != NULL && (sub(plcAd->PhECU_frame_ms , 10) == 0))
#endif
    {
        processPLCspec2shape_fx(prev_bfi, bfi, q_old_d_fx, yLen, plcAd->PhECU_oold_grp_shape_fx, plcAd->PhECU_old_grp_shape_fx);
    }
 
}


