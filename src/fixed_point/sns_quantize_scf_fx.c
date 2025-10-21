/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "functions.h"
#ifdef CR9_C_ADD_1p25MS_LRSNS
#include "constants.h"
#include <stdint.h>   /* req  for INT32_MAX in  msvc */
#endif

static Word16 stage1_base(                    /* o  :  idx                                 */
                          const Word16 *t,    /* i  :  target SCFs                         */
#ifdef ENABLE_HR_MODE
                          const Word32 *cdbk, /* i  :  SCFs cdbk                           */
#else
                          const Word16 *cdbk, /* i  :  SCFs cdbk                           */
#endif
                          const Word16  R     /* i  :  number of rows in codebook          */
)
{
    Counter row;
    Word16  k_ptr, idx;
    Word32  L_min_mse, L_mse;
    Counter col;
#ifdef ENABLE_HR_MODE
    Word32  err;
#else
    Word16  err;
#endif

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("stage1_base", sizeof(struct {
                   Counter row, col;
                   Word16  k_ptr, idx, err;
                   Word32  L_min_mse, L_mse;
               }));
#endif
    BASOP_sub_sub_start("stage1_base");

/* find first vector error energy for  */
/* loop probably works with saturation , but it should not occur anyway */
    L_min_mse = L_add(0, 0);          /*  init acc with absolute  min mse sofar */
    FOR (col = 0; col < M / 2; col++) /* fixed to 8 elements */
    {
#ifdef ENABLE_HR_MODE
        err = L_sub(cdbk[col], L_deposit_h(t[col])); /* cdbk max abs value is 2048 = 2.^11 , max nb col is 2^3  max
                                         target is approx similar (2.^14/M)*2  = +/- 2048 , errmax is 4096   */
        L_min_mse = L_add(L_min_mse, Mpy_32_32(err, err));
#else
        err = sub(cdbk[col], t[col]); /* cdbk max abs value is 2048 = 2.^11 , max nb col is 2^3  max target is approx
                                         similar (2.^14/M)*2  = +/- 2048 , errmax is 4096   */
        L_min_mse = L_mac0(L_min_mse, err, err); /*  max L_min_mse is 8*4096*4096 =2.^(3+12+12) =  2.^27  */
#endif
    }

    idx = 0; move16();

    k_ptr = M / 2; move16(); /* ptr init to second row */
    FOR (row = 1; row < R; row++)
    {
        /* loop probably works with saturation , but it should not occur anyway */

        L_mse = L_add(L_min_mse, 0);      /* init acc with min mse sofar , */
        FOR (col = 0; col < M / 2; col++) /* fixed to 8 elements */
        {
#ifdef ENABLE_HR_MODE
            err = L_sub(cdbk[k_ptr++], L_deposit_h(t[col]));
            L_mse = L_sub(L_mse, Mpy_32_32(err, err));
#else
            err   = sub(cdbk[k_ptr++], t[col]);
            L_mse = L_msu0(L_mse, err,
                           err); /* NB subtraction  from best MSE error  sofar in acc , saturation may not occur */
#endif
        }

        L_min_mse = L_sub(L_min_mse, L_max(L_mse, 0L)); /* ALWAYS update best MSE  error sofar    */

        if (L_mse > 0L) /*  if acc value  still is positive a new lower error energy vector was found in this row   */
        {
            idx = row; move16(); /* update  1-8 bits idx  */
        }

        /* this inner loop(always updating L_min_mse),          (L_msu, if )    consumes AV 19, WC  ~20 STL  cycles  ,
                                      compared to a conventional(L_mac, IF( ) )          AV 21  WC  ~23 STL  cycles per
           loop  */
    }
    ASSERT(idx >= 0 && idx < R);

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif

    return idx;
}

static void first_stage_split_search(
#ifdef ENABLE_HR_MODE      
                                     const Word32 *cbk_LF, const Word32 *cbk_HF,
#else
                                     const Word16 *cbk_LF, const Word16 *cbk_HF,
#endif
                                     const Word16 *target,
                                     const Word16 nbCbkEntries, Word16 *idxLF, Word16 *idxHF)
{
    /* find  base index for  each   SCF split  */
    *idxLF = stage1_base(target, cbk_LF, nbCbkEntries);
    *idxHF = stage1_base((&target[M / 2]), cbk_HF, nbCbkEntries);
}

static void processDeQuantize_stage1ScfDecStage1_fx(
#ifdef ENABLE_HR_MODE
                                                    const Word32 *cbk_LF, const Word32 *cbk_HF, 
#else
                                                    const Word16 *cbk_LF, const Word16 *cbk_HF,
#endif
                                                    Word16 st1_idx0, Word16 st1_idx1,
#ifdef ENABLE_HR_MODE
                                                    Word32 *st1_vector
#else
                                                    Word16 *st1_vector
#endif 
                                                    )
{
    Counter col;
    Word16 offset0, offset1;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processDeQuantize_stage1ScfDecStage1_fx", sizeof(struct {
                   Word16  offset0, offset1;
                   Counter col;
               }));
#endif

    offset0 = shl_pos(st1_idx0, 3); /* mult by M/2 */
    offset1 = shl_pos(st1_idx1, 3);
    FOR (col = 0; col < M / 2; col++)
    {
        st1_vector[col]     = cbk_LF[offset0++]; move16();
        st1_vector[col + 8] = cbk_HF[offset1++]; move16();
#ifdef ENABLE_HR_MODE
        move16(); 
        move16();
#endif
    }
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

void downshift_w32_arr(Word32 *w32_arr, Word16 *w16_arr, Word16 shft_val, Word32 len)
{
    Word32 i;
    FOR (i = 0; i < len; i++)
    {
        w16_arr[i] = extract_l(L_shr_pos(w32_arr[i], shft_val));
        move16();
    }
    return;
}

void round_w32tow16_arr(Word32 *w32_arr, Word16 *w16_arr,Word32 len)
{
    Word32 i;
    FOR (i = 0; i < len; i++)
    {
        w16_arr[i] = round_fx(w32_arr[i]);
        move16();
    }
    return;
}

static void processQuantize_stage1ScfEncStage1_fx(const Word16 *target_st1,
#ifdef ENABLE_HR_MODE
                                                  Word32 *st1_vector,
#else
                                                  Word16 *st1_vector,
#endif
                                                  Word16 *st1_idx0Ptr, Word16 *st1_idx1Ptr)

{
    BASOP_sub_sub_start("processQuantize_stage1ScfEncStage1_fx");

#ifdef ENABLE_HR_MODE
    first_stage_split_search(st1SCF0_7_base5_32x8_Q27, st1SCF8_15_base5_32x8_Q27, target_st1, SCF_STAGE1_NBCDKENTRIES, 
                             st1_idx0Ptr, st1_idx1Ptr);

    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q27, st1SCF8_15_base5_32x8_Q27, *st1_idx0Ptr,
                                            *st1_idx1Ptr, st1_vector);
#else
    first_stage_split_search(st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, target_st1, SCF_STAGE1_NBCDKENTRIES,
                             st1_idx0Ptr, st1_idx1Ptr);

    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, *st1_idx0Ptr,
                                            *st1_idx1Ptr, st1_vector);
#endif

    BASOP_sub_sub_end();

    return;
}

/* gain-shape MSE search in warped SCF-residual domain,  synthesis in SCF resiudal domain allows for easy weighting */

static void pvq_enc_find_best_submode_pre_post_fx(
#ifdef ENABLE_HR_MODE
    const Word32 *target_st2, /* this target is in the linearized  warped domain , same as input to PVQ search  */
#else
    const Word16 *target_st2, /* this target is in the linearized  warped domain , same as input to PVQ search  */
#endif
    const Word16 *enc_pulses_far, Word16 *enc_pulses_near, const Word16 *enc_pulsesA, const Word16 *enc_pulsesB,
    Word16 *sub_mode_ptr, Word16 *i_gain_ptr,
#ifdef ENABLE_HR_MODE
    Word32 *enc_adj_glob_warped_vec, 
#else
    Word16 *enc_adj_glob_warped_vec,
#endif
    Word8 *scratchBuffer) /* Size = 18 * M */
{

    Counter       L_section, idx;
#ifdef ENABLE_HR_MODE
    const Word32 *search_en1shape[N_SCF_SHAPES_ST2];
#else
    const Word16 *search_en1shape[N_SCF_SHAPES_ST2];
#endif
    const Word16 *search_gainTab[N_SCF_SHAPES_ST2];
    Word16        search_n_gains[N_SCF_SHAPES_ST2];
    Word32        L_mse, L_mse_min, L_idx;
    Word16 *      pulses_far, *pulses_near, *pulsesAB, *pulsesA;
#ifdef ENABLE_HR_MODE
    Word32 *      target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
#else
    Word16 *      target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
#endif
#ifdef ENABLE_HR_MODE
    Word32  tmp, err;
#else
    Word16  tmp, err;
#endif
    Counter i;

#ifdef DYNMEM_COUNT
#ifdef ENABLE_HR_MODE
    Dyn_Mem_In("pvq_enc_find_best_submode_pre_post_fx", sizeof(struct {
                   Counter i, L_section, idx;
                   Word32 *search_en1shape[N_SCF_SHAPES_ST2];
                   Word16 *search_gainTab[N_SCF_SHAPES_ST2];
                   Word16  search_n_gains[N_SCF_SHAPES_ST2];
                   Word32  L_mse, L_mse_min, L_idx;
                   Word16 *pulses_far, *pulses_near, *pulsesAB, *pulsesA;
                   Word32 *target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
                   Word32  tmp, err;
               }));
#else
    Dyn_Mem_In("pvq_enc_find_best_submode_pre_post_fx", sizeof(struct {
                   Counter i, L_section, idx;
                   Word16 *search_en1shape[N_SCF_SHAPES_ST2];
                   Word16 *search_gainTab[N_SCF_SHAPES_ST2];
                   Word16  search_n_gains[N_SCF_SHAPES_ST2];
                   Word32  L_mse, L_mse_min, L_idx;
                   Word16 *pulses_far, *pulses_near, *pulsesAB, *pulsesA;
                   Word16 *target_w, *shape_far, *shape_near, *shapeAB, *shapeA;
                   Word16  tmp, err;
               }));
#endif /* ENABLE_HR_MODE */
#endif /* DYNMEM_COUNT */

    pulses_near = (Word16 *)scratchAlign(scratchBuffer, 0); /* Size = 2 * M */

    pulsesAB = (Word16 *)scratchAlign(pulses_near, sizeof(*pulses_near) * M); /* Size = 2 * M */

    pulsesA = (Word16 *)scratchAlign(pulsesAB, sizeof(*pulsesAB) * M); /* Size = 2 * M */

#ifdef ENABLE_HR_MODE
    target_w = (Word32 *)scratchAlign(pulsesA, sizeof(*pulsesA) * M); /* Size = 2 * M */

    shape_near = (Word32 *)scratchAlign(target_w, sizeof(*target_w) * M); /* Size = 4 * M */

    shapeAB = (Word32 *)scratchAlign(shape_near, sizeof(*shape_near) * M); /* Size = 4 * M */

    shapeA = (Word32 *)scratchAlign(shapeAB, sizeof(*shapeAB) * M); /* Size = 4 * M */
#else
    target_w = (Word16 *)scratchAlign(pulsesA, sizeof(*pulsesA) * M); /* Size = 2 * M */

    shape_near = (Word16 *)scratchAlign(target_w, sizeof(*target_w) * M); /* Size = 2 * M */

    shapeAB = (Word16 *)scratchAlign(shape_near, sizeof(*shape_near) * M); /* Size = 2 * M */

    shapeA = (Word16 *)scratchAlign(shapeAB, sizeof(*shapeAB) * M); /* Size = 2 * M */
#endif

    pulses_far = (Word16 *)scratchAlign(shapeA, sizeof(*shapeA) * M); /* Size = 2 * M */

#ifdef ENABLE_HR_MODE
    shape_far = (Word32 *)scratchAlign(pulses_far, sizeof(*pulses_far) * M); /* Size = 2 * M */
#else
    shape_far = (Word16 *)scratchAlign(pulses_far, sizeof(*pulses_far) * M); /* Size = 2 * M */
#endif

    BASOP_sub_sub_start("pvq_enc_find_best_submode_pre_post_fx");

    /* construct pulse vectors and en1 normalized shape vectors  */ /* use shape Q in Q14 */
    basop_memmove(pulses_far, enc_pulses_far, M * sizeof(*pulses_far));
    basop_memmove(pulses_near, enc_pulses_near, M * sizeof(*pulses_near));
    basop_memmove(target_w, target_st2, M * sizeof(*target_w));

    pvq_dec_en1_normQ14_fx(shape_near, pulses_near, sns_Kval[2][0], M); /* near outlier mode  */
    pvq_dec_en1_normQ14_fx(shape_far, pulses_far, sns_Kval[3][0], M);   /* far outlier mode  */

    /* regular mode(with a split),   prepare vectors  of full length M */
    basop_memmove(pulsesAB, enc_pulsesA, N_SETA * sizeof(*pulsesAB));
    basop_memmove(pulsesA, enc_pulsesA, N_SETA * sizeof(*pulsesA));

    FOR (i = N_SETA; i < M; i++)
    {
        pulsesAB[i] = enc_pulsesB[sub(i, N_SETA)]; move16();
    }

    IF (M > N_SETA)
    {
        basop_memset(&pulsesA[N_SETA], 0, (M - N_SETA) * sizeof(*pulsesA));
    }

    pvq_dec_en1_normQ14_fx(shapeAB, pulsesAB, sns_Kval[0][0], M);
    /* regular AB , b_pulses = 1 ;*/ /* OPT: combine  with shapeA */

    pvq_dec_en1_normQ14_fx(shapeA, pulsesA, sns_Kval[1][0], M);
    /* regular A ,  b_pulses = 0 */ /* OPT:  M-> N_SETA */

    /* setup search structure */

    /* now aligned with order of  j  {regular=0, regular_lf=1, outlier_near=2, outlier far=3}  */

    search_en1shape[0] = shapeAB;
    search_gainTab[0]  = sns_gaintabPtr[0];
    search_n_gains[0]  = sns_gainSz[0]; /* assumes whole bits */

    search_en1shape[1] = shapeA;
    search_gainTab[1]  = sns_gaintabPtr[1];
    search_n_gains[1]  = sns_gainSz[1]; /* assumes whole bits */

    search_en1shape[2] = shape_near;
    search_gainTab[2]  = sns_gaintabPtr[2];
    search_n_gains[2]  = sns_gainSz[2]; /*assume whole bits */

    search_en1shape[3] = shape_far;
    search_gainTab[3]  = sns_gaintabPtr[3];
    search_n_gains[3]  = sns_gainSz[3]; /*assume whole bits */

    /* start actual search loop */

    /* basic raw MSE loop,   */
    L_mse_min = INT_MAX;         move32();
    L_idx     = L_deposit_l(-1); /* section in low 2  bits* gain idx above */

    FOR (L_section = 0; L_section < N_SCF_SHAPES_ST2; L_section++)
    {
        /* raw MSE  over gain and shape */
        FOR (idx = 0; idx < search_n_gains[L_section]; idx++)
        {
            /* MSE ( proc_target_local[i]-adjGain[i]*en1Shape[i] ) */

            L_mse = L_deposit_l(0);
            FOR (i = 0; i < M; i++)
            {
#ifdef ENABLE_HR_MODE
                tmp   = Mpy_32_16(search_en1shape[L_section][i], search_gainTab[L_section][idx]); /* Q30 + 14 - 15 = Q29 */
                err   = L_sub(target_w[i], tmp);                                                /*  both in  Q29      */
                L_mse = L_add(L_mse, L_shr_pos(Mpy_32_32(err, err), 1)); /* Q29+29-31 = Q27 */
#else
                tmp   = mult_r(search_gainTab[L_section][idx], search_en1shape[L_section][i]); /* Q15+14+1-16= Q14 */
                err   = sub(target_w[i], tmp);                                                 /*  both in  Q14      */
                L_mse = L_mac0(L_mse, err, err);                                               /* Q14+14 = Q28 */
#endif                                           /* Q14+14 = Q28 */
            }

            IF (L_sub(L_mse, L_mse_min) < 0) /* OPT: always update L_mse_min) */
            {
                L_mse_min = L_mse;                          move32();
                L_idx     = L_mac0(L_section, idx, 1 << 2); /* save both section and gain  idx */
            }
        } /* gains */
    }     /*submodes*/

    L_section = L_and(0x3L, L_idx); /* section was stored in two lowest bits */
    ASSERT(L_section >= 0 && L_section <= 3);
    *i_gain_ptr = extract_l(L_shr_pos(L_idx, 2)); /*1,2,3 bit gain */
    ASSERT(*i_gain_ptr >= 0 && *i_gain_ptr <= 7);

    /* returns a scaled and transformed vector, ___EXACTLY__ as a decoder would scale it */
    ASSERT(enc_adj_glob_warped_vec != NULL);
    {
        /* warp/rotate search result to SCF residual domain */
#ifdef ENABLE_HR_MODE
        idct32_32_fx(search_en1shape[L_section], target_w);
#else
        idct16_fx(search_en1shape[L_section], target_w); /* fwd synthesis  warping */
#endif
        /* actual synthesis gain scaling in SCF-residual domain, for easy weighting analysis  */
        pvq_dec_scale_vec_fx(target_w, search_gainTab[L_section][*i_gain_ptr], enc_adj_glob_warped_vec);
    }

    *sub_mode_ptr = extract_l(L_section);
    move16(); /* 0,1,2,3 */

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return;
}

static void processQuantize_stage2ScfEncStage2_fx(
#ifdef ENABLE_HR_MODE
                                                  const Word32 *target_st2, Word32 *st2_vector,
#else
                                                  const Word16 *target_st2, Word16 *st2_vector,
#endif
                                                  Word32 *L_prm_idx,
                                                  Word16 submodes, Word8 *scratchBuffer) /* Size = 26 * M + 48 */
{                                                                                        /*func */
#ifdef ENABLE_HR_MODE
    Word32 *proc_target;
#else
    Word16 *proc_target;
#endif
    
    Word16 *enc_pulses_far, *enc_pulses_near, *enc_pulsesA, *enc_pulsesB;

    Word16 *pulses_fin, *pulses_proj;
    Word32  L_tmp;

    Word8 *buffer_pvq_enc_find_best_submode_pre_post_fx;

    PvqEntry_fx enc_PVQ_OA, enc_PVQ_B;
    Word16      submode, i_gain, submodeMSB, submodeLSB;
    Word32 *    L_search_corr, *L_search_en;
#ifdef ENABLE_HR_MODE
    Word16 *    proc_target_lp;
#endif

#ifdef DYNMEM_COUNT
#ifdef ENABLE_HR_MODE
    Dyn_Mem_In("processQuantize_stage2ScfEncStage2_fx", sizeof(struct {
                   Word32 *proc_target;
                   Word16 *enc_pulses_far, *enc_pulses_near, *enc_pulsesA, *enc_pulsesB;
                   Word16 *pulses_fin, *pulses_proj;
                   Word32  L_tmp;
                   Word8 *buffer_pvq_enc_find_best_submode_pre_post_fx;
                   PvqEntry_fx enc_PVQ_OA, enc_PVQ_B;
                   Word16      submode, i_gain, submodeMSB, submodeLSB;
                   Word32 *    L_search_corr, *L_search_en;
                   Word16 *    proc_target_lp;
               }));
#else
    Dyn_Mem_In("processQuantize_stage2ScfEncStage2_fx", sizeof(struct {
                   Word16 *proc_target;
                   Word16 *enc_pulses_far, *enc_pulses_near, *enc_pulsesA, *enc_pulsesB;
                   Word16 *pulses_fin, *pulses_proj;
                   Word32  L_tmp;
                   Word8 *buffer_pvq_enc_find_best_submode_pre_post_fx;
                   PvqEntry_fx enc_PVQ_OA, enc_PVQ_B;
                   Word16      submode, i_gain, submodeMSB, submodeLSB;
                   Word32 *    L_search_corr, *L_search_en;
               }));
#endif /* ENABLE_HR_MODE */
#endif /* DYNMEM_COUNT */

#ifdef ENABLE_HR_MODE
    buffer_pvq_enc_find_best_submode_pre_post_fx = (Word8 *) scratchAlign(scratchBuffer, 0);
    proc_target         = (Word32 *) scratchAlign(buffer_pvq_enc_find_best_submode_pre_post_fx, sizeof(*buffer_pvq_enc_find_best_submode_pre_post_fx) * 28 * M);
#else
    buffer_pvq_enc_find_best_submode_pre_post_fx = scratchAlign(scratchBuffer, 0); /* Size = 18 * M */
    proc_target =
        (Word16 *)scratchAlign(buffer_pvq_enc_find_best_submode_pre_post_fx,
                               sizeof(*buffer_pvq_enc_find_best_submode_pre_post_fx) * 18 * M); /* Size = 2 * M */
#endif

    enc_pulses_near = (Word16 *)scratchAlign(proc_target, sizeof(*proc_target) * M);         /* Size = 2 * M */
    enc_pulsesA     = (Word16 *)scratchAlign(enc_pulses_near, sizeof(*enc_pulses_near) * M); /* Size = 2 * N_SETA */
    enc_pulsesB     = (Word16 *)scratchAlign(enc_pulsesA, sizeof(*enc_pulsesA) * N_SETA);    /* Size = 2 * N_SETB */
    pulses_fin = (Word16 *)scratchAlign(enc_pulsesB, sizeof(*enc_pulsesB) * N_SETB); /* Size = 2 * N_SCF_SHAPES_ST2 */
    pulses_proj =
        (Word16 *)scratchAlign(pulses_fin, sizeof(*pulses_fin) * N_SCF_SHAPES_ST2); /* Size = 2 * N_SCF_SHAPES_ST2 */
    L_search_corr =
        (Word32 *)scratchAlign(pulses_proj, sizeof(*pulses_proj) * N_SCF_SHAPES_ST2); /* Size = 4 * N_SCF_SHAPES_ST2 */
    L_search_en    = (Word32 *)scratchAlign(L_search_corr,
                                         sizeof(*L_search_corr) * N_SCF_SHAPES_ST2); /* Size = 4 * N_SCF_SHAPES_ST2 */
    enc_pulses_far = (Word16 *)scratchAlign(L_search_en, sizeof(*L_search_en) * N_SCF_SHAPES_ST2); /* Size = 2 * M */

#ifdef ENABLE_HR_MODE
    proc_target_lp      = (Word16 *)buffer_pvq_enc_find_best_submode_pre_post_fx; /* size = 2*M */
#endif
    
    BASOP_sub_sub_start("processQuantize_stage2ScfEncStage2_fx");

    /* fixed setup for a given  bitrate of 38 ,  no  moves needed */
    /* k_far  = sns_Kval[3][0]; */
    /* k_near = sns_Kval[2][0]; */
    /* kA     = sns_Kval[1][0]; */ /* regular, regular_lf */
                                   /* kB is always  1 */

    /* NB  these search indecese do not correspond exactly to specification shape_index j */

    pulses_fin[0] = sns_Kval[3][0]; /* far   6 */
    pulses_fin[1] = sns_Kval[2][0]; /* near  8 */
    pulses_fin[2] = sns_Kval[1][0]; /* section A     10 */
    pulses_fin[3] = sns_Kval[0][1]; /* section B     1 */

    pulses_proj[0] = sns_Kval[3][0];
    pulses_proj[1] = 0;
    pulses_proj[2] = 0;
    pulses_proj[3] = 0;

    /*  pre_process  */
#ifdef ENABLE_HR_MODE
    dct32_fx(target_st2, proc_target); /* enc analysis */
    downshift_w32_arr(proc_target, proc_target_lp, 16, M);

    /* get the initial four integer shape candidate vectors,  no normalization at this stage  */
    pvq_enc_search_fx(proc_target_lp, enc_pulses_far, enc_pulses_near, enc_pulsesA, enc_pulsesB, L_search_corr,
                      L_search_en, pulses_fin, pulses_proj, M, N_SETA);
#else
    Word32 target_st2_32[M]; Word32 proc_target_32[M]; int i;
    
    FOR (i = 0; i < M; i++)
    {
        target_st2_32[i] = L_shl_pos(target_st2[i], 16);
    }
    
    dct32_fx(target_st2_32, proc_target_32); /* enc analysis */
    
    downshift_w32_arr(proc_target_32, proc_target, 16, M);

    /* get the initial four integer shape candidate vectors,  no normalization at this stage  */
    pvq_enc_search_fx(proc_target, enc_pulses_far, enc_pulses_near, enc_pulsesA, enc_pulsesB, L_search_corr,
                      L_search_en, pulses_fin, pulses_proj, M, N_SETA);
#endif

    /* scale with gains a after a  unit energy fwd transform  */
    /* apply transform to each candidate shape vector priot  to gain-shape search loop */
    submode = submodes; /* used as input solely to debug/unit test a specific shape mode  */

    /*target should be in a  linearized residual domain target */
    /* search pre, synthesis  post*/
    pvq_enc_find_best_submode_pre_post_fx(proc_target, enc_pulses_far, enc_pulses_near, enc_pulsesA, enc_pulsesB,
                                          &submode, &i_gain, st2_vector,
                                          buffer_pvq_enc_find_best_submode_pre_post_fx); /* Q14 tr out */

    /* send parameters  to multiplexor as a series/vector  of Long Words */
    /*    0 :    0..3  submode  */
    /*    1 :    0..7  gain_ind  */
    /*    2 :    0..1  LeadSign ind */
    /*    3 :    25 bit     MPVQ index    outl_near or  A  part  */
    /*    4 :    3.7 to 21 bits  MPVQ index           B  part  OR   -2  */

    L_prm_idx[0] = L_deposit_l(submode); /*  complete submode fwd'ed  to ari_codec as  0,1,2,3  */

    submodeMSB = shr_pos(submode, 1);                       /* LSB of submode , sent as main submode bit  */
    submodeLSB = s_and(submode, 0x1); /* LSB of submode  */ /*   sent via shape param  */

    /* gain, shape indicese , incl. calls to  MPVQ indexing */
    IF (submodeMSB == 0)
    { /* regular modes::   j=0(reg=AB)  or 1(reg_lf  A)  */ /* regular mode, with two or one shape indices  */

        /* assume regular_lf part ,  shape_j == 1 */
        enc_PVQ_OA =
            mpvq_index_fx(enc_pulsesA, N_SETA, sns_Kval[submode][0]); /* o : leading_sign_index, index, size, k_val */
        L_prm_idx[2] = L_deposit_l(enc_PVQ_OA.lead_sign_ind);         /*LS set A */

        ASSERT(enc_PVQ_OA.size == (UWord32)sns_MPVQ_Sz[submode][0]);
        L_prm_idx[3] = L_add(0L, (Word32)enc_PVQ_OA.index); /* MPVQ shape index set A fractional   */

        /* section B always have low indexing dynamics and is  combined into one joint single  index */
        IF (submodeLSB == 0)
        {                                                                              /* regular   AB  , shape_j == 0*/
            L_prm_idx[1] = L_deposit_l(i_gain); /* full established gain idx fwd'ed */ /*      2  possible values */
            enc_PVQ_B    = mpvq_index_fx(enc_pulsesB, N_SETB, 1);
            ASSERT(((enc_PVQ_B.size << 1)) ==
                   (sns_MPVQ_Sz[submode][1])); /*  two lowest indeces indicate all_zero B section  */

            L_tmp        = L_shl_pos((Word32)enc_PVQ_B.index, 1);            /* 2*section B  MPVQ index */
            L_prm_idx[4] = L_add(L_tmp, enc_PVQ_B.lead_sign_ind); move32(); /* add joint section B and  LS index */

            ASSERT(L_prm_idx[4] >= 0 && L_prm_idx[4] < (Word32)sns_MPVQ_Sz[submode][0]);
        }
        ELSE
        {
            L_prm_idx[1] = L_deposit_l(i_gain);
            /* MSBs of established gain idx */ /*  2 or 4   total  possible values */
            L_prm_idx[4] = L_deposit_l(-2);
        }
    }
    ELSE
    {
        /* outlier  modes   shape_j= 2(near, LSB=0) or 3(far, LSB=1)  */

        IF (submodeLSB == 0)
        {
            L_prm_idx[1] = L_deposit_l(i_gain); /* established gain idx  */ /*   4  possible values */
            enc_PVQ_OA   = mpvq_index_fx(enc_pulses_near, M,
                                       sns_Kval[submode][0]); /* o :  leading_sign_index,  index, size, k_val        */
            ASSERT(enc_PVQ_OA.size == sns_MPVQ_Sz[submode][0]);
            L_prm_idx[3] = L_add(0L, enc_PVQ_OA.index); /* MPVQ index  fractional bits */
            L_prm_idx[4] = L_deposit_l(-1);             /* no gain LSBs  */
        }
        ELSE
        {
            L_prm_idx[1] = L_deposit_l(i_gain); /* established gain idx MSBs   */ /*   all 4 or 8   possible values */
            enc_PVQ_OA   = mpvq_index_fx(enc_pulses_far, M,
                                       sns_Kval[submode][0]); /* o :  leading_sign_index,  index, size, k_val        */
            ASSERT(enc_PVQ_OA.size == sns_MPVQ_Sz[submode][0]);
            L_prm_idx[3] = L_add(0L, enc_PVQ_OA.index); /* MPVQ index  fractional bits */
            L_prm_idx[4] = L_deposit_l(-2);             /*  */
        }
        L_prm_idx[2] = L_deposit_l(enc_PVQ_OA.lead_sign_ind); /* LS shape single bit */
    }

    BASOP_sub_sub_end();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
    return;
}

static Word16 scfdec_stage2_fx(                          /* o: ber flag */
                               const Word32 *L_prm_idx,  /* set to -1 if not used */
#ifdef ENABLE_HR_MODE
                               Word32 *      st2_vector, /*o: Q14 */
#else
                               Word16 *      st2_vector, /*o: Q14 */
#endif
                               Word8 *       scratchBuffer)
{
    /*   MPVQ deindexing, gainscaling transform and transform */
#ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(
        Word16  submode;
        Word16  submodeLSB, submodeMSB;
        Word16  gValQ13;
        Word16  idxB;
        Word16  maxK;
        Word16  BER_dec;
        Word16 *dec_pulses;
        Word32 *dec_en1_vec;
        Word32 *dec_adj_glob_vec;
    );
#else
    Dyn_Mem_Deluxe_In(
        Word16  submode;
        Word16  submodeLSB, submodeMSB;
        Word16  gValQ13;
        Word16  idxB;
        Word16  maxK;
        Word16  BER_dec;
        Word16 *dec_pulses;
        Word16 *dec_en1_vec;
        Word16 *dec_adj_glob_vec;
    );
#endif

    BASOP_sub_sub_start("scfdec_stage2_fx");

    dec_pulses       = (Word16 *)scratchAlign(scratchBuffer, 0);                      /* Size = 2 * M = 32 bytes */
#ifdef ENABLE_HR_MODE
    dec_en1_vec      = (Word32 *)scratchAlign(dec_pulses, sizeof(*dec_pulses) * M);   /* Size = 2 * M = 32 bytes */
    dec_adj_glob_vec = (Word32 *)scratchAlign(dec_en1_vec, sizeof(*dec_en1_vec) * M); /* Size = 2 * M = 32 bytes */
#else
    dec_en1_vec      = (Word16 *)scratchAlign(dec_pulses, sizeof(*dec_pulses) * M);   /* Size = 2 * M = 32 bytes */
    dec_adj_glob_vec = (Word16 *)scratchAlign(dec_en1_vec, sizeof(*dec_en1_vec) * M); /* Size = 2 * M = 32 bytes */
#endif

    /* get submode   */
    submode = extract_l(L_prm_idx[0]); /* 0..3 */

    submodeLSB = s_and(submode, 0x1);
    submodeMSB = shr_pos(submode, 1);

    /* get initial adjustment gain vector for  regular, outl_near   */
    ASSERT(L_prm_idx[1] >= 0 && L_prm_idx[1] < sns_gainSz[submode]);
    gValQ13 = sns_gaintabPtr[submode][L_prm_idx[1]];
    ASSERT(gValQ13 >= 0);

    /* gain, shape indices,  incl.calls  to MPVQ deindexing */
    IF (submodeMSB != 0)
    {
        /* outlier_near or outlier_far  mode decoding */
        maxK    = sns_Kval[submode][0]; move16();
        BER_dec = pvq_dec_deidx_fx(dec_pulses, maxK, M, extract_l(L_prm_idx[2]), (UWord32)L_prm_idx[3]);
    }
    ELSE
    { /* regular mode, with potentially two shape indices  */

        maxK    = sns_Kval[submode][0]; move16();
        BER_dec = pvq_dec_deidx_fx(dec_pulses, maxK, N_SETA, extract_l(L_prm_idx[2]), (UWord32)L_prm_idx[3]);

        IF (submodeLSB == 0)
        {
            idxB = extract_l(L_prm_idx[4]); /* 0..11 */
            ASSERT(idxB >= 0 && idxB < (Word16)sns_MPVQ_Sz[0][1]);
            BER_dec |= pvq_dec_deidx_fx(&(dec_pulses[N_SETA]), sns_Kval[submode][1], N_SETB, s_and(idxB, 0x1),
                                        (UWord32)L_deposit_l(shr_pos(idxB, 1)));
            /* maxK does not need to be increased as set B is not stacked  */
        }
        ELSE
        { /* LSB gain bit already parsed */
            ASSERT(L_prm_idx[4] < 0);
            basop_memset(&dec_pulses[N_SETA], 0, (N_SETB) * sizeof(*dec_pulses));
        }
    }

    /* normalize decoded integer vector , exactly as on encoder side !!  */
    pvq_dec_en1_normQ14_fx(dec_en1_vec, dec_pulses, maxK, M);

#ifdef ENABLE_HR_MODE
    idct32_32_fx(dec_en1_vec, dec_adj_glob_vec);
#else
    idct16_fx(dec_en1_vec, dec_adj_glob_vec); /* fwd warping  in unscaled domain */
#endif

    /* scaling aligend with encoder search  */
    pvq_dec_scale_vec_fx(dec_adj_glob_vec, gValQ13, st2_vector);

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return BER_dec;
}

void processSnsQuantizeScfEncoder_fx(Word16  scf[],        /* i: input scf M */
                                     Word32 *L_prm_idx,    /* o: indeces . negative == unused */
#ifdef ENABLE_HR_MODE
                                     Word32 *scf_q,        /* o: quantized scf M */
#else
                                     Word16 *scf_q,        /* o: quantized scf M */
#endif
                                     Word8 * scratchBuffer) /* Size = 28 * M + 52 */
{
#ifdef ENABLE_HR_MODE
    Dyn_Mem_Deluxe_In(
        Word32 *target_st2; 
        Word16 *st1_idx; /* stage 1 indices */
        Word8 * buffer_processQuantize_stage2ScfEncStage2_fx;
        Counter col;
    );
#else
    Dyn_Mem_Deluxe_In(
        Word16 *target_st2;
        Word16 *st1_idx; /* stage 1 indices */
        Word8 * buffer_processQuantize_stage2ScfEncStage2_fx;
        Counter col;
    );
#endif

#ifdef ENABLE_HR_MODE
    target_st2 = (Word32 *)scratchAlign(scratchBuffer, 0);                    /* Size = 2 * M */
#else
    target_st2 = (Word16 *)scratchAlign(scratchBuffer, 0);                    /* Size = 2 * M */
#endif
    st1_idx    = (Word16 *)scratchAlign(target_st2, sizeof(*target_st2) * M); /* Size = 2 * 2 */
    buffer_processQuantize_stage2ScfEncStage2_fx = (Word8 *)scratchAlign(st1_idx, sizeof(*st1_idx) * M);
    /* Size = 26 * M + 48 */

    /* TBD needs update  */

    /* 1st stage trained VQ   */
    processQuantize_stage1ScfEncStage1_fx(scf, scf_q, &st1_idx[0], &st1_idx[1]);
    L_prm_idx[0] = L_deposit_l(st1_idx[0]);
    L_prm_idx[1] = L_deposit_l(st1_idx[1]);

/* 2nd stage PVQ-based SCF quantizer   */
    FOR (col = 0; col < M; col++)
    {
#ifdef ENABLE_HR_MODE
        target_st2[col] = L_sub(L_deposit_h(scf[col]), scf_q[col]);
#else
        target_st2[col] = sub(scf[col], scf_q[col]);
#endif
    }

    processQuantize_stage2ScfEncStage2_fx(target_st2, scf_q, &L_prm_idx[2], VQMODES26,   /* 0xF means all submodes */
                                          buffer_processQuantize_stage2ScfEncStage2_fx); /*  PVQ  in stage 2 */
    Dyn_Mem_Deluxe_Out();
}

Word16 processSnsQuantizeScfDecoder_fx(                                      /* o: BER flag */
                                       Word32 *L_prm_idx,                    /* i: indeces */
#ifdef ENABLE_HR_MODE
                                       Word32 scf_q[],
#else 
                                       Word16 scf_q[],
#endif
                                       Word8 *scratchBuffer) /* o:  M */
{
    Dyn_Mem_Deluxe_In(
        Word16 BER_flag;
    );

    /* Decode First Stage */
#ifdef ENABLE_HR_MODE
    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q27, st1SCF8_15_base5_32x8_Q27,
                                            extract_l(L_prm_idx[0]), extract_l(L_prm_idx[1]), scf_q);
#else
    processDeQuantize_stage1ScfDecStage1_fx(st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11,
                                            extract_l(L_prm_idx[0]), extract_l(L_prm_idx[1]), scf_q);
#endif

    /* Decode Second Stage */
    BER_flag = scfdec_stage2_fx(&(L_prm_idx[2]), scf_q, scratchBuffer);

    Dyn_Mem_Deluxe_Out();
    return BER_flag;
}

#ifdef CR9_C_ADD_1p25MS_LRSNS 
Word16 snsQuantScfDecLR_fx(Word32* L_sns_vq_idx_fx,
    Word32* L_scf_q_fx, /* o: Q26 */
    Word16* scf_q_fx,   /* o: Q11 */
    Word16 pitch_rx_fx, Word16  ltpf_rx_fx, Word8 * scratch)
{

 

    Dyn_Mem_Deluxe_In(
        Counter   i;
    Word32   L_mPVQ_ind_fx;         /* can be up to 17 bits  */
    Word16   shape_idx_fx, gain_idx_fx, cb_idx_fx, aux_idx_fx, LS_ind_fx;
    Word16   env_ind_fx, shift_ind_fx, sign_ind_fx, n_signs_fx;

    Word16   *Y_shape_j_fx; /* Q0 */
    Word16   *Xq_shape_j_fx;  /* Q14 */
    Word32   *L_Xq_shape_j_fx;
    const Word16 *cb_fx;
    Word16 * st1_scf_q_fx;
    Word16  CBCmeanp_ind_fx;
    const Word16  *gainTab_fx;
    Word16 gainValQ12_fx;
    Word16 BER_dec;
    Word32 L_y_en, L_norm_factor;
    Word16 norm_factorQ, y_upshift;
    );

    st1_scf_q_fx = (Word16*)scratchAlign(scratch, 0);
    Y_shape_j_fx = (Word16*)scratchAlign(st1_scf_q_fx, sizeof(*st1_scf_q_fx) * M);                      /*1*16*/
    L_Xq_shape_j_fx = (Word32*)scratchAlign(Y_shape_j_fx, sizeof(*Y_shape_j_fx) * M);                   /*1*16*/
    Xq_shape_j_fx = (Word16*)scratchAlign(L_Xq_shape_j_fx, sizeof(*L_Xq_shape_j_fx) * M);               /*1*16*/
                                                                                                       /*1*16*/

#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
    UNUSED(ltpf_rx_fx);
#endif 
    BER_dec = 0;   /* no BER detected in PVQ indeces */
    basop_memset(Y_shape_j_fx, 0, sizeof(Word16) * M);

    gainTab_fx = lrsns_vq_gainsQ12_fx[0];       move16(); /* gcc warning  init   */
    gain_idx_fx = 0;  move16();/* gcc  warning init   */
    shift_ind_fx = 0; move16();
    env_ind_fx = 0;   move16();

    /* get indices from dec_entropy_dx */
    cb_idx_fx = extract_l(L_sns_vq_idx_fx[0]);  move16();  /* stage1 idx */
    aux_idx_fx = extract_l(L_sns_vq_idx_fx[1]);  move16();  /* mPVQ LS-bit  or the FESS s0 bit */
    shape_idx_fx = extract_l(L_sns_vq_idx_fx[2]);  move16();  /* analysis order shape idx  -9,-10    0,1,  2,3,4,5 */
    gain_idx_fx = extract_l(L_sns_vq_idx_fx[3]);  move16();  /* stage 2 gain index */

    /* Stage1 cand   */
    IF(sub(shape_idx_fx, -9) == 0)
    {
        /* minimal  2*16 SNS codebook,  no DC  */
        cb_fx = lrsns_cbA_fx;    /* ptr init */
        basop_memcpy(scf_q_fx, &(cb_fx[shl_pos(cb_idx_fx, 4)]), sizeof(Word16) * M);  /* cb_idx_fx * M */
    }
    ELSE IF(sub(shape_idx_fx, -10) == 0)
    {    /* 0..339 */ /* stage 1B or 1C only, transmitted  in 9+1= 10 bits */
        IF(sub(cb_idx_fx, 170) < 0)
        {  /* Stage 1B   */
            snslr_st1B_vector_dec_fx(cb_idx_fx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11,
                lrsns_st1B_merged170orderSortedSegmCum_fx, lrsns_st1B_merged170orderSort12bitIdx_fx, scf_q_fx);

            snslr_remove_st1_DC_fQ11_fx(scf_q_fx, M);
        }
        ELSE
        {
            cb_idx_fx = sub(cb_idx_fx, 170);
            /*  Stage 1C  harmonic outlier  CB  with   pitch dependent mean  vector   */
            /* Q11 values , so that BASOP and float  is always  BE in synthesis*/
             ASSERT(cb_idx_fx >= 0 && cb_idx_fx < 170);
       
             snslr_st1C_vector_dec_fx(cb_idx_fx, lrsns_st1C_Both_Word8_fx, lrsns_st1C_Both_scaleQ4_7p4bits_fx[1],
                 lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1], M, 170, scf_q_fx);
             /* add harmonic mean , based on pitch_info availability */

             pitch_rx_fx = extract_l(L_sns_vq_idx_fx[3]);  /* LTP active flag directly from dec_entropy */    
#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
             ltpf_rx_fx = 0;                               /* CB_C has no dependency on LTPF active flag */
#else 
             ltpf_rx_fx = extract_l(L_sns_vq_idx_fx[4]);   /* LTPF active flag directly from dec_entropy */
#endif 
             CBCmeanp_ind_fx = pitch_rx_fx; move16(); /* 0 or 1 */
#ifndef  LRSNS_CBC_NO_LTPF_DEPENDENCY
             test(); test(); test();
             if (pitch_rx_fx != 0 && ltpf_rx_fx != 0)
             {
                 CBCmeanp_ind_fx = add(CBCmeanp_ind_fx, 1);  /* high corr ltpf_rx is also active */
             }
#endif 
             const Word16 *mean_cb_fx = lrsns_st1CTrainedMapMeans_fx[CBCmeanp_ind_fx];  /* point to pitch dependent mean */
             for (i = 0; i < M; i++)
             {
                 scf_q_fx[i] = add(scf_q_fx[i], mean_cb_fx[i]);
             }
             /* remove_DC()  call is not required for section  C  */
        }
    }
    ELSE
    {     /* 0..169 */   /* st1B*   used with a stage 2  shape submode  */
        ASSERT(shape_idx_fx >= 0 && shape_idx_fx <= 5);
        snslr_st1B_vector_dec_fx(cb_idx_fx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, lrsns_st1B_merged170orderSortedSegmCum_fx, lrsns_st1B_merged170orderSort12bitIdx_fx, scf_q_fx);
        snslr_remove_st1_DC_fQ11_fx(scf_q_fx, M); /* inplace */     /* DC needs removal for st1 part B */
    }

    basop_memcpy(st1_scf_q_fx, scf_q_fx, sizeof(Word16) * M);  /* keep track of stage1 (A,B or C) contribution */

    IF(shape_idx_fx >= 0) /* stage 2 shapes 0,1,  2,3,4,5   ( negative idx used for stage1 only   */
    {
        /* stage 2 SNS VQ decoding */
        /* Decode shape_j */
        Y_shape_j_fx[0] = 0;   /* no DCT-II DC-coeff decoded */

        SWITCH(shape_idx_fx)
        {
        case 0:  /* splitLF 29 bits total  */
            LS_ind_fx = aux_idx_fx; move16();
            L_mPVQ_ind_fx = L_sns_vq_idx_fx[4];  move32(); /* mPVQ(5,6) or mPVQ(5,8) */

            test();
            IF(L_sns_vq_idx_fx[5] >= 0)
            {
                /* CFLT: MPVQdeenum_fx(5, 6, LS_ind_fx, L_mPVQ_ind_fx, &Y_shape_j_fx[1]); */
                BER_dec = pvq_dec_deidx_fx(&Y_shape_j_fx[1], PULSES_SPLIT_A_LR, NA_LR, LS_ind_fx, (UWord32)L_mPVQ_ind_fx);

                LS_ind_fx = extract_l(L_and(L_sns_vq_idx_fx[5], 0x1));
                L_mPVQ_ind_fx = L_shr_pos(L_sns_vq_idx_fx[5], 1);

                /* CFLT:    MPVQdeenum_fx(8, 2, LS_ind_fx, L_mPVQ_ind_fx, &Y_shape_j_fx[1 + 5]); */
                BER_dec = s_or(BER_dec, pvq_dec_deidx_fx(&Y_shape_j_fx[1 + NA_LR], PULSES_SPLIT_B_LR, NB_LR, LS_ind_fx, (UWord32)L_mPVQ_ind_fx));
            }
            ELSE
            {
                /* CFLT:  MPVQdeenum_fx(5, 8, LS_ind_fx, L_mPVQ_ind_fx, &Y_shape_j_fx[1]); */
                 BER_dec = pvq_dec_deidx_fx(&Y_shape_j_fx[1], (PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR), NA_LR, LS_ind_fx, (UWord32)L_mPVQ_ind_fx);
            }
            gainTab_fx = lrsns_vq_gainsQ12_fx[0];   move32(); /* 4 levels in 2 bits  */
            break;
        case 1:  /*  full (N=15,K=5)  30 bits total */
            LS_ind_fx = aux_idx_fx;  move16();
            L_mPVQ_ind_fx = L_sns_vq_idx_fx[4]; move32();

            /* CFLT: MPVQdeenum_fx(15, 5, LS_ind_fx, L_mPVQ_ind_fx, &Y_shape_j_fx[1]); */
            BER_dec = pvq_dec_deidx_fx(&Y_shape_j_fx[1], PULSES_FULL_LR, NFULL_LR, LS_ind_fx, (UWord32)L_mPVQ_ind_fx);
            gainTab_fx = lrsns_vq_gainsQ12_fx[1]; move16(); /* 8 levels in 3 bits  */

            break;
        case 2: /* fix env 0  ,  init_bell      12 signs  */
        case 3: /* fix env 1  ,  decay 12-->6   12 signs */
        case 4: /* fix env 2  ,  start bell     12 signs */
        case 5: /* fix env 3  ,  early bell     10 signs  */
            LS_ind_fx = aux_idx_fx; /* s0 */ move16();
            env_ind_fx = extract_l(L_sns_vq_idx_fx[4]); move16();
            ASSERT(env_ind_fx == (shape_idx_fx - 2));
            Word16  sign_mask_fx = 0x07ff;   move16();      /* mask for 11 remaining signs */
            n_signs_fx = NSIGNS_FIX_012;  move16(); /* number of signs including s0 */
            if (sub(env_ind_fx, 3) == 0)
            {
                n_signs_fx = NSIGNS_FIX_3;  move16();
            }
            sign_mask_fx = shr_pos_pos(sign_mask_fx, sub(NSIGNS_FIX_012, n_signs_fx));

            shift_ind_fx = extract_l(L_shr_pos(L_sns_vq_idx_fx[5], sub(n_signs_fx, 1)));
            sign_ind_fx = s_and(extract_l(L_sns_vq_idx_fx[5]), sign_mask_fx);

            /* put sign s0 , right next to s1 , to make the subsequent sign decoding loop straight fwd */
            sign_ind_fx = add(shl_pos(LS_ind_fx, sub(n_signs_fx, 1)), sign_ind_fx);   /* s0 put as MSB at 12th position   */

            /*FixEnvShiftedSigns deenumeration */
            /* part of 30b total,    4xenv,4xshifts, 10 or12 signs , spread over 15 positions,*/
            FESSdeenum_fx(NFULL_LR, N_CANDS_FIX_LR, N_SHIFT_FIX, n_signs_fx, env_ind_fx, shift_ind_fx, sign_ind_fx, &Y_shape_j_fx[1]);

            gainTab_fx = lrsns_vq_gainsQ12_fx[2];  move16(); /* 8 levels in 3 bits  */
           /* fix_envshift_nb_fx = env_ind_fx * 4 + shift_ind_fx; */ /* index for fast normalization lookup */
            break;
        default:
            ASSERT(0 && " LRSNS stage2 demux bad shape shape_idx_fx received from dec_entropy_fx ");
            break;
        }

        gainValQ12_fx = gainTab_fx[gain_idx_fx];


        /* enc/dec common Word32Q30 and Word16Q14  Unit energy normalization  of the received shapes mPVQ or FESS */
        IF(sub(shape_idx_fx, 2) < 0)
        { /* 0,1,2 ::  PVQ(N,K) shape  with up to  8  unit pulses */
            L_y_en = L_deposit_l(0);
            FOR(i = 1; i < M; i++)
            {
                L_y_en = L_mac0(L_y_en, Y_shape_j_fx[i], Y_shape_j_fx[i]);
            }
            L_norm_factor = isqrt_Q31tab[L_y_en];  move32(); /* Q31 inv_sqrt value */
            norm_factorQ = lrsns_norm_factorQ_L[shape_idx_fx]; move16();
            y_upshift = lrsns_y_up_bits[shape_idx_fx]; move16();
        }
        ELSE
        {  /* 2,3,4,5 ::  FESS shapes with up to    unit pulses */
           /* no energy calc   required, only 16 different energies are possible */
           L_norm_factor = L_lrsns_fixenv_enNormQ35[add(shl_pos(env_ind_fx , 2), shift_ind_fx)];  move32();  /* Q19+16=Q35 */
           norm_factorQ = lrsns_norm_factorQ_L[2]; move16();
           y_upshift = lrsns_y_up_bits[2];      move16();
        }

            /* reuse common encoder/decoder  unit energy normalization routine */
        ASSERT(Y_shape_j_fx[0] == 0);
        pvq_fess_dec_en1_normQ30andQ14_fx(Y_shape_j_fx, y_upshift, L_norm_factor, norm_factorQ, M, L_Xq_shape_j_fx, Xq_shape_j_fx);
        ASSERT(L_Xq_shape_j_fx[0] == 0 && Xq_shape_j_fx[0] == 0);

        /* Reconstruction of the quantized SNS scale factors */
#ifdef ENABLE_HR_MODE 
        L_Xq_shape_j_fx[0] = 0;           move32();   /* enforce no DC */
        idct32_32_fx(L_Xq_shape_j_fx, L_Xq_shape_j_fx); /* inplace idct */

        /* move stage1 W16Q11 to W32Q27 */
        FOR(i = 0; i < M; i++)
        {
            L_scf_q_fx[i] = L_deposit_h(st1_scf_q_fx[i]);  move32();/* stage1 W16Q11 to  W32Q27  */
        }
        lrsns_pvq_dec_scale_W32vec_fx(L_Xq_shape_j_fx, gainValQ12_fx, L_scf_q_fx  /* W32Q27 in, W32Q26 out */, scf_q_fx /* out Q11*/);

#else
        /* DISABLE_HR */
        Xq_shape_j_fx[0] = 0; move16(); /* no DC */
        idct16_fx(Xq_shape_j_fx, Xq_shape_j_fx);  /* inplace idct */  /* fwd in unscaled unit energy domain */

        /* move stage1B into a W16 vector for  accumulation in Q11    */
        basop_memcpy(scf_q_fx, st1_scf_q_fx, sizeof(*scf_q_fx)*M);

        /* scf_q_fx contains stage1 in Q11 */
        lrsns_pvq_dec_scale_W16vec_fx(Xq_shape_j_fx/* Q14 */, gainValQ12_fx, scf_q_fx /* Q11in, Q11out */);

        /* scf_q_fx is the final result incl. stage2  in Q11 */
#endif /* DISABLE_HR */

    }
    ELSE
    {   /* -9, -10 */
        /* LRSNS stage 1   A,B,C  */
        basop_memcpy(scf_q_fx, st1_scf_q_fx, sizeof(Word16) * M); /* Q11*/

        FOR(i = 0; i < M; i++)
        {
            L_scf_q_fx[i] = L_shr_pos(L_deposit_h(scf_q_fx[i]), 1); /* 11+16-1 => Q26 */
        }
    }
    Dyn_Mem_Deluxe_Out();

    return BER_dec;
}
 
/* split out of LRSNS stage 1 functionality */
Word32 snsQuantScfEncLRSt1ABC_fx(Word16* env, Word32* L_index, Word32 *L_min_mse_saveBCA_ptr_fx,
    Word16* ind_saveB_ptr, Word16* st1_vectors,
    Word16 pitch_rx, Word16 ltpf_rx, Word8 * scratch)
{
 
    Dyn_Mem_Deluxe_In(
        Counter i;
    Word16 *st1_vector_fx, *st1_vectorA_fx, *st1_vectorB_fx, *st1_vectorC_fx, *target_fx;
    Word32 L_min_mse_saveA_fx, L_min_mse_saveB_fx, L_min_mse_saveC_fx, L_min_mse_saveBC_fx;
    Word16  ind_saveA_fx, ind_saveC_fx;
    Word16 meanC_ind;
    Word16 stage1_mode; /*0=A, 1=B, 2=C. -1==fail*/
    Word32 L_min_mse_saveB_fxlike_fx;
    Word16 ind_saveB_fxlike_fx;
    Word16 dc_fx;
    Word32 L_tmp;
    Word16 *st1_vectorB_idx_fx;
    );

    stage1_mode = -1;   /* output */

    st1_vectorB_idx_fx = scratchAlign(scratch, 0);
    target_fx = scratchAlign(scratch, sizeof(*st1_vectorB_idx_fx) * M);
   
#ifdef  LRSNS_CBC_NO_LTPF_DEPENDENCY
    UNUSED(ltpf_rx);
#endif 
    st1_vectorA_fx = &(st1_vectors[0]);
    st1_vectorB_fx = &(st1_vectors[1 * M]);
    st1_vectorC_fx = &(st1_vectors[2 * M]);
    st1_vector_fx = &(st1_vectors[3 * M]);  /*selected winner */

    BASOP_sub_sub_start("snsQuantScfEncLRSt1ABC_fx");

    /* snslr stage1  B(170) and C(170), A(2)  evaluation */
    /*
       segm   idx9b

      ------+--------+--------+-------
        B   | 0--169 | aux==0 ,  170 entries
      ------+--------+--------
        C   | 170-339| aux==1  , 170 entries,
      ------+--------+--------+
        *   | 341-509| aux==2   , 170 entries (later decided to be B  + aux bit)
      ------+--------+---------
        A   | 510,511| 2 entries,  9 bit total
      ------+--------+--------
     */
    BASOP_sub_sub_start("MSEsearchCbA_fx");


    {   /* stage 1 section A(2),  a very small 2xM entry  cb */
        /* two crude vector MSE_calculations  needed */

        basop_memcpy(target_fx, env, sizeof(Word16)*M);

        ind_saveA_fx = MSEsearchGeneric_fx(target_fx, lrsns_cbA_fx, M, 2, &L_min_mse_saveA_fx);

        basop_memcpy(st1_vectorA_fx, &(lrsns_cbA_fx[ind_saveA_fx*M]), sizeof(Word16)*M);
    }
    BASOP_sub_sub_end();


    BASOP_sub_sub_start("MSEsearchCbB_fx");
    {  /* stage1 section B(170) MSE analysis */
        basop_memcpy(target_fx, env, sizeof(Word16)*M);
        *ind_saveB_ptr = -1;  move16();

        ind_saveB_fxlike_fx = MSEsearchCbBIdxMap_fx(target_fx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11,
            lrsns_st1B_merged170orderSortedSegmCum_fx, lrsns_st1B_merged170orderSort12bitIdx_fx, M, 170, &L_min_mse_saveB_fxlike_fx); /*st1B LF,HF idx lookup search 170 ,170 Word16s ,  0.34kB ROM */

        snslr_st1B_vector_dec_fx(ind_saveB_fxlike_fx, st1SCF0_7_base5_32x8_Q11, st1SCF8_15_base5_32x8_Q11, lrsns_st1B_merged170orderSortedSegmCum_fx, lrsns_st1B_merged170orderSort12bitIdx_fx,
            st1_vectorB_idx_fx);

        *ind_saveB_ptr = ind_saveB_fxlike_fx; move16();
        L_min_mse_saveB_fx = L_min_mse_saveB_fxlike_fx; move32();
        basop_memcpy(st1_vectorB_fx, st1_vectorB_idx_fx, sizeof(Word16)*M); /*use low ROM version */

        {
            /* remove the DC that can remain in the LF,HF index  stored stage cbB structure  */
            /* a very slight off-line decrease in perf 0.001 dB in AveSD when searching with DC above,
               but it allows much better stage1 ROM-reuse performance */

            dc_fx = snslr_remove_st1_DC_fQ11_fx(st1_vectorB_fx, M); /* in-place removal of dc in st1_vectorB */

            ASSERT(L_min_mse_saveB_fx >= 0L);

            dc_fx = shl_pos(dc_fx, 2);    /* mult dc_fx with  sqrt(M) */
            L_tmp = L_msu0(L_min_mse_saveB_fx, dc_fx, dc_fx);

            /*  force the MSE B  to be zero  or larger*/
            /* may  happen due to imperfect(truncation) DC removal */
            L_min_mse_saveB_fx = L_max(L_tmp, 0L);
        }
    } /* end of stage 1 section B , search */
    BASOP_sub_sub_end();

    BASOP_sub_sub_start("MSEsearchCbC_fx");
    meanC_ind = pitch_rx; move16();
#ifndef  LRSNS_CBC_NO_LTPF_DEPENDENCY 
    if (ltpf_rx != 0)
    {
        meanC_ind = add(meanC_ind, 1);
    }
#else 
    /* CbC only dependent on LTP transmission on/off */
#endif 
    Word32 L_mse_tmpQ22 = L_min(Mpy_32_16(L_min_mse_saveA_fx, SNSLR_A_CNST_WEIGHT), L_min_mse_saveB_fx);

    /* apply skipC threshold of 3.97 in Q22 */
    IF(L_sub(L_mse_tmpQ22, 16651387L) < 0)  /*   to skip C search if SD is already low enough < 1.5dB (msesum=3.97))  to save average WMOPS */
    {
        /* search stage1C */ /*  made adaptive on stage1 BA  MSE sofar */
        L_min_mse_saveC_fx = (23404216L << 1); move32();/* disable C selection in consecutive  logic */
        ind_saveC_fx = 0;  move16();/* a valid index that will not be used */
    }
    ELSE
    {
        /* search another set (pitch dependent section C) of 170 mean residual vectors  */
        /* pitch_info rx used to selects a   mean   and then a trained  residual 1x170W8 based harmonic outlier table */
        /* float means were based on  Word16 S16Q11 values,  so that BASOP and float  always may become BE in synthesis*/
        FOR(i = 0; i < M; i++) {
            target_fx[i] = sub(env[i], lrsns_st1CTrainedMapMeans_fx[meanC_ind][i]);
        }


        Word32 ind_saveC_ScaledW8_fx;
        Word32 L_min_mse_saveC_ScaledW8_fx;

        ind_saveC_ScaledW8_fx = MSEsearchGenericScaledW8_fx(target_fx, lrsns_st1C_Both_Word8_fx,
            lrsns_st1C_Both_scaleQ4_7p4bits_fx[1], lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1],
            M, 170, &(L_min_mse_saveC_ScaledW8_fx));

        snslr_st1C_vector_dec_fx(ind_saveC_ScaledW8_fx, lrsns_st1C_Both_Word8_fx,  lrsns_st1C_Both_scaleQ4_7p4bits_fx[1], lrsns_st1C_Both_inv_scaleQ15_7p4bits_fx[1],
            M, 170,
            st1_vectorC_fx
        );

        ind_saveC_fx = ind_saveC_ScaledW8_fx;
        L_min_mse_saveC_fx = L_min_mse_saveC_ScaledW8_fx;

        FOR(i = 0; i < M; i++)
        {
            st1_vectorC_fx[i] = add(st1_vectorC_fx[i], lrsns_st1CTrainedMapMeans_fx[meanC_ind][i]); /*  Q11 means*/
        }

#ifdef FIX_BASOP_ENC_LRSNS_CBC_MSE
        /* CB_C MSE search loop in BASOP is using lower resolution, than and the best possible calculation  */
        /* recalc CB_C MSE with Q22 precision for a best possible CB decision  */
        MSEsearchGeneric_fx(env, st1_vectorC_fx, M, 1, &L_min_mse_saveC_fx); /* reuse of search to get MSE for 1 vector */
#endif 

    }
    BASOP_sub_sub_end();


    /* BC stage1    comparison   */
    /* initially assume B as stage 1 winner */
    L_min_mse_saveBC_fx = L_min_mse_saveB_fx; move32();
    L_index[0] = L_deposit_l(*ind_saveB_ptr);  /* 0...169   */

    basop_memcpy(st1_vector_fx, st1_vectorB_fx, sizeof(Word16)*M);   /* stage 1 segmentB result without DC  copied as base for st2 */

    {
        IF(L_sub(L_min_mse_saveC_fx, L_min_mse_saveBC_fx) < 0)
        {  /* C better than B */
            L_min_mse_saveBC_fx = L_min_mse_saveC_fx;  move32();
            L_index[0] = L_add(170L, L_deposit_l(ind_saveC_fx));  /* [2x 170] (9+3)b        [ 0..169], [170-339] ,  ">=170"  is a signal to multiplexor  */
            basop_memcpy(st1_vector_fx, st1_vectorC_fx, sizeof(Word16)*M);
        }
        L_index[1] = -10;  move32();
        ASSERT(L_min_mse_saveBC_fx >= 0.0);

        /* (9(A)<->10(BC) bit weighted comparison   */
        *L_min_mse_saveBCA_ptr_fx = L_min_mse_saveBC_fx;  move32();
        IF(L_sub(Mpy_32_16(L_min_mse_saveA_fx, SNSLR_A_CNST_WEIGHT), L_min_mse_saveBC_fx) < 0)  /* a minor favouring of the 9b vector results  sqrt(0.875)  =>  approx  0.6dB level domain  */
        {
            L_index[0] = L_add(510L, L_deposit_l(ind_saveA_fx));  /*   only [510,511]  possible */
            L_index[1] = -9;

            basop_memcpy(st1_vectorA_fx, &(lrsns_cbA_fx[ind_saveA_fx*M]), sizeof(Word16)*M);
            basop_memcpy(st1_vector_fx, st1_vectorA_fx, sizeof(Word16)*M);

            *L_min_mse_saveBCA_ptr_fx = Mpy_32_16(L_min_mse_saveA_fx, SNSLR_A_CNST_WEIGHT); /* Q22 */
        }
    }

    /* index0_saveBCA = index[0];*/    /* 0...511 */
    /* index1_saveBCA = index[1];*/    /* -9 or -10 */

    stage1_mode = 2; move16(); /* cbC */
    if (L_sub(L_index[0], 510L) >= 0) {
        stage1_mode = 0; move16(); /* cbA */
    }
    if (L_sub(L_index[0], 169L) <= 0) {
        stage1_mode = 1;  move16();/* cbB */
    }

    BASOP_sub_sub_end();

    Dyn_Mem_Deluxe_Out();

    return stage1_mode; /* return best mode */
}

/*  top level  lrsns BASOP code  calling both st1 and st2 */
Word16  snsQuantScfEncLR_fx(   /* o:  bits spent on LRSNS-VQ  envelope */
    Word16  scf_fx[],            /* i: input scf M W16Q11  */
    Word32 *L_index_fx,          /* o: SNS  indeces . */
#  ifdef ENABLE_HR_MODE
    Word32 *L_scf_q_fx,          /* o: quantizefl_env scf M   ? W32Q11 or W32Q27 */
#else
    Word16 *scf_q_fx,            /* o: quantizefl_env scf M   W16Q11 */
#endif
    Word16  pitch_rx_fx,             /*i:  0 or 1 */
    Word16  ltpf_rx_fx,               /*i:  o or 1 */
    Word8 * scratch)
{

    Dyn_Mem_Deluxe_In(
        Counter    col;
 
    Word16 *st1_vectorA_fx, *st1_vectorB_fx, *st1_vectorC_fx, *st1_vector_fx;
    Word16 ind_saveB_fx;
    Word16 st1_mode_fx;
    Word16 envelope_bits_fx; /* output */
    Word32 L_min_mse_saveBCA_Q22_fx;

    Word16 gain_idx_fx; /* gain  index 0..3, 0..7,  0..7  */
    Word16 s_idx_fx; /* shape index 0 =split, 1 = full,  2= fix*/
    Word16 shape_idx_fx; /* expanded shape index 0 ..5 */
    Word32 L_mse_st1B_st2_Q22_fx;
    Word16  *y_split_fx, *y_full_fx, *y_fix_fx;
    Word16  gainValQ12_fx;   /* output from PVQ search */
    Word16  fixShapeNb;
    Word16  fixShiftIdx;
   
    /* scratch ptrs */
    Word16 *st1_vectors_fx;
    Word16 *target_st2_fx; 
    Word8  *scratch_ABC_fx;
    Word32 *L_target_st2_fx;    /*  req. for dct32 use */
    Word16 *pvq_target_fx;
    Word32 *L_pvq_target_fx;
    Word16  *y_Q0;
    Word16  *y_normQ14_fx;
    Word32  *L_y_normQ30_fx;
    Word8  *scratch_pvq_fess_fx;
    );


#ifdef   ENABLE_HR_MODE
    Word16 scf_q_fx[M];            /*   W16Q11  always in use also for HR_MODE */
    UNUSED(scf_q_fx);
#endif 
    UNUSED(st1_mode_fx);
    UNUSED(st1_vectorA_fx);
    UNUSED(st1_vectorC_fx);

    st1_mode_fx = -1;          move16();
    envelope_bits_fx = -1;     move16(); /* output information  */
    shape_idx_fx = 0;         move16();

    st1_vectors_fx    = (Word16*)scratchAlign(scratch, 0);
    scratch_ABC_fx    = (Word8*)scratchAlign(st1_vectors_fx, sizeof(*st1_vector_fx) * M * 4);
    target_st2_fx     = (Word16*)scratch_ABC_fx;
    L_target_st2_fx   = (Word32*)scratchAlign(target_st2_fx, sizeof(*target_st2_fx) * M);

    pvq_target_fx = (Word16*)scratchAlign(L_target_st2_fx, sizeof(*L_target_st2_fx) * M);
    L_pvq_target_fx = (Word32*)scratchAlign(pvq_target_fx, sizeof(*pvq_target_fx) * M);
    y_Q0 = (Word16*)scratchAlign(L_pvq_target_fx, sizeof(*L_pvq_target_fx) * M);
    y_normQ14_fx = (Word16*)scratchAlign(y_Q0, sizeof(*y_Q0) *SNSLR_MAX_PVQ_SEARCH_CAND*M);           
    L_y_normQ30_fx = (Word32*)scratchAlign(y_normQ14_fx, sizeof(*y_normQ14_fx) *SNSLR_MAX_PVQ_SEARCH_CAND*M);                               
    scratch_pvq_fess_fx = (Word8*)scratchAlign(L_y_normQ30_fx, sizeof(*L_y_normQ30_fx) * SNSLR_MAX_PVQ_SEARCH_CAND * M );  

    y_split_fx = &(y_Q0[0 * M]); /* ptr init */
    y_full_fx = &(y_Q0[1 * M]);  /* ptr init */
    y_fix_fx = &(y_Q0[2 * M]);   /* ptr init */

    st1_vectorA_fx = &(st1_vectors_fx[0]);    /* ptr init */
    st1_vectorB_fx = &(st1_vectors_fx[1 * M]);
    st1_vectorC_fx = &(st1_vectors_fx[2 * M]);
    st1_vector_fx = &(st1_vectors_fx[3 * M]);  /* best vector */

    FOR(col = 0; col < SCF_MAX_PARAM; col++)
    {
        L_index_fx[col] = L_sub(-32000, col);  move32();   /* safety init parameters to be fwd'ed to LRSNS VQ multiplexor */
    }

    /*  Stage 1 Cb's A,B,C  */
    ind_saveB_fx = -1;
    L_min_mse_saveBCA_Q22_fx = M * 32 * 32;  move32();

    st1_mode_fx = snsQuantScfEncLRSt1ABC_fx(
        scf_fx, L_index_fx, &L_min_mse_saveBCA_Q22_fx, &ind_saveB_fx, st1_vectors_fx,
        pitch_rx_fx, ltpf_rx_fx, scratch_ABC_fx);

    Word32 L_mse_lim_smooth_Q22_fx = 22691185L; /* round(5.41*pow(2.0, 22.0))*//* 1.75 dB */

    /* mse_st1B_st2_fl = 2.0* min_mse_saveBCA + 1.0;*/   /*  safety indicate that st1B+st2 is not used by setting a higher MSE than st1BCA  */
    L_mse_st1B_st2_Q22_fx = L_add(L_shl_pos(L_min_mse_saveBCA_Q22_fx, 1), (1L << 22)); /* set st1B+st2 to a safety bad MSE value */


    IF((L_sub(L_min_mse_saveBCA_Q22_fx, L_mse_lim_smooth_Q22_fx) > 0))
    {   /* stage 2 analysis */
        /* prepare stage2 W16 target */
        FOR(col = 0; col < M; col++)
        {
            target_st2_fx[col] = sub(scf_fx[col], st1_vectorB_fx[col]); /* Q11 */
        }

        /* both ENABLE_HR and DISABLE_HR runs the same analysis DCT-II(16)  */
        /* for analysis use 16 bit i/o Word16 constants,  but Word32 internal states in DCT-II(M=16) */
        basop_memcpy(pvq_target_fx, target_st2_fx, M * sizeof(*pvq_target_fx));

        dct16_W32int_fx(pvq_target_fx, pvq_target_fx);  /*  Q11 to Q11,   enc-side  analysis with  W16 i/o internally W32 precision  */

        pvq_target_fx[0] = 0; move16();

        /* PVQ FESS search, include norm , shape search and Q gain search   */
        /* NB IDCT-II  __not__ run in stage2  shape and gain search */
        pvq_fess_enc_search_fx(
            pvq_target_fx,
            y_Q0,
            y_normQ14_fx,           /*  normally  calculated for DISABLE_HR_MODE */
            L_y_normQ30_fx,         /*  calculated for both DISABLE_HR_MODE and ENABLE_HR_MODE */
            &s_idx_fx/*[0,1,2]*/,
            &gain_idx_fx/* [0...7]*/,
            &gainValQ12_fx,
            &fixShapeNb, /* 0,1,2,3, [-1] */
            &fixShiftIdx, /* 0,1,2,3, [-1] */
            &L_mse_st1B_st2_Q22_fx,
            scratch_pvq_fess_fx);


        /* update to shape idx to one of 0...5 ] */
        shape_idx_fx = s_idx_fx; move16();
        if (sub(shape_idx_fx, 2) == 0)
        {
            shape_idx_fx = add(shape_idx_fx, fixShapeNb); /* [0..5] */
        }
    }
    ELSE
    {   /* indicate complete skipping of stage2 search,  stage1 is good enough  */
        L_mse_st1B_st2_Q22_fx = INT32_MAX;   move32();
        /*copy of best stage 1 to proper output , done below */
    }

    IF(L_sub(L_mse_st1B_st2_Q22_fx, L_min_mse_saveBCA_Q22_fx) > 0)
    {
        /* skip IDCT,  as stage2 will not be used */
        //L_mse_st1B_st2_Q22_fx = L_mse_st1B_st2_Q22_fx; /* dummy op */
    }
    ELSE
    {
#ifdef ENABLE_HR_MODE  
        /* apply IDCT-II  on the Q30 unit energy normalized vector */
        Word32 L_tmp_vec_fx[M];
        Word32* L_y_norm_fx = &(L_y_normQ30_fx[s_idx_fx*M]);  move32();

        ASSERT(L_y_norm_fx[0] == 0);
        /*  32x32 bit dec IDCT-II analysis, with W32 internal constants  */
        idct32_32_fx(L_y_norm_fx, L_tmp_vec_fx);    /* currently Q30 to Q30 , can also be inplace  */
        /* idct32_32_fx::  162.712 dBSegSNR, minSNR = 157.979, WMOPS 0.26 */
#else
        Word16 tmp_vec_fx[M];
    /* apply IDCT-II  on the Q14 unit energy normalized vector */

    Word16* y_norm_fx = &(y_normQ14_fx[s_idx_fx*M]);  move32();
    ASSERT(y_norm_fx[0] == 0);
    idct16_fx(y_norm_fx, tmp_vec_fx); /*  idct16_fx::  segsnr 73.7 /minsnr 67.65  ,  0.118 WMOPS  using mult_r)  */
#endif /* HR_MODE*/

#ifdef ENABLE_HR_MODE 
        /* move stage1 W16Q11 to W32Q27 */
        FOR(col = 0; col < M; col++)
        {
            L_scf_q_fx[col] = L_deposit_h(st1_vectorB_fx[col]); /* W16Q11 to W32Q27  */
        }
        lrsns_pvq_dec_scale_W32vec_fx(L_tmp_vec_fx, gainValQ12_fx, L_scf_q_fx /* W32Q27 in, W32Q26 out */, scf_q_fx /* W16Q11 out */);
#else
        /* non HR */
       /* move stage1B into a W16 vector for  accumulation in Q11    */
        basop_memcpy(scf_q_fx, st1_vectorB_fx, sizeof(*scf_q_fx)*M);

        /* scf_q_fx contains stage1 in Q11*/
        lrsns_pvq_dec_scale_W16vec_fx(tmp_vec_fx/*Q14*/, gainValQ12_fx, scf_q_fx /* Q11 */); /* scf_q_fx is the final result incl. stage2  in Q11 */

#endif /* ENABLE_HR_MODE */
    }


        /* post-evaluate if one of (st1B, st1C, st1A) was actually better than st1B+stage2 */
   
    IF(L_sub(L_mse_st1B_st2_Q22_fx, L_min_mse_saveBCA_Q22_fx) <= 0)
    {
        /*   use stage1B + st2  at 29b or 30b bits total cost */
        L_index_fx[0] = L_deposit_l(ind_saveB_fx); move32();
        L_index_fx[1] = L_deposit_l(2930);         move32();  /* later stage2    aux  value  LS_splitLF or LS_full or s0,   put here  as a 0 or 1 */
        L_index_fx[2] = L_deposit_l(shape_idx_fx); move32();   /* 0=splitLF, 1=full,  ( 2=fixEnv0, 3=fixEnv1, 4: fixEnv2, 5: fixEnv3 )  */
        L_index_fx[3] = L_deposit_l(gain_idx_fx);  move32();   /*  gain idx  with a shape dependent number of levels  (4  or 8  levels ) */

        basop_memcpy(st1_vector_fx, st1_vectorB_fx, sizeof(*st1_vectorB_fx)*M);
        /* final  result  st1 in combination with stage 2, kept now for verification at decoder */

        envelope_bits_fx = 29;   move16(); /* 'LR_splitLF' bitrate */
        test();
        if (shape_idx_fx > 0) {
            envelope_bits_fx = add(envelope_bits_fx, 1); /*30 'LR_full/LR_fixenv' */
        }

        {
            /* DBG check values  */
            ASSERT(shape_idx_fx >= 0);
            ASSERT(envelope_bits_fx >= 29);
            ASSERT(L_index_fx[0] <= 170); /*only B allowed */
            ASSERT(L_index_fx[1] >= 0);

            ASSERT(gain_idx_fx >= 0); /*gain index*/
            ASSERT(gainValQ12_fx > 0);  /* gain value */
        }
    }
    ELSE
    {   /*stick to stage1(best of BCA)   at  9 or 10  bits */
        ASSERT(L_index_fx[1] < 0 && L_index_fx[0] >= 0 && L_index_fx[0] < 512);
        envelope_bits_fx = ((L_index_fx[0] >= 510) ? 9 : 10);
        shape_idx_fx = -envelope_bits_fx; /* signal an invalid stage2 shape number to enc-entropy */
#  ifdef ENABLE_HR_MODE
            FOR(col = 0; col < M; col++)
            {
                L_scf_q_fx[col] = L_shl_pos(st1_vector_fx[col], 16 - 1); move32(); /* W16Q11 to W32Q26  */
            }
#  else
            basop_memcpy(scf_q_fx, st1_vector_fx, M * sizeof(*scf_q_fx)); /* output */
#  endif 

            gainValQ12_fx = 0; move16();

            gain_idx_fx = -1;  move16(); /* L_index sentinel */
            L_index_fx[2] = shape_idx_fx; move32();
    }

    /******************************************************************/
    /*  signal to enc_entropy_fx for LRSNS semi-fractional multiplexing  */
    /******************************************************************/
    /* integer multiplexing   29/30 bit modes into intermediate  unmuxed integer indeces  0...7  */
    /* a bit of fractional multiplexing  for the  L_index_fx  0...7  is done later,  in function enc_entropy_fx()  */
    test();
    IF(shape_idx_fx >= 0)
    {  /* stage 2 multiplexing manipulations */
        Word16  fix_end_sign_fx;

        PvqEntry_fx enc_PVQ_A, enc_PVQ_B;

        test();
        IF(shape_idx_fx == 0)
        {   /*  splitLF  shape  */
            Word16   n5k = 0;
            FOR(col = 1; col < 6; col++)
            {
                n5k += abs_s(y_split_fx[col]);
            }

            IF(sub(n5k, PULSES_SPLIT_A_LR) == 0)
            {
                /* this (6+2) pulses over NA+NB=5+8,   is  expected to be more frequent than 8 pulses over NA=5  */
                enc_PVQ_A = mpvq_index_fx(&(y_split_fx[1]), NA_LR, PULSES_SPLIT_A_LR); /*  P(N=5,K=6) (10)=10 bit L_index    */
                L_index_fx[4] = (Word32)enc_PVQ_A.index;                 move32();
                L_index_fx[1] = L_deposit_l(enc_PVQ_A.lead_sign_ind);  /* aux bit for the  splitLF path  ,  we plant the first LS there */

                enc_PVQ_B = mpvq_index_fx(&(y_split_fx[1 + NA_LR]), NB_LR, PULSES_SPLIT_B_LR);
                L_index_fx[5] = L_mac0(L_shl_pos((Word32)enc_PVQ_B.index, 1), 1, enc_PVQ_B.lead_sign_ind); move32(); /* A full PVQ 7 bit index for the  P(N=8,K=2) B config*/

                assert(L_index_fx[5] >= 0);
            }
            ELSE
            {
                ASSERT(n5k == (PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR));/*  PVQ(N=5,K=8) (12.x   in total, i.e.  LS+ 11.x )    */
                enc_PVQ_A = mpvq_index_fx(&(y_split_fx[1]), NA_LR, PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR);
                L_index_fx[4] = (Word32)enc_PVQ_A.index;
                L_index_fx[1] = L_deposit_l(enc_PVQ_A.lead_sign_ind);   /* aux bit for the  splitLF path  ,  we plant the first LS there */
                L_index_fx[5] = L_deposit_l(-8);  /* signal LF PVQ(5,k=8)  and  zeroed HF(10,0) apart */
            }
        }
        IF(sub(shape_idx_fx, 1) == 0)
        {  /*  full (15,5),   LS kept separated  */
            enc_PVQ_A = mpvq_index_fx(&(y_full_fx[1]), NFULL_LR, PULSES_FULL_LR);/*  mPVQ 16.66 bits in index[4], and LS 1 bit in index[1] */
            L_index_fx[4] = (Word32)enc_PVQ_A.index;            move32();
            L_index_fx[1] = L_deposit_l(enc_PVQ_A.lead_sign_ind);
        }

        IF(sub(shape_idx_fx, 2) >= 0)
        {    /* fixEnv0, fixEnv1, fixEnv2, fixEnv3 */
            ASSERT(shape_idx_fx <= 5);
            /* send the fixed env subshape mode to enc_entropy  */
            L_index_fx[4] = L_deposit_l(sub(shape_idx_fx, 2));  /* env shape, 0-->"1" , 1--> "env1"  */   /* L_index[2] has original shape 0...5 */

            ASSERT(fixShapeNb == L_index_fx[4]);
            ASSERT(fixShiftIdx < (1 << 2));

            Word16 tmp1 = add(1, fixShiftIdx);

            /* aux_bit : 0 (or 1)   , will indicate the  s0 sign in the FESS fix shape */
            L_index_fx[1] = 0L; move32(); /* s0 positive */
            test();
            if (y_fix_fx[tmp1] < 0)
            {
                L_index_fx[1] = 1L;  move32(); /* s0 negative*/
            }

            fix_end_sign_fx = 12;  move16();
            if (sub(shape_idx_fx, 5) == 0)
            {
                ASSERT(L_index_fx[4] == 3L); /* final fix_envelope only has 10 signs */
                fix_end_sign_fx = 10;  move16(); /* shape 5 has 2 bits shift and a total of 10 signs =2^10*2^2 = 2^12 = 4096 */
            }

            Word16 tmp = fixShiftIdx;  move16();  /* the two shift bits will be pushed up to b11,b12, for    11 signs s1-s11  */

            /* sign loop */
            FOR(int sign_ind_fx = 1; sign_ind_fx < fix_end_sign_fx; sign_ind_fx++)  /* push the remaining  sequential signs s1-s11(or s1-s9),  into a single idx */
            {   /* s1 is in the MSB, and s11 is in the LSB*/
                tmp = shl_pos(tmp, 1);  /* shift in  a zero */
                test();
                if (y_fix_fx[add(tmp1, sign_ind_fx)] < 0) /* "1" indicates  negative, "0" means positive */
                {
                    tmp = add(tmp, 1);
                }
            }
            L_index_fx[5] = L_deposit_l(tmp);

            ASSERT(L_index_fx[5] >= 0 && L_index_fx[5] < (1 << (2 + (fix_end_sign_fx - 1))));
        }
    } /* end of stage2 premultiplexing ,  fractional packing aspects done within enc_entropy_fx()  */

    ASSERT(envelope_bits_fx == 9 || envelope_bits_fx == 10 || envelope_bits_fx == 29 || envelope_bits_fx == 30); 

    Dyn_Mem_Deluxe_Out();
    return envelope_bits_fx;
}

/* LRSNS function needed in both encoder and decoder */
Word16 snslr_remove_st1_DC_fQ11_fx(                       /* o  : dc in Q11 */
    Word16  *scfq, /* i/o: stage1B  vector in Q11 */
    Word16  len    /* i  : length Q0 */
)
{
    Dyn_Mem_Deluxe_In(
        Counter i;
    Word32  L_dcQ11;
    Word16  dcQ11;
    );
    assert(len == 16);
    BASOP_sub_sub_start("snslr_remove_st1_DC_fQ11_fx");

    L_dcQ11 = 0;   move32();

    for (i = 0; i < len; i++) {
        L_dcQ11 = L_mac0(L_dcQ11, scfq[i], 2048);
    }
    dcQ11 = extract_h(L_shl_pos(L_dcQ11, 1));        /*NB ! same truncation needed in both CFL and in basop  */

    for (i = 0; i < len; i++)
    {
        scfq[i] = sub(scfq[i], dcQ11);   move16(); /* result update */
    }

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();

    return dcQ11;  /* output used for encoder side  mse update*/
}


void snslr_st1B_vector_dec_fx(Word16 idx, const Word16* LFCB, const Word16 *HFCB, const Word16* seg_cnt_cum, const Word16* idx12b_cb, Word16 *st1B_vectorQ11)
{
    /* decompose the received 0... 169 index , into the correct integer  and float st1B vector  */
    Dyn_Mem_Deluxe_In(
        Counter i;
    const Word16 *lf_cb, *hf_cb;

    Word16 seg_fx;
    Word16 idx_12b_fx, lf_sign_fx;
    Word16 hf_sign_fx;
    Word16 idx_HF_fx;
    Word16 idx_LF_fx;
    Word16 st1B_W16Q11_fx[M];
    Word16 buf[M];
    );

    BASOP_sub_sub_start("snslr_st1B_vector_dec_fx");
    assert(idx >= 0 && idx < 170);


    seg_fx = 0; move16();
    WHILE(sub(seg_cnt_cum[add(seg_fx, 1)], idx) <= 0)
    {
        seg_fx = add(seg_fx, 1);
    }

    assert(seg_fx >= 0 && seg_fx < 4);

    idx_12b_fx = idx12b_cb[idx];  /* indirect lookup from sequential value to a coded 12b index*/

    lf_sign_fx = 1; move16(); /*   assume a 0 bit -> "+"  */
    if (s_and(0x0800, idx_12b_fx) != 0) {
        lf_sign_fx = -1; move16(); /*   assume a 1 bit -> " -"  */
    }

    hf_sign_fx = 1; move16(); /*   assume a 0 bit -> "+"  */
    if ((s_and(0x0400, idx_12b_fx)) != 0) {
        hf_sign_fx = -1; move16(); /*   assume a 1 bit -> "-"  */
    }
    idx_LF_fx = shr_pos(s_and(0x03e0, idx_12b_fx), 5);
    idx_HF_fx = s_and(0x001f, idx_12b_fx);

    /* extseg0  f,f  */
    lf_cb = &(LFCB[idx_LF_fx*(M / 2)]); L_mult(0, 0); /* mult and pointer init */
    hf_cb = &(HFCB[idx_HF_fx*(M / 2)]); L_mult(0, 0); /* mult and pointer init */

    FOR(i = 0; i < (M / 2); i++)
    {
        st1B_W16Q11_fx[i] = i_mult(lf_sign_fx, lf_cb[i]);     /* imult() to switch sign without changing dynamics*/
        st1B_W16Q11_fx[M / 2 + i] = i_mult(hf_sign_fx, hf_cb[i]);
    }

    basop_memcpy(buf, st1B_W16Q11_fx, sizeof(*buf)*M);  /* buffer cpy needed for reversal sections */

    IF(s_and(seg_fx, 0x0002) != 0)
    {   /* r,*  */     /* flip LF */
        FOR(i = 0; i < (M / 2); i++)
        {
            st1B_W16Q11_fx[i] = buf[(M / 2 - 1) - i]; move16();
        }
    }

    IF(s_and(seg_fx, 0x0001) != 0)
    {   /* *,r */ /* flip HF */
        FOR(i = 0; i < (M / 2); i++)
        {
            st1B_W16Q11_fx[(M / 2) + i] = buf[(M - 1) - i]; move16();
        }
    }

    basop_memcpy(st1B_vectorQ11, st1B_W16Q11_fx, M * sizeof(*st1B_vectorQ11));

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
}

Word16 MSEsearchCbBIdxMap_fx(const Word16 *scf, const Word16 *LFCB, const Word16 *HFCB, const Word16 *seg_cnt_cum, const Word16* idx12b_cb, Word16 v_len, Word16 cb_len, Word32* min_mse)
{
    Dyn_Mem_Deluxe_In(

        Counter  seg, i, j; /*counters */

    Word16   scfLF_Q11[M];

    Word16*  scfHF_Q11;
    const Word16 *lf_cb, *hf_cb;
    Word16  idx_12b, signbitLF, signbitHF, idx_LF, idx_HF;
    Word32 L_mse_best_fx;
    Word16 best_ind;
    );
    UNUSED(v_len);
    UNUSED(cb_len);
    Word32 L_tEnBy2;
    Word32 L_cbEnBy2;
    Word32 L_mseBy2;   /* sum( t_i^2 + cbB_i^2 - 2*t_i*cbB_i )/2    */
 
 

 
    /*   MSE  separated  into  targetEn + vectorB_en - 2 *correlation t*cbB */
    /*   vectorBenergy   is obtained from addition of two Word32 lookup  tables LFen and HFen */
 

    BASOP_sub_sub_start("MSEsearchCbBIdxMap_fx");

    L_mse_best_fx = INT_MAX; move32();
    best_ind = -1; /*for debug*/

    assert(v_len == M);

    L_tEnBy2 = L_mult0(scf[0], scf[0]);
    FOR(i = 1; i < M; i++)
    {
        L_tEnBy2 = L_mac0(L_tEnBy2, scf[i], scf[i]);   /* calc target Energy part of MSE */
    }
    L_tEnBy2 = L_shr_pos_pos(L_tEnBy2, 1);
 

    scfHF_Q11 = (&scfLF_Q11[M / 2]);  /* ptr init */

    FOR(seg = 0; seg < 4; seg++)
    {
        basop_memcpy(scfLF_Q11, scf, M * sizeof(*scfLF_Q11));
        /*seg==0: fwd, fwd *//*seg==1: fwd, rev *//*seg==2: fwd, fwd */ /*seg==3: rev, rev */

        IF(s_and(seg, 0x0002) != 0)
        {   /* {r,*}  */     /* flip LF */
            FOR(i = 0; i < (M / 2); i++)
            {
                scfLF_Q11[i] = scf[(M / 2 - 1) - i]; move16();
            }
        }
        IF(s_and(seg, 0x0001) != 0)
        {  /* {*,r} */    /* flip HF */
            FOR(i = 0; i < (M / 2); i++)
            {
                scfHF_Q11[i] = scf[(M - 1) - i]; move16();
            }
        }

        FOR(i = seg_cnt_cum[seg]; i < seg_cnt_cum[add(seg, 1)]; i++)
        {
            /* Note: these 4 subindex extractions can be in parallel */
            idx_12b = idx12b_cb[i];    /* tmp variable, indirect adressing lookup of  12b index pointing to LF and HF + individual sign swaps  */

            idx_LF = shr_pos_pos(s_and(0x03e0, idx_12b), 5); /* b9...b5              */
            idx_HF = s_and(0x001f, idx_12b);             /* b4...b0 lowest 5 bits */

            L_cbEnBy2 = L_add(lrsns_st1B_enBy2TabW32_fx[idx_LF], lrsns_st1B_enBy2TabW32_fx[32 + idx_HF]);  /* ptr init to HF tab*/
 
            signbitLF = s_and(0x0800, idx_12b);    /* b11 logical  0 or 2048   */
            signbitHF = s_and(0x0400, idx_12b);    /* b10 logical  0 or 1024    */

           /* conditional update of signed section offsets */
            if (signbitLF != 0)
            {
                idx_LF = add(idx_LF, 32);  /* point to negated part of lf_cb  */
            }
            if (signbitHF != 0)
            {
                idx_HF = add(idx_HF, 32);   /* point to negated part of hf_cb */
            }
 
            lf_cb = &(LFCB[idx_LF * M / 2]);     /* adaptive ptr init */
            hf_cb = &(HFCB[idx_HF * M / 2]);    /* adaptive ptr init */
              
            L_mseBy2 = L_add(L_tEnBy2, L_cbEnBy2);   /* in ARM one op per cbB idx [0..169] can likely be saved by adding targetEnergyBy2 last  */

            /* mse/2  =   tEn^2/2  +  cbEn^2/2    -    corr(target_v, (sign)*cb_v )   */
            FOR(j = 0; j < (M / 2); j++)
            {
                /* cycles saved by extending the LF and HF tables with negated versions  */
                L_mseBy2 = L_msu0(L_mseBy2, scfLF_Q11[j], lf_cb[j]);     /* acc with    -1 * t*c ,  cb sign set by lf_cb */
                L_mseBy2 = L_msu0(L_mseBy2, scfHF_Q11[j], hf_cb[j]);     /* acc  -1* t*c with cb sign  set by hf_cb ptr  */
            }
            ASSERT(L_mseBy2 >= 0);
            /* always update best case  */
            L_mse_best_fx = L_min(L_mseBy2, L_mse_best_fx); /* 1 cycle BASOP  update     */

            if (L_sub(L_mseBy2, L_mse_best_fx) == 0)
            {
                best_ind = ((Word16)i);    move16();  /* update winner, cond move , single BASOP */
            }
        } /* one segment seg*/
    } /* all  segms */

    *min_mse = L_add(L_mse_best_fx, L_mse_best_fx);   /* multiply halved  MSEby2 to get final MSE  */
 
    assert(best_ind >= 0 && best_ind < cb_len);

    BASOP_sub_sub_end(); /* func */
    Dyn_Mem_Deluxe_Out();

    return best_ind;
}

void snslr_st1C_vector_dec_fx(Word16 idx, const Word8* CBW8, Word16 scaleQ4, Word16 inv_scaleQ15, Word16 v_len, Word16 cb_len, Word16 *st1C_vector)
{
    /* decompose the received 0... 169 index , into the correct (integer and)  float st1C vector  */
    /* even in C-float the st1C coeffs  are put into a S16Q11 final integers domain  */
    /* Enables BE compatibility between {BASOP, float, double}  arithmetic implmentations */
    Dyn_Mem_Deluxe_In(
        Counter i;
    const Word8 *cb;
    Word32  L_tmp;
    Word16  s_tmp;
    );
    UNUSED(v_len);
    UNUSED(cb_len);
    UNUSED(inv_scaleQ15);  /* req for debugging only */

    BASOP_sub_sub_start("snslr_st1C_vector_dec_fx");

    assert(idx >= 0 && idx < cb_len);
    assert(v_len == M);

    cb = &(CBW8[idx*M]); move16();                 /* pointer init */
    FOR(i = 0; i < M; i++)
    {
        L_tmp = L_mult0((Word16)cb[i], scaleQ4); /*S8Q7 * S15Q4 */ /*sign+7bit, sign+4 bits --> sign+11bit  .lt  sign+23 bits*/
        s_tmp = extract_l(L_tmp);
        assert(L_tmp >= -32768L && L_tmp <= 32767L); /* INT16 domain check*/
        st1C_vector[i] = s_tmp;
    }

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
}

Word16 MSEsearchGenericScaledW8_fx(Word16 *scf, const Word8 *sns_CBW8, Word16 scaleQ4, Word16 inv_scaleQ15, Word16 v_len, Word16 cb_len, Word32* L_min_mse)
{
    /*  scf float input values are  typically in the range  +12.0  to -12.0.
        rom table stored in Word8 [+127,-128], format  corresponding to   ]+1.0  .. -1.0 ]
        inv_scaleQ15, [downscaling  value in Q15]  applied before search
        scaleQ12   upscaling value quantized in Q12,  used in the mse calulation and in the common float and BASOP  synthesis routines

     */
    Dyn_Mem_Deluxe_In(
        Counter i, n;
    Word32  best_ind;
    Word32 L_mse, L_mse_best;
    Word16  targetW16[M];
    const Word8  * cbW8;
    Word16 shift_mse;
    Word16 tmp;
    );

    BASOP_sub_sub_start("MSEsearchGenericScaledW8_fx");
    FOR(i = 0; i < v_len; i++)
    {
        targetW16[i] = shr_pos(mult(scf[i], inv_scaleQ15), 4); /*   target from W16Q11 to    W8Q7 domain ,  and floor */
    };

    best_ind = -1;    /* needed to avoid compiler warning  */

    /*   use   use S1.15 wise multiplication operations QSUB16 for two coeffs correlation  */
    /*   L_msu0  can use 16 bit input only e.g.  SMLSLD   */

    /*   MSE_i = sum(i=0,15,  (t_i-c_i)*(t_i-c_i) ) = sum(i=0,15,  (t_i)^2 + c_i^2 - 2*(t_i*c_i) )   */
    /*     store  EnC_i = sum(1,16) c_i^2 , in a Word32 ROM entry , 170 entries , or potentially in Word16   */
    /*   calculate           MSE_i    =  EnT+EnC   -  sum (2*t_i*c_i)   */
    /*   or even calculate   MSE_i/2  =  enT/2 + EnC/2  - sum(  t_i*c_i)   , as an approximation of MSE */

    /* commmon target contribution to MSE */

    Word32 L_targetEn = 0; move32();
    if (v_len > 0) {
        L_targetEn = L_mult0(targetW16[0], targetW16[0]);
        FOR (i = 1; i < v_len; i++) 
        {
            L_targetEn = L_mac0(L_targetEn, targetW16[i], targetW16[i]);
        }
    }
    L_mse_best = INT32_MAX; move32(); /* largest possible positive number in INT32 */

    FOR(i = 0; i < cb_len; i++)
    {
        cbW8 = &(sns_CBW8[i*M]);  /* ptr init += i*16 */
        assert(v_len == M);

        L_mse = L_deposit_l(lrsns_st1C_Both_EnBy2Tab_fx[i]);    /* scaled down CB energy by 2 from ROM to allow use of L_msu0 for best ARM ops  */
        L_mse = L_add(L_mse, L_shr_pos_pos(L_targetEn, 1));  /* scale down target energy by 2 to allow use of L_msu0 for best ARM ops  */

        FOR(n = 0; n < v_len; n++)
        {  /* correlate is fast in most ARM cases, as SIMD or as   blocks of 4 */
            L_mse = L_msu0(L_mse, targetW16[n], (Word16)(cbW8[n]));   /* no upscaling included in L_msu0()  */
        }

        assert(L_mse >= 0);  /* check that error is always positive, otherwise we might get overflow in the subtraction below */
        if (L_sub(L_mse, L_mse_best) <= 0)
        {
            best_ind = i;   move16();  /* single BASOP, conditional move,  for best_ind update */
        }
        L_mse_best = L_min(L_mse, L_mse_best); /* always update best MSE using L_min() in the idx loop, reduces WC WMOPS  */
    }
    L_mse_best = L_add(L_mse_best, L_mse_best);      /* scale up by two from the BASOP/ARM optimized search mse domain */

 
    assert(best_ind >= 0 && best_ind < cb_len);

    *L_min_mse = L_mse_best;     move32();

    tmp = i_mult(scaleQ4, scaleQ4); /* non fractional square --> Q4 to Q8 */
    tmp = shl_pos(tmp, 2);          /* a Q10 const */

    shift_mse = norm_l(*L_min_mse);

    *L_min_mse = L_shl(*L_min_mse, shift_mse);
    *L_min_mse = Mpy_32_16(*L_min_mse, tmp);     /* still use max upshifting Mpy_32_16,  here for this final Word32 L_mse  output calculation  */
    *L_min_mse = L_shr(*L_min_mse, add(shift_mse, 2 - 15));

    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();

    return best_ind;
}

Word16 MSEsearchGeneric_fx(Word16 *scf, const Word16 *sns_CB, Word16 v_len, Word16 cb_len, Word32* L_min_mse)
{
    Dyn_Mem_Deluxe_In(
        Counter   i, n;
    Word16 ind;
    Word32 L_mse_best_fx, L_mse_fx;
    Word16 tmp;
    Word16 const *cb_ptr;
    Word16 const *scf_ptr;
    );

    BASOP_sub_sub_start("MSEsearchGeneric_fx");
    ind = -1;  /* avoid compiler warning */
    L_mse_best_fx = INT_MAX; move32();

    FOR(i = 0; i < cb_len; i++)
    {
        L_mse_fx = L_deposit_l(0);
        cb_ptr = &(sns_CB[i*v_len]); /* pointer init */
        scf_ptr = &(scf[0]);         /* pointer init */
        FOR(n = 0; n < v_len; n++)
        {
            tmp = sub(*scf_ptr++, *cb_ptr++);
            L_mse_fx = L_mac0(L_mse_fx, tmp, tmp);
        }

        L_mse_best_fx = L_min(L_mse_fx, L_mse_best_fx);
        if (L_sub(L_mse_fx, L_mse_best_fx) == 0)
        {
            ind = i; move16();
        }
    }

    *L_min_mse = L_mse_best_fx; move32();

    assert(ind >= 0 && ind < cb_len);
    BASOP_sub_sub_end();
    Dyn_Mem_Deluxe_Out();
    return ind;
}

#endif
