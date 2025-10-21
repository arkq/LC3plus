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

#ifdef  CR9_C_ADD_1p25MS_LRSNS
    static void pvq_pyr_project_lrsns_adv(const Word16  dim_proj,   /*  vector dimension       */
        const Word16 *xabs,                      /* absolute vector values */
        Word32        L_xsum,                    /* absolute vector sum over dim_proj  */
        Word16        num,                       /* start/safe target number of pulses for projection to a pyramid e.g. KMAX-1   */
        Word16        num_max,                   /* max  number of pulses allowed , typically K_MAX  */
        Word16 *      y,                         /* projected integer output vector   */
        Word16 *pulse_tot_ptr, Word32 *L_xy_ptr, /* accumulated correlation  Q(in+0+1) = Qin+1 */
        Word32 *L_yy_ptr,                         /* accumulated energy  Q0  */
        Word8 * scratch
    )
    {
        Dyn_Mem_Deluxe_In(
            Counter i;
        Word32  L_tmp, L_num;
        Word16  den, shift_num, shift_den, shift_delta, proj_fac;
        Word16  *y_r, *y_r_soft;
        Word16  pulse_tot_r, pulse_tot_r_soft;
        );

#ifdef FIX_BASOP_ENC_LRSNS_ST2FULL_PROJ 
        Word16    skip_rnd_flag = 0;
#endif 

        y_r      = (Word16 *)scratchAlign(scratch, 0);                  
        y_r_soft = (Word16 *)scratchAlign(y_r, M * sizeof(*y_r));   /*  y_r      = M   */
                                                                     /* y_r_soft = M  */
 
        *pulse_tot_ptr = 0; move16();
        pulse_tot_r = 0; move16();
        pulse_tot_r_soft = 0; move16();

        shift_den = norm_l(L_xsum);                          /* x_sum input  Qin                         */
        den = extract_h(L_shl_pos(L_xsum, shift_den)); /* now in Qin+shift_den                     */

#ifdef FIX_BASOP_ENC_LRSNS_ST2FULL_PROJ 
        IF(num_max == 0 ) {
            num_max = num;    move16(); /* use projection value num as num_max  limit*/
            skip_rnd_flag = 1; move16();
        }
#endif 

        L_num = L_deposit_l(num);
        shift_num = sub(norm_l(L_num), 1);
        L_num = L_shl_pos(L_num, shift_num); /* now in Q0 +shift_num -1                  */

        proj_fac = div_l(L_num, den); /* L_num always has to be less than den<<16 , norm_l-1  makes that happen  */

        shift_delta = sub(shift_num, shift_den);

        /* we are here using a single proj factor, but different rounding offsets  0(==trunc=safe), .25(soft), .5(==round=risky)  */

        FOR(i = 0; i < dim_proj; i++)
        {
            L_tmp = L_mult(proj_fac, xabs[i]);               /* Q  shift_delta + PVQ_SEARCH_QIN */
            L_tmp = L_shr(L_tmp, shift_delta);

            y[i] = extract_h(L_tmp); move16();                /*  to  Q0 with floor   */
            *pulse_tot_ptr = add(*pulse_tot_ptr, y[i]);       /* Q0                                         */

            y_r[i] = extract_h(L_add(L_tmp, (Word32)0x00008000L));  /* regular round  i.e. +.5 and truncate */
            pulse_tot_r = add(pulse_tot_r, y_r[i]);       /* Q0                                         */

            y_r_soft[i] = extract_h(L_add(L_tmp, (Word32)0x00004000L));  /* softer round i.e. +.25 and truncate */
            pulse_tot_r_soft = add(pulse_tot_r_soft, y_r_soft[i]);       /* Q0                                         */
        }

        /* now analyze which rounding is valid and closest to num_max */

#ifdef FIX_BASOP_ENC_LRSNS_ST2FULL_PROJ 
        IF (skip_rnd_flag != 0 ) 
        {
            pulse_tot_r = 32767;       /*disable most optimistic rounding */
            pulse_tot_r_soft = 32767;  /*disable soft optimistic rounding */
        }
#endif 

        /* if y_r is successful and does not overshoot, let us use that */
        IF(sub(pulse_tot_r, num_max) <= 0)
        {
            basop_memcpy(y, y_r, dim_proj * sizeof(*y));
            *pulse_tot_ptr = pulse_tot_r;
            pulse_tot_r_soft = 32767;  /* exclude r_soft */
        }

        /* if y_r_soft is successful and does not overshoot,   use that */
        IF(sub(pulse_tot_r_soft, num_max) <= 0)
        {
            basop_memcpy(y, y_r_soft, dim_proj * sizeof(*y));
            *pulse_tot_ptr = pulse_tot_r_soft;
            pulse_tot_r = 32767;
        }

        /*sum up correlation and  energy */
        *L_xy_ptr = L_deposit_l(0);
        *L_yy_ptr = L_deposit_l(0);
        FOR(i = 0; i < dim_proj; i++)
        {
            *L_yy_ptr = L_mac0(*L_yy_ptr, y[i], y[i]);   /* Energy,   Q0 */
            *L_xy_ptr = L_mac(*L_xy_ptr, xabs[i], y[i]); /* Corr, Q11*Q0  +1  --> Q12   */
        }
        ASSERT(*pulse_tot_ptr <= num_max);

        Dyn_Mem_Deluxe_Out();
    }
#endif 

static void pvq_pyr_project(const Word16  dim_proj,                  /* end vector dimension+1       */
                            const Word16 *xabs,                      /* absolute vector values */
                            Word32        L_xsum,                    /* absolute vector sum over dim  */
                            Word16        num,                       /* target number of pulses */
                            Word16 *      y,                         /* projected output vector    */
                            Word16 *pulse_tot_ptr, Word32 *L_xy_ptr, /* accumulated correlation  Q(in+0+1) = Qin+1 */
                            Word32 *L_yy_ptr                         /* accumulated energy  Q0  */
)
{

    Dyn_Mem_Deluxe_In(
        Counter i;
        Word32  L_tmp, L_num;
        Word16  den, shift_num, shift_den, shift_delta, proj_fac;
    );

    *pulse_tot_ptr = 0; move16();
    *L_xy_ptr      = L_deposit_l(0);
    *L_yy_ptr      = L_deposit_l(0);

    shift_den = norm_l(L_xsum);                          /* x_sum input  Qin                         */
    den       = extract_h(L_shl_pos(L_xsum, shift_den)); /* now in Qin+shift_den                     */

    L_num     = L_deposit_l(num);
    shift_num = sub(norm_l(L_num), 1);
    L_num     = L_shl_pos(L_num, shift_num); /* now in Q0 +shift_num -1                  */
    proj_fac  = div_l(L_num, den); /* L_num always has to be less than den<<16 , norm_l-1  makes that happen  */

    shift_delta = sub(shift_num, shift_den);
    FOR (i = 0; i < dim_proj; i++)
    {
        L_tmp = L_mult(proj_fac, xabs[i]);                       /* Q  shift_delta + PVQ_SEARCH_QIN */
        y[i]  = extract_h(L_shr(L_tmp, shift_delta)); move16(); /*  to  Q0 with floor , and potential  sturation */
        ;

        *pulse_tot_ptr = add(*pulse_tot_ptr, y[i]);       /* Q0                                         */
        *L_yy_ptr      = L_mac0(*L_yy_ptr, y[i], y[i]);   /* Energy,   Q0 */
        *L_xy_ptr      = L_mac(*L_xy_ptr, xabs[i], y[i]); /* Corr, Q0*Q12  +1 --> Q13                   */
    }

    Dyn_Mem_Deluxe_Out();
}


static __forceinline Word16 one_pulse_search(const Word16  dim_start, /* start vector dimension       */
                                             const Word16  dim_end,   /* end vector dimension+1       */
                                             const Word16 *x_abs,     /* absolute vector values */
                                             Word16 *      y,         /* output vector    */
                                             Word16 *      pulse_tot_ptr,
                                             Word32 *      L_xy_ptr, /* accumulated correlation  Q(12+0+1) = Q13 */
                                             Word32 *      L_yy_ptr, /* accumulated energy  Q0 */
                                             Word16        max_xabs)        /* current max amplitude for target  */
{
    Dyn_Mem_Deluxe_In(
        Counter i;
        Word16  corr_tmp, corr_sq_tmp, en_max_den, cmax_num, en_tmp;
        Word32  L_tmp_en_lc, L_tmp_corr;
        Word16  corr_up_shift, imax;
    );

    /* maximize correlation precision, prior to every unit pulse addition in the vector */
    corr_up_shift = norm_l(L_mac(*L_xy_ptr, 1, max_xabs)); /* pre analyze worst case L_xy update in the dim  loop  */
    imax          = -1; /* not needed for search, only added to avoid compiler warning   */
    {
        en_max_den = 0;             move16();
        cmax_num   = -1; move16(); /* req. to force a 1st update for  n==dim_start   */

        FOR (i = dim_start; i < dim_end; i++)
        {
            L_tmp_corr = L_shl_pos(L_mac(*L_xy_ptr, 1, x_abs[i]), corr_up_shift); /*  actual in-loop target value */

            corr_tmp = round_fx_sat(L_tmp_corr);

            corr_sq_tmp = mult(corr_tmp, corr_tmp); /* CorrSq_tmp for a 16bit low complexity cross multiplication */

            L_tmp_en_lc = L_mac(*L_yy_ptr, 1, y[i]); /*Q0 x^2+ 2x ,  "+1" added once before loop , result ,  energy may
                                                        span up to ~14+1(Q1)+1(sign)=16 bits */
            /* extract_l without shift can always be used for this section as energy is guaranteed to stay in the lower
             * word*/

            en_tmp = extract_l(L_tmp_en_lc); /* L_shl + round_fx could also be used also but then adds an uphift cost */

            /* 16/32 bit comparison    WC (4 +1+1 + (1+1+1) = 9 */
            IF (L_msu(L_mult(corr_sq_tmp, en_max_den), cmax_num, en_tmp) > 0) /* use L_mult and then a L_msu */
            {
                cmax_num   = corr_sq_tmp; move16();
                en_max_den = en_tmp;      move16();
                imax       = i;           move16();
            }
        } /* dim  */
    }


    /*  finally add found unit pulse contribution to past L_xy, Lyy,  for next pulse loop    */
    *L_xy_ptr = L_mac(*L_xy_ptr, x_abs[imax], 1); /*      Qin+1 */
    *L_yy_ptr = L_mac(*L_yy_ptr, 1, y[imax]);

    y[imax]          = add(y[imax], 1); move16(); /* Q0 added pulse              */
    (*pulse_tot_ptr) = add((*pulse_tot_ptr), 1);   /* increment total pulse sum   */
    Dyn_Mem_Deluxe_Out();
    return imax;
}

#ifdef CR9_C_ADD_1p25MS_LRSNS

/* evaluate corr/sqrt(en) = corr*inv_sqrt_tab[en],    as  CorrSq/en ratio  cross-multiplication may cost more  */

static __forceinline Word16 one_pulse_search_tab_isqrt(
    const Word16  dim_start, /* start vector dimension       */
    const Word16  dim_end,   /* end vector dimension+1       */
    const Word16 *x_abs,     /* absolute vector values */
    Word16 *      y,         /* output vector    */
    Word16 *      pulse_tot_ptr,
    Word32 *      L_xy_ptr, /* accumulated correlation  Q(12+0+1) = Q13 */
    Word32 *      L_yy_ptr, /* accumulated energy  Q0 */
    Word16        max_xabs)  /* current max amplitude for target  */
{
    Dyn_Mem_Deluxe_In(
        Counter i;

    Word32  L_cand_en_lc, L_cand_corr;
    Word16  corr_up_shift, imax;
    Word32  L_cand_norm_ratio, L_max_norm_ratio;
    );

    /* maximize correlation precision, prior to every unit pulse addition in the vector */
    corr_up_shift = norm_l(L_mac0(*L_xy_ptr, 2, max_xabs)); /* pre analyze worst case L_xy update in the dim  loop  */
    imax = -1; /* init not needed for this search, only added to avoid MSVC compiler warning   */

    /* select isqrt table for this pulse addition loop , based on current energy    */
    const Word16* isqrt_QxTab = isqrt_Q16tab;   move32();   /* Q16 table valid for energies 4...64   */
    if (L_sub(*L_yy_ptr, 3L) <= 0) /* +1 for the inloop value of L_yy was preadded outside to *L_yy_ptr */
    {
        isqrt_QxTab = isqrt_Q15tab;    move32();    /* energies:  1...6   */
    }

    {
        L_max_norm_ratio = L_deposit_l(-1);  /* req. to force a 1st update for  n==dim_start   */

        FOR(i = dim_start; i < dim_end; i++)
        {
            L_cand_corr = L_shl_pos(L_mac0(*L_xy_ptr, 2, x_abs[i]), corr_up_shift); /*  actual in-loop target value */
            L_cand_en_lc = L_mac0(*L_yy_ptr, 2, y[i]); /*Q0 y^2+ 2y +1,   "+1" added once before loop  outside */

            L_cand_norm_ratio = Mpy_32_16_0_0(L_cand_corr, isqrt_QxTab[L_cand_en_lc]);  /* energy normalized ratio,  non-saturating operation  */
 
            /* 32 bit comparison    WC (1 +1 +1 = 2 */
            if (L_sub(L_cand_norm_ratio, L_max_norm_ratio) >= 0)
            {
                imax = i;           move16();  /* conditional move single basop */
            }
            L_max_norm_ratio = L_max(L_max_norm_ratio, L_cand_norm_ratio); /* always update */
        } /* dim  */
    }

    /*  finally add found unit pulse contribution to past L_xy, L_yy,  for next pulse loop    */
    *L_xy_ptr = L_mac0(*L_xy_ptr, x_abs[imax], 2); /*      Qin*Q0+1 */
    *L_yy_ptr = L_mac0(*L_yy_ptr, 2, y[imax]);     /*  en = en + 2*y + 1  ,  1 was preadded outside */

    y[imax] = add(y[imax], 1);       move16();   /* Q0 added pulse              */
    (*pulse_tot_ptr) = add((*pulse_tot_ptr), 1);   /* increment total pulse sum   */

    Dyn_Mem_Deluxe_Out();
    return imax;
}
#endif 

void pvq_enc_search_fx(
    const Word16 *x,           /* i:   target vector to quantize       Qin     */
    Word16 *      y_far,       /* o:   outl_far o, raw pulses  (non-scaled short) Q0  , length dim     */
    Word16 *      y,           /* o:   outl_near o, raw pulses  (non-scaled short) Q0  , length dim     */
    Word16 *      yA,          /* o:   A section  raw pulses  (non-scaled short) Q0 , length dimA   */
    Word16 *      yB,          /* o:   B section  raw pulses  (non-scaled short) Q0   , length dim-dimA  */
    Word32 *      L_corr,      /* o:   4 un-normalized correlation sums for outl_far, outl_near, A, AB  */
    Word32 *      L_search_en, /* o:   4  energy sums for outl_far, outl_near,  A, AB  */
    Word16 *      pulses_fin,  /* i:   number of allocated pulses  to outl_far, outl_near ,  A, AB  sections   */
    Word16 *      pulses_proj, /* i:   number of projection pulses  for outl_far, outl_near,  A, AB    */

    const Word16 dim, /* i:   Length of outlier  vector */
    const Word16 dimA /* i:   Length of vector A section */
)
{

    Dyn_Mem_Deluxe_In(
        Counter       i;
        Word16        pulse_tot_far, pulse_tot, pulse_totA, pulse_totB;
        Word16        xabs[PVQ_MAX_VEC_SIZE];
        Word16        max_xabs, max_xabsA, max_xabsB;
        Word32        L_xsum, L_xsumA;
        Word32        L_yy, L_xy;
        Word16        imax;
        Counter       k;
        Word16        dim_m1;
        Word16        dimB;
        const Word16 *xBptr;
        Word16        pulses_far, pulses, pulsesA, pulsesB;
    );

    BASOP_sub_sub_start("pvq_enc_search_fx");

    pulses_far = pulses_fin[0]; move16();
    pulses     = pulses_fin[1]; move16();
    pulsesA    = pulses_fin[2]; move16();
    pulsesB    = pulses_fin[3]; move16();

    FOR (i = 0; i < N_SCF_SHAPES_ST2; i++)
    {
        L_corr[i]      = L_deposit_l(0);
        L_search_en[i] = L_deposit_l(0);
    }

    dimB = sub(dim, dimA);

    L_xsum = L_deposit_h(0);

    max_xabs  = -1; move16();
    max_xabsA = -1; move16();
    max_xabsB = -1; move16();
    FOR (i = 0; i < dimA; i++)
    {
        xabs[i]   = abs_s(x[i]); move16();     /* Qx */
        max_xabsA = s_max(max_xabsA, xabs[i]);  /* for efficient  search correlation scaling */
        L_xsum    = L_mac0(L_xsum, 1, xabs[i]); /* stay in Qx */
    }

    basop_memset(y_far, 0, dim * sizeof(Word16));
    basop_memset(y, 0, dimA * sizeof(Word16));
    basop_memset(yA, 0, dimA * sizeof(Word16));

    L_xsumA = L_add(L_xsum, 0); /* save for section A projection */

    FOR (i = dimA; i < dim; i++)
    {
        xabs[i]   = abs_s(x[i]); move16();     /* Qx */
        max_xabsB = s_max(max_xabsB, xabs[i]);  /* for efficient  search correlation scaling */
        L_xsum    = L_mac0(L_xsum, 1, xabs[i]); /* stay in Qx */
    }

    basop_memset(&y[dimA], 0, (dim - dimA) * sizeof(Word16));

    basop_memset(yB, 0, dimB * sizeof(Word16));

    max_xabs = s_max(max_xabsA, max_xabsB); /* global max abs value */

    test();
    IF (L_xsum == 0)
    { /* no shape in any  section, projection in outl_far, outl_near, A, AB  not possible, any search meaningless  */

        dim_m1        = sub(dim, 1);
        y_far[0]      = shr_pos(pulses_far, 1);                        move16();
        y_far[dim_m1] = add(y_far[dim_m1], sub(pulses_far, y_far[0])); move16();

        dim_m1    = sub(dim, 1);
        y[0]      = shr_pos(pulses, 1);                move16();
        y[dim_m1] = add(y[dim_m1], sub(pulses, y[0])); move16();

        dim_m1     = sub(dimA, 1);
        yA[0]      = shr_pos(pulsesA, 1);                  move16();
        yA[dim_m1] = add(yA[dim_m1], sub(pulsesA, yA[0])); move16();

        dim_m1     = sub(dimB, 1);
        yB[0]      = shr_pos(pulsesB, 1);                  move16();
        yB[dim_m1] = add(yB[dim_m1], sub(pulsesB, yB[0])); move16();
    }
    ELSE
    {
        ASSERT(pulses_proj[0] > 0);
        ASSERT(L_xsum > 0);

        pvq_pyr_project(dim, xabs, L_xsum, pulses_proj[0], y_far, &pulse_tot_far, &L_xy,
                        &L_yy); /* outlier  submode projection  */

        ASSERT(pulses_far <= 127);
        FOR (k = pulse_tot_far; k < pulses_far; k++)
        {
            L_yy = L_add(L_yy, 1); /* pre add 1 in  Q0 in    L_yyQ0 = (x^2 + 2*x + 1)    */
            imax = one_pulse_search(0, dim, xabs, y_far, &pulse_tot_far, &L_xy, &L_yy, max_xabs);
        }
        ASSERT(pulse_tot_far == pulses_far);
        /* outlier far submode result vector in   y_far[0...15]  */
        L_corr[0] = L_shr_pos(L_xy, 1); /* to Qin*Q0 */

        basop_memmove(y, y_far, dim * sizeof(Word16)); /*y_far->y  */

        pulse_tot = pulse_tot_far; move16();

        ASSERT(pulses <= 127);
        FOR (k = pulse_tot; k < pulses; k++)
        {
            L_yy = L_add(L_yy, 1); /* pre add 1 in  Q0 in    L_yyQ0 = (x^2 + 2*x + 1)    */
            imax = one_pulse_search(0, dim, xabs, y, &pulse_tot, &L_xy, &L_yy, max_xabs);
        }

        /* outlier near submode result vector in   y[0...15]  */
        L_corr[1] = L_shr_pos(L_xy, 1); /* to Qin*Q0 */

        ASSERT(pulse_tot == pulses);

        IF (L_xsumA == 0)
        {
            /* no shape in A section, projection in A not possible,  search meaningless  */
            dim_m1     = sub(dimA, 1);
            yA[0]      = shr_pos(pulsesA, 1);                  move16();
            yA[dim_m1] = add(yA[dim_m1], sub(pulsesA, yA[0])); move16();
        }
        ELSE
        {
            IF (pulses_proj[2] != 0) /* fixed setup  if bitrate is fixed */
            {
                ASSERT(pulses_proj[2] > 0);
                ASSERT(L_xsumA > 0);
                pvq_pyr_project(dimA, xabs, L_xsumA, pulses_proj[2], yA, &pulse_totA, &L_xy,
                                &L_yy); /* section A  , in submode 1 projection  */
            }
            ELSE
            {
                /*  default, otherwise recalculate A   from outlier result  (to remove any section B pulses influence)
                 */
                pulse_totA = 0; move16();
                L_xy       = L_deposit_l(0);
                L_yy       = L_deposit_l(0);

                basop_memmove(yA, y, dimA * sizeof(Word16));
                FOR (i = 0; i < dimA; i++)
                {
                    pulse_totA = add(pulse_totA, yA[i]);      /* Q0                                         */
                    L_xy       = L_mac(L_xy, xabs[i], yA[i]); /* Corr, Q0*Q12  +1 --> Q13                   */
                    L_yy       = L_mac(L_yy, yA[i], yA[i]);   /* Energy, Q(0+0)+1)= Q1 */
                }
                L_yy = L_shr_pos(L_yy, 1); /* En to Q0  */
            }

            /* search remaining pulses in regular section A  */
            FOR (k = pulse_totA; k < pulsesA; k++)
            {
                L_yy = L_add(L_yy, 1); /* 1 added in Q0 */
                imax = one_pulse_search(0, dimA, xabs, yA, &pulse_totA, &L_xy, &L_yy, max_xabsA);
            }
            ASSERT(pulse_totA == pulsesA);
        } /* L_xsumA!=0 */

        /* reg Set A result vector now in  yA[0...9]  */
        L_corr[2] = L_shr_pos(L_xy, 1); /* to Qin*Q0 */

        /* search remaining pulses in regular section B,  even if energy in B is zero  */
        ASSERT(pulses_proj[3] == 0);
        pulse_totB = 0; move16();

        IF (sub(pulsesB, 1) == 0)
        { /* LC search,  sufficient to find a single max,  as pulses can not be stacked, when nb-pulses==1  */
            imax = 0; move16(); /* safety */
            FOR (i = dimA; i < dim; i++)
            {
                if (xabs[i] == max_xabsB)
                {
                    imax = sub(i, dimA);
                }
            }
            pulse_totB = 1;                                     move16();
            yB[imax]   = 1; move16();                          /* reg set B result vector in  yB[0...5]  */
            L_xy       = L_mac(L_xy, xabs[add(imax, dimA)], 1); /* calc total corr for A+B sections */
            L_yy       = L_add(L_yy, 1);
        }
        ELSE
        { /* more than one target pulse in section B */
            /* keep A pulses influence,  search  section B pulses influence */
            FOR (k = pulse_totB; k < pulsesB; k++)
            {
                L_yy = L_add(L_yy, 1); /* 1 added in Q0*/
                imax = one_pulse_search(dimA, dim, xabs, &(yB[-dimA]), &pulse_totB, &L_xy, &L_yy, max_xabsB);
            }
        }

        L_corr[3] = L_shr_pos(L_xy, 1); move32(); /* to Qin*Q0 ,  corr of combined A and B */

        ASSERT(pulse_totB == pulsesB);
        /* reg set B result vector now in  yB[0...5]  */
    } /* L_xsum != 0 */

/* apply sign of (x) to  first orthant result */
    FOR (i = 0; i < dim; i++)
    {
        if (x[i] < 0)
        {
            y_far[i] = negate(y_far[i]); /* apply sign for outlier far */
        }
    }

    FOR (i = 0; i < dim; i++)
    {
        if (x[i] < 0)
        {
            y[i] = negate(y[i]); /* apply sign for outliers near */
        }
    }

    xBptr = &(x[dimA]); move32(); /* ptr init to B target section */
    FOR (i = 0; i < dimA; i++)
    {
        if (x[i] < 0)
        {
            yA[i] = negate(yA[i]); /* apply sign  in N_SETA */
        }
    }

    FOR (i = 0; i < (dimB); i++)
    {
        if (xBptr[i] < 0)
        {
            yB[i] = negate(yB[i]); /* apply sign  in N_SETB */
        }
    }

    Dyn_Mem_Deluxe_Out();
    BASOP_sub_sub_end();
}

#ifdef CR9_C_ADD_1p25MS_LRSNS

/* 29/30 bits  optimized search functions for PVQ and FESS  */
/* stage 2 submode shape  0:  "splitLF"  (N=5,K=6)(N=8,K=2)(N=2,K=0),  or (N=5,K=8)(N=10,K=0)  , 4 gains */
/* stage 2 submode shape  1:  "full" FB (N=15,K=5),        , 8 gains */
/* stage 2 submode shape  2-5: "fixed env "  (N=13-15,K=10-12), 4xfixed  , 8 gains */

#define SC 1  /* starting coeff for LRSNS VQ searches, as DC is   zero in the DCT-II domain target */
void pvq_fess_enc_search_fx(
    const Word16 *x_in,          /* i:   target vector to quantize       Qin (=Q11)   0...M-1    */

    Word16 *      y_Q0,          /* o:   raw integer pulses  (non-scaled short) Q0  , length 3*M          */
    Word16 *      y_normQ14,     /* o:   normalized  integer pulses  (non-scaled short) Q14 , length 3*M       */
    Word32 *      L_y_normQ30,   /* o:   normalized  integer pulses  (non-scaled short) Q30 , length 3*M       */
    Word16 *      s_idxPtr,      /* o:   quantized shape index,   0,1,2   */
    Word16 *      g_idxPtr,      /* o:    quantized gain index,   0..3, or 0..7 */
    Word16 *      g_qvalQ12Ptr,  /* o:    quantized gain value,   in Q12   -32767 == 7.99975  for best shape     */
    Word16 *      fixShapeNbPtr, /* o:  idx for the selected fix shape y_fix  0...3, only relevant in the case  s_idx==2  */
    Word16 *      fixShiftIdxPtr, /* o:  idx for the selected fix shift  0...3,   only relevant in the case  s_idx==2  */
    Word32 *      L_MSEQ22Ptr,   /* o:  1   Q11+Q0+1 -->  Q22  */
    Word8 * scratch    
)
{
 
    Dyn_Mem_Deluxe_In(
        Counter       i, k, n;
    Word16        pulse_tot;
    Word16        max_xabs, max_xabsA, max_xabsB, max_xabsAB, max_xabsC;
    Word32        L_xsum, L_xsumA, L_xsumB, L_xsumAB;
    Word32        L_yy, L_xy;

    Word16        imax;
    Word16        *y_split, *y_full, *y_fix;
    Word16        *y_splitAB;

    Word16        *y_splitA0;
    Word16        best_env_ind;
    Word16        best_shift_ind;
    Word32        L_targetEnNeg;

    Word16 shift_ind, fix_ind;
    Word16* xabs1;
    const Word16* envPtr;
    Word32 L_corr;
 

    Word16 *y, *y_norm, tmp;
    Word32 *L_y_norm;
    Word32 L_tmp;
    Word32 L_min_mse_opt;
 
    Word16 gain_idx_opt, shape_idx_opt;
    Word16 best_ind;
    Word16 gidx;
    Word32 L_g_q_tmp;
    const Word16 *gTabPtr;
    Word32 L_mse;
    Word32 L_xy_a6_mem;
    Word32  L_MSEQ22_recalc;  
    Word32  L_normcorrQy;     
    Word16 norm_factors[N_SCF_SEARCH_SHAPES_ST2_LR];
    Word32 L_norm_factors[N_SCF_SEARCH_SHAPES_ST2_LR];

    Word32        *L_search_corr;
    Word32        *L_search_en;
    Word16        *xabs;
    Word8* scratch_top_proj;
    Word32 *L_corr_fixenv;
    Word32 *L_normcorr_fixenv;
    Word16 *gain_idx_opt_save;
    Word32 *L_min_mse_opt_save;
    );
    UNUSED(imax);  /* avoid gcc compiler warning */
    UNUSED(tmp);

    L_search_corr      = (Word32 *)scratchAlign(scratch, 0);
    L_search_en        = (Word32 *)scratchAlign(L_search_corr    , N_SCF_SEARCH_SHAPES_ST2_LR * sizeof(*L_search_corr));
    xabs               = (Word16 *)scratchAlign(L_search_en      , N_SCF_SEARCH_SHAPES_ST2_LR * sizeof(*L_search_en));
    scratch_top_proj   = (Word8 *)scratchAlign(xabs, PVQ_MAX_VEC_SIZE * sizeof(*xabs) );

    L_corr_fixenv      = (Word32 *)scratch_top_proj ;
    L_normcorr_fixenv  = (Word32 *)scratchAlign(L_corr_fixenv    , SNSLR_N_FIXENV*SNSLR_N_FIXENV_SHIFTS * sizeof(*L_corr_fixenv));
    gain_idx_opt_save  = (Word16 *)scratchAlign(L_normcorr_fixenv, SNSLR_N_FIXENV*SNSLR_N_FIXENV_SHIFTS * sizeof(*L_normcorr_fixenv));
    L_min_mse_opt_save = (Word32 *)scratchAlign(gain_idx_opt_save, N_SCF_SEARCH_SHAPES_ST2_LR * sizeof(*gain_idx_opt_save));
    /* scratch_top     = (Word8 *)scratchAlign(L_min_mse_opt_save, N_SCF_SEARCH_SHAPES_ST2_LR *sizeof(L_min_mse_opt_save)); */
   
    BASOP_sub_sub_start("pvq_fess_enc_search_fx");

    best_env_ind = -1;
    best_shift_ind = -1;
    gain_idx_opt = 0;

    y_split = &(y_Q0[0]); move32();
    y_full = &(y_Q0[1 * M]); move32();
    y_fix = &(y_Q0[2 * M]); move32();


    /* init */
    basop_memset(L_search_corr, 0, N_SCF_SEARCH_SHAPES_ST2_LR * sizeof(Word32));
    basop_memset(L_search_en, 0, N_SCF_SEARCH_SHAPES_ST2_LR * sizeof(Word32));

    max_xabs = -1; move16();
    max_xabsC = -1; move16();

    basop_memset(y_split, 0, PVQ_MAX_VEC_SIZE * sizeof(Word16));
    basop_memset(y_full, 0, PVQ_MAX_VEC_SIZE * sizeof(Word16));
    basop_memset(y_fix, 0, PVQ_MAX_VEC_SIZE * sizeof(Word16));

    y_splitAB = y_split;   /* ptr init */
    y_splitA0 = y_fix;    /* ptr init */  /*   y_fix temporarlily used for   splitA0 */

    xabs[0] = 0;  move16();  /* always no DC in LRSNS  scheme */
    FOR(i = SC; i < PVQ_MAX_VEC_SIZE; i++)
    {
        xabs[i] = abs_s(x_in[i]); move16();      /* input Qx */
    }

    /* A section */
    max_xabsA = -1; move16();
    L_xsumA = L_deposit_h(0);
    FOR(i = SC; i < SC + NA_LR; i++)
    {
        max_xabsA = s_max(max_xabsA, xabs[i]);     /* max value used for near optimal search correlation scaling */
        L_xsumA = L_mac0(L_xsumA, 1, xabs[i]);   /* stay in Qx */
    }
    /* B section */
    max_xabsB = -1; move16();
    L_xsumB = L_deposit_h(0);
    FOR(i = SC + NA_LR; i < SC + NA_LR + NB_LR; i++)
    {
        max_xabsB = s_max(max_xabsB, xabs[i]);       /* later used for  near optimal   search correlation scaling */
        L_xsumB = L_mac0(L_xsumB, 1, xabs[i]);      /* stay in Qx */
    }
    L_xsumAB = L_add(L_xsumA, L_xsumB);

    /* "C" section, 2 coeffs, not coded by split shapes {A,B}, {A,0} ,  but potentially coded by full and fix shapes */
    max_xabsC = s_max(xabs[PVQ_MAX_VEC_SIZE - 2], xabs[PVQ_MAX_VEC_SIZE - 1]);
    L_xsum = L_mac0(L_xsumAB, 1, xabs[PVQ_MAX_VEC_SIZE - 2]);
    L_xsum = L_mac0(L_xsum, 1, xabs[PVQ_MAX_VEC_SIZE - 1]);

    /*globalMax abs , for split and full search  and projection */
    max_xabsAB = s_max(max_xabsA, max_xabsB);
    max_xabs = s_max(max_xabsAB, max_xabsC); /* global max abs value over all M-1 coeffs */

    L_targetEnNeg = L_deposit_l(0);
    FOR(i = SC; i < PVQ_MAX_VEC_SIZE; i++)
    {
        /*target MSE part for final  winning stage 2 MSE shape */
        L_targetEnNeg = L_msu0(L_targetEnNeg, x_in[i], x_in[i]);  /* for final MSE output calculation  */
    }
 
    IF(L_xsum == 0)
    { /* no shape in any section, projection, A, AB, ABC   not possible::   -->  search is meaningless  */
        y_split[SC] = PULSES_SPLIT_A_LR;       move16();
        y_split[SC + NA_LR] = PULSES_SPLIT_B_LR; move16();

        y_full[SC] = PULSES_FULL_LR;     move16();

        FOR(i = SC;  i < SC + NSIGNS_FIX_3; i++)
        {
            y_fix[i] = 1;                 move16();
        }
        *fixShapeNbPtr = 3;    move16();

        /* L_search_corr[],   [0,1,2]  all stays zero valued  */
        ASSERT(L_search_corr[0] == 0 && L_search_corr[1] == 0 && L_search_corr[2] == 0);
        L_search_en[0] = PULSES_SPLIT_A_LR * PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR * PULSES_SPLIT_B_LR; move32();
        L_search_en[1] = PULSES_FULL_LR * PULSES_FULL_LR;  move32();
        L_search_en[2] = NSIGNS_FIX_3;  move32(); /* read from table as shape is not all ones */
    }
    ELSE
    {  /*   projection, to a valid  pyramid or sub-pyramid  */
       /* use  the most optimistic projection  */
#ifdef FIX_BASOP_ENC_LRSNS_ST2FULL_PROJ 
       /* always use at least one loop of single pulse optimization */
      pvq_pyr_project_lrsns_adv(SC + NFULL_LR, xabs, L_xsum, (PULSES_FULL_LR - 1), 0, y_full, &pulse_tot, &L_xy, &L_yy, scratch_top_proj);   /* 0 --> old floor( 1/( K-1) ) */
#else 
       pvq_pyr_project_lrsns_adv(SC + NFULL_LR, xabs, L_xsum, (PULSES_FULL_LR), PULSES_FULL_LR, y_full, &pulse_tot, &L_xy, &L_yy, scratch_top_proj);
#endif

        FOR(k = pulse_tot; k < PULSES_FULL_LR; k++)
        {
           L_yy = L_add(L_yy, 1);    /* pre add 1 in  Q0 in    L_yyQ0 = (x^2 + 2*x + 1)    */
           imax = one_pulse_search_tab_isqrt(SC, PVQ_MAX_VEC_SIZE, xabs, y_full, &pulse_tot, &L_xy, &L_yy, max_xabs);
        }
        ASSERT(pulse_tot == PULSES_FULL_LR);

        /* full  subshape   result vector in    y_full[0...15]  */
        /* L_mac  result from splitAB  search loop */
        L_search_corr[1] = L_shr_pos(L_xy,1);     /*  Q11*Q0   */
        L_search_en[1] = L_yy;  move32();        /*Q0*/

   
        IF( L_xsumA == 0 )
        {   /* no shape in A section,  search is meaningless  */
            y_split[SC] = PULSES_SPLIT_A_LR;         move16();   /*splitAB A-shape set  */
            y_split[SC + NA_LR] = PULSES_SPLIT_B_LR; move16();   /*splitAB B-shape set  */
            L_search_corr[0] = L_deposit_h(0);
            L_search_en[0] = (PULSES_SPLIT_A_LR*PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR + PULSES_SPLIT_B_LR); move32(); /*Q0*/
        }
        ELSE
        {   /* recalculate A section from the  partially available full result ( to remove section B pulses influence in  corr^2/en  )  */
            pulse_tot = 0; move16();
            L_xy = L_deposit_l(0);
            L_yy = L_deposit_l(0);

            basop_memmove(y_splitAB, y_full, (SC + NA_LR) * sizeof(Word16));     /* copy full(15) -> y_splitAB(5+8) (5 first here only)  */

            FOR(i = SC; i < (SC + NA_LR); i++)
            {
                      pulse_tot = add(pulse_tot, y_splitAB[i]);         /* Q0                                         */
                      L_xy = L_mac(L_xy, xabs[i], y_splitAB[i]);        /* Corr, Q0*Qx  +1 --> Qx+1                 */
                      L_yy = L_mac0(L_yy, y_splitAB[i], y_splitAB[i]);  /* Energy, Q(0+0)+1)= Q0 */
            }

            /* search any remaining unit pulses in   section A  */
            /* (continue up to 6 LF section A pulses)  */

            FOR(k = pulse_tot; k < PULSES_SPLIT_A_LR; k++)
            {
                L_yy = L_add(L_yy, 1); /* 1 added in Q0 */
                imax = one_pulse_search_tab_isqrt(SC, SC + NA_LR, xabs, y_splitAB, &pulse_tot, &L_xy, &L_yy, max_xabsA); 
            }

            /* L_mac  result from splitAB  search loop */
            L_xy_a6_mem = L_xy; move32();
            L_search_corr[0] = L_shr_pos(L_xy,1);   /* A-part of AB temporarily saved */
            L_search_en[0] = L_yy; move32();   /* A-part of AB  temporarily saved */

            ASSERT(pulse_tot == PULSES_SPLIT_A_LR);

            basop_memmove(y_splitA0, y_splitAB, (SC + NA_LR) * sizeof(Word16));        /* copy_splitAB(5 first)   -> y_splitA (5) */

            /* search any remaining unit pulses in  sectionA  for the PVQ(n=5,k=8)+PVQ(n=9,k=0)  alternative */
            /* (i.e. add 2 LF pulses, up to 8)  */

            FOR(k = pulse_tot; k < (PULSES_SPLIT_A_LR + PULSES_SPLIT_B_LR); k++)
            {
                L_yy = L_add(L_yy, 1); /* 1 added in Q0 */
                imax = one_pulse_search_tab_isqrt(SC, SC + NA_LR, xabs, y_splitA0, &pulse_tot, &L_xy, &L_yy, max_xabsA);
            }


            /* note: a local temporary  use of the fixed env location "2"  for corr and energy  memory*/
            L_search_corr[2] = L_xy; move32(); /* A-part of A,0 saved temporarily */
            L_search_en[2] = L_yy; move32();   /* A-part of A,0 saved */

            /* search B section of splitAB  (add 2 HF pulses) */
            L_xy = L_xy_a6_mem;            move32();
            L_yy = L_search_en[0];        move32();
            pulse_tot = PULSES_SPLIT_A_LR;  move16();

            IF(L_xsumB != 0)
            {
                FOR(k = 0; k < PULSES_SPLIT_B_LR; k++)
                {
                    L_yy = L_add(L_yy, 1); /* 1 added in Q0 */
                    imax = one_pulse_search_tab_isqrt(SC + NA_LR, SC + NA_LR + NB_LR, xabs, y_splitAB, &pulse_tot, &L_xy, &L_yy, max_xabsB);
                }
            }
            ELSE
            {
                 y_splitAB[SC + NA_LR] = PULSES_SPLIT_B_LR;   /* correlation L_xy not increased as x is zero */
                 L_yy = L_add(L_yy, PULSES_SPLIT_B_LR * PULSES_SPLIT_B_LR); /* y_splitAB energy increased by  2^2   */
            }

                /*   compare the best split shape out of (splitAB), vs  (splitA,0) */
                /* if   corrAB^2/enAB > corrA^2/enA,  then choose AB else A,0 */
                /* energy in L_yy can at most be 64 , and low value is  8     */
                /*  however   as energy can only take on ~64 values  [8 ... 64] */
                /*  we use :.  if   corrAB*sqrt(1/enAB) > corrA* sqrt(1/enA),  then choose AB else A,0  */
            {
                Word32 L_tmpA0, L_tmpAB;
                Word16 tmp_shift;

                /* maximize precision in corr values */
                tmp_shift = s_min(norm_l(L_search_corr[2]), norm_l(L_xy));

                ASSERT(tmp_shift >= 0);

                ASSERT(L_yy >= 4 && L_search_en[2] > 4 && "isqrt_Q16tab requires energy >=4 ");
                L_tmpA0 = Mpy_32_16_0_0(L_shl_pos(L_search_corr[2], tmp_shift), isqrt_Q16tab[L_search_en[2]]);  /*  nonsaturating Mpy32_16_0_0 used for ARM-SMULWB reuse */
                L_tmpAB = Mpy_32_16_0_0(L_shl_pos(L_xy, tmp_shift), isqrt_Q16tab[L_yy]);

                L_search_corr[0] = L_shr_pos_pos(L_xy, 1);   /* assume AB corr as winner */
                L_search_en[0] = L_yy;    move32();         /* assume AB energy as winner */

                IF(L_sub(L_tmpA0, L_tmpAB) > 0)
                {
                    basop_memcpy(y_split, y_splitA0, M * sizeof(Word16));   /* cpy complete  LFonly A,0  as the final splitLF   output   */
                    L_search_corr[0] = L_shr_pos_pos(L_search_corr[2], 1); /* set  A,0 as winner */
                    L_search_en[0] = L_search_en[2];   move32();     /* set  A,0 as winner */
                }
            }
        } /* L_xsumA!=0 */
    } /* L_xsum != 0 */

    {
        /* Fixed envelope signband coding,  including  sign coding  of a shifted  block  */
        /* subshape  idx 2,3,4     , "1"/env    1(s0)+ 2(shift)+ 11bits(s1..s11)   +   always a 3 bits gain  */
        /* subshape  idx 5         , "1"/env    1(s0)+ 2(shift)+  9bits(s1..s9)    +   always a 3 bits gain  */
        /*
         2,  init_bell_12signs ,         [ 8,8,8, 7,7 ... ]
         3,  decaying envelope 12 signs, [ 12,12,11,11, ... ]
         4,  start_bell_12signs ,        [ 7,7,8,8,8, 7,... ]
         5,  early_bell_10signs,         [ 6,6, 7,7,8,8,8,7,... ]*/

         /* find the minimum shape error across all possible 4 envelopes  and all 4 shifts */
         /* maximise normalized cross correlation  target*y /sqrt(en(y))   to minimize shape error */

        basop_memset(y_fix, 0, PVQ_MAX_VEC_SIZE * sizeof(Word16)); /* y_fix buffer was used for y_splitA0 above, reset needed */

        FOR(fix_ind = 0; fix_ind < SNSLR_N_FIXENV; fix_ind++)
        {
            FOR(shift_ind = 0; shift_ind < SNSLR_N_FIXENV_SHIFTS; shift_ind++)
            {
                L_corr = L_deposit_l(0);
                envPtr = &(lrsns_fix_env_fx[fix_ind][shift_ind]);     /* ptr init to start of the selected Q0 unit pulses  */
                xabs1 = &(xabs[SC + shift_ind]);                      /* ptr init */
                FOR(n = 0; n < lrsns_signs_fix_fx[fix_ind]; n++)
                {
                    ASSERT(envPtr[n] >= 0 && xabs1[n] >= 0);
                    L_corr = L_mac0(L_corr, xabs1[n], envPtr[n]);   /* runtime first octant correlation calc  */
                }

                tmp = add(shl_pos(fix_ind, 2), shift_ind);  /*  fix_ind * 4 + shift_ind */
                L_corr_fixenv[tmp] = L_corr;  move32();
                L_corr = L_shl_pos(L_corr, 9);  /* always 100% safe to shift up by 9 ,  tested with input  M* 32767*/
                /*    energy normalize, for the specific shift and env,
                      otherwise it is not possible to compare among the fixed shape envelopes  */
                L_normcorr_fixenv[tmp] = Mpy_32_16_0_0(L_corr, lrsns_fixenv_enNormQ19[tmp]);
            }
        }

        /* now actually pick the fixenv[0,1,2,3] and env shift[0,1,2,3] option which maximizes the normcorr
          . (also minimizes shape MSE as all fixed envelopes later use the very same 8 level gain quantizer) */

        Word32 *pL = &(L_normcorr_fixenv[0]); /* ptr init */

        best_ind = 0;  move16();
        FOR(n = 1; n < (SNSLR_N_FIXENV * SNSLR_N_FIXENV_SHIFTS); n++)
        {
            if (L_sub(pL[n], pL[best_ind]) > 0)
            {
                best_ind = n; move16();  /* conditional move single BASOP */
            }
        }

        /* decompose index */
        best_env_ind = shr_pos_pos(best_ind, 2);
        best_shift_ind = sub(best_ind, i_mult(best_env_ind, 4));

        ASSERT(best_env_ind >= 0 && best_env_ind < SNSLR_N_FIXENV);
        ASSERT(best_shift_ind >= 0 && best_shift_ind < SNSLR_N_FIXENV_SHIFTS);


        *fixShapeNbPtr = best_env_ind; move16();
        *fixShiftIdxPtr = best_shift_ind; move16();

        basop_memcpy(&(y_fix[add(SC, best_shift_ind)]), &(lrsns_fix_env_fx[best_env_ind][best_shift_ind]), lrsns_signs_fix_fx[best_env_ind] * sizeof(Word16));   /* cpy complete fixed envelope to  output   */

        tmp = add(shl_pos(best_env_ind, 2), best_shift_ind);
        L_search_corr[2] = L_corr_fixenv[tmp];  move32();         /* L_mac0 result from fixenv normcorr search loop*/
        L_search_en[2] = L_deposit_l(lrsns_fixenv_enQ0[tmp]);   /* tabled energy */
    }

    /* sign application for all three shape options    */
    /* apply sign of (x_in) to  first orthant result in Q0 */

    FOR(i = SC; i < (SC + NA_LR + NB_LR); i++)
    {
        if (x_in[i] < 0)
        {
            y_split[i] = negate(y_split[i]);  move16();  /* apply sign for split (AB or A0) */
        }
    }

    FOR(i = SC; i < (SC + NFULL_LR); i++)
    {
        if (x_in[i] < 0)
        {
            y_full[i] = negate(y_full[i]); move16(); /* apply sign for full  */
        }
    }

    FOR(i = SC; i < (SC + NFULL_LR); i++)
    {
        if (x_in[i] < 0)
        {
            y_fix[i] = negate(y_fix[i]);  move16(); /* apply sign for fixed env  */
        }
    }

    /*  prepare Q14 normalized values and run the gain quantization with all remaining 3 shape options  */
    /*  and prepare Q30 normalized values  for HR mode   */
    {

        /* calculate each normalized envelope  y_Q0[0...M-1]/sqrt(en) and provide result in Q14 */

        /* make sure Word16 inv_sqrt table lookup is ok */
        ASSERT(L_search_en[0] <= SQRT_EN_MAX_FX && L_search_en[1] <= SQRT_EN_MAX_FX);
        ASSERT(L_search_en[0] > 4 && L_search_en[1] > 4); /* inv_sqrt table only valid from  energy value 5  */


        norm_factors[0] = isqrt_Q16tab[L_search_en[0]];                            move16(); move16(); /* Q16 */
        norm_factors[1] = isqrt_Q16tab[L_search_en[1]];                            move16(); move16(); /* Q16 */
        tmp = add(shl_pos(best_env_ind, 2), best_shift_ind);
        norm_factors[2] = lrsns_fixenv_enNormQ19[tmp]; move16();                          /* Q19 */

        L_norm_factors[0] = isqrt_Q31tab[L_search_en[0]];                            move32(); move32(); /* Q31 */
        L_norm_factors[1] = isqrt_Q31tab[L_search_en[1]];                            move32(); move32(); /* Q31 */
        L_norm_factors[2] = L_lrsns_fixenv_enNormQ35[tmp];  move32();

        FOR(n = 0; n < (N_SCF_SEARCH_SHAPES_ST2_LR); n++)
        {
            /* PVQ and FESS unit energy normalization  */
            /*   some precomputed norm factor and pre-known normalization shift  */
            y = &(y_Q0[n*M]);                /* ptr init */
            L_y_norm = &(L_y_normQ30[n*M]);  /* ptr init */
            y_norm = &(y_normQ14[n*M]);      /* ptr init to signed seq Q14 */

            ASSERT(y[0] == 0);
            pvq_fess_dec_en1_normQ30andQ14_fx(
                y /*Q0*/,
                lrsns_y_up_bits[n],
                L_norm_factors[n],
                lrsns_norm_factorQ_L[n],
                M,           /*M==16 always used to enable simple unroll into 4 */
                L_y_norm,  /* for MSEQ22 estimation, and for ENABLE_HR  */
                y_norm);   /* for gain Q*/
            ASSERT(L_y_norm[0] == 0 && y_norm[0] == 0);

            L_y_norm[0] = 0; move32(); /* DC is always zero in LRSNS stage2 */
            y_norm[0] = 0; move16();


            /*  establish normalized  correlation to  target  */

            Word32 L_normcorr_pre_fx = Mpy_32_16(L_search_corr[n], norm_factors[n]);
            /* no need to recalc correlation,   simply scale with factor  */
            /*here Mpy_32_16  shifts left by 1 and maintains maximum precision */
            /*For ARM the somewhat better/native 32_16_0_0() without saturation,  however not used here    */

            Word16 g_shiftright_opt[3] = { 0 , 0 , 3 };  /*  down_shift(right) to get  result as  "g_opt*2"  in Q11 */
            Word32 L_normcorrX2_opt_fx = L_shr_pos_pos(L_normcorr_pre_fx, g_shiftright_opt[n]);   /* -2*g_opt for MSE calc , in Q11   */

            /*   MSE= t_2 + g^2*1.0 -(2*g_opt)*g , with g_opt = corr(t_2,x_n)  */

            gTabPtr = (const Word16 *)lrsns_vq_gainsQ12_fx[n];  move32(); /* q gains  in Q12,  n-adaptive ptr init    */
            L_mse = L_deposit_l(0);
            L_min_mse_opt = L_add(MAX_32, 0);   /* min_MSE w/o target energy gain  loop */

            FOR(gidx = 0; gidx < lrsns_vq_gain_lvls_fx[n]; gidx++)
            {
                /* adding at the same Q required, we use Q22 in search loop  for now */
                /* g_q*(g_q - 2*g_opt)/ */
                 /* (-2*g_opt)*g_q + g_q*g_q    =  g_q*(g_q-2*g_opt)  ,   Q11*Q11->Q22  */

                L_g_q_tmp = L_shr_pos(L_deposit_h(gTabPtr[gidx]), 2); /* from Q12 to  12+16-2-> Q26  */
                L_tmp = L_shl_pos(L_normcorrX2_opt_fx, 15);  /*11+15 -> 26 */
                ASSERT(L_normcorrX2_opt_fx <= 32767 && "L_normcorrX2_opt_fx <= 32767 "); /* should   fit in Word16  */

                L_tmp = L_sub(L_g_q_tmp, L_tmp); /* Q26-Q26  --> Q26 */

                L_mse = Mpy_32_16_0_0(L_tmp, gTabPtr[gidx]);  /* Q26+Q12 - 16 => 22  */

                L_min_mse_opt = L_min(L_min_mse_opt, L_mse);
                if (L_sub(L_mse, L_min_mse_opt) == 0)
                {
                    gain_idx_opt = gidx; move16();  /* single BASOP, conditional move  */
                }
            }
            gain_idx_opt_save[n] = gain_idx_opt;  move16();
            L_min_mse_opt_save[n] = L_min_mse_opt; move32();

        } /* shape n */

        /* get the best min(MSE) among the 3 shapes */
        shape_idx_opt = 0;  move16();
        L_min_mse_opt = L_min_mse_opt_save[shape_idx_opt]; move32();
        FOR(n = 1; n < (N_SCF_SEARCH_SHAPES_ST2_LR); n++)
        {
            L_min_mse_opt = L_min(L_min_mse_opt_save[n], L_min_mse_opt);
            if (L_sub(L_min_mse_opt, L_min_mse_opt_save[n]) == 0)
            {
                shape_idx_opt = n;  move16(); /*   single basop  */
            }
        }
        /* provide function output */
        *s_idxPtr = shape_idx_opt;                  move16();
        *g_idxPtr = gain_idx_opt_save[*s_idxPtr];  move16();
        *L_MSEQ22Ptr = L_sub(L_min_mse_opt, L_targetEnNeg);        /* search based MSE, add target energy to gain part of MSE (sum t*2) */
        *g_qvalQ12Ptr = lrsns_vq_gainsQ12_fx[*s_idxPtr][*g_idxPtr];

        ASSERT(L_min_mse_opt == L_min_mse_opt_save[*s_idxPtr]);

 

     /** Recalculate  winning stage2 MSE with high precision    **/
        L_MSEQ22_recalc = 0;  
        L_normcorrQy = 0;     /* INT32_MIN; */

        /* shift up  input signal  to the max with a  margin of 4 bits. */
        Word16 x_in_upshift = norm_s(max_xabs);

        /*  accumulate 15 values, would typically require ~4 bits margin , however we know only 12  can have a non-zero value  */
        x_in_upshift = s_max(sub(x_in_upshift, 4), 0);

        FOR(i = 1; i < M; i++)
        {
            L_tmp = Mpy_32_16_0_0(L_y_normQ30[*s_idxPtr*M + i], shl_pos(x_in[i], x_in_upshift)); /* signed*signed */
            L_normcorrQy = L_add(L_normcorrQy, L_tmp);
        }
        L_normcorrQy = L_shr_pos(L_normcorrQy, x_in_upshift); /* back to Q25 */ 


        /*  x_in^2  +  gq^2*1  - 2*gq*gopt  */
        L_MSEQ22_recalc = L_negate(L_targetEnNeg);                    
        L_tmp = L_shr_pos(L_mult0(*g_qvalQ12Ptr, *g_qvalQ12Ptr), 2);  
        L_MSEQ22_recalc = L_add(L_MSEQ22_recalc, L_tmp);
        L_tmp = L_shl_pos(Mpy_32_16_0_0(L_normcorrQy, *g_qvalQ12Ptr), (26 + 1) - 25);  /*  Qy+12-16+lshift = 22,  -->  lshift = 26-Qy , 2.0 factor  -> add 1 leftshift*/
                                                                                   
        L_MSEQ22_recalc = L_sub(L_MSEQ22_recalc, L_tmp);              

        *L_MSEQ22Ptr = L_max(L_MSEQ22_recalc, 0L);  move32();  /*use high precision result recalculated MSE for winner */

 

        ASSERT(*L_MSEQ22Ptr >= 0 && " warning negative total MSE ");


    } /*  shape norm and gain calc ,  and  shape decision */

    BASOP_sub_sub_end();

    Dyn_Mem_Deluxe_Out();
}
#endif

