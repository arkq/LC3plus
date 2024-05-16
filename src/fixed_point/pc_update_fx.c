/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "defines.h"
#include "constants.h"
#include "functions.h"

#ifdef ENABLE_HR_MODE

void Copy_Scale_sig_32(const Word32 x[], /* i  : signal to scale input           Qx        */
    Word16       y[], /* o  : scaled signal output            Qx        */
    const Word16 lg,  /* i  : size of x[]                     Q0        */
    const Word16 exp0 /* i  : exponent: x = round(x << exp)   Qx ?exp  */
)
{
    Counter i;
    Word16 tmp;

    if (exp0 == 0)
    {
        for (i = 0; i < lg; i++)
        {
            y[i] = extract_h(x[i]);
        }

        return;
    }

    tmp = s_max(exp0, -31);
    tmp = s_min(tmp, 31);

    for (i = 0; i < lg; i++)
    {
        y[i] = extract_h(L_shr_r(x[i], -tmp));
    }
}
#endif

void processPCupdate_fx(Word16 bfi, Word16 yLen, Word16 q_old_res_fx[], Word16 *q_old_res_fx_exp,
#ifdef ENABLE_HR_MODE
                        Word32 q_res_fx[],
#else
                        Word16 q_res_fx[],
#endif
                        Word16 spec_inv_idx, Word16 gg_idx, Word16 gg_idx_off,
                        Word16 *prev_gg, Word16 *prev_gg_e, Word16 rframe, Word16 *BW_cutoff_idx_nf,
                        Word16 *prev_BW_cutoff_idx_nf, Word16 fac_ns_idx, Word16 *prev_fac_ns_fx, Word16 fac,
                        Word16 fac_e)
{
    Word16  global_gain, global_gain_e, s, s2, s3, tmp16;
    Word32  tmp32;

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Word16  global_gain, global_gain_e, s, s2, s3, tmp16;
        Word32  tmp32;
    };
    Dyn_Mem_In("processPCupdate_fx", sizeof(struct _dynmem));
#endif

    tmp32         = L_shl_pos(L_mult0(add(gg_idx, gg_idx_off), 0x797D), 7);
    global_gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1);
    global_gain = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000)));

    *prev_gg   = global_gain;  move16();
    *prev_gg_e = global_gain_e;  move16();

#ifdef ENABLE_HR_MODE
    s = getScaleFactor32(q_res_fx, spec_inv_idx); /* exp = 0 */
#else
    s = getScaleFactor16(q_res_fx, spec_inv_idx); /* exp = 0 */
#endif
    
    IF (bfi == 0)
    {
        *q_old_res_fx_exp = negate(s);
#ifdef ENABLE_HR_MODE
        Copy_Scale_sig_32(q_res_fx, q_old_res_fx, yLen, s);
        *q_old_res_fx_exp = *q_old_res_fx_exp + 16;
#else
        Copy_Scale_sig(q_res_fx, q_old_res_fx, yLen, s);
#endif
    }
    ELSE
    {
#ifdef ENABLE_HR_MODE
        s2 = getScaleFactor32(&q_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx)); /* exp = q_old_res_fx_exp */
        IF (s2 == 0)
            s2 = 16;
        s3 = add(s, *q_old_res_fx_exp);
        IF (sub(s3, s2) > 0)
        {
            tmp16 = sub(s3, s2);
            s     = sub(s, tmp16);
        }
        s2                = add(s, *q_old_res_fx_exp);
        *q_old_res_fx_exp = negate(s) + 16;


        if (s2 > -32)
        {
            Copy_Scale_sig_32(&q_res_fx[spec_inv_idx], &q_old_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx), s2);
        }
#else
        s2 = getScaleFactor16(&q_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx)); /* exp = q_old_res_fx_exp */
        s3 = add(s, *q_old_res_fx_exp);
        IF (sub(s3, s2) > 0)
        {
            tmp16 = sub(s3, s2);
            s = sub(s, tmp16);
        }
        s2 = add(s, *q_old_res_fx_exp);
        *q_old_res_fx_exp = negate(s);

        s = s_max(s_min(s, 15), -15);
        s2 = s_max(s_min(s2, 15), -15);
        Copy_Scale_sig(q_res_fx, q_old_res_fx, spec_inv_idx, s);
        Copy_Scale_sig(&q_res_fx[spec_inv_idx], &q_old_res_fx[spec_inv_idx], sub(yLen, spec_inv_idx), s2);
#endif
    }

    IF (rframe == 0)
    {
        *prev_BW_cutoff_idx_nf = *BW_cutoff_idx_nf;
        *prev_fac_ns_fx = shl_pos(sub(8, fac_ns_idx), 11);
    }
    ELSE IF(sub(bfi, 2) == 0 && sub(*BW_cutoff_idx_nf, *prev_BW_cutoff_idx_nf) != 0
            && sub(spec_inv_idx, yLen) < 0)
    {
        *BW_cutoff_idx_nf = *prev_BW_cutoff_idx_nf;
        *prev_fac_ns_fx = shl_sat(mult(*prev_fac_ns_fx, fac), fac_e);
        *prev_fac_ns_fx = s_max(*prev_fac_ns_fx, 2048);
        *prev_fac_ns_fx = s_min(*prev_fac_ns_fx, 16384);
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


