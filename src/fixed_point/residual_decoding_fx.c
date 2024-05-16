/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef ENABLE_HR_MODE /* HRMODE enables packing of residual bits */

void processResidualDecoding_fx(Word32 x[], Word16 x_e, Word16 L_spec, UWord8 prm[], Word16 resQBits
#ifdef ENABLE_HR_MODE
                                , Word16 hrmode
#endif
)
{

    Counter i;
    Word32  fac_m, fac_p;
    Word16  s, bits;
    Word32  tmp;
    Counter idx;
    Word16 N_nz = 0;
    Word16 iter;
    Word32 fac_hr;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualDecoding_fx", sizeof(struct {
                   Counter i;
                   Word32  fac_m, fac_p;
                   Word16  s, bits;
                   Word32  tmp;
                   Counter idx;
                   Word16 N_nz;
                   Word16 *nz_idx;
                   Word16 iter;
               }));
#endif

    tmp = 0;
    s   = sub(x_e, 1);
    s = s_min(s, 31);

    IF (hrmode)
    {
        fac_hr = L_shr(0x10000000, s); /* 0.25 in 1Q30 */
    }
    ELSE
    {
        fac_m = L_shr(0xC000000, s);  /* 0.1875 in 1Q30 */
        fac_p = L_shr(0x14000000, s); /* 0.3125 in 1Q30 */
    }

    bits = 0;
    move16();

     Word16 nz_idx[MAX_LEN];

    IF (hrmode)
    {
        FOR (i = 0; i < L_spec; i++)
        {
            IF (x[i])
            {
                nz_idx[N_nz] = i; move16();
                N_nz = add(N_nz, 1);
            }
        }
        FOR (iter = 0; iter < EXT_RES_ITER_MAX; iter++)
        {
            IF (sub(bits, resQBits) >= 0)
            {
                BREAK;
            }
            FOR (i = 0; i < N_nz; i++)
            {
                idx = nz_idx[i]; move16();

                IF (sub(bits, resQBits) >= 0)
                {
                    BREAK;
                }

                IF (! (s_and(prm[shr(bits, RESBITS_PACK_SHIFT)], shl(1, s_and(bits, RESBITS_PACK_MASK)))))
                {
                    tmp = L_sub_sat(x[idx], fac_hr);
                }
                ELSE
                {
                    tmp = L_add_sat(x[idx], fac_hr);
                }
                x[idx] = tmp;
                move32();
                bits = add(bits, 1);
            }
            fac_hr = L_shr(fac_hr, 1);
        }
    }
    ELSE
    {
        FOR (i = 0; i < L_spec; i++)
        {
            IF (sub(bits, resQBits) >= 0)
            {
                BREAK;
            }

            IF (x[i] != 0)
            {
                IF (! (s_and(prm[shr(bits, RESBITS_PACK_SHIFT)], shl(1, s_and(bits, RESBITS_PACK_MASK)))))
                {
                    if (x[i] > 0)
                        tmp = L_sub(x[i], fac_m);
                    if (x[i] < 0)
                        tmp = L_sub(x[i], fac_p);
                }
                ELSE
                {
                    if (x[i] > 0)
                        tmp = L_add(x[i], fac_p);
                    if (x[i] < 0)
                        tmp = L_add(x[i], fac_m);
                }
                x[i] = tmp;
                move32();
                bits = add(bits, 1);
            }
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#else

void processResidualDecoding_fx(Word32 x[], Word16 x_e, Word16 L_spec, UWord8 prm[], Word16 resQBits
)
{

    Counter i;
    Word32  fac_m, fac_p;
    Word16  s, bits;
    Word32  tmp;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualDecoding_fx", sizeof(struct {
                   Counter i;
                   Word32  fac_m, fac_p;
                   Word16  s, bits;
                   Word32  tmp;
               }));
#endif

    tmp = 0;
    s   = sub(x_e, 1);
    s = s_min(s, 31);

    {
        fac_m = L_shr(0xC000000, s);  /* 0.1875 in 1Q30 */
        fac_p = L_shr(0x14000000, s); /* 0.3125 in 1Q30 */
    }

    bits = 0;
    move16();

    {
        FOR (i = 0; i < L_spec; i++)
        {
            IF (sub(bits, resQBits) >= 0)
            {
                BREAK;
            }

            IF (x[i] != 0)
            {
                IF (prm[bits] == 0)
                {
                    if (x[i] > 0)
                        tmp = L_sub(x[i], fac_m);
                    if (x[i] < 0)
                        tmp = L_sub(x[i], fac_p);
                }
                ELSE
                {
                    if (x[i] > 0)
                        tmp = L_add(x[i], fac_p);
                    if (x[i] < 0)
                        tmp = L_add(x[i], fac_m);
                }
                x[i] = tmp;
                move32();
                bits = add(bits, 1);
            }
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#endif
