/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef ENABLE_HR_MODE

void processResidualCoding_fx(Word16 x_e, Word32 x[], 
                              Word32 xq[], Word32 gain,
                              Word16 gain_e, Word16 L_spec,
                              Word16 targetBits, Word16 nBits, UWord8 *resBits, Word16 *numResBits
                              , Word16 hrmode
)
{

    Counter i;
    Word16  s, n, m;
    Word32  L_tmp;
    Word16  iter = 0;
    Counter idx;
    Word16 N_nz = 0;
    Word16 n1, n2;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualCoding_fx", sizeof(struct {
                   Counter i;
                   Word16  s, n, m;
                   Word32  L_tmp;
                   Word16  iter;
                   Counter idx;
                   Word16 N_nz;
                   Word16 *nz_idx;
                   Word16 n1, n2;
               }));
#endif

    n = 0;
    move16();
    
    IF (hrmode)
    {
        m = add(sub(targetBits, nBits), 14);
        assert(m <= (MAX_RESBITS_LEN << RESBITS_PACK_SHIFT));
    }
    ELSE
    {
        m = add(sub(targetBits, nBits), 4);
        if (m > L_spec)
        {
            m = L_spec;
        }
    }

    basop_memset(resBits, 0, sizeof(*resBits)*((m + RESBITS_PACK_MASK) >> RESBITS_PACK_SHIFT));

    s = sub(add(15, gain_e), x_e);

    IF(hrmode)
    {
        /* 
        x = xval * 2^(31 - x_e)
        gain = gainval * 2^(15 - gain_e)

        To bring gain to the same q of x :

        gain (in_Qx_e) = gain * 2^(16 + gain_e - x_e)

        gain*offset (in Qx_e) = gain (in_Qx_e) * 0.25 = gain (in_Qx_e) / 4 = gain * 2^(16 - 2 + gain_e - x_e)
        */
        Word16 shift_val;
            shift_val = sub(sub(gain_e, x_e), 2);

        Word32 gain_offset;

        IF (shift_val <= -32)
            gain_offset = 0;
        ELSE
        {
            gain_offset = L_shl_sat(gain, shift_val);
        }
        
        Word16 exit_iter   = 0;
        
        Word16 nz_idx[MAX_LEN];

        /* enumerate non-zero coefficients */
        FOR (i = 0; i < L_spec; i++)
        {
            IF (xq[i])
            {
                nz_idx[N_nz] = i; move16();
                N_nz = add(N_nz, 1);
            }
        }

        s = sub(31, sub(x_e, gain_e));
        FOR (iter = 0; iter < EXT_RES_ITER_MAX; iter++)
        {
            FOR (i = 0; i < N_nz; i++)
            {
                idx = nz_idx[i];

                L_tmp = L_sub(x[idx], Mpy_32_32(L_shl(xq[idx], s), gain));

                IF (L_tmp >= 0)
                {
                    n1          = shr(n, RESBITS_PACK_SHIFT);
                    n2          = s_and(n, RESBITS_PACK_MASK);
                    resBits[n1] = (UWord8)s_or(resBits[n1], shl(1, n2));
                    move16();
                    x[idx] = L_sub(x[idx], gain_offset);
                    move16();
                }
                ELSE
                {
                    move16();
                    x[idx] = L_add(x[idx], gain_offset);
                    move16();
                }
                n = add(n, 1);
                IF (sub(n, m) == 0)
                {
                    exit_iter = 1;
                    BREAK;
                }
            }
            gain_offset = L_shr(gain_offset, 1); /* offset *= 0.5 */

            IF (exit_iter)
            {
                BREAK;
            }
        }
    }
    ELSE
    {
        s = add(s, 16);
        gain = round_fx(gain);

        FOR (i = 0; i < L_spec; i++)
        {
            IF (xq[i] != 0)
            {
                L_tmp = L_sub(x[i], Mpy_32_16(L_shl(xq[i], s), gain));
                if (L_tmp >= 0)
                {
                    n1 = shr(n, RESBITS_PACK_SHIFT);
                    n2 = s_and(n, RESBITS_PACK_MASK);
                    resBits[n1] = (UWord8) s_or(resBits[n1], shl(1, n2));
                    move16();
                }
                n = add(n, 1);
                IF (sub(n, m) == 0)
                {
                    BREAK;
                }
            }
        }
    }
    *numResBits = n;
    move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#else

void processResidualCoding_fx(Word16 x_e, Word32 x[], 
                              Word16 xq[], 
                              Word16 gain, Word16 gain_e, Word16 L_spec,
                              Word16 targetBits, Word16 nBits, UWord8 *resBits, Word16 *numResBits
)
{

    Counter i;
    Word16  s, n, m;
    Word32  L_tmp;
#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processResidualCoding_fx", sizeof(struct {
                   Counter i;
                   Word16  s, n, m;
                   Word32  L_tmp;
               }));
#endif

    n = 0;
    move16();
    {
        m = add(sub(targetBits, nBits), 4);
        if (m > L_spec)
        {
            m = L_spec;
        }
    }


    s = sub(add(15, gain_e), x_e);
    {
        FOR (i = 0; i < L_spec; i++)
        {
            IF (xq[i] != 0)
            {
                L_tmp = L_sub(x[i], L_shl(L_mult(xq[i], gain), s));
                if (L_tmp < 0)
                {
                    resBits[n] = 0;
                    move16();
                }
                if (L_tmp >= 0)
                {
                    resBits[n] = 1;
                    move16();
                }
                n = add(n, 1);
                IF (sub(n, m) == 0)
                {
                    BREAK;
                }
            }
        }
    }
    *numResBits = n;
    move16();

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

#endif
