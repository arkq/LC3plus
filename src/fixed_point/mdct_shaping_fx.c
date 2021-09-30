/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"



void processMdctShaping_fx(Word32 x[],
#ifdef ENABLE_HR_MODE
                           Word32 scf[], 
#else
                           Word16 scf[],
#endif
                           Word16 scf_exp[], const Word16 bands_offset[], Word16 fdns_npts)
{
    Counter i, j;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processMdctShaping_fx", sizeof(struct { Counter i, j; }));
#endif

    j = 0; move16();
    FOR (i = 0; i < fdns_npts; i++)
    {
        FOR (; j < bands_offset[i + 1]; j++)
        {
#ifdef ENABLE_HR_MODE
            x[j] = L_shl(Mpy_32_32(x[j], scf[i]), scf_exp[i]); move32();
#else
            x[j] = L_shl(Mpy_32_16(x[j], scf[i]), scf_exp[i]); move32();
#endif
        }
    }

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}


void processScfScaling(Word16 scf_exp[], Word16 fdns_npts, Word16 *x_e)
{
    Counter i;
    Word16  scf_exp_max;

#ifdef DYNMEM_COUNT
    Dyn_Mem_In("processLpcGainScaling", sizeof(struct {
                   Counter i;
                   Word16  scf_exp_max;
               }));
#endif

    scf_exp_max = scf_exp[0]; move16();

    FOR (i = 1; i < fdns_npts; i++)
    {
        scf_exp_max = s_max(scf_exp_max, scf_exp[i]);
    }

    FOR (i = 0; i < fdns_npts; i++)
    {
        scf_exp[i] = sub(scf_exp[i], scf_exp_max); move16();
    }

    *x_e = add(*x_e, scf_exp_max); move16();
#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

