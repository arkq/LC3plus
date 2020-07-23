/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"

void processPlcMain_fl(LC3_FLOAT *q_d_fl_c, LC3_FLOAT *syntM_fl_c, LC3_Dec* decoder, DecSetup* h_DecSetup, LC3_INT bfi
, PlcSetup *PlcSetup, PlcNsSetup *PlcNsSetup, LC3_INT plcMeth, LC3_INT ltpf_pitch_int, LC3_INT ltpf_pitch_fr, LC3_INT tilt, const LC3_INT *bands_offset, LC3_INT bands_number, const LC3_INT *bands_offsetPLC, LC3_INT n_bandsPLC
)
{
    
    
    UNUSED(n_bandsPLC);
    UNUSED(bands_offsetPLC);
    UNUSED(bands_offset);
    UNUSED(tilt);
    UNUSED(ltpf_pitch_fr);
    UNUSED(ltpf_pitch_int);
    UNUSED(plcMeth);
    UNUSED(syntM_fl_c);
    UNUSED(bands_number);

    if (bfi == 1)
    {
        PlcSetup->nbLostCmpt = PlcSetup->nbLostCmpt + 1;
    }
    
    if (bfi == 1)
    {
        
    
        switch (h_DecSetup->concealMethod)
        {
            case 0:
                processNoiseSubstitution0_fl(PlcSetup->q_d_prev, decoder->yLen, &PlcSetup->nbLostCmpt, &PlcNsSetup->cum_alpha, &PlcNsSetup->seed, q_d_fl_c);
                break;
            default:
                assert("Invalid PLC method!"); 
        }
    }
    
    if (bfi == 0)
    {
        processPlcUpdateSpec_fl(PlcSetup->q_d_prev, q_d_fl_c, decoder->yLen);
    }
}

