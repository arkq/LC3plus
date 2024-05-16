/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPlcDampingScramblingMain_fl(LC3_INT32 *ns_seed,
                                        LC3_INT32 *pc_seed, LC3_INT32 ns_nbLostCmpt_pc,
                                        LC3_INT32 ns_nbLostCmpt, LC3_FLOAT *stabFac, LC3_FLOAT *cum_fading_slow, LC3_FLOAT *cum_fading_fast,
                                        LC3_FLOAT *spec_prev, LC3_FLOAT *spec, LC3_INT32 spec_inv_idx, LC3_INT32 yLen, LC3_INT32 bfi,
                                        LC3_INT32 frame_dms, LC3_INT32 concealMethod, LC3_INT32 pitch_present_bfi1, LC3_INT32 pitch_present_bfi2,
                                        LC3_FLOAT *cum_fflcAtten
                                        , LC3_UINT8 plc_fadeout_type
                                        )
{

        LC3_INT32 processDampScramb;

    processDampScramb = 0;
    

    if ( bfi != 0 )
    {
        if (concealMethod == 4 || bfi == 2)
        {
            processDampScramb = 1;
        }
        
        if (ns_nbLostCmpt == 1)
        {
            *cum_fading_slow = 1;
            *cum_fading_fast = 1;
            *cum_fflcAtten   = 1;
        }     
        
        if ( bfi == 1 )
        {
            processPlcDampingScrambling_fl(spec, yLen, ns_nbLostCmpt, stabFac, processDampScramb, cum_fflcAtten,
                                           pitch_present_bfi1, frame_dms, cum_fading_slow, cum_fading_fast, ns_seed, 0
                                           , plc_fadeout_type
                                          );
        }
        else /* bfi == 2 */
        {
            processPlcDampingScrambling_fl(spec, yLen, ns_nbLostCmpt_pc, stabFac, processDampScramb, cum_fflcAtten,
                                           pitch_present_bfi2, frame_dms, cum_fading_slow, cum_fading_fast, pc_seed, spec_inv_idx
                                           , plc_fadeout_type
                                          );
            processPlcUpdateSpec_fl(spec_prev, spec, yLen);
        }
    }
}

void processPlcDampingScrambling_fl(LC3_FLOAT *spec, LC3_INT32 yLen, LC3_INT32 nbLostCmpt, LC3_FLOAT *stabFac, LC3_INT32 processDampScramb,
                            LC3_FLOAT *cum_fflcAtten, LC3_INT32 pitch_present, LC3_INT32 frame_dms, LC3_FLOAT *cum_fading_slow,
                            LC3_FLOAT *cum_fading_fast, LC3_INT32 *seed, LC3_INT32 spec_inv_idx
                            , LC3_UINT8 plc_fadeout_type
                            )
{
    LC3_INT32 plc_start_inFrames, plc_end_inFrames, plc_duration_inFrames, x, b, i, ad_ThreshFac_start;
    LC3_FLOAT slow, fast, linFuncStartStop, randThreshold, ad_ThreshFac_end, ad_threshFac, frame_energy, mean_energy, energThreshold, fac, m, n, fflcAtten, cum_fading_slow_local, cum_fading_fast_local;
    
    frame_energy = 0;

    slow = 0.8 + 0.2 * (*stabFac);
    fast = 0.3 + 0.2 * (*stabFac);

    switch (frame_dms)
    {
    case 25:
        slow = LC3_SQRT(LC3_SQRT(slow));
        fast = LC3_SQRT(LC3_SQRT(fast));
        break;
    case 50:
        slow = LC3_SQRT(slow);
        fast = LC3_SQRT(fast);
        break;
    case 75:
        slow = LC3_SQRT(LC3_SQRT(slow*slow*slow));
        fast = LC3_SQRT(LC3_SQRT(fast*fast*fast));
        break;
    }

    if (plc_fadeout_type == 0)
    {
        *cum_fading_slow = *cum_fading_slow * slow;
        *cum_fading_fast = *cum_fading_fast * fast;
    }
    
	if (processDampScramb == 1)
	{
		if (plc_fadeout_type != 0)
		{
			if (nbLostCmpt < (4 * (100.0 / (LC3_FLOAT)frame_dms))) {
				cum_fading_slow_local = 1.0;
			}
			else if (nbLostCmpt < (8 * (100.0 / (LC3_FLOAT)frame_dms))) {
				cum_fading_slow_local = 0.9;
			}
			else {
				cum_fading_slow_local = 0.85;
			}

			*cum_fading_slow = *cum_fading_slow * cum_fading_slow_local;
			cum_fading_slow_local = *cum_fading_slow;		 
		}
		else {
			fflcAtten = 1;
			cum_fading_slow_local = *cum_fading_slow;
			cum_fading_fast_local = *cum_fading_fast;

			if (spec_inv_idx == 0)
			{
				if (nbLostCmpt * frame_dms > PLC_FADEOUT_IN_MS * 10)
				{
					fflcAtten = 0;
					*cum_fflcAtten = 0;
				}
				else if (nbLostCmpt * frame_dms > 200)
				{
					switch (frame_dms)
					{
					case  25: fflcAtten = PLC34_ATTEN_FAC_025; break;
					case  50: fflcAtten = PLC34_ATTEN_FAC_050; break;
                    case  75: fflcAtten = PLC34_ATTEN_FAC_075; break;
					case 100: fflcAtten = PLC34_ATTEN_FAC_100; break;
					}
				}


				*cum_fflcAtten = *cum_fflcAtten * fflcAtten;
				cum_fading_slow_local = *cum_fading_slow * *cum_fflcAtten;
				cum_fading_fast_local = *cum_fading_fast * *cum_fflcAtten;
			}

			if (pitch_present == 0)
			{
				plc_start_inFrames = 1;
			}
			else {
				plc_start_inFrames = floor(PLC4_TRANSIT_START_IN_MS / (frame_dms / 10.0));
			}

			plc_end_inFrames = floor(PLC4_TRANSIT_END_IN_MS / (frame_dms / 10.0));
			plc_duration_inFrames = plc_end_inFrames - plc_start_inFrames;

			if (nbLostCmpt <= plc_start_inFrames)
			{
				linFuncStartStop = 1;
			}
			else if (nbLostCmpt >= plc_end_inFrames)
			{
				linFuncStartStop = 0;
			}
			else {
				x = nbLostCmpt;
				m = -1.0 / plc_duration_inFrames;
				b = -plc_end_inFrames;
				linFuncStartStop = m * (x + b);
			}

			randThreshold = -32768 * linFuncStartStop;
		}
        
        for (i = spec_inv_idx; i < yLen; i++)
        {
            *seed = 16831 + *seed * 12821;

            *seed = (LC3_INT16)(*seed);

            if (*seed < 0)
            {
                if (plc_fadeout_type != 0 || pitch_present == 0 || *seed < randThreshold )
                {
                    spec[i] = -spec[i];
                }
            }
        }

        if (plc_fadeout_type == 0)
        {
            ad_ThreshFac_start = 10;
            ad_ThreshFac_end   = 1.2;
            ad_threshFac       = (ad_ThreshFac_start - ad_ThreshFac_end) * linFuncStartStop + ad_ThreshFac_end;

            if (spec_inv_idx < yLen)
            {
                for (i = spec_inv_idx; i < yLen; i++)
                {
                    frame_energy = frame_energy + (spec[i] * spec[i]);
                }

                mean_energy = frame_energy * 1 / (yLen - spec_inv_idx);
            }
            else
            {
                mean_energy = 0;
            }

            energThreshold = LC3_SQRT(ad_threshFac * mean_energy);
            fac = (cum_fading_slow_local - cum_fading_fast_local) * energThreshold;
        }

        for (i = spec_inv_idx; i < yLen; i++)
        {
            if (plc_fadeout_type != 0 || LC3_FABS(spec[i]) < energThreshold  )
            {
                m = cum_fading_slow_local;
                n = 0;
            }
            else
            {
                m = cum_fading_fast_local;

                if (spec[i] > 0)
                {
                    n = fac;
                }
                else if (spec[i] == 0)
                {
                    n = 0;
                }
                else
                {
                    n = -fac;
                }
            }

            spec[i] = m * spec[i] + n;
        }
    }
}

