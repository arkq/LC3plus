/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static LC3_INT32 change_bit_at_position(LC3_INT32 value, LC3_UINT8 bit_position, LC3_INT8 bit)
{
    LC3_INT32 helper_mask = ~(1 << bit_position);
    LC3_INT32 tmp = value & helper_mask;
    tmp = tmp | (bit << bit_position);
    return tmp;
}

static void update_bit_and_byte_positions(LC3_INT16 longterm_analysis_counter_max_bytebuffer, LC3_INT8 *byte_position, LC3_INT8 *bit_position)
{
    if (*bit_position == 29)
    {
        *bit_position = 0;
        
        if ((*byte_position - longterm_analysis_counter_max_bytebuffer) < -1)
        {
            *byte_position = *byte_position + 1;
        } else {
            *byte_position = 0;
        }
    } else {
        *bit_position = *bit_position + 1;
    }
}

static void array_insert_and_shift(LC3_INT32 *array, LC3_UINT8 value, LC3_INT16 longterm_analysis_counter_max, LC3_INT16 *overall_counter, LC3_INT8 *byte_position, LC3_INT8 *bit_position)
{
    LC3_INT32 current_byte = 0;
    
    if (overall_counter != NULL) 
    {
        *overall_counter = MIN(*overall_counter + 1, longterm_analysis_counter_max);
    }
    
    current_byte = array[*byte_position];
    
    current_byte = change_bit_at_position(current_byte, *bit_position, value);
    
    array[*byte_position] = current_byte;
}

static void array_calculate(LC3_INT32 *array_tdc, LC3_INT32 *array_ns, int length, LC3_INT16 *counter_tdc, LC3_INT16 *counter_ns, LC3_INT16 longterm_analysis_counter_max)
{
    int i, k;
    LC3_INT32 current_byte_tdc = 0, current_byte_ns = 0;
    LC3_INT16 counter_loc_tdc = 0, counter_loc_ns = 0, counter_tmp = 0;

    for (i = length - 1; i >= 0; i--)
    {
        current_byte_tdc = array_tdc[i];
        current_byte_ns = array_ns[i];
        
        for (k = 0; k < 30; k++)
        {
            counter_loc_tdc += ((current_byte_tdc >> k) & 1);
            counter_loc_ns += ((current_byte_ns >> k) & 1);
            counter_tmp++;
            
            /* Break from both loops if full 2s buffer has been evaluated */
            if (counter_tmp >= longterm_analysis_counter_max)
            {
                i = -1;
                k = 30;
                break;
            }
        }
    }
    
    *counter_tdc = counter_loc_tdc;
    *counter_ns = counter_loc_ns;
}

static void plc_xcorr_lc(LC3_FLOAT *pcmbufHist, LC3_INT32 max_len_pcm_plc, LC3_INT32 pitch_int, LC3_INT32 framelength, LC3_INT32 frame_dms, LC3_INT32 fs, LC3_FLOAT *xcorr);
static void spectral_centroid_lc(LC3_FLOAT *gains, LC3_INT32 tilt, const LC3_INT *bands_offset, LC3_INT32 bands_number, LC3_INT32 framelength, LC3_INT32 fs, LC3_FLOAT *sc);

void processPlcClassify_fl(LC3_INT plcMeth, LC3_INT *concealMethod, LC3_INT32 *nbLostCmpt, LC3_INT32 bfi,
                           LC3_FLOAT *xcorr, LC3_INT32 framelength, LC3_INT32 frame_dms, LC3_INT32 pitch_int,
                           LC3_INT32 fs, const LC3_INT *band_offsets, LC3_INT32 bands_number, LC3_INT32 tilt, PlcAdvSetup *plcAd
                           , LC3_INT32 hrmode
)
{
        LC3_FLOAT sc, class;
        int fs_idx_tmp;

    if (plcAd)
    {
        *xcorr = 0;
    }
    
    fs_idx_tmp = FS2FS_IDX(fs);
    /* Save statistics for 24 kHz, 48 kHz and 96 kHz */
    if ((bfi == 1) || ((bfi >= 0) && (bfi <= 2) && ((fs_idx_tmp == 2) || (fs_idx_tmp == 4) || (fs_idx_tmp == 5))))   /* Partial Concealment PC(bfi==2) requires allowing value  2 to pass thru  as well */
    {
        if (bfi == 1)
        {
            *nbLostCmpt = *nbLostCmpt + 1;
        }
        
        /* Use pitch correlation at ltpf integer lag if available */
        if ((*nbLostCmpt == 1) || (bfi != 1) )/* PC(bfi==2) requires allowing 2 to pass thru  as well */
        {
            *concealMethod = 4;  /* Noise Substitution */
            UNUSED(plcMeth);

            /* Advanced PLC */
            if (pitch_int > 0)
            {
                *concealMethod = 3; /* Timedomain PLC assumed */
                plc_xcorr_lc(plcAd->pcmbufHist, plcAd->max_len_pcm_plc, pitch_int, framelength, frame_dms, fs, xcorr);

                spectral_centroid_lc(plcAd->scf_q_old, tilt, band_offsets, bands_number, framelength, fs, &sc);
                class = *xcorr * 7640.0 / 32768.0 - sc - 5112.0 / 32768.0;
                
                if (class <= 0)
                {
                    if (frame_dms == 100 && hrmode == 0)
                    {
                        *concealMethod = 2; /* PhaseEcu selected */
                        array_insert_and_shift(plcAd->plc_longterm_advc_tdc, 0, plcAd->longterm_analysis_counter_max, &plcAd->overall_counter, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                        array_insert_and_shift(plcAd->plc_longterm_advc_ns, 0, plcAd->longterm_analysis_counter_max, NULL, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                    }
                    else
                    {
                        array_insert_and_shift(plcAd->plc_longterm_advc_tdc, 0, plcAd->longterm_analysis_counter_max, &plcAd->overall_counter, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                        array_insert_and_shift(plcAd->plc_longterm_advc_ns, 0, plcAd->longterm_analysis_counter_max, NULL, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                    }
                }
                else {
                    array_insert_and_shift(plcAd->plc_longterm_advc_tdc, 1, plcAd->longterm_analysis_counter_max, &plcAd->overall_counter, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                    array_insert_and_shift(plcAd->plc_longterm_advc_ns, 0, plcAd->longterm_analysis_counter_max, NULL, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                }
            }
            else
            {
                *concealMethod = 4; /* Noise Substitution */
                array_insert_and_shift(plcAd->plc_longterm_advc_tdc, 0, plcAd->longterm_analysis_counter_max, &plcAd->overall_counter, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
                array_insert_and_shift(plcAd->plc_longterm_advc_ns, 1, plcAd->longterm_analysis_counter_max, NULL, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
            }
            
            array_calculate(plcAd->plc_longterm_advc_tdc, plcAd->plc_longterm_advc_ns, plcAd->longterm_analysis_counter_max_bytebuffer, &plcAd->longterm_counter_plcTdc, &plcAd->longterm_counter_plcNsAdv, plcAd->longterm_analysis_counter_max);
            update_bit_and_byte_positions(plcAd->longterm_analysis_counter_max_bytebuffer, &plcAd->longterm_counter_byte_position, &plcAd->longterm_counter_bit_position);
        }
    }
}

static void spectral_centroid_lc(LC3_FLOAT *gains, LC3_INT32 tilt, const LC3_INT *band_offsets, LC3_INT32 bands_number, LC3_INT32 framelength, LC3_INT32 fs, LC3_FLOAT *sc)
{
        LC3_FLOAT gains_lin[M], gains_dee[M], numerator, denumerator;
        LC3_INT32 i, j, sum, len, start, stop;
        LC3_INT band_offsets_local[MAX_BANDS_NUMBER + 1];
    
    numerator = 0;

    for (i = 0; i < M; i++)
    {
        gains_lin[i] = LC3_POW(2, gains[i]);
    }
    
    for (i = 0; i < M; i++)
    {
        gains_dee[i] = gains_lin[i] / LC3_POW(10, i * (LC3_FLOAT) tilt / (LC3_FLOAT) (M - 1) / 10.0);
    }
    
    if (bands_number == 64)
    {
        memmove(band_offsets_local, band_offsets, (bands_number + 1) * sizeof(LC3_INT));
    }

    if (bands_number < 32)
    {
        band_offsets_local[0] = 0;
        j = 32 - bands_number;
        for (i = bands_number - 1; i >= j; i--)
        {
            band_offsets_local[(i + j) * 2 + 1 + 1] = band_offsets[i + 1];
            band_offsets_local[(i + j) * 2 + 0 + 1] = band_offsets[i + 1];
        }
        for (i = j - 1; i >= 0; i--)
        {
            band_offsets_local[i * 4 + 3 + 1] = band_offsets[i + 1];
            band_offsets_local[i * 4 + 2 + 1] = band_offsets[i + 1];
            band_offsets_local[i * 4 + 1 + 1] = band_offsets[i + 1];
            band_offsets_local[i * 4 + 0 + 1] = band_offsets[i + 1];
        }
    }
    else
    if (bands_number < 64)
    {
        band_offsets_local[0] = 0;
        j = 64 - bands_number;
        for (i = bands_number - 1; i >= j; i--)
        {
            band_offsets_local[i + j + 1] = band_offsets[i + 1];
        }
        for (i = j - 1; i >= 0; i--)
        {
            band_offsets_local[i * 2 + 1 + 1] = band_offsets[i + 1];
            band_offsets_local[i * 2 + 0 + 1] = band_offsets[i + 1];
        }
    }

    denumerator = 0.001;
    
    for (i = 0; i < M; i++)
    {
        sum = 0; len = 0;
        start = band_offsets_local[i * 4] + 1; stop = band_offsets_local[i * 4 + 4];
        
        for (j = stop; j >= start; j--)
        {
            sum += j;
            len++;
        }
        
        numerator += gains_dee[i] * ((LC3_FLOAT) sum / (LC3_FLOAT) framelength);
        denumerator += gains_dee[i] * len;
    }
    
    *sc = numerator / denumerator;
    *sc = *sc * (LC3_FLOAT) fs / 48000.0; /* scaling, because training is for 48kHz */
}

static void plc_xcorr_lc(LC3_FLOAT *pcmbufHist, LC3_INT32 max_len_pcm_plc, LC3_INT32 pitch_int, LC3_INT32 framelength,
                         LC3_INT32 frame_dms, LC3_INT32 fs, LC3_FLOAT *xcorr)
{
        LC3_INT32 max_corr_len, pitch_min, corr_len, min_corr_len, pcm_max_corr_len, range1Start, range2Start, i;
        LC3_FLOAT norm_w, norm_w_t;
    
    norm_w_t = 0; norm_w = 0;

    assert(pitch_int >= 0);
    assert(pitch_int <= MAX_LEN*100*MAX_PITCH_12K8/12800);
    
    *xcorr = 0;

    if (pitch_int > 0)
    {
        pitch_min = fs * MIN_PITCH_12K8/12800;
        pcm_max_corr_len = max_len_pcm_plc - pitch_int;

        min_corr_len = 2 * pitch_min;              /* at least 5 ms (=2*pitchmin*) corr length */
        max_corr_len = framelength*100/frame_dms;  /* maximum 10 ms */
        max_corr_len = MIN( max_corr_len, pcm_max_corr_len );

        corr_len = MIN( max_corr_len, pitch_int ); /* pitch_int is prefered, but maximum 10ms or left pcm buf size */
        corr_len = MAX( min_corr_len, corr_len );

        range1Start = max_len_pcm_plc - corr_len;
        range2Start = range1Start - pitch_int;

        assert( corr_len >= min_corr_len );
        assert( corr_len <= max_corr_len );
        assert( range2Start >= 0 );

        for (i = 0; i < corr_len; i++)
        {
            norm_w += pcmbufHist[range1Start + i] * pcmbufHist[range1Start + i];
        }

        for (i = 0; i < corr_len; i++)
        {
            norm_w_t += pcmbufHist[range2Start + i] * pcmbufHist[range2Start + i];
        }

        for (i = 0; i < corr_len; i++)
        {
            *xcorr = *xcorr + pcmbufHist[range1Start + i] * pcmbufHist[range2Start + i];
        }

        *xcorr = *xcorr / sqrt(norm_w * norm_w_t + 0.1);
        *xcorr = MAX(0, *xcorr);
    } else {
        *xcorr = 0;
    }
}

