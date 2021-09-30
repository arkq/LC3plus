/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "lc3.h"
#include "defines.h"
#include "functions.h"
#include <stdio.h>

#include "setup_dec_lc3.h"
#include "setup_enc_lc3.h"

#define RETURN_IF(cond, error) \
    if (cond)                  \
    return (error)

/* ensure api header constants are up to date */
STATIC_ASSERT(LC3PLUS_MAX_SAMPLES >= MAX_LEN);
STATIC_ASSERT(LC3PLUS_MAX_CHANNELS >= MAX_CHANNELS);
STATIC_ASSERT(LC3PLUS_MAX_BYTES >= BYTESBUFSIZE);

/* misc functions ************************************************************/

int lc3plus_version(void)
{
    return LC3PLUS_VERSION;
}

int lc3plus_channels_supported(int channels)
{
    return channels >= 1 && channels <= MAX_CHANNELS;
}

int lc3plus_samplerate_supported(int samplerate)
{
    switch (samplerate) {
    case 8000:
        return 1;
    case 16000:
        return 1;
    case 24000:
        return 1;
    case 32000:
        return 1;
    case 44100:
        return 1;
    case 48000:
        return 1;
#ifdef ENABLE_HR_MODE_FL_FLAG
    case 96000:
        return 1;
#endif
    default:
        return 0;
    }
}

static int lc3plus_plc_mode_supported(LC3PLUS_PlcMode plc_mode)
{
    switch ((int)plc_mode)
    {
    case LC3PLUS_PLC_ADVANCED: return 1;
    default: return 0;
    }
}

static int lc3plus_frame_size_supported(float frame_ms)
{
    switch ((int)(ceil(frame_ms * 10)))
    {
    case 25: /* fallthru */
    case 50: /* fallthru */
    case 100: return 1;
    default: return 0;
    }
}

static int null_in_list(void **list, int n)
{
    while (--n >= 0)
        RETURN_IF(list[n] == NULL, 1);
    return 0;
}


/* encoder functions *********************************************************/

LC3PLUS_Error lc3plus_enc_init(LC3PLUS_Enc *encoder, int samplerate, int channels)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((uintptr_t)encoder % 4 != 0, LC3PLUS_ALIGN_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
    return FillEncSetup(encoder, samplerate, channels); /* real bitrate check happens here */
}

int lc3plus_enc_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    return alloc_encoder(NULL, channels);
}

int lc3plus_enc_get_input_samples(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length;
}

int lc3plus_enc_get_num_bytes(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    
    return encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
}

int lc3plus_enc_get_real_bitrate(const LC3PLUS_Enc *encoder)
{
    int ch = 0, totalBytes = 0;
    RETURN_IF(encoder == NULL, 0);
    RETURN_IF(!encoder->lc3_br_set, LC3PLUS_BITRATE_UNSET_ERROR);
    for (ch = 0; ch < encoder->channels; ch++)
    {
        totalBytes += encoder->channel_setup[ch]->targetBytes;
    }
    int bitrate = (totalBytes * 80000)/ encoder->frame_dms;
    if (encoder->fs_in == 44100)
    {
        int rem = bitrate % 480;
        bitrate = ((bitrate - rem) / 480)* 441 + (rem * 441) / 480;
    }
    return bitrate;
}

LC3PLUS_Error lc3plus_enc_set_bitrate(LC3PLUS_Enc *encoder, int bitrate)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(bitrate <= 0, LC3PLUS_BITRATE_ERROR);
#ifndef STRIP_HR_MODE_API
    RETURN_IF(encoder->fs_idx == 5 && encoder->hrmode == 0, LC3PLUS_HRMODE_ERROR);
#endif
    return update_enc_bitrate(encoder, bitrate);
}

int lc3plus_enc_get_delay(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length - 2 * encoder->la_zeroes;
}

LC3PLUS_Error lc3plus_enc_set_frame_ms(LC3PLUS_Enc *encoder, float frame_ms)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_ms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(encoder->lc3_br_set, LC3PLUS_BITRATE_SET_ERROR);
    encoder->frame_dms = (int)(frame_ms * 10);
    encoder->frame_ms = frame_ms;
    set_enc_frame_params(encoder);
    return LC3PLUS_OK;
}

#ifndef STRIP_HR_MODE_API
LC3PLUS_Error lc3plus_enc_set_hrmode(LC3PLUS_Enc* encoder, int hrmode)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(encoder->fs_in < 48000 && hrmode != 0, LC3PLUS_SAMPLERATE_ERROR);
    encoder->hrmode = hrmode > 0;
    set_enc_frame_params(encoder);
    return LC3PLUS_OK;
}
#endif

LC3PLUS_Error lc3plus_enc_set_bandwidth(LC3PLUS_Enc *encoder, int bandwidth)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
#ifdef ENABLE_HR_MODE_FL_FLAG
    RETURN_IF(encoder->hrmode == 1, LC3PLUS_HRMODE_BW_ERROR);
#endif
    LC3_INT effective_fs = encoder->fs_in;
    if (encoder->bandwidth != bandwidth) {
        if (encoder->fs_in > 40000) {
            effective_fs = 40000;
        }
        if ((bandwidth * 2) > effective_fs) {
            return LC3PLUS_BW_WARNING;
        }
        else {
            encoder->bandwidth = bandwidth;
            encoder->bandwidth_preset = bandwidth;
            encoder->bw_ctrl_active = 1;     
            update_enc_bitrate(encoder, encoder->bitrate);
        }
    }
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc_set_lfe(LC3PLUS_Enc *encoder, int lfe)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    encoder->lfe = (lfe != 0);
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc16(LC3PLUS_Enc* encoder, int16_t** input_samples, void* output_bytes, int* num_bytes)
{
    return lc3plus_enc_fl(encoder, (void**)input_samples, 16, output_bytes, num_bytes);
}

LC3PLUS_Error lc3plus_enc24(LC3PLUS_Enc* encoder, int32_t** input_samples, void* output_bytes, int* num_bytes)
{
    return lc3plus_enc_fl(encoder, (void**)input_samples, 24, output_bytes, num_bytes);
}


LC3PLUS_Error lc3plus_enc_fl(LC3PLUS_Enc* encoder, void** input_samples, int bitdepth, void* output_bytes, int* num_bytes)
{
    RETURN_IF(!encoder || !input_samples || !output_bytes || !num_bytes, LC3PLUS_NULL_ERROR);
    RETURN_IF(null_in_list(input_samples, encoder->channels), LC3PLUS_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24 && bitdepth != 32, LC3PLUS_ERROR);
    *num_bytes = Enc_LC3PLUS_fl(encoder, input_samples, output_bytes, bitdepth
                            , *num_bytes == -1
                           );
    assert(*num_bytes == lc3plus_enc_get_num_bytes(encoder));
    return LC3PLUS_OK;
}

/* decoder functions *********************************************************/

LC3PLUS_Error lc3plus_dec_init(LC3PLUS_Dec* decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
    RETURN_IF(!lc3plus_plc_mode_supported(plc_mode), LC3PLUS_PLCMODE_ERROR);
    return FillDecSetup(decoder, samplerate, channels, plc_mode);
}

int lc3plus_dec_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    return alloc_decoder(NULL, samplerate, channels);
}

LC3PLUS_Error lc3plus_dec_set_frame_ms(LC3PLUS_Dec* decoder, float frame_ms)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_ms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(decoder->plcMeth == 2 && frame_ms != 10, LC3PLUS_FRAMEMS_ERROR);
    decoder->frame_ms = frame_ms;
    decoder->frame_dms = (LC3_INT) (frame_ms * 10);
    set_dec_frame_params(decoder);
    return LC3PLUS_OK;
}

int lc3plus_dec_get_output_samples(const LC3PLUS_Dec* decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length;
}

int lc3plus_dec_get_delay(const LC3PLUS_Dec* decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length - 2 * decoder->la_zeroes;
}

LC3PLUS_Error lc3plus_dec_fl(LC3PLUS_Dec* decoder, void* input_bytes, int num_bytes, void** output_samples, int bps, int bfi_ext)
{
    RETURN_IF(!decoder || !input_bytes || !output_samples, LC3PLUS_NULL_ERROR);
    RETURN_IF(null_in_list((void**)output_samples, decoder->channels), LC3PLUS_NULL_ERROR);
    return Dec_LC3PLUS_fl(decoder, input_bytes, num_bytes, output_samples, bps, bfi_ext);
}

LC3PLUS_Error lc3plus_dec16(LC3PLUS_Dec* decoder, void* input_bytes, int num_bytes, int16_t** output_samples, int bfi_ext)
{
    return lc3plus_dec_fl(decoder, input_bytes, num_bytes, (void**)output_samples, 16, bfi_ext);
}

LC3PLUS_Error lc3plus_dec24(LC3PLUS_Dec* decoder, void* input_bytes, int num_bytes, int32_t** output_samples, int bfi_ext)
{
    return lc3plus_dec_fl(decoder, input_bytes, num_bytes, (void**)output_samples, 24, bfi_ext);
}

/* memory functions *********************************************************/

LC3PLUS_Error lc3plus_enc_free_memory(LC3PLUS_Enc* encoder)
{
    RETURN_IF(!encoder, LC3PLUS_NULL_ERROR);

    lc3plus_free_encoder_structs(encoder);
    free(encoder);

    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_dec_free_memory(LC3PLUS_Dec* decoder)
{
    RETURN_IF(!decoder, LC3PLUS_NULL_ERROR);

    lc3plus_free_decoder_structs(decoder);
    free(decoder);

    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_free_encoder_structs(LC3PLUS_Enc* encoder)
{
    RETURN_IF(!encoder, LC3PLUS_NULL_ERROR);

    int ch = 0;
    for (ch = 0; ch < encoder->channels; ch++) {
        mdct_free(&encoder->channel_setup[ch]->mdctStruct);
        dct2_free(&encoder->channel_setup[ch]->dct2StructSNS);
    }

    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_free_decoder_structs(LC3PLUS_Dec* decoder)
{
    RETURN_IF(!decoder, LC3PLUS_NULL_ERROR);

    int ch = 0;
    for (ch = 0; ch < decoder->channels; ch++) {
        dct4_free(&decoder->channel_setup[ch]->dct4structImdct);
        real_fft_free(&decoder->channel_setup[ch]->PlcAdvSetup->PlcPhEcuSetup.PhEcu_Fft);
        real_fft_free(&decoder->channel_setup[ch]->PlcAdvSetup->PlcPhEcuSetup.PhEcu_Ifft);
    }

    return LC3PLUS_OK;
}

#ifndef STRIP_HR_MODE_API_MODE_API
LC3PLUS_Error lc3plus_dec_set_hrmode(LC3PLUS_Dec* decoder, int hrmode)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(decoder->fs_idx < 4 && hrmode != 0, LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF((decoder->fs_idx == 5) && (hrmode == 0), LC3PLUS_HRMODE_ERROR);
    decoder->hrmode = hrmode > 0;
    set_dec_frame_params(decoder);
    return LC3PLUS_OK;
}
#endif

LC3PLUS_Error lc3plus_dec_get_ep_mode_request(const LC3PLUS_Dec* decoder, LC3PLUS_EpModeRequest* const epmr)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    *epmr = (LC3PLUS_EpModeRequest)decoder->epmr;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_dec_get_error_report(const LC3PLUS_Dec *decoder, int32_t* const errorReport)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    *errorReport = decoder->error_report == 2047 ? -1 : decoder->error_report & 0x07FF;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc_set_ep_mode(LC3PLUS_Enc *encoder, LC3PLUS_EpMode epmode)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((unsigned)epmode > LC3PLUS_EP_HIGH, LC3PLUS_EPMODE_ERROR);
    LC3PLUS_EpMode oldEpmode = encoder->epmode;
    encoder->epmode = epmode;
    LC3PLUS_Error error = encoder->lc3_br_set ? update_enc_bitrate(encoder, encoder->bitrate) : LC3PLUS_OK;
    if (error != LC3PLUS_OK)
    {
        encoder->epmode = oldEpmode;  // preserve old epmode in case of failure
    }
    return error;
}

LC3PLUS_Error lc3plus_enc_set_ep_mode_request(LC3PLUS_Enc *encoder, LC3PLUS_EpModeRequest epmr)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((unsigned)epmr > LC3PLUS_EPMR_HIGH, LC3PLUS_EPMODE_ERROR);
    encoder->epmr = epmr;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_dec_set_ep_enabled(LC3PLUS_Dec *decoder, int32_t ep_enabled)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    decoder->ep_enabled = ep_enabled != 0;
    decoder->epmr       = LC3PLUS_EPMR_ZERO;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_dec_get_epok_flags(LC3PLUS_Dec* const decoder, int32_t* const epokFlags)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    *epokFlags = decoder->error_report >> 11;
    return LC3PLUS_OK;
}

#ifndef STRIP_ERROR_PROTECTION_API_FL
#endif /* STRIP_ERROR_PROTECTION_API_FL */

#ifndef STRIP_ERROR_PROTECTION_API_FL
#endif /* STRIP_ERROR_PROTECTION_API_FL */

