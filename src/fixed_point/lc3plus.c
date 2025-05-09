/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                              
#include "defines.h"
#include "functions.h"
#include "lc3plus.h"
#include "setup_dec_lc3plus.h"
#include "setup_enc_lc3plus.h"

#define RETURN_IF(cond, error)                                                                                         \
    if (cond)                                                                                                          \
    return (error)

#ifdef SUBSET_NB
#pragma message("- SUBSET_NB")
#endif
#ifdef SUBSET_WB
#pragma message("- SUBSET_WB")
#endif
#ifdef SUBSET_SSWB
#pragma message("- SUBSET_SSWB")
#endif
#ifdef SUBSET_SWB
#pragma message("- SUBSET_SWB")
#endif
#ifdef SUBSET_FB
#pragma message("- SUBSET_FB")
#endif
#ifdef SUBSET_UB
#pragma message("- SUBSET_UB")
#endif

/* ensure api header constants are up to date */
STATIC_ASSERT(LC3PLUS_MAX_SAMPLES >= MAX_LEN);
STATIC_ASSERT(LC3PLUS_MAX_CHANNELS >= MAX_CHANNELS);
STATIC_ASSERT(LC3PLUS_MAX_BYTES >= BYTESBUFSIZE);
STATIC_ASSERT(LC3PLUS_ENC_MAX_SIZE >= ENC_MAX_SIZE);
STATIC_ASSERT(LC3PLUS_DEC_MAX_SIZE >= DEC_MAX_SIZE);
STATIC_ASSERT(LC3PLUS_ENC_MAX_SCRATCH_SIZE >= SCRATCH_BUF_LEN_ENC_TOT);
STATIC_ASSERT(LC3PLUS_DEC_MAX_SCRATCH_SIZE >= SCRATCH_BUF_LEN_DEC_TOT);
STATIC_ASSERT(PLC_FADEOUT_IN_MS >= 20);


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
    switch (samplerate)
    {
#ifdef SUBSET_NB
    case 8000: return 1;
#endif
#ifdef SUBSET_WB
    case 16000: return 1;
#endif
#ifdef SUBSET_SSWB
    case 24000: return 1;
#endif
#ifdef SUBSET_SWB
    case 32000: return 1;
#endif
#ifdef SUBSET_FB
    case 44100: return 1;
    case 48000: return 1;
#endif
#ifdef ENABLE_HR_MODE
    case 96000: return 1;
#endif
    default: return 0;
    }
}

static int lc3plus_plc_mode_supported(LC3PLUS_PlcMode plc_mode)
{
    switch ((int)plc_mode)
    {
    case LC3PLUS_PLC_ADVANCED: /* fallthru */
        return 1;
    default: return 0;
    }
}

static int lc3plus_frame_size_supported(int frame_dms)
{
    switch (frame_dms)
    {
    case 25: /* fallthru */
    case 50: /* fallthru */
    case 75: /* fallthru */
    case 100:
            return 1;
    default: return 0;
    }
}

static int null_in_list(void **list, int n)
{
    while (--n >= 0)
        RETURN_IF(list[n] == NULL, 1);
    return 0;
}

/* return pointer to aligned base + base_size, *base_size += size + 4 bytes align */
void *balloc(void *base, size_t *base_size, size_t size)
{
    uintptr_t ptr = ((uintptr_t)base + *base_size + 3) & ~3;
    assert((uintptr_t)base % 4 == 0); /* base must be 4-byte aligned */
    *base_size = (*base_size + size + 3) & ~3;
    return (void *)ptr;
}

int32_t lc3_enc_supported_lfe(void)
{
    return 1;
}

/* encoder functions *********************************************************/

LC3PLUS_Error lc3plus_enc_init(LC3PLUS_Enc *encoder, int samplerate, int channels
#ifdef ENABLE_HR_MODE
                               , int hrmode
#endif
                               , int32_t lfe_channel_array[]
                              )
{
    int ch = 0;

    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((uintptr_t)encoder % 4 != 0, LC3PLUS_ALIGN_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
#ifdef ENABLE_HR_MODE
    RETURN_IF(samplerate==96000 && hrmode == 0, LC3PLUS_HRMODE_ERROR);
#endif

    if (lfe_channel_array != NULL)
    {
        for (ch = 0; ch < channels; ch++)
        {
            RETURN_IF(!lc3_enc_supported_lfe() && lfe_channel_array[ch], LC3PLUS_LFE_MODE_NOT_SUPPORTED);
        }
    }

#ifdef ENABLE_HR_MODE
    return FillEncSetup(encoder, samplerate, channels, hrmode
                        , lfe_channel_array
    ); /* real bitrate check happens here */
#else
    return FillEncSetup(encoder, samplerate, channels
                        , lfe_channel_array
    ); /* real bitrate check happens here */
#endif
}

int lc3plus_enc_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    return alloc_encoder(NULL, samplerate, channels);
}

int lc3plus_enc_get_scratch_size(const LC3PLUS_Enc *encoder)
{
    int size = 0;
     RETURN_IF(encoder == NULL, 0);

#ifdef ENABLE_HR_MODE
    size = 47 * MAX(encoder->frame_length, 160) + 64;
#else
    size = 14 * MAX(encoder->frame_length, 160) + 64;
#endif
    assert(size <= LC3PLUS_ENC_MAX_SCRATCH_SIZE);
    return size;
}

int lc3plus_enc_get_input_samples(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length;
}

int lc3plus_enc_get_num_bytes(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return (Word32)encoder->bitrate * encoder->frame_length / (8 * encoder->fs_in);
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
    
    int bitrate = (totalBytes * 80000.0 + encoder->frame_dms - 1) / encoder->frame_dms;
    
    if (encoder->fs_in == 44100)
    {
        int rem = bitrate % 480;
        bitrate = ((bitrate - rem) / 480) * 441 + (rem * 441) / 480;
    }
    
    return bitrate;
}

LC3PLUS_Error lc3plus_enc_set_bitrate(LC3PLUS_Enc *encoder, int bitrate)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(bitrate <= 0, LC3PLUS_BITRATE_ERROR);
    return update_enc_bitrate(encoder, bitrate);
}

int lc3plus_enc_get_delay(const LC3PLUS_Enc *encoder)
{
    RETURN_IF(encoder == NULL, 0);
    return encoder->frame_length - 2 * encoder->la_zeroes;
}

LC3PLUS_Error lc3plus_enc_set_ep_mode(LC3PLUS_Enc *encoder, LC3PLUS_EpMode epmode)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((unsigned)epmode > LC3PLUS_EP_HIGH, LC3PLUS_EPMODE_ERROR);
    encoder->epmode = epmode;
    return encoder->lc3_br_set ? update_enc_bitrate(encoder, encoder->bitrate) : LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc_set_ep_mode_request(LC3PLUS_Enc *encoder, LC3PLUS_EpModeRequest epmr)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((unsigned)epmr > LC3PLUS_EPMR_HIGH, LC3PLUS_EPMR_ERROR);
    encoder->epmr = epmr;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc_set_frame_dms(LC3PLUS_Enc *encoder, int frame_dms)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_dms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(encoder->lc3_br_set, LC3PLUS_BITRATE_SET_ERROR);
    encoder->frame_dms = frame_dms;
    set_enc_frame_params(encoder);
    return LC3PLUS_OK;
}


LC3PLUS_Error lc3plus_enc_set_bandwidth(LC3PLUS_Enc *encoder, int bandwidth)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    Word32 effective_fs = encoder->fs_in;
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
            encoder->bw_ctrl_active   = 1;
            update_enc_bitrate(encoder, encoder->bitrate);
        }
    }
    return LC3PLUS_OK;
}


static LC3PLUS_Error lc3plus_enc(LC3PLUS_Enc *encoder, void **input_samples, int bitdepth, void *output_bytes, int *num_bytes,
                         void *scratch)
{
    RETURN_IF(!encoder || !input_samples || !output_bytes || !num_bytes || !scratch, LC3PLUS_NULL_ERROR);
    RETURN_IF(null_in_list(input_samples, encoder->channels), LC3PLUS_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3PLUS_ERROR);
    RETURN_IF(!encoder->lc3_br_set, LC3PLUS_BITRATE_UNSET_ERROR);
    *num_bytes = Enc_LC3PLUS(encoder, input_samples, bitdepth, output_bytes, scratch, *num_bytes == -1);
    
    assert(*num_bytes == lc3plus_enc_get_num_bytes(encoder));
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc16(LC3PLUS_Enc *encoder, int16_t **input_samples, void *output_bytes, int *num_bytes, void *scratch)
{
    return lc3plus_enc(encoder, (void **)input_samples, 16, output_bytes, num_bytes, scratch);
}

LC3PLUS_Error lc3plus_enc24(LC3PLUS_Enc *encoder, int32_t **input_samples, void *output_bytes, int *num_bytes, void *scratch)
{
    return lc3plus_enc(encoder, (void **)input_samples, 24, output_bytes, num_bytes, scratch);
}

/* decoder functions *********************************************************/

LC3PLUS_Error lc3plus_dec_init(LC3PLUS_Dec *decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode
#ifdef ENABLE_HR_MODE
                            , int hrmode
#endif
)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
    RETURN_IF(!lc3plus_plc_mode_supported(plc_mode), LC3PLUS_PLCMODE_ERROR);
#ifdef ENABLE_HR_MODE
    RETURN_IF(samplerate==96000 && hrmode == 0, LC3PLUS_HRMODE_ERROR);
#endif

    return FillDecSetup(decoder, samplerate, channels, plc_mode
#ifdef ENABLE_HR_MODE
                        , hrmode
#endif
                       );
} 

int lc3plus_dec_get_size(int samplerate, int channels, LC3PLUS_PlcMode plc_mode)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    RETURN_IF(!lc3plus_plc_mode_supported(plc_mode), 0);
    return alloc_decoder(NULL, samplerate, channels);
}

int lc3plus_dec_get_scratch_size(const LC3PLUS_Dec *decoder)
{
    int size = 0;
    RETURN_IF(decoder == NULL, 0);

#ifdef ENABLE_HR_MODE
    size = 30 * DYN_MAX_LEN(decoder->fs) + 2866;
    size += 4 * MAX_LGW + 8 * DYN_MAX_LPROT(decoder->fs) + 16 * DYN_MAX_LEN(decoder->fs);
#else
    size = 12 * DYN_MAX_LEN(decoder->fs) + 752;
    size += 2 * MAX_LGW + 8 * DYN_MAX_LPROT(decoder->fs) + 8 * DYN_MAX_LEN(decoder->fs);
    size += 3720;
#endif

    assert(size <= LC3PLUS_DEC_MAX_SCRATCH_SIZE);
    return size;
}

LC3PLUS_Error lc3plus_dec_set_ep_enabled(LC3PLUS_Dec *decoder, int ep_enabled)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    decoder->ep_enabled = ep_enabled != 0;
    decoder->epmr       = LC3PLUS_EPMR_ZERO;
    return LC3PLUS_OK;
}

int lc3plus_dec_get_error_report(const LC3PLUS_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->error_report == 2047 ? -1 : decoder->error_report & 0x07FF;
}

int lc3plus_dec_get_epok_flags(const LC3PLUS_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->error_report >> 11;
}

LC3PLUS_EpModeRequest lc3plus_dec_get_ep_mode_request(const LC3PLUS_Dec *decoder)
{
    RETURN_IF(decoder == NULL, LC3PLUS_EPMR_ZERO);
    return (LC3PLUS_EpModeRequest)decoder->epmr;
}

LC3PLUS_Error lc3plus_dec_set_frame_dms(LC3PLUS_Dec *decoder, int frame_dms)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_dms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(decoder->plcMeth == 2 && frame_dms != 100, LC3PLUS_FRAMEMS_ERROR);

    decoder->frame_dms = frame_dms;
    set_dec_frame_params(decoder);
    return LC3PLUS_OK;
}


int lc3plus_dec_get_output_samples(const LC3PLUS_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length;
}

int lc3plus_dec_get_delay(const LC3PLUS_Dec *decoder)
{
    RETURN_IF(decoder == NULL, 0);
    return decoder->frame_length - 2 * decoder->la_zeroes;
}

static LC3PLUS_Error lc3plus_dec(LC3PLUS_Dec *decoder, void *input_bytes, int num_bytes, void **output_samples, int bitdepth,
                         void *scratch, int bfi_ext)
{
    if (bfi_ext == 1)
    {
        RETURN_IF(!decoder || !output_samples || !scratch, LC3PLUS_NULL_ERROR);
    } else {
        RETURN_IF(!decoder || !input_bytes || !output_samples || !scratch, LC3PLUS_NULL_ERROR);
    }

    RETURN_IF(null_in_list(output_samples, decoder->channels), LC3PLUS_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3PLUS_ERROR);
    return Dec_LC3PLUS(decoder, input_bytes, num_bytes, output_samples, bitdepth, scratch, bfi_ext);
}

LC3PLUS_Error lc3plus_dec16(LC3PLUS_Dec *decoder, void *input_bytes, int num_bytes, int16_t **output_samples, void *scratch, int bfi_ext)
{
    return lc3plus_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 16, scratch, bfi_ext);
}

LC3PLUS_Error lc3plus_dec24(LC3PLUS_Dec *decoder, void *input_bytes, int num_bytes, int32_t **output_samples, void *scratch, int bfi_ext)
{
    return lc3plus_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 24, scratch, bfi_ext);
}
