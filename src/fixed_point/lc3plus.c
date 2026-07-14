/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
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
#if defined(ENABLE_HR_MODE) || defined(SUBSET_UB)
    case 96000: return 1;
#endif
#ifdef SUBSET_UUB
    case 192000: return 1;
#endif
    default: return 0;
    }
    return 0;
}

static int lc3plus_plc_mode_supported(LC3PLUS_PlcMode plc_mode)
{
    switch ((int)plc_mode)
    {
    case LC3PLUS_PLC_ADVANCED: /* fallthru */
        return 1;
    default: return 0;
    }
    return 0;
}

static int lc3plus_frame_size_supported(LC3PLUS_FrameDuration frame_dms)
{
    switch (frame_dms)
    {
#ifdef CR9_C_ADD_1p25MS
    case LC3PLUS_FRAME_DURATION_1p25MS: /* fallthru */
#endif
    case LC3PLUS_FRAME_DURATION_2p5MS: /* fallthru */
    case LC3PLUS_FRAME_DURATION_5MS: /* fallthru */
    case LC3PLUS_FRAME_DURATION_7p5MS: /* fallthru */
    case LC3PLUS_FRAME_DURATION_10MS:
            return 1;
    default: return 0;
    }
    return 0;
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
                               , int32_t lfe_channel_array[], int32_t* const scratchSize
#ifdef CR14_A_ADD_LOSSLESS_MODE
                               , int wavFormat, int padding
#endif
                              )
{
    int ch = 0;
    LC3PLUS_Error err = LC3PLUS_OK;

    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF((uintptr_t)encoder % 4 != 0, LC3PLUS_ALIGN_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
#ifdef ENABLE_HR_MODE
    RETURN_IF(samplerate == 96000 && hrmode == 0, LC3PLUS_HRMODE_ERROR);
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    RETURN_IF( (samplerate < 44100 && hrmode != 0), LC3PLUS_SAMPLERATE_ERROR );
    RETURN_IF( samplerate == 192000 && hrmode <= 1, LC3PLUS_HRMODE_ERROR );
#endif

    if (lfe_channel_array != NULL)
    {
        for (ch = 0; ch < channels; ch++)
        {
            RETURN_IF(!lc3_enc_supported_lfe() && lfe_channel_array[ch], LC3PLUS_LFE_MODE_NOT_SUPPORTED);
        }
    }

#ifdef ENABLE_HR_MODE
    err = FillEncSetup(encoder, samplerate, channels, hrmode
                        , lfe_channel_array
#ifdef CR14_A_ADD_LOSSLESS_MODE
                        , wavFormat, padding
#endif
    ); /* real bitrate check happens here */
#else
    err = FillEncSetup(encoder, samplerate, channels
                        , lfe_channel_array
    ); /* real bitrate check happens here */
#endif

    err = lc3plus_enc_get_scratch_size( encoder, scratchSize );  // dummy-encoder run
    if ( err != LC3PLUS_OK )
    {
        return err;
    }

    memset(encoder, 0, lc3plus_enc_get_size(samplerate, channels)); // clean-up after dummy-encoder run
#ifdef ENABLE_HR_MODE
    err = FillEncSetup(encoder, samplerate, channels, hrmode
                        , lfe_channel_array
#ifdef CR14_A_ADD_LOSSLESS_MODE
                        , wavFormat, padding
#endif
    ); /* real bitrate check happens here */
#else
    err = FillEncSetup(encoder, samplerate, channels
                        , lfe_channel_array
    ); /* real bitrate check happens here */
#endif
    
    encoder->scratch_max_size = *scratchSize;
    encoder->lc3_scratch_initialized = 1;
    
    return err;
}

int lc3plus_enc_get_size(int samplerate, int channels)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    return alloc_encoder(NULL, samplerate, channels);
}

int lc3plus_enc_get_scratch_size(LC3PLUS_Enc *encoder, int* const scratch_size)
{
    size_t size = 0;
    LC3PLUS_Error err = LC3PLUS_OK;
    Word16* input16[1];
#ifdef CR14_A_ADD_LOSSLESS_MODE
    int nBytes[LC3PLUS_MAX_CHANNELS] = { 0 };
#endif
    Word32 buf_32[LC3PLUS_MAX_SAMPLES] = { 0 };
    UWord8 bytes[LC3PLUS_MAX_BYTES] = { 0 };
    UWord8 tmp_scratch[LC3PLUS_ENC_MAX_USER_SYSTEM_SCRATCH_SIZE]; /* temp/dummy scratch buffer needed to run the encoder with a zero-frame to calculate the max. used scratch size */

    Word16* buf_16 = (int16_t*) (void*) buf_32;
    Word32 bitrate;

    lc3_scratch_t enc_scratch = lc3_scratch_init( (void*) tmp_scratch, LC3PLUS_ENC_MAX_SCRATCH_SIZE, SCRATCH_BUFFER_ALIGNMENT_BITS, SCRATCH_ALLOCATOR_CALCULATE_MAX );

    input16[0] = buf_16;

#  ifdef ENABLE_HR_MODE
    IF( encoder->hrmode )
    {
        SWITCH( encoder->frame_dms )
        {
        case LC3PLUS_FRAME_DURATION_1p25MS:
            bitrate = 672000;
            BREAK;
        case LC3PLUS_FRAME_DURATION_2p5MS:
            bitrate = 672000;
            BREAK;
        case LC3PLUS_FRAME_DURATION_5MS:
            bitrate = 600000;
            BREAK;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            bitrate = 500000;
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS:
            bitrate = 500000;
            BREAK;
        default:
            return LC3PLUS_HRMODE_ERROR;
        }
    }
    ELSE
#  endif
    {
        SWITCH( encoder->frame_dms )
        {
        case LC3PLUS_FRAME_DURATION_1p25MS:
            bitrate = MAX_BR;
            BREAK;
        case LC3PLUS_FRAME_DURATION_2p5MS:
            bitrate = MAX_BR;
            BREAK;
        case LC3PLUS_FRAME_DURATION_5MS:
            bitrate = MAX_BR;
            SWITCH( encoder->fs_in )
            {
            case 8000:
                bitrate = MAX_BR_050DMS_NB;
                BREAK;
            default:
                BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_7p5MS:
            bitrate = MAX_BR_075DMS;
            SWITCH( encoder->fs_in )
            {
            case 8000:
                bitrate = MAX_BR_075DMS_NB;
                BREAK;
            case 16000:
                bitrate = MAX_BR_075DMS_WB;
                BREAK;
            case 24000:
                bitrate = MAX_BR_075DMS_SSWB;
                BREAK;
            default:
                BREAK;
            }
            BREAK;
        case LC3PLUS_FRAME_DURATION_10MS:
            bitrate = MAX_BR;
            SWITCH( encoder->fs_in )
            {
            case 8000:
                bitrate = MAX_BR_100DMS_NB;
                BREAK;
            case 16000:
                bitrate = MAX_BR_100DMS_WB;
                BREAK;
            case 24000:
                bitrate = MAX_BR_100DMS_SSWB;
                BREAK;
            default:
                BREAK;
            }
            BREAK;
        default:
            return LC3PLUS_FRAMEMS_ERROR;
        }
        if ( encoder->fs_in == 44100 )
        {
            bitrate = Mpy_32_32( bitrate, 1973000602 );
        }
    }
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    /* Required for correct scratch calculation */
    if (encoder->fs_in == 44100 && encoder->hrmode)
    {
        bitrate = 294000;
    }
#endif
  
    SWITCH( encoder->frame_dms )
    {
    case LC3PLUS_FRAME_DURATION_10MS:
        bitrate = MIN( bitrate, 320000 );
        BREAK;
    case LC3PLUS_FRAME_DURATION_7p5MS:
        bitrate = MIN( bitrate, 320000 );
        BREAK;
    case LC3PLUS_FRAME_DURATION_5MS:
        bitrate = MIN( bitrate, 480000 );
        BREAK;
    case LC3PLUS_FRAME_DURATION_2p5MS:
        bitrate = MIN( bitrate, 960000 );
        BREAK;
    case LC3PLUS_FRAME_DURATION_1p25MS:
        bitrate = MIN( bitrate, 960000 );
        BREAK;
    default:
        BREAK;
    }

    err = lc3plus_enc_set_bitrate( encoder, bitrate );
    if ( err != LC3PLUS_OK )
    {
        return LC3PLUS_SCRATCH_INVALID_ERROR;
    }

    err = lc3plus_enc_set_ep_mode( encoder, LC3PLUS_EP_HIGH );
    if ( err != LC3PLUS_OK )
    {
        return LC3PLUS_SCRATCH_INVALID_ERROR;
    }

#ifndef CR14_A_ADD_LOSSLESS_MODE
    int nbytes_ret =     Enc_LC3PLUS( encoder, (void**) input16, 
                16, 
                bytes, enc_scratch, 0);
    UNUSED(nbytes_ret);
#else
    Enc_LC3PLUS( encoder, (void**) input16, 
                encoder->wavFormat, 
                bytes, enc_scratch, 0, 
                nBytes);
#endif
    size = enc_scratch->max_alloc_size;

    assert( size <= LC3PLUS_ENC_MAX_USER_SYSTEM_SCRATCH_SIZE );
    assert( size <= LC3PLUS_ENC_MAX_SCRATCH_SIZE );

#  ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "determined scratch size during init: %zu\n", size );
#  endif

    *scratch_size = size;
    return LC3PLUS_OK;
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
#ifdef CR9_C_ADD_1p25MS
    int frame_ns = (int)(1250L*(encoder->frame_dms));
    int bitrate = ((long long int)totalBytes * 8L * 1000000L + (frame_ns - 1L)) / frame_ns;
#else
    int bitrate = (totalBytes * 80000.0 + encoder->frame_dms - 1) / encoder->frame_dms;
#endif

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
#ifdef CR14_A_ADD_LOSSLESS_MODE
    RETURN_IF(encoder->hrmode < 2 && bitrate <= 0, LC3PLUS_BITRATE_ERROR);
#ifdef PADDING_CBR
    RETURN_IF(encoder->padding == 1 && (encoder->hrmode != 2 || bitrate == 0), LC3PLUS_PADDING_ERROR);
#endif
#else
    RETURN_IF(bitrate <= 0, LC3PLUS_BITRATE_ERROR);
#endif
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

#ifdef CR14_A_ADD_LOSSLESS_MODE
LC3PLUS_Error lc3plus_enc_set_relative_priority(LC3PLUS_Enc *encoder, int enable)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(enable != 0 && enable != 1, LC3PLUS_EPMODE_ERROR);
    encoder->b_relative = (Word16) enable;
    return LC3PLUS_OK;
}
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
LC3PLUS_Error lc3plus_enc_set_ll_shift(LC3PLUS_Enc *encoder, int shift)
{
    RETURN_IF(encoder == NULL, LC3PLUS_SHIFT_ERROR );
    RETURN_IF(shift < 0 || shift > 8, LC3PLUS_SHIFT_ERROR );
    RETURN_IF(encoder->wavFormat == 16, LC3PLUS_SHIFT_ERROR );
    encoder->ll_shift = (Word16) shift;

    if ( encoder->lossless
         && shift > 0
#  ifdef ENABLE_HR_MODE
         && encoder->hrmode
#  endif
         && (encoder->fs == 48000 || encoder->fs == 96000 || encoder->fs == 192000)
         && encoder->wavFormat == 24 )
    {
        Word32 sh_minBR = 0;
        if ( encoder->fs == 48000 )
        {
            SWITCH( encoder->frame_dms )
            {
            case LC3PLUS_FRAME_DURATION_10MS:   sh_minBR = 124800; BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:  sh_minBR = 124800; BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:    sh_minBR = 148800; BREAK;
            case LC3PLUS_FRAME_DURATION_2p5MS:  sh_minBR = 172800; BREAK;
            case LC3PLUS_FRAME_DURATION_1p25MS: sh_minBR = 204800; BREAK;
            default:                                               BREAK;
            }
        }
        else if ( encoder->fs == 96000 )
        {
            SWITCH( encoder->frame_dms )
            {
            case LC3PLUS_FRAME_DURATION_10MS:   sh_minBR = 149600; BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:  sh_minBR = 149600; BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:    sh_minBR = 174400; BREAK;
            case LC3PLUS_FRAME_DURATION_2p5MS:  sh_minBR = 198400; BREAK;
            case LC3PLUS_FRAME_DURATION_1p25MS: sh_minBR = 230400; BREAK;
            default:                                               BREAK;
            }
        }
        else /* fs == 192000 */
        {
            SWITCH( encoder->frame_dms )
            {
            case LC3PLUS_FRAME_DURATION_10MS:   sh_minBR = 199200; BREAK;
            case LC3PLUS_FRAME_DURATION_7p5MS:  sh_minBR = 199200; BREAK;
            case LC3PLUS_FRAME_DURATION_5MS:    sh_minBR = 225600; BREAK;
            case LC3PLUS_FRAME_DURATION_2p5MS:  sh_minBR = 249600; BREAK;
            case LC3PLUS_FRAME_DURATION_1p25MS: sh_minBR = 249600; BREAK;
            default:                                               BREAK;
            }
        }

        if(encoder->channels > 1)
        {
            sh_minBR *= encoder->channels;
        }

        if ( sh_minBR > 0 && encoder->bitrate < sh_minBR )
        {
            return LC3PLUS_SHIFT_ERROR;
        }
    }
    return LC3PLUS_OK;
}
#endif

LC3PLUS_Error lc3plus_enc_set_frame_dms(LC3PLUS_Enc *encoder, LC3PLUS_FrameDuration frame_dms)
{
    RETURN_IF(encoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_dms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(encoder->lc3_br_set, LC3PLUS_BITRATE_SET_ERROR);
#ifdef CR9_C_ADD_1p25MS
    RETURN_IF(encoder->fs == 8000 && frame_dms == LC3PLUS_FRAME_DURATION_1p25MS, LC3PLUS_SAMPLERATE_ERROR);
#endif
#ifndef CR15_A_LOSSLESS_1p25MS
    RETURN_IF(encoder->hrmode == 2 && frame_dms == LC3PLUS_FRAME_DURATION_1p25MS, LC3PLUS_FRAMEMS_ERROR);
#endif
  
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
                         lc3_scratch_t scratch)
{
    RETURN_IF(!encoder || !input_samples || !output_bytes || !num_bytes || !scratch, LC3PLUS_NULL_ERROR);
    RETURN_IF(null_in_list(input_samples, encoder->channels), LC3PLUS_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3PLUS_ERROR);
    RETURN_IF(!encoder->lc3_br_set, LC3PLUS_BITRATE_UNSET_ERROR);
  
    lc3_scratch_t enc_scratch = lc3_scratch_init( scratch, encoder->scratch_max_size, SCRATCH_BUFFER_ALIGNMENT_BITS, SCRATCH_ALLOCATOR_NORMAL_OPERATION );
  
#ifdef CR14_A_ADD_LOSSLESS_MODE
    Enc_LC3PLUS(encoder, input_samples, bitdepth, output_bytes, enc_scratch, *num_bytes == -1, num_bytes);
#else
    *num_bytes = Enc_LC3PLUS(encoder, input_samples, bitdepth, output_bytes, enc_scratch, *num_bytes == -1);
#endif
  
#  ifdef DEBUG
    max_enc_scratch = MAX( (size_t) max_enc_scratch, enc_scratch->max_alloc_size );
    max_enc_stack_index = MAX( max_enc_stack_index, enc_scratch->max_stack_index );
    if ( enc_scratch->stack_index != 0 )
    {
        fprintf( stderr, "warning: scratch not free at encoder exit:\n" );
        lc3_scratch_print_status( enc_scratch );
    }

#    ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "Encoder scratch status at exit (%d):\n", encoder->scratch_max_size );
    lc3_scratch_print_status( enc_scratch );
#    endif
#  endif
    
#ifndef CR14_A_ADD_LOSSLESS_MODE
    assert(*num_bytes == lc3plus_enc_get_num_bytes(encoder));
#endif
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc16(LC3PLUS_Enc *encoder, int16_t **input_samples, void *output_bytes, int *num_bytes, lc3_scratch_t scratch)
{
    return lc3plus_enc(encoder, (void **)input_samples, 16, output_bytes, num_bytes, scratch);
}

LC3PLUS_Error lc3plus_enc24(LC3PLUS_Enc *encoder, int32_t **input_samples, void *output_bytes, int *num_bytes, lc3_scratch_t scratch)
{
    return lc3plus_enc(encoder, (void **)input_samples, 24, output_bytes, num_bytes, scratch);
}

/* decoder functions *********************************************************/

LC3PLUS_Error lc3plus_dec_init(LC3PLUS_Dec *decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode
#ifdef ENABLE_HR_MODE
                            , int hrmode
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            , int wavFormat
#endif
                            , int32_t* const scratchSize
)
{
    LC3PLUS_Error err = LC3PLUS_OK;
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), LC3PLUS_SAMPLERATE_ERROR);
    RETURN_IF(!lc3plus_channels_supported(channels), LC3PLUS_CHANNELS_ERROR);
    RETURN_IF(!lc3plus_plc_mode_supported(plc_mode), LC3PLUS_PLCMODE_ERROR);
#ifdef ENABLE_HR_MODE
    RETURN_IF(samplerate == 96000 && hrmode == 0, LC3PLUS_HRMODE_ERROR);
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    RETURN_IF( (samplerate < 44100 && hrmode != 0), LC3PLUS_SAMPLERATE_ERROR );
    RETURN_IF( samplerate == 192000 && hrmode <= 1, LC3PLUS_HRMODE_ERROR );
#endif

    err = FillDecSetup(decoder, samplerate, channels, plc_mode
#ifdef ENABLE_HR_MODE
                        , hrmode
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                        , wavFormat
#endif
                       );
                       
    err = lc3plus_dec_get_scratch_size( decoder, scratchSize );
    if ( err != LC3PLUS_OK )
    {
        return err;
    }
    
    memset(decoder, 0, lc3plus_dec_get_size(samplerate, channels, plc_mode));  // clean-up after dummy-decoder run
    
    err = FillDecSetup(decoder, samplerate, channels, plc_mode
#ifdef ENABLE_HR_MODE
                        , hrmode
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                        , wavFormat
#endif
                       );
                       
    decoder->scratch_max_size = *scratchSize;
    decoder->lc3_scratch_initialized = 1;

    return err;
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
LC3PLUS_Error lc3plus_dec_get_lossless_status( LC3PLUS_Dec* const decoder, int16_t* const is_lossless )
{
    RETURN_IF( decoder == NULL, LC3PLUS_NULL_ERROR );
    *is_lossless = decoder->lossless ? decoder->ll_adap_flag : 0;
    return LC3PLUS_OK;
}
#endif

int lc3plus_dec_get_size(int samplerate, int channels, LC3PLUS_PlcMode plc_mode)
{
    RETURN_IF(!lc3plus_samplerate_supported(samplerate), 0);
    RETURN_IF(!lc3plus_channels_supported(channels), 0);
    RETURN_IF(!lc3plus_plc_mode_supported(plc_mode), 0);
    return alloc_decoder(NULL, samplerate, channels);
}

#  ifdef DEBUG
extern int32_t max_enc_scratch;
extern int32_t max_enc_stack_index;
#  endif

#  ifdef DEBUG
extern int32_t max_dec_scratch;
extern int32_t max_dec_stack_index;
#  endif

int lc3plus_dec_get_scratch_size(LC3PLUS_Dec *decoder, int* const scratch_size)
{
    size_t size = 0;
    Word32 minBytes;
    Word32 maxBytes;  // dummy, not used
    LC3PLUS_Error err = LC3PLUS_OK;

    Word16* output16[1];
    Word32 buf_32[LC3PLUS_MAX_SAMPLES] = { 0 };
    UWord8 bytes[LC3PLUS_MAX_BYTES] = { 0 };
    UWord8 tmp_scratch[LC3PLUS_DEC_MAX_USER_SYSTEM_SCRATCH_SIZE]; /* temp/dummy scratch buffer needed to run the encoder with a zero-frame to calculate the max. used scratch size */

    Word16* buf_16 = (Word16*) (void*) buf_32;

    lc3_scratch_t dec_scratch = lc3_scratch_init( (void*) tmp_scratch, LC3PLUS_DEC_MAX_SCRATCH_SIZE, SCRATCH_BUFFER_ALIGNMENT_BITS, SCRATCH_ALLOCATOR_CALCULATE_MAX );

    output16[0] = buf_16;

    err = lc3plus_get_decoder_min_max_bytes( decoder, &minBytes, &maxBytes );
    if ( err )
    {
        return LC3PLUS_SCRATCH_INVALID_ERROR;
    }

    if ( decoder->ep_enabled )
    {
        minBytes = MAX( minBytes, FEC_SLOT_BYTES_MIN );
        minBytes = MAX( maxBytes, FEC_SLOT_BYTES_MAX );
    }

    err = Dec_LC3PLUS( decoder, bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                      &maxBytes, 
#else
                      maxBytes, 
#endif
                      (void**) output16, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                      decoder->wavFormat,
#else
                      16,
#endif
                      dec_scratch, 0 );
    if ( err != LC3PLUS_OK )
    {
        return LC3PLUS_SCRATCH_INVALID_ERROR;
    }
    size = dec_scratch->max_alloc_size;

    err = Dec_LC3PLUS( decoder, bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                      &maxBytes, 
#else
                      maxBytes, 
#endif
                      (void**) output16, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                      decoder->wavFormat,
#else
                      16,
#endif
                      dec_scratch, 1 );
    if ( err != LC3PLUS_DECODE_ERROR )
    {
        return LC3PLUS_SCRATCH_INVALID_ERROR;
    }
    size = MAX( size, dec_scratch->max_alloc_size );

    assert( size <= LC3PLUS_DEC_MAX_USER_SYSTEM_SCRATCH_SIZE );
    assert( size <= LC3PLUS_DEC_MAX_SCRATCH_SIZE );

#  ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "determined scratch size during init: %zu\n", size );
#  endif

    *scratch_size = size;
    return LC3PLUS_OK;
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

LC3PLUS_Error lc3plus_dec_set_frame_dms(LC3PLUS_Dec *decoder, LC3PLUS_FrameDuration frame_dms)
{
    RETURN_IF(decoder == NULL, LC3PLUS_NULL_ERROR);
    RETURN_IF(!lc3plus_frame_size_supported(frame_dms), LC3PLUS_FRAMEMS_ERROR);
    RETURN_IF(decoder->plcMeth == 2 && frame_dms != LC3PLUS_FRAME_DURATION_10MS, LC3PLUS_FRAMEMS_ERROR);
#ifdef CR9_C_ADD_1p25MS
    RETURN_IF(decoder->fs == 8000 && frame_dms == LC3PLUS_FRAME_DURATION_1p25MS, LC3PLUS_SAMPLERATE_ERROR);
#endif
  
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

static LC3PLUS_Error lc3plus_dec(LC3PLUS_Dec *decoder, void *input_bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                 int *num_bytes,                 
#else
                                 int num_bytes, 
#endif
                                 void **output_samples, int bitdepth,
                         lc3_scratch_t scratch, int bfi_ext)
{
    LC3PLUS_Error err = LC3PLUS_OK;
  
    if (bfi_ext == 1)
    {
        RETURN_IF(!decoder || !output_samples || !scratch, LC3PLUS_NULL_ERROR);
    } else {
        RETURN_IF(!decoder || !input_bytes || !output_samples || !scratch, LC3PLUS_NULL_ERROR);
    }

    RETURN_IF(null_in_list(output_samples, decoder->channels), LC3PLUS_NULL_ERROR);
    RETURN_IF(bitdepth != 16 && bitdepth != 24, LC3PLUS_ERROR);
  
    lc3_scratch_t dec_scratch = lc3_scratch_init( scratch, decoder->scratch_max_size, SCRATCH_BUFFER_ALIGNMENT_BITS, SCRATCH_ALLOCATOR_NORMAL_OPERATION );
  
    err = Dec_LC3PLUS(decoder, input_bytes, num_bytes, output_samples, bitdepth, dec_scratch, bfi_ext);
  
#  ifdef DEBUG
    max_dec_scratch = MAX( (size_t) max_dec_scratch, dec_scratch->max_alloc_size );
    max_dec_stack_index = MAX( max_dec_stack_index, dec_scratch->max_stack_index );
    if ( dec_scratch->stack_index != 0 )
    {
        fprintf( stderr, "warning: scratch not free at decoder exit:\n" );
        lc3_scratch_print_status( dec_scratch );
    }
#    ifdef VERBOSE_SCRATCH_ALLOC
    fprintf( stderr, "Decoder scratch status at exit (%d):\n", decoder->scratch_max_size );
    lc3_scratch_print_status( dec_scratch );
#    endif /* VERBOSE_SCRATCH_ALLOC */
#  endif   /* DEBUG */
  
    return err;
}

LC3PLUS_Error lc3plus_dec16(LC3PLUS_Dec *decoder, void *input_bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            int *num_bytes, 
#else
                            int num_bytes, 
#endif
                            int16_t **output_samples, lc3_scratch_t scratch, int bfi_ext)
{
    return lc3plus_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 16, scratch, bfi_ext);
}

LC3PLUS_Error lc3plus_dec24(LC3PLUS_Dec *decoder, void *input_bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                            int *num_bytes, 
#else
                            int num_bytes, 
#endif
                            int32_t **output_samples, lc3_scratch_t scratch, int bfi_ext)
{
    return lc3plus_dec(decoder, input_bytes, num_bytes, (void **)output_samples, 24, scratch, bfi_ext);
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
LC3PLUS_Error lc3_enc_get_max_frame_bytes( int32_t const samplerate, int32_t const channels, int32_t const bitsPerSample, LC3PLUS_FrameDuration const frame_dms, int32_t* const maxBytes )
{
    int32_t fs = 0;
    RETURN_IF( !lc3plus_samplerate_supported( samplerate ) || samplerate < 44100, LC3PLUS_SAMPLERATE_ERROR );
    RETURN_IF( !lc3plus_channels_supported( channels ), LC3PLUS_CHANNELS_ERROR );
    RETURN_IF( bitsPerSample != 16 && bitsPerSample != 24, LC3PLUS_HRMODE_ERROR );
    RETURN_IF( !lc3plus_frame_size_supported( frame_dms ), LC3PLUS_FRAMEMS_ERROR );

    SWITCH(samplerate)
    {
        case 8000:
        fs = 80;
        BREAK;
        case 16000:
        fs = 160;
        BREAK;
        case 24000:
        fs = 240;
        BREAK;
        case 32000:
        fs = 320;
        BREAK;
        case 44100:
        case 48000:
        fs = 480;
        BREAK;
        case 96000:
        fs = 960;
        BREAK;
        #ifdef SUBSET_UUB
        case 192000:
        fs = 1920;
        BREAK;
        #endif 
        default:
        UNUSED(fs);
        assert(0);
    }
    *maxBytes = ((channels * (bitsPerSample + 1) * fs * frame_dms)>>6) + 1;
    return LC3PLUS_OK;
}

LC3PLUS_Error lc3plus_enc_get_lossless_status( LC3PLUS_Enc* const encoder, int16_t* const is_lossless )
{
    RETURN_IF( encoder == NULL, LC3PLUS_NULL_ERROR );
    *is_lossless = encoder->lossless ? encoder->ll_adap_flag : 0;
    return LC3PLUS_OK;
}

#endif
