/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

/*! \file lc3.h
 *  This header provides the API for LC3.
 *
 *  This library is targeting devices with extreme memory limitations, so memory management
 *  must be handeled by the user. This includes allocating memory for the structs and scratch
 *  memory. The structs are persistent between function calls. The scratch memory is working
 *  memory that does not persist between function calls.
 *
 *  The amount of memory needed for various configurations can be obtained from the lc3_*_get_size
 *  and lc3_*_get_scratch_size functions. If memory usage is not a concern the LC3_*_MAX_SIZE
 *  LC3_*_MAX_SCRATCH_SIZE macros can be used for all configurations.
 *
 *  Depending on the build configuration some functions might not be available.
 */

#ifndef LC3PLUS_H
#define LC3PLUS_H



#ifndef _MSC_VER
#include <stdint.h>
#else
typedef __int16 int16_t;
typedef __int32 int32_t;
#endif

/*! Construct version number from major/minor/micro values. */
#define LC3PLUS_VERSION_INT(major, minor, micro) (((major) << 16) | ((minor) << 8) | (micro))

/*! Version number to ensure header and binary are matching. */
#define LC3PLUS_VERSION LC3PLUS_VERSION_INT(1, 6, 3)

/*! Maximum number of supported channels. The actual binary might support
 *  less, use lc3_channels_supported() to check. */
#define LC3PLUS_MAX_CHANNELS 2

/*! Maximum number of samples per channel that can be stored in one LC3 frame. */
#ifdef ENABLE_HR_MODE
#define LC3PLUS_MAX_SAMPLES 960
#else
#define LC3PLUS_MAX_SAMPLES 480
#endif

/*! Maximum number of bytes of one LC3 frame. */
#ifdef ENABLE_HR_MODE
#define LC3PLUS_MAX_BYTES   (625 * LC3PLUS_MAX_CHANNELS)
#else
#define LC3PLUS_MAX_BYTES 870
#endif

/*! Maximum size needed to store encoder state. */
#ifdef ENABLE_HR_MODE
#define LC3PLUS_ENC_MAX_SIZE 12628
#else
#define LC3PLUS_ENC_MAX_SIZE 7226
#endif

/*! Maximum size needed to store decoder state. */
#ifdef ENABLE_HR_MODE
#define LC3PLUS_DEC_MAX_SIZE 42488
#else
#define LC3PLUS_DEC_MAX_SIZE 28446
#endif

/*! Maximum scratch size needed by lc3_enc16() or lc3_enc24().*/
#ifdef ENABLE_HR_MODE
#  define LC3PLUS_ENC_MAX_SCRATCH_SIZE 45624
#else
#  define LC3PLUS_ENC_MAX_SCRATCH_SIZE 6784
#endif

/*! Maximum scratch size needed by LC3PLUS_Dec16() or LC3PLUS_Dec24(). */
#ifdef ENABLE_HR_MODE
#define LC3PLUS_DEC_MAX_SCRATCH_SIZE 59768
#else
#define LC3PLUS_DEC_MAX_SCRATCH_SIZE 27474
#endif
/*! Decoder packet loss concealment mode */
typedef enum
{
    LC3PLUS_PLC_ADVANCED = 1  /*!< Enhanced concealment method */
} LC3PLUS_PlcMode;

/*! Error protection mode. LC3_EP_ZERO differs to LC3PLUS_EP_OFF in that
 *  errors can be detected but not corrected. */
typedef enum
{
    LC3PLUS_EP_OFF    = 0, /*!< Error protection is disabled */
    LC3PLUS_EP_ZERO   = 1, /*!< Error protection with 0 bit correction */
    LC3PLUS_EP_LOW    = 2, /*!< Error protection correcting one symbol per codeword */
    LC3PLUS_EP_MEDIUM = 3, /*!< Error protection correcting two symbols per codeword */
    LC3PLUS_EP_HIGH   = 4  /*!< Error protection correcting three symbols per codeword */
} LC3PLUS_EpMode;

/*! Error protection mode request. On the encoder sidem, LC3PLUS_EPMR_ZERO to LC3PLUS_EPMR_HIGH
 *  can be set. The decoder returns mode requests with different confidences. */
typedef enum
{
    LC3PLUS_EPMR_ZERO      = 0,  /*!< Request no error correction. High confidence if returned by decoder. */
    LC3PLUS_EPMR_LOW       = 1,  /*!< Request low error correction. High confidence if returned by decoder. */
    LC3PLUS_EPMR_MEDIUM    = 2,  /*!< Request medium error correction. High confidence if returned by decoder. */
    LC3PLUS_EPMR_HIGH      = 3,  /*!< Request high error correction. High confidence if returned by decoder. */
    LC3PLUS_EPMR_ZERO_MC   = 4,  /*!< No error correction requested, medium confidence. */
    LC3PLUS_EPMR_LOW_MC    = 5,  /*!< Low error correction requested, medium confidence. */
    LC3PLUS_EPMR_MEDIUM_MC = 6,  /*!< Medium error correction requested, medium confidence. */
    LC3PLUS_EPMR_HIGH_MC   = 7,  /*!< High error correction requested, medium confidence. */
    LC3PLUS_EPMR_ZERO_NC   = 8,  /*!< No error correction requested, unvalidated. */
    LC3PLUS_EPMR_LOW_NC    = 9,  /*!< Low error correction requested, unvalidated. */
    LC3PLUS_EPMR_MEDIUM_NC = 10, /*!< Medium error correction requested, unvalidated. */
    LC3PLUS_EPMR_HIGH_NC   = 11  /*!< High error correction requested, unvalidated. */
} LC3PLUS_EpModeRequest;

/*! Error codes returned by functions. */
typedef enum
{
    LC3PLUS_OK                  = 0,  /*!< No error occurred */
    LC3PLUS_ERROR               = 1,  /*!< Function call failed */
    LC3PLUS_DECODE_ERROR        = 2,  /*!< Frame failed to decode and was concealed */
    LC3PLUS_NULL_ERROR          = 3,  /*!< Pointer argument is null */
    LC3PLUS_SAMPLERATE_ERROR    = 4,  /*!< Invalid samplerate value */
    LC3PLUS_CHANNELS_ERROR      = 5,  /*!< Invalid channels value */
    LC3PLUS_BITRATE_ERROR       = 6,  /*!< Invalid bitrate value */
    LC3PLUS_NUMBYTES_ERROR      = 7,  /*!< Invalid num_bytes value */
    LC3PLUS_PLCMODE_ERROR       = 8,  /*!< Invalid plc_method value */
    LC3PLUS_EPMODE_ERROR        = 9,  /*!< Invalid epmode value */
    LC3PLUS_FRAMEMS_ERROR       = 10, /*!< Invalid frame_ms value */
    LC3PLUS_ALIGN_ERROR         = 11, /*!< Unaligned pointer */
    LC3PLUS_EPMR_ERROR          = 12, /*!< Invalid epmr value */
    LC3PLUS_HRMODE_ERROR        = 13, /*!< Invalid usage of hrmode, sampling rate and frame size */
    LC3PLUS_BITRATE_UNSET_ERROR = 14, /*!< Function called before bitrate has been set */
    LC3PLUS_BITRATE_SET_ERROR   = 15, /*!< Function called after bitrate has been set */
    LC3PLUS_BFI_EXT_ERROR       = 16, /*!< Invalid external bad frame index */

    /* START WARNING */
    LC3PLUS_WARNING    = 17,
    LC3PLUS_BW_WARNING = 18 /*!< Invalid bandwidth cutoff frequency */

} LC3PLUS_Error;

typedef struct LC3PLUS_Enc LC3PLUS_Enc; /*!< Opaque encoder struct. */
typedef struct LC3PLUS_Dec LC3PLUS_Dec; /*!< Opaque decoder struct. */

/*! \addtogroup Misc
 *  \{ */

/*! Return library version number. It should match LC3PLUS_VERSION. */
int lc3plus_version(void);

/*! Tests if the library supports number of channels.
 *
 *  \param[in]  channels    Number of channels.
 *  \return                 1 for true, 0 for false.
 */
int lc3plus_channels_supported(int channels);

/*! Tests if the library supports a sampling rate.
 *
 *  \param[in]  samplerate  Sampling rate
 *  \return                 1 for true, 0 for false
 */
int lc3plus_samplerate_supported(int samplerate);

/*! \}
 *  \addtogroup Encoder
 *  \{ */

/*!
 *  Initialize LC3 encoder.
 *
 *  This function is used to fill a user-allocated encoder struct. This is typically
 *  called once for a samplerate / channel configuration. After init and before encoding
 *  the first frame you must call lc3_enc_set_bitrate().
 *
 *  \param[out] encoder     Pointer to allocated encoder memory. It must have a size provided
 *                          by lc3_enc_get_size() for matching samplerate / channels
 *                          configuration or LC3_ENC_MAX_SIZE.
 *  \param[in]  channels    Number of channels.
 *  \param[in]  samplerate  Input sampling rate. Allowed sampling rates are:
 *                          8000, 16000, 24000, 32000, 44100, 48000
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_init(LC3PLUS_Enc *encoder, int samplerate, int channels
#ifdef ENABLE_HR_MODE
                               , int hrmode
#endif
                              );

/*!
 *  Encode LC3 frame with 16 bit input.
 *
 *  Each call consumes a fixed number of samples. The number of input samples
 *  can be obtained from lc3_enc_get_input_samples().
 *
 *  \param[in]  encoder         Encoder handle initialized by lc3_enc_init().
 *  \param[in]  input_samples   Input samples. The left channel is stored in input_samples[0],
 *                              the right channel in input_samples[1]. The input is not changed
 *                              by the encoder.
 *  \param[out] output_bytes    Output buffer. It must have a at least lc3_enc_get_num_bytes()
 *                              or at most LC3PLUS_MAX_BYTES.
 *  \param[out] num_bytes       Number of bytes written to output_bytes.
 *  \param      scratch         A pointer to an allocated work buffer of at least
 *                              lc3_enc_get_scratch_size() or at most LC3_ENC_SCRATCH_SIZE bytes.
 *                              The buffer does not have to persist, so it can be used for other
 *                              purposes in between calls.
 *  \return                     LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc16(LC3PLUS_Enc *encoder, int16_t **input_samples, void *output_bytes, int *num_bytes, void *scratch);

/*! Encode LC3 frame with 24 bit input.
 *
 *  The input samples are expected to be 24-bit values, sign-extended to 32-bit.
 *  See lc3_enc16() for parameter documentation.
 */
LC3PLUS_Error lc3plus_enc24(LC3PLUS_Enc *encoder, int32_t **input_samples, void *output_bytes, int *num_bytes, void *scratch);

/*! Get the size of the LC3 encoder struct for a samplerate / channel configuration.
 *  If memory is not restricted LC3_ENC_MAX_SIZE can be used for all configurations.
 *
 *  \param[in]  samplerate  Sampling rate.
 *  \param[in]  channels    Number of channels.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3plus_enc_get_size(int samplerate, int channels);

/*! Get the size of the scratch buffer required by lc3_enc16() or lc3_enc24() for the current
 *  encoder configuration. If memory is not restricted, LC3_ENC_MAX_SCRATCH_SIZE can be used for
 *  all configurations.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3plus_enc_get_scratch_size(const LC3PLUS_Enc *encoder);

/*! Get number of samples per channel expected by lc3_enc16() or lc3_enc24().
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Number of samples or 0 on error.
 */
int lc3plus_enc_get_input_samples(const LC3PLUS_Enc *encoder);

/*! Get real internal bitrate of the encoder. It might differ from the requested bitrate due
 *  to error protection or 44.1 kHz input.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Bitrate in bits per second or 0 on error.
 */
int lc3plus_enc_get_real_bitrate(const LC3PLUS_Enc *encoder);

/*! Get the maximum number of bytes produced by lc3_enc16() or lc3_enc24() for the current
 *  bitrate. It should be equal to the num_bytes output of lc3_enc16().
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3plus_enc_get_num_bytes(const LC3PLUS_Enc *encoder);

/*! Set encoder bitrate for all channels.
 *  This function must be called at least once before encoding the first frame, but
 *  after other configuration functions such as lc3_enc_set_frame_dms().
 *
 *  Recommended bitrates for input sampling rates with 10 ms framing:
 *  kHz     | kbps
 *  --------|-----
 *  8       | 24
 *  16      | 32
 *  24      | 48
 *  32      | 64
 *  44.1/48 | 80(voice) 128(music)
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  bitrate     Bitrate in bits per second.
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_bitrate(LC3PLUS_Enc *encoder, int bitrate);

/*! Set encoder Low-frequency effect moded. deactivates LTPF, TNS, NF
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  lfe         LFE mode flag
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_lfe(LC3PLUS_Enc* encoder, int lfe);

/*! Get the encoder delay in number of samples.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \return                 Encoder in samples or 0 on error.
 */
int lc3plus_enc_get_delay(const LC3PLUS_Enc *encoder);

/*! Set the frame length for LC3 encoder in deci milliseconds.
 *  Not all lengths may be enabled, in that case LC3PLUS_FRAMEMS_ERROR is returned.
 *  This function must be called before lc3_enc_set_bitrate(). The decoder must be
 *  configured with LC3PLUS_Dec_set_frame_dms() with the same value.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  frame_ms    Frame length in ms.
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_frame_dms(LC3PLUS_Enc *encoder, int frame_ms);

/*! Set error protection mode. The default is LC3PLUS_EP_OFF. It is possible to switch between
 *  different modees during encoding. Dynamic switching is only allowed between LC3_EP_ZERO,
 *  LC3_EP_LOW, LC3_EP_MEDIUM, and LC3PLUS_EP_HIGH. The the decoder must be notified with
 *  LC3PLUS_Dec_set_ep_enabled() to expect protected data if epmode is other than LC3PLUS_EP_OFF.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  epmode      Error protection mode.
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_ep_mode(LC3PLUS_Enc *encoder, LC3PLUS_EpMode epmode);

/*! Sets error protection mode request transmitted in each channel encoded frame.
 *  The channel coder includes an error protection mode request (EPMR) in every frame.
 *  The EPMR takes value 0, 1, 2, and 3 which request ep modes 1, 2, 3, and 4 from the
 *  decoding device. The EPMR can be retrieved from the channel decoder via the interface
 *  routine LC3PLUS_Dec_get_ep_mode_request().
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  epmr        Error Protection Mode Request
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_ep_mode_request(LC3PLUS_Enc *encoder, LC3PLUS_EpModeRequest epmr);

/*! Set encoder bandwidth to a different value. All frequency bins above the cutoff
 *  frequency are cut off. Allowed frequencies are: 4 kHz, 8 kHz, 12 kHz, 16 kHz and 24 kHz.
 *
 *  \param[in]  encoder     Encoder handle.
 *  \param[in]  bandwidth   Cutoff Frequency in Hz
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_enc_set_bandwidth(LC3PLUS_Enc *encoder, int bandwidth);

/*! \}
 *  \addtogroup Decoder
 *  \{ */


/*!
 *  Initialize LC3 decoder.
 *
 *  This function is used to fill a user-allocated decoder struct. This is typically
 *  called once for a samplerate / channel / plc_mode configuration.
 *
 *  The samplerate and channel arguments must have the same values that were used for encoding.
 *  LC3 does not provide a signalling scheme, transporting these values is the responsibility
 *  of the application.
 *
 *  \param[out] decoder         Pointer to decoder memory. It must have as size of least
 *                              LC3PLUS_Dec_get_size() or at most LC3PLUS_Dec_MAX_SIZE.
 *  \param[in]  samplerate      Bitstream sampling rate.
 *  \param[in]  channels        Bitstream number of channels.
 *  \param[in]  plc_mode        Packet loss concealment mode.
 *
 *  \return                     LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_dec_init(LC3PLUS_Dec *decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode
#ifdef ENABLE_HR_MODE
                               , int hrmode
#endif
                              );

/*!
 *  Decode compressed LC3 frame to 16 bit PCM output.
 *
 *  Each call decodes a fixed number of samples. Use LC3PLUS_Dec_get_output_samples() to obtain this
 *  number. When the input is corrupted and can not be decoded, LC3PLUS_DECODE_ERROR is returned and
 *  packet loss concealment is applied, so the output is still usable.
 *  If error protection is enabled and the errors can be corrected the frame is corrected and
 *  normally decoded. Use LC3PLUS_Dec_get_error_report() to check if errors were corrected.
 *
 *  \param[in]  decoder         Decoder initialized by LC3PLUS_Dec_init().
 *  \param[in]  input_bytes     Input bytes. If error protection is enabled the input bytes can be
 *                              altered when error correction is applied. This is why this buffer
 *                              must be writable.
 *  \param[in]  num_bytes       Number of valid bytes in input_bytes. To signal a lost frame and
 *                              generate concealment output this value must be set to 0.
 *  \param[out] output_samples  Array of pointers to output channel buffers. Each channel buffer
 *                              should provide enough space to hold at most LC3PLUS_MAX_SAMPLES. The
 *                              left channel is stored in output_samples[0], the right channel in
 *                              output_samples[1].
 *  \param      scratch         A pointer to an allocated work buffer of at least
 *                              LC3PLUS_Dec_get_scratch_size() or at most LC3PLUS_Dec_MAX_SCRATCH_SIZE
 *                              bytes. The scratch buffer does not have to persist, so it can be
 *                              used for other purposes in between calls.
 *  \return                     Returns LC3PLUS_OK on success or appropriate error code. Note there is
 *                              a special case for LC3PLUS_DECODE_ERROR where the output is still valid.
 */
LC3PLUS_Error lc3plus_dec16(LC3PLUS_Dec *decoder, void *input_bytes, int num_bytes, int16_t **output_samples, void *scratch,
                    int bfi_ext);

/*! Decode compressed LC3 frame to 24 bit PCM output.
 *
 *  The output samples are 24-bit values, sign-extended to 32-bit.
 *  See LC3PLUS_Dec16() for parameter documentation.
 */
LC3PLUS_Error lc3plus_dec24(LC3PLUS_Dec *decoder, void *input_bytes, int num_bytes, int32_t **output_samples, void *scratch,
                    int bfi_ext);

/*! Get the size of the LC3 decoder struct for a samplerate / channel / plc_mode configuration.
 *  If memory is not restricted LC3PLUS_Dec_MAX_SIZE can be used for all configurations.
 *
 *  \param[in]  channels    Number of channels.
 *  \param[in]  samplerate  Sampling rate.
 *  \param[in]  plc_mode    Packet loss concealment mode.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3plus_dec_get_size(int samplerate, int channels, LC3PLUS_PlcMode plc_mode);

/*! Get the size of the scratch buffer required by LC3PLUS_Dec16() or LC3PLUS_Dec24() for the current
 *  decoder configuration. If memory is not restricted LC3PLUS_Dec_MAX_SCRATCH_SIZE can be used for
 *  all configurations.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Size in bytes or 0 on error.
 */
int lc3plus_dec_get_scratch_size(const LC3PLUS_Dec *decoder);

/*! Get the number of samples per channel produced by LC3PLUS_Dec16() or LC3PLUS_Dec24.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Number of samples or 0 on error.
 */

int lc3plus_dec_get_output_samples(const LC3PLUS_Dec *decoder);

/*! Get the decoder delay in number of samples.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Delay in samples or 0 on error.
 */
int lc3plus_dec_get_delay(const LC3PLUS_Dec *decoder);

/*! Set the frame length for LC3 decoder in deci milliseconds.
 *  Not all lengths may be enabled, in that case LC3PLUS_FRAMEMS_ERROR is returned.
 *  This only works correcly if the encoder was configured with the same vale.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \param[in]  frame_ms    Frame length in ms.
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_dec_set_frame_dms(LC3PLUS_Dec *decoder, int frame_ms);

/*! Enable or disable error protection. Default value is 0 (disabled). If error protection is
 *  enabled, the decoder expects that the frames were encoded with error protection mode
 *  LC3_EP_ZERO or higher.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \param[in]  ep_enabled  1 (or any nonzero) for true, 0 for false.
 *  \return                 LC3PLUS_OK on success or appropriate error code.
 */
LC3PLUS_Error lc3plus_dec_set_ep_enabled(LC3PLUS_Dec *decoder, int ep_enabled);

/*! Retrieves the error protection mode reqeust from channel decoder.
 *
 *  The return value encodes both the error protection mode request (EPMR)
 *  and the confidence of the method by which it was retrieved.
 *
 *  The requested error protection mode is (epmr % 4) + 1, where epmr is the
 *  function's return value. The confidence is specified as follows.
 *
 *  Confidence | Range
 *  -----------|-------------
 *  high       | 0 <= epmr < 4
 *  medium     | 4 <= epmr < 8
 *  no         | 8 <= epmr < 12
 *
 *  When receiving stereo content of separately channel encoded audio frames the
 *  return value is the minimum of two values retrieved from the individual channels.
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Error protection mode reqeust.
 */
LC3PLUS_EpModeRequest lc3plus_dec_get_ep_mode_request(const LC3PLUS_Dec *decoder);

/*! Get the number of corrected bit errors in the last decoded frame. This only works if
 *  error protection is active. If the number of errors is greater than the current error
 *  protection mode can correct, -1 is returned. If the last frame had no errors or the
 *  decoder handle is NULL, 0 is returned,
 *
 *  \param[in]  decoder     Decoder handle.
 *  \return                 Number of corrected bits or -1. See description for details.
 */
int lc3plus_dec_get_error_report(const LC3PLUS_Dec *decoder);
/*! This function returns an set of flags indicating whether the last frame
 *  would have been channel decodable in epmode m, m ranging from 1 to 4. Note that
 *  this information is not available in case the last frame was not channel
 *  decodable in which case the return value is 0. If the last frame would have
 *  been decodable in epmode m,  m-1th of the return value will be 1.
 *  Otherwise, if the frame would not have been decodable or if this information
 *  cannot be retrieved, the m-1th bit of the return value will be 0.
 */
int lc3plus_dec_get_epok_flags(const LC3PLUS_Dec *decoder);

/*! \} */
#endif /* LC3 */
