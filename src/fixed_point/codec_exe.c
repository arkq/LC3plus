/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h" /* needed for basop instrumentation */
#include "lc3plus.h"
#include "tinywavein_c.h"
#include "tinywaveout_c.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

/* struct to hold command line arguments */
typedef struct
{
    char *inputFilename;
    char *outputFilename;
    int   bitrate;
    char *bitrate_file;
    int   encoder_only;
    int   decoder_only;
    int   bipsOut;
    int   formatG192;
    char *configFilenameG192;
    float frame_ms;
    int   hide_counter;
    int   verbose;
    int   plcMeth;
    char *epf;
    int   epmode;
    char *epmode_file;
    char *edf;
    int   ept;
    int   hrmode;
    int   dc;
    char *bandwidth;
    char *channel_coder_vars_file;
    int32_t   startFrame;
    int32_t   stopFrame;
    int32_t   lfe[LC3PLUS_MAX_CHANNELS];
    int32_t   lfeChanCnt;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    int padding;
    int   rel_prio;
    char *ll_shift;
#endif
    int   bipsOut_set;
} Arguments;

/* local helper functions */
static void  parseCmdl(int ac, char **av, Arguments *arg);
static FILE *open_bitstream_reader(const char *file, uint32_t *samplerate, int *bitrate, short *channels,
                                   uint32_t *signal_len, float *frame_ms, int *epmode, int *hrmode, int g192,
                                   const char *file_cfg
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    , int *bitsPerSample
#endif
                                   );
static FILE *open_bitstream_writer(const char *file, uint32_t samplerate, int bitrate, short channels,
                                   uint32_t signal_len, float frame_ms, int epmode, int32_t hrmode, int g192, const char *file_cfg
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    , int bitsPerSample
#endif
);
static void    write_bitstream_frame(FILE *bitstream_file, uint8_t *bytes,
                                     int size,
                                     int g192
#ifdef CR14_A_ADD_LOSSLESS_MODE
, int *nBytes, int nChannels, int lossless
#endif
);

#ifdef CR14_A_ADD_LOSSLESS_MODE
static void
#else
static int
#endif
read_bitstream_frame(FILE *bitstream_file, uint8_t *bytes,
                                    int size,
                                    int g192, int *bfi_ext
#ifdef CR14_A_ADD_LOSSLESS_MODE
, int *nBytes, int nChannels, int lossless
#endif
);

static FILE *  fopen_with_ext(const char *file, const char *ext, const char *mode);
static void    cleanup(void);
static void    safe_fclose(FILE *f);
static int16_t loopy_read16(FILE *f);
static int64_t loopy_read64(FILE *f);
static void    exit_if(int condition, const char *message);
static void    scale_24_to_16(const int32_t *in, int16_t *out, int n);

static void    deinterleave(int32_t *in, int32_t **out, int n, int channels);

static void    interleave_short(int16_t** in, int16_t* out, int32_t n, int32_t channels);
static void    interleave_int  (int32_t** in, int32_t* out, int32_t n, int32_t channels);

/* needed by cleanup function */
static WAVEFILEIN * input_wav;
static WAVEFILEOUT *output_wav;
static FILE *       output_bitstream;
static FILE *       input_bitstream;
#ifdef G192_BITSTREAM_SPLIT
static FILE *       output_bitstream_second_channel;
static FILE *       input_bitstream_secondchannel;
#endif
static FILE *       error_pattern_file;
static FILE *       error_detection_file;
static FILE *       bitrate_switching_file;
static FILE *       epmode_switching_file;
static FILE *bandwidth_switching_file;
static FILE *channel_decoder_debug_file_bfi;
static FILE *channel_decoder_debug_file_epmr;
static FILE *channel_decoder_debug_file_error_report;

#include "license.h" /* provides LICENSE string */

static const char *const USAGE_MESSAGE =
/* Lines must not be longer than this! --------------------------------------->| */
    "Usage: LC3plus [OPTIONS] INPUT OUTPUT BITRATE\n"
    "\n"
    "  INPUT and OUTPUT are wav files, unless another mode is selected in OPTIONS.\n"
    "  BITRATE is specified in bits per second. Alternatively a switching file can\n"
    "  be provided.\n"
#ifdef CR14_A_ADD_LOSSLESS_MODE
    " If no bitrate is given in lossless mode, the codec will determine the bitrate\n"
    " required to code the frame lossless.\n"
#endif
    "\nGeneral options:\n"
    "  -E                      Encode mode. INPUT is a wav file, OUTPUT is a binary file.\n"
    "  -D                      Decode mode. INPUT is a binary file, OUTPUT is a wav file.\n"
    "                          In decode mode the BITRATE parameter is ignored.\n"
    "  -bps NUM                Output bits per sample. NUM must be 16 (default) or 24.\n"
    "  -swf FILE               Use a bitrate switching file instead of fixed bitrate.\n"
    "  -dc NUM                 0: Don't use delay compensation\n"
    "                          1: Compensate delay in decoder (default)\n"
    "                          2: Split delay equally between encoder and decoder\n"
    "  -frame_ms               NUM Frame length in ms. NUM must be 10 (default), 7.5, 5, 2.5 or 1.25.\n"
    "  -bandwidth NUM|FILE     Select audio bandwidth limitation via value in Hz or switching file.\n"
    "                          NUM can be any integer value describing the bandwidth; max NUM=20000 Hz\n"
    "  -q                      Disable frame counter printout\n"
    "  -v                      Verbose switching commands\n"
    "  -y                      StartFrame: frame number where encoding/decoding shall start\n"
    "  -z                      StopFrame: frame number where encoding/decoding shall stop\n"
    "\nFormat options:\n"
    "  -formatG192             Activate G192 bitstream format. A filename.cfg will be used to\n"
    "                          store/load decoder info.\n"
    "  -cfgG192 FILE           Specify a configuration file for G192 bitstream format.\n"
    "\nPLC options:\n"
    "  -epf FILE               Enable packet loss simulation using error pattern from FILE.\n"
    "  -ept                    Use together with -E -epf FILE to signal lost frames within\n"
    "                          the LC3plus bitstream.\n"
    "  -edf FILE               Write error detection pattern to FILE.\n"
    "\nChannel coder options:\n"
    "  -epmode NUM|FILE        Error protection mode. NUM must be one of the following:\n"
    "                          0: Error protection disabled\n"
    "                          1: Minimum error protection, detection only\n"
    "                          2: Moderate error protection\n"
    "                          3: Strong error protection\n"
    "                          4: Maximum error protection\n"
    "  -ep_dbg FILE            Save variables bfi, epmr and error report to binary files\n"
    "                          FILE.bfi, FILE.epmr and FILE.error_report\n"
    "\nLow-frequency effects options:\n"
    "  -lfe NUM[,NUM[,NUM]...] Set lfe flags for all audio channels of the input waveform\n"
    "                          NUM is interpreted as follows:\n"
    "                          0: normal channel\n"
    "                          1: low frequeny enhancement channel\n"
    "                          All channels are treated as normal, if this parameter is omitted.\n"
#ifdef ENABLE_HR_MODE
    "\nHigh resolution mode options:\n"
    "  -hrmode                 Enable high resolution mode.\n"
    "\nLossless mode options:\n"
    "  -lossless               Lossless Mode.\n"
    "  -padding                Enable padding in lossless mode when bitrate is specified.\n"
    "                          Lossless frames are padded to maintain a constant frame size\n"
    "                          according to the specified bitrate.\n"
    "  -rel_prio <0|1>         Use relative priority for residual LSBs. Encoder option,\n"
    "                          requires -lossless. Default: 0.\n"
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    "  -ll_shift NUM           Signal shift option in lossless mode: shift every input sample right by N bits before\n"
    "                          and shift the decoder output left by N bits. NUM = value from 0 to 8. \n"
    "                          The signal shift works only for 24 bit input signals and is limited to 8 bits.\n"
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
    " If no bitrate is given in lossless mode, the codec will determine the bitrate\n"
    " required to code the frame lossless.\n"
#endif
    ;

static const char *const MISSING_ARGUMENT_MESSAGE = "Not enough parameters! Use -h to show help.";

static const char* ERROR_MESSAGE[] = {
    "",                                                                     /* LC3PLUS_OK                  */
    "Function call failed!",                                                /* LC3PLUS_ERROR               */
    "Frame failed to decode and was concealed!",                            /* LC3PLUS_DECODE_ERROR        */
    "Pointer argument is null!",                                            /* LC3PLUS_NULL_ERROR          */
    "Invalid sampling rate!",                                               /* LC3PLUS_SAMPLERATE_ERROR    */
    "Invalid number of channels!",                                          /* LC3PLUS_CHANNELS_ERROR      */
    "Invalid bitrate!",                                                     /* LC3PLUS_BITRATE_ERROR       */
    "Invalid number of bytes!",                                             /* LC3PLUS_NUMBYTES_ERROR      */
    "Invalid ep mode!",                                                     /* LC3PLUS_EPMODE_ERROR        */
    "Invalid frame ms value!",                                              /* LC3PLUS_FRAMEMS_ERROR       */
    "Unaligned pointer!",                                                   /* LC3PLUS_ALIGN_ERROR         */
    "96 kHz sampling rate cannot be used without -hrmode option!",          /* LC3PLUS_HRMODE_ERROR        */
    "Bitrate has not been set!",                                            /* LC3PLUS_BITRATE_UNSET_ERROR */
    "Function can't be called after bitrate was set!",                      /* LC3PLUS_BITRATE_SET_ERROR   */
    "High resolution mode and bandwidth switching are exclusive!",          /* LC3PLUS_HRMODE_BW_ERROR     */
    "Invalid PLC method!",                                                  /* LC3PLUS_PLCMODE_ERROR       */
    "Invalid EPMR value!",                                                  /* LC3PLUS_EPMR_ERROR          */
    "Incorrect padding!",                                                   /* LC3PLUS_PADDING_ERROR       */
    "Incorrect signal shift configuration!",                                /* LC3PLUS_SHIFT_ERROR       */
    "Incorrect frame size during decoding!",                                /* FRAMESIZE_ERROR             */
    "LFE support not available!",                                           /* LC3PLUS_LFE_MODE_NOT_SUPPORTED             */
    "Scratch not allocated!",                                           /* LC3PLUS_SCRATCH_INVALID_ERROR             */
    "Generic Warning",                                                      /* LC3PLUS_WARNING             */
    "Invalid bandwidth frequency!"                                          /* LC3PLUS_BW_WARNING          */
};

#ifdef G192_BITSTREAM_SPLIT
static char *insert_before_ext(const char *filename, const char *suffix)
{
    unsigned int n        = (unsigned int)strlen(filename);
    char *       new_name = (char *)malloc(n + strlen(suffix) + 2);
    const char * dot      = strrchr(filename, '.');
    if (dot && dot != filename)
    {
        int len = (int)(dot - filename);
        snprintf(new_name, n + strlen(suffix) + 2, "%.*s_%s%s", len, filename, suffix, dot);
    }
    else
    {
        snprintf(new_name, n + strlen(suffix) + 2, "%s_%s", filename, suffix);
    }
    return new_name;
}

static char *file_with_ext_split(const char *file, const char *ext)
{
    char *tmp = (char *)malloc(strlen(file) + strlen(ext) + 1);
    sprintf(tmp, "%s%s", file, ext);
    return tmp;
}
#endif

int main(int ac, char **av)
{
    Arguments arg;
    uint32_t  nSamples = 0, nSamplesRead = 0, nSamplesFile = 0xffffffff, sampleRate = 0;
    short     nChannels = 0, bipsIn = 0;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    int nBytes[LC3PLUS_MAX_CHANNELS] = {0};
#else
    int       nBytes = 0;
#endif
    int       real_bitrate = 0, frame = 1, delay = 0;
    int       encoder_size = 0, decoder_size = 0, scratch_size = 0;
    int       bfi_ext = 0;
    LC3PLUS_Enc * encoder = NULL;
    LC3PLUS_Dec * decoder = NULL;
    void *    scratch = NULL;
    LC3PLUS_Error err     = LC3PLUS_OK;
    int32_t   sample_buf[LC3PLUS_MAX_CHANNELS * LC3PLUS_MAX_SAMPLES];
    int32_t   buf_24[LC3PLUS_MAX_CHANNELS * LC3PLUS_MAX_SAMPLES];
    int16_t   buf_16[LC3PLUS_MAX_CHANNELS * LC3PLUS_MAX_SAMPLES];
    uint8_t   bytes[LC3PLUS_MAX_BYTES];
    int       dc2_extra_frame = 0;
    int32_t scratch_size_enc = 0;
    int32_t scratch_size_dec = 0;
#ifdef CR14_A_ADD_LOSSLESS_MODE
    int16_t lossless_status_enc = 0;
    int16_t lossless_status_dec = 0;
#endif
    
    int32_t sample_buf_int[LC3PLUS_MAX_CHANNELS * LC3PLUS_MAX_SAMPLES] = {0};
    int16_t* sample_buf_short = (int16_t*)(void*)sample_buf_int;
    int i;


    /* Parse Command-line */
    printf(LICENSE, LC3PLUS_VERSION >> 16, (LC3PLUS_VERSION >> 8) & 255, LC3PLUS_VERSION & 255);
    parseCmdl(ac, av, &arg);

#ifdef STAMEM_COUNT
    Sta_Mem_Init();
#endif
#ifdef DYNMEM_COUNT
    Dyn_Mem_Init();
#endif

    /* exit handler to clean up resources */
    atexit(cleanup);

    LC3PLUS_FrameDuration frameDuration = LC3PLUS_FRAME_DURATION_10MS;
    switch ((int) (arg.frame_ms*100))
    {
#ifdef CR9_C_ADD_1p25MS
          case 125:
            frameDuration = LC3PLUS_FRAME_DURATION_1p25MS; break;
#endif
          case 250:
            frameDuration = LC3PLUS_FRAME_DURATION_2p5MS; break;
          case 500:
            frameDuration = LC3PLUS_FRAME_DURATION_5MS; break;
          case 750:
            frameDuration = LC3PLUS_FRAME_DURATION_7p5MS; break;
          case 1000:
            frameDuration = LC3PLUS_FRAME_DURATION_10MS; break;
          case LC3PLUS_FRAME_DURATION_UNDEFINED:
            assert(0);
    }

    if (!arg.decoder_only)
    {
        /* Open Input Wav File */
        input_wav = OpenWav(arg.inputFilename, &sampleRate, &nChannels, &nSamplesFile, &bipsIn);
        exit_if(!input_wav, "Error opening wav file!");
        exit_if(bipsIn != 16 && bipsIn != 24, "Only 16 or 24bits per sample are supported!");

#ifdef CR14_A_ADD_LOSSLESS_MODE
        exit_if(arg.hrmode == 2 && arg.bipsOut_set, "-bps option is not supported for lossless mode!");
        if (arg.hrmode == 2)
        {
            arg.bipsOut = bipsIn;
        }
#endif

        /* Check if LFE flag was set for each channel */
        if (arg.lfeChanCnt != 0 && arg.lfeChanCnt != nChannels)
        {
            fprintf(stderr, "%" PRIu32 " channel(s) with -lfe configured, but waveform consists of %i channel(s)\n", arg.lfeChanCnt, nChannels);
            exit_if(1, "Inadequate LFE channel config given!");
        }

        /* Setup Encoder */
        encoder_size = lc3plus_enc_get_size(sampleRate, nChannels);
        encoder      = malloc(encoder_size);
        err          = lc3plus_enc_init(encoder, sampleRate, nChannels
#ifdef ENABLE_HR_MODE
                                        , arg.hrmode
#endif
                                        , arg.lfe, &scratch_size_enc
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                        , bipsIn, arg.padding
#endif
            );
        exit_if(err, ERROR_MESSAGE[err]);

        err = lc3plus_enc_set_frame_dms(encoder, frameDuration);
        exit_if(err, ERROR_MESSAGE[err]);

        err = lc3plus_enc_set_ep_mode(encoder, (LC3PLUS_EpMode)arg.epmode);
        exit_if(err, ERROR_MESSAGE[err]);

        if ( !arg.bitrate_file ) {
            err = lc3plus_enc_set_bitrate(encoder, arg.bitrate);
            exit_if(err, ERROR_MESSAGE[err]);
        }

#ifdef CR14_A_ADD_LOSSLESS_MODE
        /* get to know the max bitrate using lc3_enc_get_max_frame_bytes()
           when fully lossless mode and arg.bitrate is set to 0 */
        /* Code duplication from update_enc_bitrate() @ setup_enc_lc3plus.c */
        if (arg.hrmode == 2 && arg.bitrate==0) {
            int bitsPerSample = (bipsIn == 24) ? 24 : 16;
            int maxBytes = 0;
            /*arg.bitrate = nChannels * (bitsPerSample+3) * sampleRate;*/
            err = lc3_enc_get_max_frame_bytes( sampleRate, nChannels, bitsPerSample, frameDuration, &maxBytes );
            exit_if( err, ERROR_MESSAGE[err] );
            arg.bitrate = maxBytes * 8 * 800 / frameDuration + 1; /* add 1 to avoid rounding problems for some configs */
            if ( sampleRate == 44100 )
            {
                /* scale back to real bitrate for 44.1 kHz */
                arg.bitrate = (arg.bitrate * 441) / 480 + 1;
            }
        }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
        if (arg.rel_prio)
        {
            lc3plus_enc_set_relative_priority(encoder, 1);
        }
#endif

#ifdef CR14_A_ADD_LOSSLESS_MODE
        if ( arg.ll_shift )
        {
            int initial_shift = atoi( arg.ll_shift );
            err = lc3plus_enc_set_ll_shift( encoder, initial_shift );
            exit_if( err, ERROR_MESSAGE[err] );
        }
#endif

        delay        = arg.dc ? lc3plus_enc_get_delay(encoder) / arg.dc : 0;
        nSamples     = lc3plus_enc_get_input_samples(encoder);
        real_bitrate = lc3plus_enc_get_real_bitrate(encoder);

        if (arg.bandwidth && atoi(arg.bandwidth) == 0)
        {
            bandwidth_switching_file = fopen(arg.bandwidth, "rb");
            exit_if(bandwidth_switching_file == NULL, "Error opening bandwidth switching file!");
            puts("Using bandwidth switching file!");
        }
    }
    else /* !arg->decoder_only */
    {
        /* Open Input Bitstream File */
#ifdef CR14_A_ADD_LOSSLESS_MODE
        int bipsOutBitstream = 0;
#endif
                                    
        input_bitstream =
            open_bitstream_reader(arg.inputFilename, &sampleRate, &arg.bitrate, &nChannels, &nSamplesFile,
                                  &arg.frame_ms, &arg.epmode, &arg.hrmode, arg.formatG192, arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                  , &bipsOutBitstream
#endif
                                 );
#ifndef ENABLE_HR_MODE
        exit_if(arg.hrmode, "HR bitstreams not supported!");
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
        exit_if(arg.hrmode == 2 && arg.bipsOut_set, "-bps option is not supported for lossless mode!");
        if (arg.hrmode == 2 && bipsOutBitstream)
        {
            arg.bipsOut = bipsOutBitstream;
        }
#endif
#ifdef G192_BITSTREAM_SPLIT
        if (arg.formatG192 && arg.hrmode == 2 && nChannels > 1)
        {
            if (!arg.configFilenameG192)
            {
                arg.configFilenameG192 = file_with_ext_split(arg.inputFilename, ".cfg");
            }
            if (input_bitstream)
            {
                fclose(input_bitstream);
                input_bitstream = NULL;
            }
            const char *inputFilename_L = insert_before_ext(arg.inputFilename, "L");
            const char *inputFilename_R = insert_before_ext(arg.inputFilename, "R");
            input_bitstream = open_bitstream_reader(inputFilename_L, &sampleRate, &arg.bitrate, &nChannels,
                                                    &nSamplesFile, &arg.frame_ms, &arg.epmode, &arg.hrmode,
                                                    arg.formatG192, arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                                    , &bipsOutBitstream
#endif
                                                   );
            exit_if(!input_bitstream, "Error opening L channel bitstream file!");
            input_bitstream_secondchannel = open_bitstream_reader(inputFilename_R, &sampleRate, &arg.bitrate,
                                                                  &nChannels, &nSamplesFile, &arg.frame_ms,
                                                                  &arg.epmode, &arg.hrmode, arg.formatG192,
                                                                  arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                                                  , &bipsOutBitstream
#endif
                                                                 );
            exit_if(!input_bitstream_secondchannel, "Error opening R channel bitstream file!");
        }
        else
#endif
        {
            exit_if(!input_bitstream, "Error opening bitstream file!");
        }
    }

    if (!arg.encoder_only)
    {
        /* Setup Decoder */
        decoder_size = lc3plus_dec_get_size(sampleRate, nChannels, (LC3PLUS_PlcMode)arg.plcMeth);
        decoder      = malloc(decoder_size);
        err          = lc3plus_dec_init(decoder, sampleRate, nChannels, (LC3PLUS_PlcMode)arg.plcMeth
#ifdef ENABLE_HR_MODE
                                        , arg.hrmode
#endif
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                        , arg.bipsOut
#endif
                                        , &scratch_size_dec
            );
        exit_if(err, ERROR_MESSAGE[err]);
      
        switch ((int) (arg.frame_ms*100))
        {
#ifdef CR9_C_ADD_1p25MS
          case 125:
            frameDuration = LC3PLUS_FRAME_DURATION_1p25MS; break;
#endif
          case 250:
            frameDuration = LC3PLUS_FRAME_DURATION_2p5MS; break;
          case 500:
            frameDuration = LC3PLUS_FRAME_DURATION_5MS; break;
          case 750:
            frameDuration = LC3PLUS_FRAME_DURATION_7p5MS; break;
          case 1000:
            frameDuration = LC3PLUS_FRAME_DURATION_10MS; break;
          case LC3PLUS_FRAME_DURATION_UNDEFINED:
            assert(0);
        }

        err = lc3plus_dec_set_frame_dms(decoder, frameDuration);
        exit_if(err, ERROR_MESSAGE[err]);

        err = lc3plus_dec_set_ep_enabled(decoder, arg.epmode != 0);
        exit_if(err, ERROR_MESSAGE[err]);

        delay    = arg.dc ? lc3plus_dec_get_delay(decoder) / arg.dc : 0;
        nSamples = lc3plus_dec_get_output_samples(decoder);

        /* Open Output Wav File */
        output_wav = CreateWav(arg.outputFilename, sampleRate, nChannels, arg.bipsOut);
        exit_if(!output_wav, "Error creating wav file!");
    }
    else /* !arg->encoder_only */
    {
        /* Open Output Bitstream File */
#ifdef G192_BITSTREAM_SPLIT
        if (arg.formatG192 && arg.hrmode == 2 && nChannels > 1)
        {
            if (!arg.configFilenameG192)
            {
                arg.configFilenameG192 = file_with_ext_split(arg.outputFilename, ".cfg");
            }
            const char *outputFilename_L = insert_before_ext(arg.outputFilename, "L");
            const char *outputFilename_R = insert_before_ext(arg.outputFilename, "R");
            output_bitstream = open_bitstream_writer(outputFilename_L, sampleRate, arg.bitrate,
                                                     nChannels, nSamplesFile, arg.frame_ms, arg.epmode,
                                                     arg.hrmode, arg.formatG192, arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                                     , bipsIn
#endif
                                                    );
            exit_if(!output_bitstream, "Error creating L channel bitstream file!");
            output_bitstream_second_channel = open_bitstream_writer(outputFilename_R, sampleRate, arg.bitrate,
                                                                    nChannels, nSamplesFile, arg.frame_ms,
                                                                    arg.epmode, arg.hrmode, arg.formatG192,
                                                                    arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                                                    , bipsIn
#endif
                                                                   );
            exit_if(!output_bitstream_second_channel, "Error creating R channel bitstream file!");
        }
        else
#endif
        {
            output_bitstream = open_bitstream_writer(arg.outputFilename, sampleRate,
                                                                         arg.bitrate
                                                     , nChannels, nSamplesFile,
                                                     arg.frame_ms, arg.epmode, arg.hrmode, arg.formatG192, arg.configFilenameG192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                                    , bipsIn
#endif
            );
            exit_if(!output_bitstream, "Error creating bitstream file!");
        }
    }

    /* open auxillary files */
    if (arg.epf)
    {
        error_pattern_file = fopen(arg.epf, "rb");
        exit_if(!error_pattern_file, "Error opening error pattern file!");
    }
    if (arg.bitrate_file)
    {
        bitrate_switching_file = fopen(arg.bitrate_file, "rb");
        exit_if(!bitrate_switching_file, "Error opening bitrate switching file!");
    }
    if (arg.epmode_file)
    {
        epmode_switching_file = fopen(arg.epmode_file, "rb");
        exit_if(epmode_switching_file == NULL, "Error opening epmode switching file!");
    }
    if (arg.edf)
    {
        error_detection_file = fopen(arg.edf, "wb");
        exit_if(!error_detection_file, "Error creating error detection file!");
    }
    if (arg.channel_coder_vars_file)
    {
        channel_decoder_debug_file_bfi          = fopen_with_ext(arg.channel_coder_vars_file, ".bfi", "wb");
        channel_decoder_debug_file_epmr         = fopen_with_ext(arg.channel_coder_vars_file, ".epmr", "wb");
        channel_decoder_debug_file_error_report = fopen_with_ext(arg.channel_coder_vars_file, ".error_report", "wb");
        exit_if(!channel_decoder_debug_file_bfi || !channel_decoder_debug_file_epmr ||
                    !channel_decoder_debug_file_error_report,
                "Error creating channel decoder debug files!");
    }
    
    scratch_size = MAX( scratch_size_dec, scratch_size_enc );

    scratch = malloc( scratch_size );
    exit_if( !scratch, "Failed to allocate scratch memory!" );

#ifndef NO_SCRATCH_STATS
    UWord32 *sc32 = (UWord32*)scratch;
    for( i = 0; i < scratch_size>>2; i++)
    {
        sc32[i] = 0xDEADCAFE;
    }
#endif

#ifdef STAMEM_COUNT
    Sta_Mem_Add("Encoder", encoder_size);
    Sta_Mem_Add("Decoder", decoder_size);
#endif

    /* Print info */
    printf("Encoder size:       %i\n", encoder_size);
    printf("Decoder size:       %i\n", decoder_size);
    printf("Scratch size:       %i\n", scratch_size);
    printf("Sample rate:        %i\n", sampleRate);
    printf("Channels:           %i\n", nChannels);
    printf("Signal length:      %u\n", nSamplesFile);
    printf("Frame length:       %i\n", nSamples);
    printf("Output format:      %i bits\n", arg.bipsOut);
    printf("Target bitrate:     %i\n", arg.bitrate);
    if (!arg.decoder_only)
    {
        printf("Real bitrate:       %i\n\n", real_bitrate);
    }
    printf("Bandwidth cutoff:       %s\n", arg.bandwidth ? arg.bandwidth : "-");
    printf("PLC mode:               %i\n", arg.plcMeth);
#ifdef ENABLE_HR_MODE
    printf("\n");
    printf("Improved precision coding is active. (ENABLE_HR_MODE preprocessor define active)\n");
    printf("High resolution mode:  %8s\n", arg.hrmode ? "on" : "off");
#else
    printf("\n");
    printf("Normal precision coding is active. (ENABLE_HR_MODE preprocessor define inactive)\n");
#endif
    printf("\n");

    /* delay compensation */
    if (arg.dc == 2 && !arg.decoder_only)
    {
        ReadWavInt(input_wav, sample_buf, nChannels * delay, &nSamplesRead);
    }

    setFrameRate(sampleRate, nSamples);
    Init_WMOPS_counter();

    /* Encoder + Decoder loop */
    while (1)
    {
        if (!arg.hide_counter)
        {
            if ((frame >= arg.startFrame && frame <= arg.stopFrame && arg.stopFrame != 0) || (arg.stopFrame == 0 && arg.startFrame == 0))
            {
                printf("\rProcessing frame %i", frame);
                fflush(stdout);
            }
        }
        if (!arg.decoder_only)
        {
            /* Encoder */
            int32_t *input24[] = {buf_24, buf_24 + nSamples};

            /* read bitrate switching file and set new bitrate */
            if (bitrate_switching_file)
            {
                int32_t new_bitrate = (int32_t)(loopy_read64(bitrate_switching_file));
                if (arg.verbose && encoder->bitrate != new_bitrate * nChannels)
                {
                    printf("Switching rate from %d to %d\n", encoder->bitrate, new_bitrate * nChannels);
                }
                err = lc3plus_enc_set_bitrate(encoder, new_bitrate * nChannels);
                exit_if(err, ERROR_MESSAGE[err]);
            }
            /* read epmode switching file and set error protection */
            if (epmode_switching_file)
            {
                /* gen_rate_profile tool can only write values starting from 100 */
                int16_t epmode = (int16_t)(loopy_read64(epmode_switching_file) / 100);
                assert(epmode > 0 && epmode < 5); /* epmode 0 not allowed when switching */
                if (arg.verbose && encoder->epmode != epmode)
                {
                    printf("Switching epmode from %d to %d\n", encoder->epmode, epmode);
                }
                err = lc3plus_enc_set_ep_mode(encoder, (LC3PLUS_EpMode)epmode);
                exit_if(err, ERROR_MESSAGE[err]);
            }
            /* read bandwidth switching file and set bandwidth */
            if (arg.bandwidth || bandwidth_switching_file)
            {
                int32_t bw =
                    bandwidth_switching_file ? (int32_t)loopy_read64(bandwidth_switching_file) : atoi(arg.bandwidth);
                int32_t bw_old = encoder->bandwidth;
                err            = lc3plus_enc_set_bandwidth(encoder, bw);
                if (arg.verbose && bw_old != bw && err == LC3PLUS_OK)
                {
                    printf("Switching bandwidth from %i to %i\n", bw_old, bw);
                }
                exit_if(err, ERROR_MESSAGE[err]);
            }

            /* read audio data */
            ReadWavInt(input_wav, sample_buf, nSamples * nChannels, &nSamplesRead);
            /* zero out rest of last frame */
            memset(sample_buf + nSamplesRead, 0, (nSamples * nChannels - nSamplesRead) * sizeof(sample_buf[0]));
            if (frame < arg.startFrame)
                goto while_end;

            if ( (arg.hrmode == 0 && arg.dc != 2) || arg.dc == 0)
            {
                if (nSamplesRead == 0)
                {
                    break;
                }
            }
            else
            {
                if ( arg.hrmode > 1 )
                {
                    if ( nSamplesRead == 0 )
                    {
                        if ( dc2_extra_frame == 1 )
                        {
                            break;
                        }
                        dc2_extra_frame = 1;
                    }
                }
                else
                {
                    if (nSamplesRead != (nSamples * nChannels))
                    {
                        Word16 padded_samples = ((nSamples * nChannels) - nSamplesRead) / nChannels;
                        Word16 delay_samples  = lc3plus_enc_get_delay(encoder) / 2;

                        if (padded_samples >= delay_samples)
                        {
                            if (dc2_extra_frame == 1)
                            {
                                break;
                            }
                            dc2_extra_frame = 1;
                        }
                    }
                }
            }

            if (arg.ept && loopy_read16(error_pattern_file))
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                for ( i = 0; i < nChannels; i++ )
                {
                    nBytes[i] = -1;
                }
#else
                nBytes = -1; /* tell encoder packet is lost and trigger PLC */
#endif
            }

            /* deinterleave channels */
            deinterleave(sample_buf, input24, nSamples, nChannels);

            /* encode */
            if (bipsIn == 24)
            {
                err = lc3plus_enc24(encoder, input24, bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    nBytes, 
#else
                                    &nBytes, 
#endif
                                    scratch);
            }
            else
            {
                int16_t *input16[] = {buf_16, buf_16 + nSamples};
                scale_24_to_16(buf_24, buf_16, nSamples * nChannels);
                err = lc3plus_enc16(encoder, input16, bytes, 
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    nBytes, 
#else
                                    &nBytes, 
#endif
                                    scratch);
            }

            exit_if(err, ERROR_MESSAGE[err]);
#ifdef CR14_A_ADD_LOSSLESS_MODE
            lc3plus_enc_get_lossless_status( encoder, &lossless_status_enc );
            UNUSED(lossless_status_enc);
#endif
        }
        else /* !arg.decoder_only */
        {
            /* Read bitstream */
#ifdef CR14_A_ADD_LOSSLESS_MODE
#  ifdef G192_BITSTREAM_SPLIT
            if (arg.formatG192 && arg.hrmode == 2 && nChannels > 1)
            {
                read_bitstream_frame(input_bitstream, bytes, sizeof(bytes), arg.formatG192, &bfi_ext, &nBytes[0], 1, 1);
                if (nBytes[0] >= 0)
                {
                    read_bitstream_frame(input_bitstream_secondchannel, bytes + nBytes[0],
                                         (int)sizeof(bytes) - nBytes[0], arg.formatG192, &bfi_ext, &nBytes[1], 1, 1);
                }
                else
                {
                    nBytes[1] = -1;
                }
            }
            else
#  endif
            read_bitstream_frame(input_bitstream, bytes, sizeof(bytes), arg.formatG192, &bfi_ext, nBytes, nChannels, arg.hrmode == 2);

            int exit_loop = 0;
            for ( i = 0; i < nChannels; i++ )
            {
                if ( nBytes[i] < 0 )
                {
                    exit_loop = 1;
                }
            }
            
            if (exit_loop)
            {
                break;
            }
#else
            nBytes = read_bitstream_frame(input_bitstream, bytes, sizeof(bytes), arg.formatG192, &bfi_ext);
            if (nBytes < 0)
            {
                break;
            }
#endif

        }

        if (!arg.encoder_only)
        {
            /* Decoder */
            /* read error pattern */
            if (error_pattern_file && loopy_read16(error_pattern_file))
            {
#ifdef CR14_A_ADD_LOSSLESS_MODE
                for ( i = 0; i < nChannels; i++ )
                {
                    nBytes[i] = 0;
                }
#else
                nBytes = 0; /* tell decoder packet is lost and needs to be concealed */
#endif
            }
            
            int16_t* output16[LC3PLUS_MAX_CHANNELS];
            int32_t* output24[LC3PLUS_MAX_CHANNELS];
            
            for (i = 0; i < nChannels; i++)
            {
                output16[i] = buf_16 + i * nSamples;
                output24[i] = buf_24 + i * nSamples;
            }

            /* Run Decoder */
            if (arg.bipsOut == 24)
            {
                if (bfi_ext == 1)
                {
                    err = lc3plus_dec24(decoder, NULL, nBytes, output24, scratch, bfi_ext);
                } else {
                    err = lc3plus_dec24(decoder, bytes, nBytes, output24, scratch, bfi_ext);
                }
            }
            else
            {
                if (bfi_ext == 1)
                {
                    err = lc3plus_dec16(decoder, NULL, nBytes, output16, scratch, bfi_ext); 
                } else {
                    err = lc3plus_dec16(decoder, bytes, nBytes, output16, scratch, bfi_ext);
                }
            }
            exit_if(err && err != LC3PLUS_DECODE_ERROR, ERROR_MESSAGE[err]);

            /* write error detection to file */
            if (error_detection_file != NULL)
            {
                int16_t tmp = (err == LC3PLUS_DECODE_ERROR);
                fwrite(&tmp, 2, 1, error_detection_file);
            }

            /* write bfi, empr and error report to files */
            if (arg.channel_coder_vars_file)
            {
                int16_t tmp = (err == LC3PLUS_DECODE_ERROR) ? 1 : LC3PLUS_OK;
                fwrite(&tmp, 2, 1, channel_decoder_debug_file_bfi);
                tmp = lc3plus_dec_get_ep_mode_request(decoder);
                fwrite(&tmp, 2, 1, channel_decoder_debug_file_epmr);
                tmp = lc3plus_dec_get_error_report(decoder);
                fwrite(&tmp, 2, 1, channel_decoder_debug_file_error_report);
            }

#ifdef CR14_A_ADD_LOSSLESS_MODE
            lc3plus_dec_get_lossless_status( decoder, &lossless_status_dec );
            UNUSED(lossless_status_dec);
#endif

                uint32_t out_samples = MIN((uint32_t)nSamples - delay, nSamplesFile);
            
                switch (arg.bipsOut)
                {
                case 16:
                    interleave_short(output16, sample_buf_short, nSamples, nChannels);
                    WriteWavShort(output_wav, sample_buf_short + delay * nChannels, out_samples * nChannels);
                    break;
                case 24:
                    interleave_int  (output24, sample_buf_int  , nSamples, nChannels);
                    WriteWavLong (output_wav, sample_buf_int   + delay * nChannels, out_samples * nChannels);
                    break;
                }
            
                nSamplesFile -= out_samples;
            
                delay = 0;
        }
        else /* !arg.encoder_only */
        {
#ifdef CR14_A_ADD_LOSSLESS_MODE
#  ifdef G192_BITSTREAM_SPLIT
            if (arg.formatG192 && arg.hrmode == 2 && nChannels > 1)
            {
                write_bitstream_frame(output_bitstream, bytes, nBytes[0],
                                      arg.formatG192, &nBytes[0], 1, 1);
                write_bitstream_frame(output_bitstream_second_channel, bytes + nBytes[0], nBytes[1],
                                      arg.formatG192, &nBytes[1], 1, 1);
            }
            else
#  endif
            {
                int numBytesTotal = nChannels == 2 ? nBytes[0] + nBytes[1] : nBytes[0];
                write_bitstream_frame(output_bitstream, bytes, numBytesTotal,
                                      arg.formatG192, nBytes, nChannels, arg.hrmode == 2);
            }
#else
            write_bitstream_frame(output_bitstream, bytes, nBytes, arg.formatG192);
#endif
        }

 while_end:
        frame++;

        if (frame > arg.stopFrame && arg.stopFrame != 0)
        {
            break;
        }

        BASOP_frame_update();
    }

    if (!arg.encoder_only && nSamplesFile > 0 && nSamplesFile <= (uint32_t)nSamples)
    {
        memset(sample_buf, 0, (nSamplesFile * nChannels) * sizeof(sample_buf[0]));
        WriteWavLong(output_wav, sample_buf, nSamplesFile * nChannels);
    }

    puts("\nProcessing done!");
    if (output_wav)
    {
        printf("%i samples clipped!\n", (int) output_wav->clipCount);
    }

    free(encoder);
    free(decoder);

#ifndef NO_SCRATCH_STATS
    Word32 usedScratch = 0;
    Word32 unusedScratch = 0;
    Word32 fragmentedScratch = 0;
    Word32 maxUsedScratchIndex = 0;

    // initial continuous block:
    while(usedScratch < (scratch_size >> 2) && 0xDEADCAFE != sc32[usedScratch])
    {
        usedScratch++;
    }
    maxUsedScratchIndex = usedScratch;
    // later scratch usage:
    for( i = usedScratch; i < scratch_size>>2; i++)
    {
        if(0xDEADCAFE == sc32[i])
        {
            unusedScratch++;
        }
        else
        {
            fragmentedScratch++;
            maxUsedScratchIndex = i;
        }
    }
    printf("\nScratch statistics (granularity: 4 bytes):\n");
    printf("Scratch allocated: %d bytes\n", scratch_size);
    printf("Scratch used (initial continuous block): %d bytes\n", usedScratch*4);
    printf("Scratch used (fragmented later usage): %d bytes\n", fragmentedScratch*4);
    printf("Scratch used (initial continuous block + fragmented later usage): %d bytes\n", (usedScratch+fragmentedScratch)*4);
    printf("Scratch gap: %d bytes\n", ((maxUsedScratchIndex+1-usedScratch-fragmentedScratch)*4));
    printf("Scratch unused: %d bytes\n", unusedScratch*4);
    printf("Scratch occupied: %d bytes\n\n", (maxUsedScratchIndex+1)*4);
#endif

    free(scratch);

#ifdef WMOPS
    BASOP_end;
#else
    BASOP_end_noprint;
#endif
#ifdef STAMEM_COUNT
    Sta_Mem_Exit();
#endif
#ifdef DYNMEM_COUNT
    Dyn_Mem_Exit();
#endif
}

/* open file with extra extension */
static FILE *fopen_with_ext(const char *file, const char *ext, const char *mode)
{
    FILE *f   = NULL;
    char *tmp = malloc(strlen(file) + strlen(ext) + 1);
    sprintf(tmp, "%s%s", file, ext);
    f = fopen(tmp, mode);
    free(tmp);
    return f;
}

/* close file ignoring NULL pointer */
static void safe_fclose(FILE *f)
{
    if (f != NULL)
        fclose(f);
}

/* ensure clean exit so valgrind & co. don't complain */
void cleanup(void)
{
    CloseWavIn(input_wav);
    CloseWav(output_wav);
    safe_fclose(output_bitstream);
    safe_fclose(input_bitstream);
#ifdef G192_BITSTREAM_SPLIT
    safe_fclose(output_bitstream_second_channel);
    safe_fclose(input_bitstream_secondchannel);
#endif
    safe_fclose(error_pattern_file);
    safe_fclose(error_detection_file);
    safe_fclose(bitrate_switching_file);
    safe_fclose(epmode_switching_file);
    safe_fclose(bandwidth_switching_file);
    safe_fclose(channel_decoder_debug_file_bfi);
    safe_fclose(channel_decoder_debug_file_epmr);
    safe_fclose(channel_decoder_debug_file_error_report);
}

static void parseCmdl(int ac, char **av, Arguments *arg)
{
    int pos = 1;
    memset(arg, 0, sizeof(*arg));
    arg->bipsOut  = 16;
    arg->frame_ms = 10;
    arg->dc       = 1;
    arg->plcMeth = LC3PLUS_PLC_ADVANCED;
    exit_if(ac <= 1, USAGE_MESSAGE);

    /* parse options in any order */
    for (; pos < ac && av[pos][0] == '-'; pos++)
    {
        if (!strcmp(av[pos], "-h"))
        {
            puts(USAGE_MESSAGE);
            exit(0);
        }
        if (!strcmp(av[pos], "-q"))
        {
            arg->hide_counter = 1;
        }
        if (!strcmp(av[pos], "-v"))
        {
            arg->verbose = 1;
        }
        if (!strcmp(av[pos], "-y") && pos + 1 < ac)
        {
            arg->startFrame = atoi(av[++pos]);
            printf("Start Frame: %d!\n", arg->startFrame);
            continue;
        }
        if (!strcmp(av[pos], "-z") && pos + 1 < ac)
        {
            arg->stopFrame = atoi(av[++pos]);
            printf("Stop Frame: %d!\n", arg->stopFrame);
            continue;
        }
        if (!strcmp(av[pos], "-E"))
        {
            arg->encoder_only = 1;
            puts("Using only encoder!");
        }
        if (!strcmp(av[pos], "-D"))
        {
            arg->decoder_only = 1;
            puts("Using only decoder!");
        }
        if (!strcmp(av[pos], "-formatG192"))
        {
            arg->formatG192 = 1;
            puts("Reading/writing bitstream in G192 format!");
        }
        if (!strcmp(av[pos], "-cfgG192") && pos + 1 < ac)
        {
            arg->configFilenameG192 = av[++pos];
            puts("Using user defined configuration file for G192 bitstream format!");
        }
        /* error pattern */
        if (!strcmp(av[pos], "-epf") && pos + 1 < ac)
        {
            arg->epf = av[++pos];
            puts("Using error pattern file for frame loss simulation!");
        }
        /* trigger PLC with special decoder modes */
        if (!strcmp(av[pos], "-ept"))
        {
            arg->ept = 1;
            puts("Simulating frame loss by writing reserved values into the LC3plus bitstream!");
        }
        /* Bits per sample */
        if (!strcmp(av[pos], "-bps") && pos + 1 < ac)
        {
            arg->bipsOut = atoi(av[++pos]);
            exit_if(arg->bipsOut != 16 && arg->bipsOut != 24,
                    "Only 16, or 24 bits per sample are supported!");
            arg->bipsOut_set = 1;
        }

        /* delay compensation */
        if (!strcmp(av[pos], "-dc") && pos + 1 < ac)
        {
            arg->dc = atoi(av[++pos]);
            exit_if(arg->dc < 0 || arg->dc > 2, "dc must be 0, 1 or 2!");
        }
        /* select bandwidth */
        if (!strcmp(av[pos], "-bandwidth") && pos + 1 < ac)
        {
            arg->bandwidth = av[++pos];
        }
#ifdef CR14_A_ADD_LOSSLESS_MODE

        if (!strcmp(av[pos], "-ll_shift") && pos + 1 < ac)
        {
            arg->ll_shift = av[++pos];
        }
#endif
        /* frame length in ms */
        if (!strcmp(av[pos], "-frame_ms") && pos + 1 < ac)
        {
            arg->frame_ms = (float)atof(av[++pos]);
        }
        
#ifdef ENABLE_HR_MODE
        if (!strcmp(av[pos], "-hrmode"))
        {
            arg->hrmode = 1;
            printf("Enabling hrmode!\n");
        }
#endif
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if ( !strcmp( av[pos], "-lossless" ) )
        {
            arg->hrmode = 2;
            printf( "Enabling lossless mode; also enables hrmode!\n" );
        }
      
        /* padding */
        if ( !strcmp( av[pos], "-padding" ) )
        {
            arg->padding = 1;
            printf( "Enable padding for lossless mode!\n" );
        }

        if ( !strcmp( av[pos], "-rel_prio" ) )
        {
            arg->rel_prio = atoi( av[++pos] );
            printf( "Relative priority for residual LSBs: %s\n", arg->rel_prio ? "on" : "off" );
        }
#endif

            /* lfe mode */
            if (lc3_enc_supported_lfe())
            {
                if (!strcmp(av[pos], "-lfe"))
                {
                    char *lfe_string = av[++pos];
                    char *some_pointer;
                    char  lfeChans[LC3PLUS_MAX_CHANNELS] = {0};
                    int   lfeChansCnt = 0;
                    int   i;
                    some_pointer = strtok (lfe_string,",");

#ifdef DEBUG
/* Allow old-style usage of -lfe flag for ETSI vs Trunk tests */
                    if ((strstr(some_pointer, ".wav") != NULL) || (strstr(some_pointer, ".bin") != NULL)) {
                        for (i = 0; i < LC3PLUS_MAX_CHANNELS; i ++)
                        {
                            arg->lfe[i] = 1;
                        }

                        pos--;
                        continue;
                    }
#endif

                    while (some_pointer != NULL) {
                        if (arg->lfeChanCnt >= LC3PLUS_MAX_CHANNELS)
                        {
                            fprintf(stderr, "at least %" PRId32 "channels with -lfe configured, but only %i channels are supported\n", arg->lfeChanCnt+1, LC3PLUS_MAX_CHANNELS);
                            exit_if(1, "More LFE parameters given than supported channels exist!");
                        }
                        arg->lfe[arg->lfeChanCnt] = atoi(some_pointer);

                        if (arg->lfe[arg->lfeChanCnt]) lfeChans[lfeChansCnt++] = arg->lfeChanCnt;

                        some_pointer = strtok (NULL, ",");
                        arg->lfeChanCnt++;
                    }

                    for (i = 0; i < lfeChansCnt; i ++)
                    {
                        printf("Enabling LFE mode for channel %i\n", lfeChans[i]);
                    }

                    continue;
                }
            }

        /* Bitrate switching file */
        if (!strcmp(av[pos], "-swf") && pos + 1 < ac)
        {
#ifdef ENABLE_HR_MODE
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if(arg->hrmode == 2)
        {
            if(pos >= ac)
            {
                puts("Lossless mode: Using maximum bitrate for lossless coding of all frames");
                arg->bitrate = 0; /* dummy value for encoder init */
            }
            else
            {
                puts( "Lossless mode with fixed bitrate." );
                arg->bitrate = atoi( av[pos] );
                if ( arg->bitrate == 0 )
                {
                    arg->bitrate = 64000 * (1 + 3 * arg->hrmode); /* dummy value for encoder init */
                    arg->bitrate_file = av[pos];
                    printf( "Using bitrate switching file!\n" );
                }
            }
        }
        else
#endif
        {
            arg->bitrate = 64000 * (1 + 3 * arg->hrmode); /* dummy value for encoder init */
        }
#else
            arg->bitrate = 64000; /* dummy value for encoder init */
#endif
            arg->bitrate_file = av[++pos];
            puts("Using bitrate switching file!");
        }
        /* Error protection mode */
        if (!strcmp(av[pos], "-epmode") && pos + 1 < ac)
        {
            arg->epmode = atoi(av[++pos]);
            exit_if((unsigned)arg->epmode > 5, "EP mode must be in range [0-5]");
            if (arg->epmode == 0 && strcmp(av[pos], "0"))
            {
                arg->epmode      = 1;
                arg->epmode_file = av[pos];
                puts("Using epmode switching file!");
            }
            else
            {
                printf("Error protection %sabled (%i). ", arg->epmode ? "en" : "dis", arg->epmode);
            }
        }
        /* Error detection pattern */
        if (!strcmp(av[pos], "-edf") && pos + 1 < ac)
        {
            arg->edf = av[++pos];
            puts("Writing error detection file!");
        }

        /* error pattern */
        if (!strcmp(av[pos], "-ep_dbg") && pos + 1 < ac)
        {
            arg->channel_coder_vars_file = av[++pos];
            puts("Saving channel decoder debug information to files!");
        }
    }

    exit_if(arg->encoder_only && arg->decoder_only, "Encoder and decoder modes are exclusive!");
    exit_if(arg->ept && (!arg->epf || !arg->encoder_only), "Use -ept only with -E -epf FILE!");
#ifdef CR14_A_ADD_LOSSLESS_MODE
    exit_if(arg->rel_prio && (arg->decoder_only || arg->hrmode != 2), "-rel_prio requires encoder mode with -lossless!");
#endif
    exit_if(pos + 1 >= ac, MISSING_ARGUMENT_MESSAGE);

    arg->inputFilename  = av[pos++];
    arg->outputFilename = av[pos++];

    /* Bitrate */
    if (!arg->decoder_only)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if(arg->hrmode == 2)
        {
            if(pos >= ac)
            {
                puts("Lossless mode: Using maximum bitrate for lossless coding of all frames");
                arg->bitrate = 0; /* dummy value for encoder init */
            }
            else
            {
                puts( "Lossless mode with fixed bitrate." );
                arg->bitrate = atoi( av[pos] );
                if ( arg->bitrate == 0 )
                {
                    arg->bitrate = 64000 * (1 + 3 * arg->hrmode); /* dummy value for encoder init */
                    arg->bitrate_file = av[pos];
                    printf( "Using bitrate switching file!\n" );
                }
            }
        }
        else
#endif
        {
            exit_if(pos >= ac, MISSING_ARGUMENT_MESSAGE);
            arg->bitrate = atoi(av[pos]);
            if (arg->bitrate == 0)
            {
#ifdef ENABLE_HR_MODE
                arg->bitrate = 64000 * (1 + 3 * arg->hrmode); /* dummy value for encoder init */  
#else
                arg->bitrate      = 64000; /* dummy value */
#endif
                arg->bitrate_file = av[pos];
                puts("Using bitrate switching file!");
            }
        }
    }
    putchar('\n');
}

/* check condition and if it fails, exit with error message */
static void exit_if(int condition, const char *message)
{
    if (condition)
    {
        puts(message);
        if (condition < LC3PLUS_WARNING)
        {
            exit(1);
        }
    }
}

/* open file with .cfg suffix if file_cfg is null */
static FILE *fopen_cfg(const char *file, const char *file_cfg, const char *mode)
{
    return file_cfg ? fopen(file_cfg, mode) : fopen_with_ext(file, ".cfg", mode);
}

static FILE *open_bitstream_writer(const char *file, uint32_t samplerate, int bitrate, short channels,
                                   uint32_t signal_len, float frame_ms, int epmode, int32_t hrmode, int g192, const char *file_cfg
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    , int bitsPerSample
#endif
)
{
    FILE *f     = fopen(file, "wb");
    FILE *f_use = f;
    FILE *f_cfg = NULL;

    if (g192)
    {
        f_cfg = fopen_cfg(file, file_cfg, "wb");
        exit_if(f_cfg == NULL, "Error opening G192 configuration-file!");
        f_use = f_cfg;
    }

    if (f_use)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        uint16_t bps = 0;
        if (bitsPerSample == 16)
        {
            bps = 1;
        } else if (bitsPerSample == 24)
        {
            bps = 2;
        }
        
        uint16_t header[11] = 
#else
        uint16_t header[10] = 
#endif
        {0xcc1c,        sizeof(header), samplerate / 100,
                              bitrate / 100, channels,       
                              (uint16_t)(frame_ms * 100),
                              epmode > 0 ? 1 : 0,   signal_len,     signal_len >> 16, hrmode
#ifdef CR14_A_ADD_LOSSLESS_MODE
                              , bps
#endif
                              };
        fwrite(&header, sizeof(header), 1, f_use);
    }

    safe_fclose(f_cfg);
    return f;
}

static FILE *open_bitstream_reader(const char *file, unsigned int *samplerate, int *bitrate, short *channels,
                                   uint32_t *signal_len, float *frame_ms, int *epmode, int *hrmode, int g192,
                                   const char *file_cfg
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                    , int *bitsPerSample
#endif
                                   )
{
    FILE *f     = fopen(file, "rb");
    FILE *f_use = f;
    FILE *f_cfg = NULL;
    int32_t tmp_return_val;

    if (g192)
    {
        f_cfg = fopen_cfg(file, file_cfg, "rb");
        exit_if(f_cfg == NULL, "Error opening G192 configuration-file!");
        f_use = f_cfg;
    }

    if (f_use)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        uint16_t header[11] = {0};
#else
        uint16_t header[10] = {0};
#endif
        tmp_return_val = fread(header, sizeof(header), 1, f_use);
        
        if (header[0] != 0xcc1c)
        {
            /* Old style bitstream header */
            *samplerate = header[0] * 100;
            *bitrate    = header[1] * 100;
            *channels   = header[2];
            fseek(f_use, 6, SEEK_SET);
        }
        else
        {         
#ifdef CR14_A_ADD_LOSSLESS_MODE
            if (header[1] > 18 && header[9] > 1) 
            {
                assert(header[1] >= 20);
            }
#else
            assert(header[1] >= 18);
#endif
            *samplerate = header[2] * 100;
            *bitrate    = header[3] * 100;
            *channels   = header[4];
            *frame_ms   = (float)(header[5] / 100.0);
            *epmode     = header[6];
            *signal_len = (uint32_t)header[7] | ((uint32_t)header[8] << 16);
            *hrmode     = header[1] > 18 ? header[9] : 0;
#ifdef CR14_A_ADD_LOSSLESS_MODE
            *bitsPerSample = header[1] > 20 ? header[10] : 0;
#endif
            fseek(f_use, header[1], SEEK_SET);
        }
    }
    
#ifdef CR14_A_ADD_LOSSLESS_MODE
            if (*bitsPerSample == 1)
            {
                *bitsPerSample = 16;
            } else if (*bitsPerSample == 2) {
                *bitsPerSample = 24;
            } else {
                *bitsPerSample = 0;
            }
#endif

    (void) tmp_return_val;
    safe_fclose(f_cfg);
    return f;
}

static void write_bitstream_frame_G192(FILE *bitstream_file, uint8_t *bytes, int size)
{
    int      i           = 0;
    Word16   currentByte = 0;
    Word16   bit = 0, bitNumber = 0, syncWord = 0;
    uint16_t nbits = size * 8; /* G192 expects number of bits */

    /* Write good/bad frame info -> encoder writes only good frames */
    syncWord = G192_GOOD_FRAME;
    fwrite(&syncWord, sizeof(Word16), 1, bitstream_file);

    /* Write length info */
    fwrite(&nbits, sizeof(nbits), 1, bitstream_file);

    for (i = 0; i < size; i++)
    {
        currentByte = bytes[i];

        /* Start with LSB */
        for (bitNumber = 1; bitNumber < 9; bitNumber++)
        {
            bit = (currentByte & (1 << (bitNumber - 1))) != 0;
            bit = bit ? G192_ONE : G192_ZERO;
            fwrite(&bit, sizeof(bit), 1, bitstream_file);
        }
    }
}

static void write_bitstream_frame(FILE *bitstream_file, uint8_t *bytes,
                                  int size,
                                  int g192
#ifdef CR14_A_ADD_LOSSLESS_MODE
                                  , int *nBytes, int nChannels, int lossless
#endif
                                  )
{
    if (g192)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if (lossless && nChannels > 1)
        {
            /* one G192 frame per channel for lossless 2ch */
            int32_t offset = 0;
            for (int ch = 0; ch < nChannels; ch++)
            {
                write_bitstream_frame_G192(bitstream_file, bytes + offset, nBytes[ch]);
                offset += nBytes[ch];
            }
        }
        else
        {
            write_bitstream_frame_G192(bitstream_file, bytes, size);
        }
#else
        write_bitstream_frame_G192(bitstream_file, bytes, size);
#endif
    }
    else
    {
        int      i      = 0;
        uint16_t nbytes = size;

#ifdef CR14_A_ADD_LOSSLESS_MODE
        uint16_t total_numbytes = 0;
        for (int k = 0; k < nChannels; k++)
        {
            nbytes = nBytes[k];
            total_numbytes += nbytes;

            fwrite( &nbytes, sizeof( nbytes ), 1, bitstream_file );

            for ( ; i < total_numbytes; i++ )
            {
                putc( bytes[i], bitstream_file );
            }
        }
#else
        fwrite(&nbytes, sizeof(nbytes), 1, bitstream_file);
        for (i = 0; i < size; i++)
        {
            putc(bytes[i], bitstream_file);
        }
#endif
    }
}

static int read_bitstream_frame_G192(FILE *bitstream_file, 
                                     int size, 
                                     uint8_t *bytes, int *bfi_ext)
{
    int      i = 0, j = 0, read = 0;
    uint16_t nbits      = 0;
    int16_t  currentBit = 0, frameIndicator = 0, nbytes = 0;

    /* Read frame indicator info -> good/bad/redundancy frame */
    read = (int)fread(&frameIndicator, sizeof(frameIndicator), 1, bitstream_file);
    if (read != 1)
    {
        return -1;
    }

    /* Read length info */
    read = (int)fread(&nbits, sizeof(nbits), 1, bitstream_file);

    nbytes = nbits / 8;

    exit_if(frameIndicator != G192_GOOD_FRAME && frameIndicator != G192_BAD_FRAME &&
                frameIndicator != G192_REDUNDANCY_FRAME,
            "Wrong G192 format detected in bitstream file! The sync word could not be recognized!");

    for (i = 0; i < nbytes && i < size; i++)
    {
        int byte = 0;
        for (j = 0; j < 8; j++)
        {
            read = (int)fread(&currentBit, sizeof(currentBit), 1, bitstream_file);
            if (currentBit == G192_ONE)
            {
                byte |= 1UL << j;
            }
        }
        bytes[i] = (uint8_t)byte;
    }
    if (frameIndicator == G192_GOOD_FRAME)
    {
        *bfi_ext = 0;
    }
    else if (frameIndicator == G192_BAD_FRAME)
    {
        nbytes   = 0;
        *bfi_ext = 1;
    }
    else if (frameIndicator == G192_REDUNDANCY_FRAME)
    {
        *bfi_ext = 3;
    }

    return nbytes;
}

#ifdef CR14_A_ADD_LOSSLESS_MODE
static void 
#else
static int 
#endif
read_bitstream_frame(FILE *bitstream_file, uint8_t *bytes,
                     int size,
                     int g192, int *bfi_ext
#ifdef CR14_A_ADD_LOSSLESS_MODE
, int *nBytes, int nChannels, int lossless
#endif
)
{
    if (g192)
    {
#ifdef CR14_A_ADD_LOSSLESS_MODE
        if (lossless && nChannels > 1)
        {
            int32_t offset = 0;
            for (int ch = 0; ch < nChannels; ch++)
            {
                int32_t got = read_bitstream_frame_G192(bitstream_file, size - offset, bytes + offset, bfi_ext);
                if (got < 0)
                {
                    for (int j = ch; j < nChannels; j++) { nBytes[j] = -1; }
                    return;
                }
                nBytes[ch] = got;
                offset += got;
            }
        }
        else
        {
            nBytes[0] = read_bitstream_frame_G192(bitstream_file, size, bytes, bfi_ext);
        }
#else
        return read_bitstream_frame_G192(bitstream_file, size, bytes, bfi_ext);
#endif
    }
    else
    {
        int      i      = 0;
        uint16_t nbytes = 0;
      
#ifdef CR14_A_ADD_LOSSLESS_MODE
        uint16_t total_numbytes = 0;
        for (int k = 0; k < nChannels; k++)
        {
            if ( fread( &nBytes[k], sizeof( nbytes ), 1, bitstream_file ) != 1 )
            {
                nBytes[k] = -1;
                return; /* End of file reached */
            }
            
            total_numbytes += nBytes[k];

            for ( ; i < total_numbytes && i < size; i++ )
            {
                bytes[i] = (uint8_t) getc( bitstream_file );
            }

            if ( total_numbytes != i )
            {
                nBytes[k] = -1;
            }
        }
#else
        if (fread(&nbytes, sizeof(nbytes), 1, bitstream_file) != 1)
        {
            return -1; /* End of file reached */
        }
        for (i = 0; i < nbytes && i < size; i++)
        {
            bytes[i] = getc(bitstream_file);
        }
        return nbytes;
#endif
    }
}

/* read value from file and rewind if end is reached */
static int16_t loopy_read16(FILE *f)
{
#ifdef READ_G192_FER_BYTE
    int8_t tmp8 = -8;
#endif
    int16_t tmp = 0;
    int32_t tmp_return_val;
#ifdef READ_G192FER
    static int16_t format_start_check = -1;
#endif
#ifdef READ_G192_FER_BYTE
    if (format_start_check == -1) {
        /* first time always read an int16  (two bytes) */
        if (fread(&tmp, sizeof(tmp), 1, f) != 1) {
            printf("\n Warning !!!   loopy_read16 requires at least two initial bytes (one  int16) in FER file  \n ");
            fflush(stdout);
        }

        if (tmp == 0x2021 || tmp == 0x2020 || tmp == 0x2120 || tmp == 0x2121)
        {   /* G192  BYTE format determined on  the first two bytes    */
            format_start_check = 8;
        }
      
        fseek(f, 0, SEEK_SET); /* rewind */
    }

    /* restart loopy read with knowledge of G192_byte  ( or int16_t )  */
    if (format_start_check == 8)
    {
        if (fread(&tmp8, sizeof(tmp8), 1, f) != 1) {
            fseek(f, 0, SEEK_SET); /* rewind */
             fread(&tmp8, sizeof(tmp8), 1, f);
        }
    }
    else
    {
        if (fread(&tmp, sizeof(tmp), 1, f) != 1)
        {
            fseek(f, 0, SEEK_SET);
            tmp_return_val = fread(&tmp, sizeof(tmp), 1, f);
        }
    }
#else 
    if (fread(&tmp, sizeof(tmp), 1, f) != 1)
    {
        fseek(f, 0, SEEK_SET);
        tmp_return_val = fread(&tmp, sizeof(tmp), 1, f);
    }
#endif 

  

#ifdef READ_G192FER
    if (format_start_check < 0)
    {
        format_start_check = tmp; /*save first 16 bit  FER value  */
    }

    if (format_start_check >= 0 && format_start_check <= 1)
    {
        if (tmp != 0 && tmp != 1)
        {
            printf("\n Warning !!! assumed [0, 1] strange FER file values %d  %d \n ", format_start_check, tmp);
            fflush(stdout);
        }
    }

    if (format_start_check == G192_BAD_FRAME || format_start_check == G192_GOOD_FRAME)
    {


        if ((tmp != G192_BAD_FRAME && tmp != G192_GOOD_FRAME))
        {
            printf("\n Warning !!! assumed g.192 [0x6b21, 0x6b20,] , strange FER file values %d  %d \n ",
                   format_start_check, tmp);
            fflush(stdout);
        }
        else
        {
            tmp = (G192_GOOD_FRAME - tmp); /* convert g192 synch word to   1 and 0 , note PC byte order assumed */
        }
    }

#ifdef READ_G192_FER_BYTE
    if (format_start_check == 8)
    {   /* G192  BYTE format reading   */
        if ((tmp8 != 0x21 && tmp8 != 0x20))
        {
            printf("\n Warning !!! assumed g.192 byte  [0x21, 0x20,] , strange byte FER file values %d  %d \n ", format_start_check, tmp8);
            fflush(stdout);
        }
        tmp = (int16_t)(0x21 - tmp8); /*convert to   bfi (0 ==good, 1 = bad)*/
    }
#endif 

    ASSERT(tmp == 1 || tmp == 0);
#endif 
    
    (void) tmp_return_val;
    return tmp;
}

static int64_t loopy_read64(FILE *f)
{
    int64_t tmp = 0;
    int32_t tmp_return_val;

    if (fread(&tmp, sizeof(tmp), 1, f) != 1)
    {
        fseek(f, 0, SEEK_SET);
        tmp_return_val = fread(&tmp, sizeof(tmp), 1, f);
    }
    
    (void) tmp_return_val;
    return tmp;
}

static void scale_24_to_16(const int32_t *in, int16_t *out, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        out[i] = in[i];
    }
}


static void deinterleave(int32_t *in, int32_t **out, int n, int channels)
{
    int ch, i;
    for (ch = 0; ch < channels; ch++)
    {
        for (i = 0; i < n; i++)
        {
            out[ch][i] = in[i * channels + ch];
        }
    }
}

static void interleave_short(int16_t** in, int16_t* out, int32_t n, int32_t channels)
{
    int32_t ch, i;
    for (ch = 0; ch < channels; ch++) {
        for (i = 0; i < n; i++) {
            out[i * channels + ch] = in[ch][i];
        }
    }
}

static void interleave_int(int32_t** in, int32_t* out, int32_t n, int32_t channels)
{
    int32_t ch, i;
    for (ch = 0; ch < channels; ch++) {
        for (i = 0; i < n; i++) {
            out[i * channels + ch] = in[ch][i];
        }
    }
}
