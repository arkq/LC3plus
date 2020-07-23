/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "functions.h"
#include "lc3.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PACK 1
#define UNPACK 2

/* struct to hold command line arguments */
typedef struct
{
    char *inputFilename;
    char *outputFilename;
    int   bitrate;
    char *bitrate_file;
    int   encoder_only;
    int   decoder_only;
    int   formatG192;
    char *configFilenameG192;
    float frame_ms;
    int   plcMeth;
    char *epf;
    int   epmode;
    int   hrmode;
    int   mode; /* unpack / pack*/
    int   gross_bytes;
} Arguments;

/* local functions */
static void  parseCmdl(int ac, char **av, Arguments *arg);
static FILE *open_bitstream_reader(const char *file, uint32_t *samplerate, int *bitrate, short *channels,
                                   uint32_t *signal_len, float *frame_ms, int *epmode, int *hrmode, int g192,
                                   const char *file_cfg);
static FILE *open_bitstream_writer(const char *file, uint32_t samplerate, int bitrate, short channels,
                                   uint32_t signal_len, float frame_ms, int epmode, int g192, const char *file_cfg);
static void  write_bitstream_frame(FILE *bitstream_file, uint8_t *bytes, int size, int g192, Word16 syncWord);
static int   read_bitstream_frame(FILE *bitstream_file, uint8_t *bytes, int size, int g192, int *bfi_ext);
static FILE *fopen_with_ext(const char *file, const char *ext, const char *mode);
static void  cleanup(void);
static void  exit_if(int condition, const char *message);
LC3_Error    channel_coder_pack(LC3_Dec *decoder, UWord8 *input_bytes, int num_bytes_in, void *scratch, int bfi_ext,
                                int *lc3_num_bytes, int gross_bytes, int epmode);
LC3_Error    channel_coder_unpack(LC3_Dec *decoder, UWord8 *input_bytes, int num_bytes_in, void *scratch, int bfi_ext,
                                  int *lc3_num_bytes);
void processReorderBitstream_dec_fx(UWord8 *bytes, Word16 n_pccw, Word16 n_pc, Word16 b_left, Word8 *scratchBuffer);

static FILE *output_bitstream;
static FILE *input_bitstream;

#include "license.h" /* provides LICENSE string */

static const char *const USAGE_MESSAGE =
    /* Lines must not be longer than this! --------------------------------------->| */
    "Channel Coder Converter\n\n"
    "ccConvert unpacks protected lc3 bitstream to unprotected lc3 bitstream and vice versa.\n"
    "\n"
    "Usage:\n"
    "   pack mode:   ./ccConvert -pack gross_bytes ep_mode in.lc3 out.lc3\n"
    "   unpack mode: ./ccConvert -unpack in.lc3 out.lc3\n"
    "\n"
    "   example:     ./ccConvert -pack 80 3 in_ep0.lc3 out_ep3.lc3\n"
    "   For G.192 format use filename extension \".g192\"\n"
    "   example:     ./ccConvert -unpack in_ep3.g192 out_ep0.g192\n"
    "\n"
    "\nChannel coder options:\n"
    "  epmode must be one of the following:\n"
    "    0: Error protection disabled\n"
    "    1: Minimum error protection, detection only\n"
    "    2: Moderate error protection\n"
    "    3: Strong error protection\n"
    "    4: Maximum error protection\n";

static const char *ERROR_MESSAGE[18] = {
    "",                                                /* LC3_OK                  */
    "Function call failed!",                           /* LC3_ERROR               */
    "Frame failed to decode and was concealed!",       /* LC3_DECODE_ERROR        */
    "Pointer argument is null!",                       /* LC3_NULL_ERROR          */
    "Invalid sampling rate!",                          /* LC3_SAMPLERATE_ERROR    */
    "Invalid number of channels!",                     /* LC3_CHANNELS_ERROR      */
    "Invalid bitrate!",                                /* LC3_BITRATE_ERROR       */
    "Invalid number of bytes!",                        /* LC3_NUMBYTES_ERROR      */
    "Invalid PLC method!",                             /* LC3_PLCMODE_ERROR       */
    "Invalid EP mode!",                                /* LC3_EPCLASS_ERROR       */
    "Invalid frame ms value!",                         /* LC3_FRAMEMS_ERROR       */
    "Unaligned pointer!",                              /* LC3_ALIGN_ERROR         */
    "Invalid channel mode request!",                   /* LC3_CMR_ERROR           */
    "Bitrate has not been set!",                       /* LC3_BITRATE_UNSET_ERROR */
    "Function can't be called after bitrate was set!", /* LC3_BITRATE_SET_ERROR   */
    "Invalid external bad frame index!",               /* LC3_BFI_EXT_ERROR       */
    "Generic Warning",                                 /* LC3_WARNING             */
    "Invalid bandwidth frequency!"                     /* LC3_BW_WARNING          */
};

int main(int ac, char **av)
{
    Arguments arg;
    uint32_t  nSamplesFile = 0xffffffff, sampleRate = 0;
    short     nChannels    = 0;
    int       num_bytes_in = 0, nSamples = 0;
    int       encoder_size = 0, decoder_size = 0, scratch_size = 0;
    int       bfi_ext = 0;
    LC3_Dec * decoder = NULL;
    void *    scratch = NULL;
    LC3_Error err     = LC3_OK;
    uint8_t   bytes[LC3_MAX_BYTES];
    int       lc3_num_bytes = 0;
    Word16    syncWord      = G192_GOOD_FRAME;
    int       bitrate_in, epmode_in = 0, epmode_out = 0;
    
    UNUSED(encoder_size);

    /* Parse Command-line */
    parseCmdl(ac, av, &arg);

#ifdef STAMEM_COUNT
    Sta_Mem_Init();
#endif
#ifdef DYNMEM_COUNT
    Dyn_Mem_Init();
#endif

    /* exit handler to clean up resources */
    atexit(cleanup);

    /* Open Input Bitstream File */
    input_bitstream =
        open_bitstream_reader(arg.inputFilename, &sampleRate, &bitrate_in, &nChannels, &nSamplesFile, &arg.frame_ms,
                              &epmode_in, &arg.hrmode, arg.formatG192, arg.configFilenameG192);
    exit_if(!input_bitstream, "Error opening bitstream file!");
    exit_if(arg.hrmode, "HR bitstreams not supported!");

    /* Setup Decoder */
    decoder_size = lc3_dec_get_size(sampleRate, nChannels, (LC3_PlcMode)arg.plcMeth);
    decoder      = malloc(decoder_size);
    err          = lc3_dec_init(decoder, sampleRate, nChannels, (LC3_PlcMode)arg.plcMeth);
    exit_if(err, ERROR_MESSAGE[err]);

    err = lc3_dec_set_frame_dms(decoder, (int)(arg.frame_ms * 10));
    exit_if(err, ERROR_MESSAGE[err]);

    err = lc3_dec_set_ep_enabled(decoder, arg.epmode != 0);
    exit_if(err, ERROR_MESSAGE[err]);

    nSamples = lc3_dec_get_output_samples(decoder);

    /* Open Output Bitstream File */
    if (arg.mode == PACK)
    {
        printf("Converting unprotected LC3 bitstream %s to protected LC3 bitstream %s...\n", arg.inputFilename,
               arg.outputFilename);
        exit_if(epmode_in != 0, USAGE_MESSAGE);
        epmode_out  = arg.epmode;
    }
    if (arg.mode == UNPACK)
    {
        printf("Converting protected LC3 bitstream %s to unprotected LC3 bitstream %s...\n", arg.inputFilename,
               arg.outputFilename);
        exit_if(epmode_in == 0, USAGE_MESSAGE);
        epmode_out  = 0;
    }
    output_bitstream = open_bitstream_writer(arg.outputFilename, sampleRate, bitrate_in, nChannels, nSamplesFile,
                                             arg.frame_ms, epmode_out, arg.formatG192, arg.configFilenameG192);
    exit_if(!output_bitstream, "Error creating bitstream file!");

    scratch_size = lc3_dec_get_scratch_size(decoder);
    scratch      = malloc(scratch_size);
    exit_if(!scratch, "Failed to allocate scratch memory!");

#ifdef STAMEM_COUNT
    Sta_Mem_Add("Encoder", encoder_size);
    Sta_Mem_Add("Decoder", decoder_size);
#endif

    setFrameRate(sampleRate, nSamples);
    Init_WMOPS_counter();

    /* Encoder + Decoder loop */
    while (1)
    {
        syncWord = G192_GOOD_FRAME;
        /* Read bitstream */
        num_bytes_in = read_bitstream_frame(input_bitstream, bytes, sizeof(bytes), arg.formatG192, &bfi_ext);
        if (num_bytes_in < 0)
        {
            break;
        }

        if (bfi_ext == 0)
        {
            lc3_num_bytes = 0;

            if (arg.mode == PACK)
            {
                channel_coder_pack(decoder, bytes, num_bytes_in, scratch, bfi_ext, &lc3_num_bytes, arg.gross_bytes,
                                   arg.epmode);
            }
            if (arg.mode == UNPACK)
            {
                channel_coder_unpack(decoder, bytes, num_bytes_in, scratch, bfi_ext, &lc3_num_bytes);
            }

            exit_if(err && err != LC3_DECODE_ERROR, ERROR_MESSAGE[err]);

            /* write output bitstream */

            if (arg.formatG192)
            {
                if (err == LC3_DECODE_ERROR || bfi_ext != 0)
                {
                    syncWord      = G192_BAD_FRAME;
                    lc3_num_bytes = 0;
                }
            }
        }
        else /* if bfi_ext == 0 */
        {
            syncWord      = G192_BAD_FRAME;
            lc3_num_bytes = 0;
        }

        write_bitstream_frame(output_bitstream, bytes, lc3_num_bytes, arg.formatG192, syncWord);

        BASOP_frame_update();
    }
    puts("\nProcessing done!");

    free(decoder);
    free(scratch);

#if WMOPS
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
    safe_fclose(output_bitstream);
    safe_fclose(input_bitstream);
}

static void parseCmdl(int ac, char **av, Arguments *arg)
{
    exit_if(ac <= 1, USAGE_MESSAGE);
    exit_if(!strcmp(av[1], "-h"), USAGE_MESSAGE);
    int pos = 1;
    memset(arg, 0, sizeof(*arg));
    arg->frame_ms     = 10;
    arg->decoder_only = 1;
    arg->bitrate      = 0;
    if (!strcmp(av[pos], "-unpack"))
    {
        exit_if(ac != 4, USAGE_MESSAGE);
        arg->mode = UNPACK;
        pos++;
    }
    else if (!strcmp(av[pos], "-pack"))
    {
        exit_if(ac != 6, USAGE_MESSAGE);
        arg->mode = PACK;
        pos++;
        arg->gross_bytes = atoi(av[pos++]);
        arg->epmode      = atoi(av[pos++]);
        exit_if(arg->epmode == 0 || arg->epmode > 4, USAGE_MESSAGE);
    }
    else
    {
        arg->mode = -1;
    }
    exit_if(arg->mode == -1, USAGE_MESSAGE);
    arg->inputFilename  = av[pos++];
    arg->outputFilename = av[pos++];

    int in_len = strlen(arg->inputFilename);
    if (!strcmp(&arg->inputFilename[in_len - 4], "g192"))
    {
        arg->formatG192 = 1;
    }
    else
    {
        arg->formatG192 = 0;
    }
    putchar('\n');
}

/* check condition and if it fails, exit with error message */
static void exit_if(int condition, const char *message)
{
    if (condition)
    {
        puts(message);
        if (condition < LC3_WARNING)
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
                                   uint32_t signal_len, float frame_ms, int epmode, int g192, const char *file_cfg)
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
        uint16_t header[9] = {0xcc1c,        sizeof(header), samplerate / 100,
                              bitrate / 100, channels,       (uint16_t)(frame_ms * 100),
                              epmode,        signal_len,     signal_len >> 16};
        fwrite(&header, sizeof(header), 1, f_use);
    }

    safe_fclose(f_cfg);
    return f;
}

static FILE *open_bitstream_reader(const char *file, unsigned int *samplerate, int *bitrate, short *channels,
                                   uint32_t *signal_len, float *frame_ms, int *epmode, int *hrmode, int g192,
                                   const char *file_cfg)
{
    FILE *f     = fopen(file, "rb");
    FILE *f_use = f;
    FILE *f_cfg = NULL;

    if (g192)
    {
        f_cfg = fopen_cfg(file, file_cfg, "rb");
        exit_if(f_cfg == NULL, "Error opening G192 configuration-file!");
        f_use = f_cfg;
    }

    if (f_use)
    {
        uint16_t header[10] = {0};
        fread(header, sizeof(header), 1, f_use);
        {
            assert(header[1] >= 18);
            *samplerate = header[2] * 100;
            *bitrate    = header[3] * 100;
            *channels   = header[4];
            *frame_ms   = (float)(header[5] / 100.0);
            *epmode     = header[6];
            *signal_len = (uint32_t)header[7] | ((uint32_t)header[8] << 16);
            *hrmode     = header[1] > 18 ? header[9] : 0;
            fseek(f_use, header[1], SEEK_SET);
        }
    }

    safe_fclose(f_cfg);
    return f;
}

static void write_bitstream_frame_G192(FILE *bitstream_file, uint8_t *bytes, int size, Word16 syncWord)
{
    int      i           = 0;
    Word16   currentByte = 0;
    Word16   bit = 0, bitNumber = 0;
    uint16_t nbits = size * 8; /* G192 expects number of bits */

    /* Write good/bad frame info */
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

static void write_bitstream_frame(FILE *bitstream_file, uint8_t *bytes, int size, int g192, Word16 syncWord)
{
    if (g192)
    {
        write_bitstream_frame_G192(bitstream_file, bytes, size, syncWord);
    }
    else
    {
        int      i      = 0;
        uint16_t nbytes = size;
        fwrite(&nbytes, sizeof(nbytes), 1, bitstream_file);
        for (i = 0; i < size; i++)
        {
            putc(bytes[i], bitstream_file);
        }
    }
}

static int read_bitstream_frame_G192(FILE *bitstream_file, int size, uint8_t *bytes, int *bfi_ext)
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

static int read_bitstream_frame(FILE *bitstream_file, uint8_t *bytes, int size, int g192, int *bfi_ext)
{
    if (g192)
    {
        return read_bitstream_frame_G192(bitstream_file, size, bytes, bfi_ext);
    }
    else
    {
        int      i      = 0;
        uint16_t nbytes = 0;
        if (fread(&nbytes, sizeof(nbytes), 1, bitstream_file) != 1)
        {
            return -1; /* End of file reached */
        }
        for (i = 0; i < nbytes && i < size; i++)
        {
            bytes[i] = getc(bitstream_file);
        }
        return nbytes;
    }
}

LC3_Error channel_coder_pack(LC3_Dec *decoder, UWord8 *bytes, int num_bytes_in, void *scratch, int bfi_ext,
                             int *lc3_num_bytes, int gross_bytes, int epmode)
{
    int       ch = 0, bfi = bfi_ext;
    LC3_Error err         = LC3_OK;
    int       channel_bfi = 0;

    BASOP_sub_start("Pack");

    Word16 b_left;
    Word16 gain;
    Word16 gg_idx, fac_ns_idx;
    Word16 bp_side, mask_side;
    Word16 tns_numfilters, lsbMode, lastnz;
    Word16 ltpf_idx[3];
    Word16 nbbits = 0;

    /* Buffers */
    Word16 *int_scf_fx_exp, tns_order[TNS_NUMFILTERS_MAX];
    UWord8 *resBitBuf;
    Word16 *sqQdec;
    Word16 *int_scf_fx, /* *x_fx,*/ *indexes, *scf_q;
    Word32 *L_scf_idx;
    Word32 *q_d_fx;
    Word8 * currentScratch;
    // DecSetup *h_DecSetup = decoder->channel_setup[channel];

    /* BUFFER INITIALISATION. Some buffers may overlap since they are not used in the whole decoding process */
    q_d_fx = scratchAlign(scratch, 0); /* Size = 4 * MAX_LEN bytes */
    resBitBuf =
        scratchAlign(q_d_fx, sizeof(*q_d_fx) * decoder->frame_length); /* Size = 2 * NPRM_RESQ = 2 * MAX_LEN bytes */
    indexes = scratchAlign(
        resBitBuf, sizeof(*resBitBuf) * decoder->frame_length); /* Size = 2 * TNS_NUMFILTERS_MAX * MAXLAG = 32 bytes */
    L_scf_idx      = scratchAlign(indexes, sizeof(*indexes) * TNS_NUMFILTERS_MAX *
                                          MAXLAG); /* Size = 4 * SCF_MAX_PARAM = 28 bytes -> aligned to 32 bytes */
    sqQdec         = scratchAlign(L_scf_idx, sizeof(*L_scf_idx) * (SCF_MAX_PARAM));   /* Size = 2 * MAX_LEN bytes */
    scf_q          = scratchAlign(sqQdec, sizeof(*sqQdec) * (decoder->frame_length)); /* Size = 2 * M = 32 bytes */
    int_scf_fx_exp = scratchAlign(scf_q, sizeof(*scf_q) * M); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    int_scf_fx     = scratchAlign(int_scf_fx_exp,
                              sizeof(*int_scf_fx_exp) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    currentScratch = scratchAlign(int_scf_fx, sizeof(*int_scf_fx) * MAX_BANDS_NUMBER); /* Size = 4 * MAX_LEN */

    if (bfi == 0)
    {
        bfi = !num_bytes_in;
    }

    // add padding enc -> from where ?
    for (ch = 0; ch < decoder->channels; ch++)
    {
        DecSetup *h_DecSetup = decoder->channel_setup[ch];
        // setup ep parameters:
        h_DecSetup->targetBytes = fec_get_data_size(epmode, 0, gross_bytes); // combined channel_coding = 0
        assert(h_DecSetup->targetBytes == num_bytes_in);
        decoder->n_pccw = fec_get_n_pccw(gross_bytes, epmode, 0);
        decoder->n_pc   = fec_get_n_pc(epmode, decoder->n_pccw, gross_bytes);

        int data_bytes = h_DecSetup->targetBytes;

        err = update_dec_bitrate(decoder, ch, data_bytes);
        if (err)
            return err;

        decoder->channel_setup[ch]->last_size = data_bytes;

        Word16 ch_bfi = channel_bfi;
        nbbits        = shl_pos(h_DecSetup->targetBytes, 3);
        processDecoderEntropy_fx(bytes, &bp_side, &mask_side, nbbits, decoder->yLen, decoder->fs_idx,
                                 decoder->BW_cutoff_bits, &tns_numfilters, &lsbMode, &lastnz, &ch_bfi, tns_order,
                                 &fac_ns_idx, &gg_idx, &gg_idx, ltpf_idx, L_scf_idx, decoder->frame_dms);
        channel_bfi = ch_bfi;

        processAriDecoder_fx(bytes, &bp_side, &mask_side, nbbits, decoder->yLen, decoder->fs_idx,
                             h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                             fac_ns_idx, gg_idx, decoder->frame_dms, decoder->n_pc, 0, shr_pos(nbbits, 3), 1, &gain,
                             &b_left, &gain, sqQdec, &gain, resBitBuf, indexes, &gain, currentScratch
        );

        IF (b_left > 0)
        {
            processReorderBitstream_fx(bytes, decoder->n_pccw, decoder->n_pc, b_left, currentScratch);
        }
        // end reordering

        if (epmode)
        { // target bytes = data bytes , gross_bytes = slot bytes
            fec_encoder(epmode, decoder->epmr, bytes, data_bytes, gross_bytes, decoder->n_pccw, scratch);

            *lc3_num_bytes += gross_bytes;
            // output_size += gross_bytes;
        }
        else
        {
            *lc3_num_bytes += h_DecSetup->targetBytes;
            // output_size += h_DecSetup->targetBytes;
        }
    }

    BASOP_sub_end();

    return bfi == 1 ? LC3_DECODE_ERROR : LC3_OK;
}

LC3_Error channel_coder_unpack(LC3_Dec *decoder, UWord8 *input_bytes, int num_bytes_in, void *scratch, int bfi_ext,
                               int *lc3_num_bytes)
{
    int       ch = 0, bfi = bfi_ext;
    LC3_Error err = LC3_OK;
    int       fec_num_bytes;
    int       channel_bfi, out_bfi;
    Word16    channel_epmr;

    BASOP_sub_start("Unpack");

    Word16 b_left;
    Word16 gain;
    Word16 gg_idx, fac_ns_idx;
    Word16 bp_side, mask_side;
    Word16 tns_numfilters, lsbMode, lastnz;
    Word16 ltpf_idx[3];
    Word16 nbbits;

    /* Buffers */
    Word16 *int_scf_fx_exp, tns_order[TNS_NUMFILTERS_MAX];
    UWord8 *resBitBuf;
    Word16 *sqQdec;
    Word16 *int_scf_fx, *indexes, *scf_q;
    Word32 *L_scf_idx;
    Word32 *q_d_fx;
    Word8 * currentScratch;

    /* BUFFER INITIALISATION. Some buffers may overlap since they are not used in the whole decoding process */
    q_d_fx = scratchAlign(scratch, 0); /* Size = 4 * MAX_LEN bytes */
    resBitBuf =
        scratchAlign(q_d_fx, sizeof(*q_d_fx) * decoder->frame_length); /* Size = 2 * NPRM_RESQ = 2 * MAX_LEN bytes */
    indexes = scratchAlign(
        resBitBuf, sizeof(*resBitBuf) * decoder->frame_length); /* Size = 2 * TNS_NUMFILTERS_MAX * MAXLAG = 32 bytes */
    L_scf_idx      = scratchAlign(indexes, sizeof(*indexes) * TNS_NUMFILTERS_MAX *
                                          MAXLAG); /* Size = 4 * SCF_MAX_PARAM = 28 bytes -> aligned to 32 bytes */
    sqQdec         = scratchAlign(L_scf_idx, sizeof(*L_scf_idx) * (SCF_MAX_PARAM));   /* Size = 2 * MAX_LEN bytes */
    scf_q          = scratchAlign(sqQdec, sizeof(*sqQdec) * (decoder->frame_length)); /* Size = 2 * M = 32 bytes */
    int_scf_fx_exp = scratchAlign(scf_q, sizeof(*scf_q) * M); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    int_scf_fx     = scratchAlign(int_scf_fx_exp,
                              sizeof(*int_scf_fx_exp) * MAX_BANDS_NUMBER); /* Size = 2 * MAX_BANDS_NUMBER = 128 bytes */
    currentScratch = scratchAlign(int_scf_fx, sizeof(*int_scf_fx) * MAX_BANDS_NUMBER); /* Size = 4 * MAX_LEN */

    if (bfi == 0)
    {
        bfi = !num_bytes_in;
    }

    decoder->epmr = 12;
    out_bfi       = 0;

    for (ch = 0; ch < decoder->channels; ch++)
    {
        fec_num_bytes = num_bytes_in / decoder->channels + (ch < (num_bytes_in % decoder->channels));

        channel_bfi = bfi;
        decoder->error_report =
            fec_decoder(input_bytes, fec_num_bytes, lc3_num_bytes, &channel_epmr, decoder->combined_channel_coding,
                        &decoder->n_pccw, &channel_bfi, &decoder->be_bp_left, &decoder->be_bp_right, &decoder->n_pc,
                        &decoder->m_fec, scratch);

        decoder->epmr = MIN(decoder->epmr, channel_epmr);

        if (channel_bfi != 1 && *lc3_num_bytes != decoder->channel_setup[ch]->last_size)
        {
            err = update_dec_bitrate(decoder, ch, *lc3_num_bytes);
            if (err)
                return err;

            decoder->channel_setup[ch]->last_size = *lc3_num_bytes;
        }

        channel_bfi = channel_bfi > 0;

        //#ifdef ENABLE_PADDING
        //        if (channel_bfi != 1)
        //        {
        //            Word16 padding_len, np_zero;
        //
        //            if (paddingDec_fx(input_bytes, shl(*lc3_num_bytes, 3), decoder->yLen, decoder->BW_cutoff_bits,
        //            decoder->ep_enabled, &padding_len, &np_zero))
        //            {
        //                channel_bfi = 1;
        //            }
        //
        //            input_bytes         = input_bytes + np_zero;
        //            decoder->n_pc = s_max(decoder->n_pc - (2 * np_zero), 0);
        //
        //            if (channel_bfi == 2)
        //            {
        //                if (decoder->be_bp_right < (8*np_zero))
        //                {
        //                    channel_bfi = 0;
        //                    decoder->be_bp_left = -1;
        //                    decoder->be_bp_right = -1;
        //                }
        //                else
        //                {
        //                    decoder->be_bp_right = decoder->be_bp_right - (8 * np_zero);
        //                    decoder->be_bp_left  = s_max(decoder->be_bp_left - (8 * np_zero), 0);
        //                }
        //            }
        //
        //            *lc3_num_bytes = *lc3_num_bytes - padding_len;
        //        }
        //#endif

        // reordering for channelCodecConverter
        DecSetup *h_DecSetup = decoder->channel_setup[ch];

        IF (sub(channel_bfi, 1) != 0)
        {
            Word16 ch_bfi = channel_bfi;
            nbbits        = shl_pos(*lc3_num_bytes, 3);
            processDecoderEntropy_fx(input_bytes, &bp_side, &mask_side, nbbits, decoder->yLen, decoder->fs_idx,
                                     decoder->BW_cutoff_bits, &tns_numfilters, &lsbMode, &lastnz, &ch_bfi, tns_order,
                                     &fac_ns_idx, &gg_idx, &gg_idx, ltpf_idx, L_scf_idx, decoder->frame_dms);
            // BW_cutoff_idx_nf = BW_cutoff_idx;  move16();
            channel_bfi = ch_bfi;
        }

        IF (decoder->combined_channel_coding == 0 && decoder->n_pc > 0)
        {
            processAriDecoder_fx(input_bytes, &bp_side, &mask_side, nbbits, decoder->yLen, decoder->fs_idx,
                                 h_DecSetup->enable_lpc_weighting, tns_numfilters, lsbMode, lastnz, &gain, tns_order,
                                 fac_ns_idx, gg_idx, decoder->frame_dms, decoder->n_pc, 0, shr_pos(nbbits, 3), 2, &gain,
                                 &b_left, &gain, sqQdec, &gain, resBitBuf, indexes, &gain, currentScratch
            );

            IF (b_left > 0)
            {
                processReorderBitstream_dec_fx(input_bytes, decoder->n_pccw, decoder->n_pc,
                                               b_left - shr_sat(add(decoder->n_pc, 1), 1), currentScratch);
            }
        }
        // end reordering

        out_bfi |= channel_bfi;
        input_bytes += fec_num_bytes;
    }
    bfi = out_bfi & 1;

    BASOP_sub_end();

    return bfi == 1 ? LC3_DECODE_ERROR : LC3_OK;
}

void processReorderBitstream_dec_fx(UWord8 *bytes, Word16 n_pccw, Word16 n_pc, Word16 b_left, Word8 *scratchBuffer)
{
    Word16  block_bytes;
    UWord8 *bytes_tmp;

    bytes_tmp = (UWord8 *)scratchAlign(scratchBuffer, 0); /* Size = LC3_MAX_BYTES */

    if (n_pccw == 0)
    {
        return;
    }

    assert(b_left >= 0);

    /* set block size in bits and full bytes */
    block_bytes = shr_sat(add(n_pc, 1), 1);

    /* de-rearrange bitstream */
    basop_memmove(&bytes_tmp[0], &bytes[block_bytes], b_left * sizeof(UWord8));
    basop_memmove(&bytes_tmp[b_left], &bytes[0], block_bytes * sizeof(UWord8));
    basop_memmove(&bytes[0], &bytes_tmp[0], add(block_bytes, b_left) * sizeof(UWord8));
}
