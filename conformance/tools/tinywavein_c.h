/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __TINYWAVEIN_C_H__
#define __TINYWAVEIN_C_H__

/*#define TWI_SUPPORT_BWF*/
/* #define PRINT_HDR */     /* debug functionality to print header info */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***** Interface *********************************************************/

#ifndef TWI_UINT64
  #if ( __STDC_VERSION__ >= 199901L || ( defined(_MSC_VER) &&_MSC_VER >= 1600 ) )
    #include <stdint.h>
    #define TWI_UINT64 uint64_t
  #elif ( defined(_MSC_VER) &&_MSC_VER < 1600 )
    typedef unsigned __int64 TWI_UINT64;
  #elif ( defined(__arm__) || defined(__xtensa__) )
    #define TWI_UINT64 uint64_t
  #else
    #error "C99 or later required for 64-bit unsigned integer!"
  #endif
#endif

typedef struct WAVEFILEIN WAVEFILEIN;
#ifdef TWI_SUPPORT_BWF
typedef struct WAVEIN_LOUDNESSINFO WAVEIN_LOUDNESSINFO;
#endif

#define __TWI_SUCCESS  (0)
#define __TWI_ERROR    (-1)


/*!
 *  \brief Read header from a WAVEfile. Host endianess is handled accordingly.
 *  \return WAVEFILEIN handle on success and NULL on failure.
 *
 *  This open can cope with RF64 files, it returns samplesInFile == UINT_MAX when fileSize > 4 GB
 */
static WAVEFILEIN* OpenWav(
                           const char*     filename,
                           unsigned int*   samplerate,
                           short*          channels,
                           unsigned int*   samplesInFile,
                           short*          bps
                           );

/* read normalized 16-bit values in the range +32767..-32768 */
static int ReadWavShort(
                        WAVEFILEIN*   self,
                        short         sampleBuffer[],
                        unsigned int  nSamplesToRead,
                        unsigned int* nSamplesRead
                        );

/* read normalized 24-bit values in the range +8388607..-8388608 */
static int ReadWavInt(
                      WAVEFILEIN*   self,
                      int           sampleBuffer[],
                      unsigned int  nSamplesToRead,
                      unsigned int* nSamplesRead
                      );

/* read normalized single-precision values in the range +-1.0f */
static int ReadWavFloat(
                        WAVEFILEIN*   self,
                        float         sampleBuffer[],
                        unsigned int  nSamplesToRead,
                        unsigned int* nSamplesRead
                        );

#ifdef TWI_SUPPORT_BWF
/* read loudness information from BWF file */
static void ReadBWF(
                    WAVEFILEIN*           self,
                    WAVEIN_LOUDNESSINFO** wavInLoudness
                   );
#endif

/* close reader - always call when done reading to free resources */
static int CloseWavIn(WAVEFILEIN* self);

/* reset read pointer to first audio sample */
static int ResetWavIn(WAVEFILEIN* self);

/***** Implementation *********************************************************/

#if defined(__i386__) || defined(_M_IX86) || defined(__x86_64__) ||            \
    defined(_M_X64)   || defined(__arm__) || defined(__xtensa__) || defined(__aarch64__) || defined(__EMSCRIPTEN__) || defined(__hexagon__)
#define __TWI_LE    /* _T_iny _W_ave _I_n _L_ittle _E_ndian */
#endif

#if defined(__POWERPC__)
#define __TWI_BE    /* _T_iny _W_ave _I_n _B_ig _E_ndian */
#endif

#if !defined(__TWI_LE) && !defined(__TWI_BE)
#error unknown processor
#endif

/*--- local types/structs ----------------------------------*/

#ifdef TWI_SUPPORT_BWF
struct WAVEIN_LOUDNESSINFO {
  float loudnessVal;
  float loudnessRange;
  float maxTruePeakLevel;
  float maxMomentaryLoudnes;
  float maxShortTermLoudness;
};
#endif

struct WAVEFILEIN {
  FILE*             theFile;
  fpos_t            dataChunkPos;
  TWI_UINT64        position;      /* in pcm samples */
  TWI_UINT64        length;        /* in pcm samples */
  unsigned int      bps;
#ifdef TWI_SUPPORT_BWF
  WAVEIN_LOUDNESSINFO* loudnessInfo;
#endif
};

typedef struct
{
  short         formatType;     /* WAVE_FORMAT_PCM = 0x0001, etc. */
  short         channelCount;   /* 1 = mono, 2 = stereo, etc. */
  unsigned int  sampleRate;     /* 32000, 44100, 48000, etc. */
  unsigned int  bytesPerSecond; /* only important for compressed formats */
  short         blockAlignment; /* container size (in bytes) of one set of samples */
  short         bitsPerSample;  /* valid bits per sample 16, 20 or 24 */
  /* short extraFormatBytes ; */
} SWavInfo;

#ifdef TWI_SUPPORT_BWF
typedef struct {
  unsigned char description[256];
  unsigned char originator[32];
  unsigned char originatorReference[32];
  unsigned char originatorDate[10]; /* ASCII: <<yyyy:mm:dd>> */
  unsigned char originationTime[8]; /* ASCII: <<hh:mm:ss>> */
  unsigned int timeReferenceLow;
  unsigned int timeReferenceHigh;
  unsigned short version;
  unsigned char UMID[64]; /* Binary Bytes of SMPTE UMID */

  signed short loudnessVal;
  signed short loudnessRange;
  signed short maxTruePeakLevel;
  signed short maxMomentaryLoudnes;
  signed short maxShortTermLoudness;

  unsigned char Reserved[180];

  unsigned char codingHistory; /* ASCII: <<History coding>> */
} SBwfWav;
#endif

typedef struct {
  char chunkID[4];
  unsigned int chunkSize;
  /* long dataOffset ; */ /* never used */
} SChunk;

/* local wrapper, always returns correct endian */
static size_t fread_LE(void *ptr, size_t size, size_t nmemb, FILE *stream);
#ifdef TWI_SUPPORT_BWF
static int __ReadBextChunk( WAVEFILEIN* self, unsigned int* chunkSize );
#endif
static int __ReadDs64Chunk( WAVEFILEIN* self, unsigned int* chunkSize );

#ifdef __TWI_BE
static short BigEndian16(short v);
static int BigEndian32(int v);
/* static TWI_UINT64 BigEndian64(TWI_UINT64); */
#endif

/*!
 *  \brief Read header from a WAVEfile. Host endianess is handled accordingly.
 *  \return handle on success and NULL on failure.
 *
 *  This open can cope with RF64 files, it returns samplesInFile == UINT_MAX when fileSize > 4 GB
 */
static WAVEFILEIN* OpenWav(
                           const char*     filename,
                           unsigned int*   samplerate,
                           short*          channels,
                           unsigned int*   samplesInFile,
                           short*          bps
                           )
{
  WAVEFILEIN* self = NULL;
  size_t tmp_return_val;
  SChunk fmt_chunk, data_chunk;
  char tmpchunkID[4];
  int offset;
  unsigned int dmy;
  char tmpFormat[4];
  SWavInfo wavinfo = {0, 0, 0, 0, 0, 0};

  /* pseudo use to avoid unused symbols */
  (void)ReadWavShort;
  (void)ReadWavInt;
  (void)ReadWavFloat;
  #ifdef TWI_SUPPORT_BWF
  (void)ReadBWF;
  #endif
  (void)ResetWavIn;
  (void)(tmp_return_val);

  /* param check */
  if (!filename)      goto bail;
  if (!samplerate)    goto bail;
  if (!channels)      goto bail;
  if (!samplesInFile) goto bail;
  if (!bps)           goto bail;

  self = (WAVEFILEIN*)calloc(1, sizeof(WAVEFILEIN));
  if (!self)          goto bail; /* return NULL; */

  self->theFile = fopen(filename, "rb");
  if (!self->theFile) goto bail;


  /* read RIFF-chunk */
  if (fread(tmpFormat, 1, 4, self->theFile) != 4) {
    goto bail;
  }

  if (strncmp("RIFF", tmpFormat, 4)) {
    if (strncmp("RF64", tmpFormat, 4)) {
      goto bail;
    }
  }

  /* Read RIFF size. Ignored. */
  fread_LE(&dmy, 4, 1, self->theFile);

  /* read WAVE-chunk */
  if (fread(tmpFormat, 1, 4, self->theFile) != 4) {
    goto bail;
  }

  if (strncmp("WAVE", tmpFormat, 4)) {
    goto bail;
  }

  /* read format/bext-chunk */
  if (fread(tmpchunkID, 1, 4, self->theFile) != 4) {
    goto bail;
  }

  /* read Chunk loop */
  /* skip some potential chunks up to fmt chunk, but read ds64 and bext */
  while ( strncmp("fmt ", tmpchunkID, 4) != 0 ) {
    unsigned int chunkSize = 0;

    /* stop when we can't read the chunk length */
    if (fread_LE(&chunkSize, 1, 4, self->theFile) != 4) {
      goto bail;
    }

#ifdef TWI_SUPPORT_BWF
    if ( strncmp("bext", tmpchunkID, 4) == 0 ) {
      int err = __ReadBextChunk( self, &chunkSize );
      if (err != __TWI_SUCCESS) goto bail;
    }
#endif

    if ( strncmp("ds64", tmpchunkID, 4) == 0 ) {
      int err = __ReadDs64Chunk( self, &chunkSize );
      if (err != __TWI_SUCCESS) goto bail;
    }

    /* skip remaining chunk data */
    while (chunkSize > 0) {
      int nulbuf;
      if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1) {
        goto bail;
      }
      chunkSize -= 1;
    }

    /* read next chunk header */
    if (fread(tmpchunkID, 1, 4, self->theFile) != 4) {
      goto bail;
    }
  }

  /* after the above while() we should now be at the fmt-chunk */
  if (fread_LE(&fmt_chunk.chunkSize, 4, 1, self->theFile) != 1) {
    goto bail;
  }

  /* read fmt info */
  fread_LE(&(wavinfo.formatType),     2, 1, self->theFile);
  fread_LE(&(wavinfo.channelCount),   2, 1, self->theFile);
  fread_LE(&(wavinfo.sampleRate),     4, 1, self->theFile);
  fread_LE(&(wavinfo.bytesPerSecond), 4, 1, self->theFile);
  fread_LE(&(wavinfo.blockAlignment), 2, 1, self->theFile);
  fread_LE(&(wavinfo.bitsPerSample),  2, 1, self->theFile);

  if (wavinfo.formatType == -2) {     // WAVE_FORMAT_EXTENSIBLE
    fseek(self->theFile, 8, SEEK_CUR);                         // skip channel mask
    fread_LE(&(wavinfo.formatType), 2, 1, self->theFile);      // part of GUID
    fseek(self->theFile, 14, SEEK_CUR);                        // skip rest of GUID
    offset = fmt_chunk.chunkSize - 40;
  }
  else {
    offset = fmt_chunk.chunkSize - 16;
  }

  if (wavinfo.formatType == 0x0001) {     // WAVE_FORMAT_PCM
    if((wavinfo.bitsPerSample != 16) && (wavinfo.bitsPerSample != 24)) {
      /* we do only support 16 and 24 bit PCM audio */
      goto bail;
    }
  }
  else if(wavinfo.formatType == 0x0003) { // WAVE_FORMAT_IEEE_FLOAT
    if(wavinfo.bitsPerSample != 32) {
      /* we do only support 32 bit IEEE float audio */
      goto bail;
    }
  }
  else {
    /* we support only formatType 0x01 and 0x03 */
    goto bail;
  }

  /* Skip rest of fmt header if any. */
  for (; offset > 0; offset--) {
    tmp_return_val = fread(&dmy, 1, 1, self->theFile);
  }

  /* endless chunk reading loop, this loop exits when we reach the "data" chunk or eof */
  do {
    int tmpChunkSize = 0;

    /* Read data chunk ID */
    if (fread(data_chunk.chunkID, 1, 4, self->theFile) != 4) {
      goto bail;
    }

    /* Read chunk length */
    if (fread_LE(&tmpChunkSize, 4, 1, self->theFile) != 1) {
      goto bail;
    }

    /* Check for data chunk signature. */
    if (strncmp("data", data_chunk.chunkID, 4) == 0) {
      data_chunk.chunkSize = tmpChunkSize;
      break;
    }

    /* unused 1 byte present, if size is odd */
    /* see https://www.daubnet.com/en/file-format-riff */
    if ( tmpChunkSize%2 ){
      tmpChunkSize++;
    }

    /* Jump over non data chunk. */
    for (;tmpChunkSize > 0; tmpChunkSize--) {
      tmp_return_val = fread(&dmy, 1, 1, self->theFile);
    }

  } while (!feof(self->theFile));


  fgetpos(self->theFile, &self->dataChunkPos);

  /* we should now be at the data chunk at the start of PCM data */
  *samplerate     = wavinfo.sampleRate;
  *channels       = wavinfo.channelCount;
  if (data_chunk.chunkSize == 0xffffffff) {
    *samplesInFile = 0xffffffff;  /* for RF64 we return "-1" */
  } else {
    *samplesInFile  = data_chunk.chunkSize / wavinfo.channelCount;
    *samplesInFile /= ((wavinfo.bitsPerSample + 7) / 8);
  }
  *bps = wavinfo.bitsPerSample;

  self->position = 0;
  self->bps      = wavinfo.bitsPerSample;

  /* hack for broken lengths */
  if (data_chunk.chunkSize == 0xffffffff) {
    if(self->length == 0) {
      long dataChunkStart = ftell(self->theFile);
      fseek(self->theFile, 0, SEEK_END);
      long dataChunkEnd = ftell(self->theFile);
      fseek(self->theFile,dataChunkStart, SEEK_SET);
      self->length = (dataChunkEnd - dataChunkStart);
    } else {
      self->length *= wavinfo.channelCount;
    }
  } else {
    self->length = *samplesInFile * wavinfo.channelCount;
  }

  return self;

bail:
  if ( NULL != self ) {
      free(self);
  }
  
  return NULL;
}

#ifdef TWI_SUPPORT_BWF
static int __ReadBextChunk( WAVEFILEIN* self, unsigned int* chunkSize )
{
  unsigned int bextSize = *chunkSize;

  self->loudnessInfo = (WAVEIN_LOUDNESSINFO *) calloc(1, sizeof(WAVEIN_LOUDNESSINFO));
  if (self->loudnessInfo == NULL) return __TWI_ERROR;

  if (bextSize>=602) {       /* minimum size bext-data, w/o 'CodingHistory' */
    int i;
    signed short readBuf=0;
    signed int nulbuf=0;

    /* first skip all descriptive data */
    for(i=0; i<412; i++) {
      if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
        return __TWI_ERROR;
      }
      bextSize -=1;
    }
    /* second, read loudness data */
    fread_LE(&readBuf, 2, 1, self->theFile);
    bextSize -=2;
    self->loudnessInfo->loudnessVal = (float)readBuf * 0.01f;

    fread_LE(&readBuf, 2, 1, self->theFile);
    bextSize -=2;
    self->loudnessInfo->loudnessRange = (float)readBuf * 0.01f;

    fread_LE(&readBuf, 2, 1, self->theFile);
    bextSize -=2;
    self->loudnessInfo->maxTruePeakLevel = (float)readBuf * 0.01f;

    fread_LE(&readBuf, 2, 1, self->theFile);
    bextSize -=2;
    self->loudnessInfo->maxMomentaryLoudnes = (float)readBuf * 0.01f;

    fread_LE(&readBuf, 2, 1, self->theFile);
    bextSize -=2;
    self->loudnessInfo->maxShortTermLoudness = (float)readBuf * 0.01f;

    /* skip reserved data */
    for(i=0; i<180; i++) {
      if (fread_LE(&nulbuf, 1, 1, self->theFile) != 1){
        return __TWI_ERROR;
      }
      bextSize -= 1;
    }
  }

  *chunkSize = bextSize;
  return __TWI_SUCCESS;
}


static void ReadBWF(
                    WAVEFILEIN*           self,
                    WAVEIN_LOUDNESSINFO** wavInLoudness
                   )
{
  *wavInLoudness   = self->loudnessInfo;
}
#endif


static int __ReadDs64Chunk( WAVEFILEIN* self, unsigned int* chunkSize )
{
  unsigned int nulbuf_hi = 0;
  unsigned int nulbuf_lo = 0;
  int ds64Size = *chunkSize;

  if (ds64Size < 28) return __TWI_ERROR;

  /* skip RIFF size low+high */
  if (fread_LE(&nulbuf_lo, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  if (fread_LE(&nulbuf_hi, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  ds64Size -= 8;
#ifdef PRINT_HDR
  printf("ds64:riff size %li\n", ((TWI_UINT64)nulbuf_hi << 32) + (TWI_UINT64)nulbuf_lo);
#endif

  /* skip datasize size low+high */
  if (fread_LE(&nulbuf_lo, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  if (fread_LE(&nulbuf_hi, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  ds64Size -= 8;
#ifdef PRINT_HDR
  printf("ds64:data size %li\n", ((TWI_UINT64)nulbuf_hi << 32) + (TWI_UINT64)nulbuf_lo);
#endif

  /* read sampleCount */
  if (fread_LE(&nulbuf_lo, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  if (fread_LE(&nulbuf_hi, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  self->length = ((TWI_UINT64)nulbuf_hi << 32) + (TWI_UINT64)nulbuf_lo;
  ds64Size -= 8;
#ifdef PRINT_HDR
  printf("ds64:sample count %li\n", self->length);
#endif

  /* skip tablesize */
  if (fread_LE(&nulbuf_lo, 4, 1, self->theFile) != 1) return __TWI_ERROR;
  ds64Size -= 4;
#ifdef PRINT_HDR
  printf("ds64:table size %i\n", nulbuf_lo);
#endif

  /* any table entries are implicitly skipped outside */

  if (ds64Size < 0) return __TWI_ERROR;

  *chunkSize = ds64Size;
  return __TWI_SUCCESS;
}


static int __ReadSample16(
                          WAVEFILEIN* self, 
                          int*        sample,
                          int         scale
                          )
{
  size_t cnt;
  short v = 0;

  cnt = fread(&v, 2, 1, self->theFile);

  if (cnt != 1) {
    return __TWI_ERROR;
  }

  self->position += 1;

#ifdef __TWI_BE
  v = BigEndian16(v);
#endif

  if ((scale - 16) > 0)
    *sample = v << (scale - 16);
  else
    *sample = v >> (16 - scale);

  return __TWI_SUCCESS;
}


static int __ReadSample24(
                          WAVEFILEIN* self, 
                          int*        sample,
                          int         scale
                          )
{
  size_t cnt;
  int v = 0;

  cnt = fread(&v, 3, 1, self->theFile);

  if (cnt != 1) {
    return __TWI_ERROR;
  }

  self->position += 1;

#ifdef __TWI_BE
  v = BigEndian32(v);
#endif

  if (v >= 0x800000) {
    v |= 0xff000000;
  }

  if ((scale - 24) > 0)
    *sample = v << (scale - 24);
  else
    *sample = v >> (24 - scale);

  return __TWI_SUCCESS;
}

static int __ReadSample32(
                          WAVEFILEIN* self,
                          float*      sample
                          )
{
  size_t cnt;
  union fl_int {
    float v_float;
    int   v_int;
  };
  union fl_int v;
      
  cnt = fread(&v, 4, 1, self->theFile);
  if (cnt != 1) {
    return __TWI_ERROR;
  }

  self->position += 1;
#ifdef __TWI_BE
  v.v_int = BigEndian32(v.v_int);
#endif

  *sample = v.v_float;

  return __TWI_SUCCESS;
}



static int __ReadSampleInternal(
                                WAVEFILEIN* self, 
                                int*        sample,
                                int         scale
                                )
{
  int err;

  if (!self) {
    return __TWI_ERROR;
  }

  switch (self->bps) {

  case 16:
    err = __ReadSample16(self, sample, scale);
    break;

  case 24:
    err = __ReadSample24(self, sample, scale);
    break;

  default:
    err = __TWI_ERROR;
    break;
  }

  return err;
}

/* not fully tested */
/* this function returns normalized values in the range +32767..-32768 */
static int ReadWavShort(
                        WAVEFILEIN*   self, 
                        short         sampleBuffer[], 
                        unsigned int  nSamplesToRead,
                        unsigned int* nSamplesRead
                        )
{
  unsigned int i;
  int err = __TWI_SUCCESS;

  if (!sampleBuffer) return __TWI_ERROR;

  /* check if we have enough samples left, if not,
     set nSamplesToRead to number of samples left. */
  if (self->position + nSamplesToRead > self->length) {
    nSamplesToRead = self->length - self->position;
  }

  for (i=0; i< nSamplesToRead; i++) {
    if(self->bps == 32)
    {
      float tmp;
      err = __ReadSample32(self, &tmp);
      if (err != __TWI_SUCCESS) return err;

      sampleBuffer[i] = (int)(tmp * 32768.0f);
    }
    else
    {
      int tmp;
      err = __ReadSampleInternal(self, &tmp, 16);
      if (err != __TWI_SUCCESS) return err;

      sampleBuffer[i] = (short)tmp;
    }
    *nSamplesRead += 1;
  }

  return __TWI_SUCCESS;
}

/* not fully tested */
/* this function returns normalized values in the range +8388607..-8388608 */
static int ReadWavInt(
                      WAVEFILEIN*   self,
                      int           sampleBuffer[],
                      unsigned int  nSamplesToRead,
                      unsigned int* nSamplesRead
                      )
{
  unsigned int i;
  int err = __TWI_SUCCESS;

  if (!sampleBuffer) return __TWI_ERROR;

  /* check if we have enough samples left, if not,
     set nSamplesToRead to number of samples left. */
  if (self->position + nSamplesToRead > self->length) {
    nSamplesToRead = self->length - self->position;
  }

  for (i = 0; i < nSamplesToRead; i++) {
    if (self->bps == 32)
    {
      float tmp;
      err = __ReadSample32(self, &tmp);
      if (err != __TWI_SUCCESS) return err;

      sampleBuffer[i] = (int)(tmp * 8388608.0f);
    }
    else
    {
      int tmp;
      err = __ReadSampleInternal(self, &tmp, 24);
      if (err != __TWI_SUCCESS) return err;

      sampleBuffer[i] = tmp;
    }
    *nSamplesRead += 1;
  }

  return __TWI_SUCCESS;
}

/* this function returns normalized values in the range +-1.0 */
static int ReadWavFloat(
                        WAVEFILEIN*   self, 
                        float         sampleBuffer[], 
                        unsigned int  nSamplesToRead,
                        unsigned int* nSamplesRead
                        )
{
  unsigned  i;
  int err = __TWI_SUCCESS;

  if (!sampleBuffer) return __TWI_ERROR;
	
  /* check if we have enough samples left, if not,
     set nSamplesToRead to number of samples left. */
  if (self->position + nSamplesToRead > self->length) {
    nSamplesToRead = self->length - self->position;
  }

  for (i=0; i< nSamplesToRead; i++) {
    int tmp;

    if(self->bps == 32)
    {
      err = __ReadSample32(self, &sampleBuffer[i]);
      if (err != __TWI_SUCCESS) return err;
    }
    else
    {
      err = __ReadSampleInternal(self, &tmp, 24);
      if (err != __TWI_SUCCESS) return err;

      sampleBuffer[i] = (float)tmp / 8388608.0f;
    }
    *nSamplesRead += 1;
  }

  return __TWI_SUCCESS;
}


static int CloseWavIn(WAVEFILEIN* self)
{
  if (self) {
    if (self->theFile) {
      fclose(self->theFile);
    }
  }
  free(self);

  return __TWI_SUCCESS;
}

static int ResetWavIn(WAVEFILEIN* self)
{
  if (self) {
    if (self->theFile) {
      fsetpos(self->theFile, &self->dataChunkPos);
      self->position = 0;
    }
  }
  return __TWI_SUCCESS;
}

/*------------- local subs ----------------*/

static size_t fread_LE(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
#ifdef __TWI_LE
  return fread(ptr, size, nmemb, stream);
#endif
#ifdef __TWI_BE
  unsigned char x[sizeof(int)];
  unsigned char* y = (unsigned char*)ptr;
  int i;
  int len;

  len = fread(x, size, nmemb, stream);

  for (i = 0; i < size * nmemb; i++) {
    *y++ = x[size * nmemb - i - 1];
  }

  return len;
#endif
}

#ifdef __TWI_BE
static short BigEndian16(short v)
{
  short a = (v & 0x0ff);
  short b = (v & 0x0ff00) >> 8;

  return a << 8 | b;
}


static int BigEndian32(int v)
{
  int a = (v & 0x0ff);
  int b = (v & 0x0ff00) >> 8;
  int c = (v & 0x0ff0000) >> 16;
  int d = (v & 0xff000000) >> 24;

  return a << 24 | b << 16 | c << 8 | d;
}
#endif

#endif  /* __TINYWAVEIN_C_H__ */
