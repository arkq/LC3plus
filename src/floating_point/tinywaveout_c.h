/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __TINYWAVEOUT_C_H__
#define __TINYWAVEOUT_C_H__

/*#define TWO_SUPPORT_BWF*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#ifdef TWO_SUPPORT_BWF
#include <string.h>
#endif

/***** Interface *********************************************************/

#ifndef TWO_UINT64
  #if !(defined(WIN32))
    #include <stdint.h>
    #define TWO_UINT64 uint64_t
  #else
    #define TWO_UINT64 unsigned __int64
  #endif
#endif

typedef struct WAVEFILEOUT WAVEFILEOUT;
#ifdef TWO_SUPPORT_BWF
typedef struct WAVEOUT_LOUDNESSINFO WAVEOUT_LOUDNESSINFO;
#endif

#define __TWO_SUCCESS  (0)
#define __TWO_ERROR    (-1)

static WAVEFILEOUT* CreateWav(
                              const char *fileName, 
                              const unsigned int sampleRate, 
                              const unsigned int numChannels, 
                              const unsigned int bps
                              );

#ifdef TWO_SUPPORT_BWF
static WAVEFILEOUT* CreateWavBWF(
                                 const char                 *fileName, 
                                 const unsigned int          sampleRate, 
                                 const unsigned int          numChannels, 
                                 const unsigned int          bps,
                                 const WAVEOUT_LOUDNESSINFO *hBwfData
                                 );
#endif

/* this function expects values in the 16 bit range +32767..-32768 */
static int WriteWavShort(
                         WAVEFILEOUT* self, 
                         short        sampleBuffer[], 
                         unsigned int nSamples
                         );

/* this function expects values in the 24 bit range +8388607..-8388608 */
static int WriteWavLong(
                        WAVEFILEOUT* self, 
                        int          sampleBuffer[], 
                        unsigned int nSamples
                        );

/* this function expects normalized values in the range +-1.0f */
static int WriteWavFloat(
                         WAVEFILEOUT* self, 
                         float        sampleBuffer[], 
                         unsigned int nSamples
                         );

static int CloseWav(WAVEFILEOUT* self);
#ifdef TWO_SUPPORT_BWF
static int CloseWavBWF(WAVEFILEOUT* self, WAVEOUT_LOUDNESSINFO bextData);
#endif

/***** Implementation *********************************************************/

#if defined (__i386__) || defined (_M_IX86) || defined (_M_X64) || defined (__x86_64__) || defined (__arm__) || defined (__xtensa__) || defined (__aarch64__) || defined (__EMSCRIPTEN__)
#define __TWO_LE    /* _T_iny _W_ave _O_ut _L_ittle _E_ndian */
#endif

#if defined (__POWERPC__)
#define __TWO_BE    /* _T_iny _W_ave _O_ut _B_ig _E_ndian */
#endif

#if defined (__sparc__)
#define __TWO_BE    /* _T_iny _W_ave _O_ut _B_ig _E_ndian */
#endif

#if ! defined (__TWO_LE) && ! defined (__TWO_BE)
#error unknown processor
#endif

/*--- local types/structs ----------------------------------*/

#if defined(_MSC_VER)
  #pragma pack(push, 1)
#else
  #pragma pack(1)
#endif


#ifdef TWO_SUPPORT_BWF
struct WAVEOUT_LOUDNESSINFO {
  float loudnessVal;
  float loudnessRange;
  float maxTruePeakLevel;
  float maxMomentaryLoudnes;
  float maxShortTermLoudness;
};
#endif


typedef struct __tinyWaveOutHeader
{
  unsigned int   riffType;       /* 'RIFF/RF64' */
  unsigned int   riffSize;       /* file size/-1 */
  unsigned int   waveType;       /* 'WAVE' */
} __tinyWaveOutHeader;

typedef struct __tinyWaveOutDs64Chunk
{
  unsigned int   formatType;  /* = 'JUNK/ds64' */
  unsigned int   formatSize;  /* size info */
  TWO_UINT64     riffSize64;
  TWO_UINT64     dataSize64;
  TWO_UINT64     sampleCount64;
  unsigned int   tableLength;   /* optional tables, always 0 for tinywaveout */
  /* here: optional tables */
} __tinyWaveOutDs64Chunk;

#ifdef TWO_SUPPORT_BWF
typedef struct __tinyWaveOutBextChunk
{
  unsigned int   formatType;  /* = 'bext' */
  unsigned int   formatSize;  /* size info */

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

  unsigned char codingHistory; /* ASCII: <<History coding>> - undefined length! */
                               /* for variable length, mve this out of this struct */
} __tinyWaveOutBextChunk;
#endif

typedef struct __tinyWaveOutFmtChunk
{
  unsigned int   formatType;
  unsigned int   formatSize;

  unsigned short formatTag;
  unsigned short numChannels;
  unsigned int   sampleRate;
  unsigned int   bytesPerSecond;
  unsigned short blockAlignment;
  unsigned short bitsPerSample;

  /* wav fmt ext hdr here */
} __tinyWaveOutFmtChunk;

typedef struct __tinyWaveOutDataChunk
{
  unsigned int   dataType;
  unsigned int   dataSize;

} __tinyWaveOutDataChunk;


struct WAVEFILEOUT {
  /* for reasons of memory alignment, have 64 bit data types first */
  TWO_UINT64        dataSize;
  TWO_UINT64        dataSizeLimit;    /* maximum size for data chunk for 4 GB files */              
  TWO_UINT64        dataSizeLimit64;  /* maximum size for data chunk for 2^64 B addressable files */
  TWO_UINT64        clipCount;
  FILE*             theFile;
  unsigned int      junkChunkOffset;
  unsigned int      fmtChunkOffset;
#ifdef TWO_SUPPORT_BWF
  unsigned int      bextChunkOffset;
#endif
  unsigned int      dataChunkOffset;
  unsigned int      bps;

  /* only needed for RF64: */
  unsigned int      blockAlign;     /* only needed for close() of ds64 chunk */
};


/*--- local protos --------------------------------------------------*/
static __inline unsigned int BigEndian32(char, char, char, char);
static __inline unsigned int LittleEndian32(unsigned int);
static __inline unsigned int LittleEndian32s(int);
static __inline short LittleEndian16(short);
static __inline TWO_UINT64 LittleEndian64(TWO_UINT64);
#ifdef TWO_SUPPORT_BWF
static unsigned int EncodeLoudness(float);
#endif
static __inline int __dataSizeChk( WAVEFILEOUT* self, int newbytes );

#if defined(_MSC_VER)
  #pragma pack(pop)
#else
  #pragma pack()
#endif


#ifdef TWO_SUPPORT_BWF
static void setDefaultLoudness(WAVEOUT_LOUDNESSINFO *x)
{
  x->loudnessVal          = 1.0f;
  x->loudnessRange        = 2.0f;
  x->maxTruePeakLevel     = 3.0f;
  x->maxMomentaryLoudnes  = 4.0f;
  x->maxShortTermLoudness = 5.0f;
}
#endif

static WAVEFILEOUT* __CreateWavInternal(
                                        const char                *fileName, 
                                        const unsigned int         sampleRate, 
                                        const unsigned int         numChannels, 
                                        const unsigned int         bps
#ifdef TWO_SUPPORT_BWF
                                        ,
                                        const WAVEOUT_LOUDNESSINFO *hBwfData
#endif
                                        )
{
  WAVEFILEOUT* self = NULL;
  __tinyWaveOutHeader whdr;
  __tinyWaveOutDs64Chunk ds64ch;
#ifdef TWO_SUPPORT_BWF
  __tinyWaveOutBextChunk wbextch;
#endif
  __tinyWaveOutFmtChunk wfch;
  __tinyWaveOutDataChunk wdch;
  unsigned int blockAlignment = 0;
  unsigned int ByteCnt = 0; /* Byte counter for fwrite */

  /* pseudo use to avoid unused symbols */
  (void)WriteWavShort;
  (void)WriteWavLong;
  (void)WriteWavFloat;
#ifdef TWI_SUPPORT_BWF
  (void)CreateWavBWF;
  (void)CloseWavBWF;
  (void)setDefaultLoudness;
#endif

  /* param check */
  if (!fileName)                           goto bail;
  if (sampleRate == 0)                     goto bail;
  if (sampleRate > 768000)                 goto bail;
  if (numChannels == 0)                    goto bail;
  if (numChannels > 64)                    goto bail;
  if (bps != 16 && bps != 24 && bps != 32) goto bail;

  self = (WAVEFILEOUT*)calloc(1, sizeof(WAVEFILEOUT));
  if (!self) goto bail; /* return NULL; */

  self->theFile = fopen(fileName, "wb+");
  if (!self->theFile)                      goto bail;

  /* WAV-Header */
  whdr.riffType = BigEndian32('R','I','F','F');
  whdr.riffSize = LittleEndian32(0xffffffff);  /* set to maximum, if fseek() doesn't work later */
  whdr.waveType = BigEndian32('W','A','V','E');
  /* write to file */
  ByteCnt = 0;
  ByteCnt += fwrite(&whdr, 1, sizeof(whdr), self->theFile);

  /* ds64/JUNK-Chunk */
  ds64ch.formatType = BigEndian32('J','U','N','K');
  ds64ch.formatSize = LittleEndian32(sizeof(__tinyWaveOutDs64Chunk) - 8);
  ds64ch.riffSize64    = (TWO_UINT64) -1;
  ds64ch.dataSize64    = (TWO_UINT64) -1;
  ds64ch.sampleCount64 = (TWO_UINT64) -1;
  ds64ch.tableLength   = 0;
  self->junkChunkOffset = ByteCnt;
  ByteCnt += fwrite(&ds64ch, 1, sizeof(ds64ch), self->theFile);

#ifdef TWO_SUPPORT_BWF
  /* BEXT-Chunk */
  if (hBwfData) {
    memset(&wbextch, 0, sizeof(__tinyWaveOutBextChunk));
    wbextch.formatType = BigEndian32('b','e','x','t');
    wbextch.formatSize = LittleEndian32(sizeof(__tinyWaveOutBextChunk) - 8);
    wbextch.version    = 0x0002;
    wbextch.loudnessVal          = EncodeLoudness(hBwfData->loudnessVal);
    wbextch.loudnessRange        = EncodeLoudness(hBwfData->loudnessRange);
    wbextch.maxTruePeakLevel     = EncodeLoudness(hBwfData->maxTruePeakLevel);
    wbextch.maxMomentaryLoudnes  = EncodeLoudness(hBwfData->maxMomentaryLoudnes);
    wbextch.maxShortTermLoudness = EncodeLoudness(hBwfData->maxShortTermLoudness);
    /* t.b.d.:  more values */

    /* write to file */
    self->bextChunkOffset = ByteCnt;
    ByteCnt += fwrite(&wbextch, 1, sizeof(wbextch), self->theFile);
  }
#endif

  /* FMT-Chunk */
  wfch.formatType = BigEndian32('f','m','t',' ');
  wfch.formatSize = LittleEndian32(16);
  switch (bps) {
  case 16:
  case 24:
    wfch.formatTag    = LittleEndian16(0x0001);  /* WAVE_FORMAT_PCM */
    break;
  case 32:
    wfch.formatTag    = LittleEndian16(0x0003);  /* WAVE_FORMAT_IEEE_FLOAT */
    break;
  default:
    goto bail;
  }
  self->bps           = bps;
  wfch.bitsPerSample  = LittleEndian16(bps);
  wfch.numChannels    = LittleEndian16(numChannels);
  blockAlignment      = numChannels * (bps >> 3);
  wfch.blockAlignment = LittleEndian16(blockAlignment);
  wfch.sampleRate     = LittleEndian32(sampleRate);
  wfch.bytesPerSecond = LittleEndian32(sampleRate * blockAlignment);
  /* tbd: wavfmt ext hdr here */
  /* write to file */
  self->fmtChunkOffset = ByteCnt;
  ByteCnt += fwrite(&wfch, 1, sizeof(wfch), self->theFile);

  /* DATA-Chunk */
  self->dataChunkOffset = ByteCnt;
  wdch.dataType = BigEndian32('d','a','t','a');
  wdch.dataSize = LittleEndian32(0xffffffff - ByteCnt);  /* yet unknown. set to maximum of 4 GB file */
  /* write to file */
  ByteCnt += fwrite(&wdch, 1, sizeof(wdch), self->theFile);

  self->dataSize = 0;

  /* self->dataSizeLimit = 0x7fffffff - ByteCnt; */            /* maximum size for data chunk for 2 GB files */
  self->dataSizeLimit = 0xffffffff - ByteCnt;                  /* maximum size for data chunk for 4 GB files */
  self->dataSizeLimit64 =  0xffffffffffffffff - ByteCnt;       /* maximum size for data chunk for 64 bit addressable files */

  self->clipCount = 0;
  self->blockAlign = blockAlignment;

  return self;

 bail:
  if ( NULL != self) {
      free(self);
  }
  
  return NULL;
}

static WAVEFILEOUT* CreateWav(
                              const char*        fileName, 
                              const unsigned int sampleRate, 
                              const unsigned int numChannels, 
                              const unsigned int bps
                              )
{
#ifdef TWO_SUPPORT_BWF
  return __CreateWavInternal(fileName, sampleRate, numChannels, bps, NULL);
#else
  return __CreateWavInternal(fileName, sampleRate, numChannels, bps);
#endif
}

#ifdef TWO_SUPPORT_BWF
static WAVEFILEOUT* CreateWavBWF(
                                 const char                 *fileName, 
                                 const unsigned int          sampleRate, 
                                 const unsigned int          numChannels, 
                                 const unsigned int          bps,
                                 const WAVEOUT_LOUDNESSINFO *hBwfData
                                 )
{
  return __CreateWavInternal(fileName, sampleRate, numChannels, bps, hBwfData);
}
#endif

#define MAX_PCM16 (+32767)
#define MIN_PCM16 (-32768)
static __inline int CLIP_PCM16(int sample, TWO_UINT64* clipcount)
{
  int tmp = sample;

  if (sample >= MAX_PCM16) {
    tmp = MAX_PCM16;
    (*clipcount)++;
  } else {
    if (sample <= MIN_PCM16) {
      tmp = MIN_PCM16;
      (*clipcount)++;
    }
  }

  return tmp;
}

#define MAX_PCM24 (+8388607)
#define MIN_PCM24 (-8388608)
static __inline int CLIP_PCM24(int sample, TWO_UINT64* clipcount)
{
  int tmp = sample;

  if (sample >= MAX_PCM24) {
    tmp = MAX_PCM24;
    (*clipcount)++;
  } else {
    if (sample <= MIN_PCM24) {
      tmp = MIN_PCM24;
      (*clipcount)++;
    }
  }

  return tmp;
}

#define MAX_FLOAT32 (+1.0f)
#define MIN_FLOAT32 (-1.0f)
static __inline float CLIP_FLOAT32(float sample, TWO_UINT64* clipcount)
{
  float tmp = sample;

  if (sample >= MAX_FLOAT32) {
    tmp = MAX_FLOAT32;
    (*clipcount)++;
  } else {
    if (sample <= MIN_FLOAT32) {
      tmp = MIN_FLOAT32;
      (*clipcount)++;
    }
  }

  return tmp;
}



static int __WriteSample16(
                           WAVEFILEOUT* self, 
                           int          sample,
                           int          scale
                           )
{
  size_t cnt;
  short v;

  if ((scale - 16) > 0)
    sample = sample >> (scale - 16);
  else
    sample = (int) ((uint32_t) sample << (16 - scale));

  v = (short)CLIP_PCM16(sample, &(self->clipCount));
#ifdef __TWO_BE
  v = LittleEndian16(v);
#endif

  cnt = fwrite(&v, sizeof(short), 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 2;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}


static int __WriteSample24(
                           WAVEFILEOUT* self, 
                           int          sample,
                           int          scale
                           )
{
  size_t cnt;
  int v;

  if ((scale - 24) > 0)
    sample = sample >> (scale - 24);
  else
    sample = (int) (((unsigned int)sample) << (24 - scale));

  v = (int)CLIP_PCM24(sample, &(self->clipCount));
#ifdef __TWO_BE
  v = LittleEndian32s(v);
#endif
  cnt = fwrite(&v, 3, 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 3;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}


static int __WriteSample32(
                           WAVEFILEOUT* self, 
                           float        sample
                           )
{
  size_t cnt;
  union fl_int {
    float v_float;
    int   v_int;
  };
  union fl_int v;

#if CLIP_FLOAT
  v.v_float = CLIP_FLOAT32(sample, &(self->clipCount));
#else
  v.v_float = sample;
  if((sample > 1.0f) || (sample <-1.0f))
    self->clipCount++;
#endif

#ifdef __TWO_BE
  v.v_int = LittleEndian32s(v.v_int);
#endif
  cnt = fwrite(&v, 4, 1, self->theFile);

  if (cnt == 1) {
    self->dataSize += 4;
    return __TWO_SUCCESS;
  }

  return __TWO_ERROR;
}


static int __WriteSampleInt(
                            WAVEFILEOUT* self, 
                            int          sample,
                            int          scale
                            )
{
  int err;

  if (!self) return __TWO_ERROR;


  switch (self->bps) {

  case 16:
    err = __WriteSample16(self, sample, scale);
    break;

  case 24:
    err = __WriteSample24(self, sample, scale);
    break;

  default:
    err = __TWO_ERROR;
    break;
  }

  return err;
}


/* this function expects values in the 16 bit range +-32767/8 */
static int WriteWavShort(
                         WAVEFILEOUT* self, 
                         short        sampleBuffer[], 
                         unsigned int nSamples
                         )
{
  unsigned long i;
  int err = __TWO_SUCCESS;

  if (!self)         return __TWO_ERROR;
  if (!sampleBuffer) return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(short))) return __TWO_ERROR;

  for (i=0; i< nSamples; i++) {
    if (self->bps == 32)
    {
      err = __WriteSample32(self, sampleBuffer[i] / 32768.0f);
    }
    else
    {
      err = __WriteSampleInt(self, (int)sampleBuffer[i], 16);
    }
    if (err != __TWO_SUCCESS) return err;
  }

  return __TWO_SUCCESS;
}


/* this function expects values in the 24 bit range +-8388607/8 */
static int WriteWavLong(
                        WAVEFILEOUT* self, 
                        int          sampleBuffer[], 
                        unsigned int nSamples
                        )
{
  unsigned long i;
  int err = __TWO_SUCCESS;

  if (!self)                                       return __TWO_ERROR;
  if (!sampleBuffer)                               return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(int))) return __TWO_ERROR;

  for (i = 0; i < nSamples; i++) {
    if (self->bps == 32)
    {
      err = __WriteSample32(self, sampleBuffer[i] / 8388608.0f);
    }
    else
    {
      err = __WriteSampleInt(self, sampleBuffer[i], 24);
    }
    if (err != __TWO_SUCCESS) return err;
  }

  return __TWO_SUCCESS;
}


/* this function expects normalized values in the range +-1.0 */
#define MAX_FL (+2.0f * 8388608.0f )
#define MIN_FL (-2.0f * 8388608.0f )
#define CLIP_FL(x)   ( ((x) >= MAX_FL) ? MAX_FL : (((x) <= MIN_FL) ? MIN_FL : (x)) )
static int WriteWavFloat(
                         WAVEFILEOUT* self, 
                         float        sampleBuffer[], 
                         unsigned int nSamples
                         )
{
  unsigned int i;
  int err = __TWO_SUCCESS;

  if (!self)         return __TWO_ERROR;
  if (!sampleBuffer) return __TWO_ERROR;
  if (__dataSizeChk(self, nSamples * sizeof(float))) return __TWO_ERROR;

  for (i=0; i<nSamples; i++) {
    if(self->bps == 32)
    {
      err = __WriteSample32(self, sampleBuffer[i]);
    }
    else
    {
      float tmp = CLIP_FL(sampleBuffer[i] * 8388608.0f);   /* CLIP_FL is just to avoid an INT overrun before the actual cast to int, real clipping and counting is done below */
      err = __WriteSampleInt(self, (int) tmp, 24);
    }
    if (err != __TWO_SUCCESS) return err;
  }

  return __TWO_SUCCESS;
}


static int CloseWav(WAVEFILEOUT* self)
{
  unsigned int riffSize_le = 0;
  unsigned int dataSize_le = 0;
  TWO_UINT64 riffSize64    = 0;
  TWO_UINT64 dataSize64    = 0;
  TWO_UINT64 sampleCount64 = 0;    /* nr of samples in the WAVE sense: 1 sample are all pcm samples of all channels of one time slice */
  int mustWriteRF64        = 0;

  if (!self) return __TWO_ERROR;

  /* check for 4 GB (switch to RF64) */
  if ( self->dataSize > self->dataSizeLimit ) {
    /* when we exceed 4 GB: switch from std wave header to RF64 header */
    mustWriteRF64 = 1;
  }

  /* calc header values */
  if (mustWriteRF64 == 0) {
    /* write padding byte if dataSize is uneven */
    int pad = 0;
    if (self->dataSize % 2 > 0) {
      char tmp= 0x00;
      fwrite(&tmp, sizeof(char), 1, self->theFile);
      pad = 1;
    }
    riffSize_le = LittleEndian32(self->dataChunkOffset - 8 + 8 + (unsigned int)self->dataSize + pad);  /* sizeof(hdr) - (8 bytes of riff chunk header) + (8 bytes data chunk header) + sizeof(raw-pcm-data) + padding Byte */
    dataSize_le = LittleEndian32((unsigned int)self->dataSize);
  } else {
    riffSize_le   = 0xffffffff;
    dataSize_le   = 0xffffffff;
    riffSize64    = LittleEndian64( self->dataSize + (TWO_UINT64)self->dataChunkOffset );
    dataSize64    = LittleEndian64( self->dataSize );
    sampleCount64 = LittleEndian64( self->dataSize / (TWO_UINT64)self->blockAlign );
  }

  /* now overwrite length/size values in header with the actual/real ones */
  if (mustWriteRF64 == 0) {
    /* riffsize32 */
    fseek(self->theFile, 4, SEEK_SET);
    fwrite(&riffSize_le, sizeof(riffSize_le), 1, self->theFile);
    /* datasize32 */
    fseek(self->theFile, self->dataChunkOffset + 4, SEEK_SET);
    fwrite(&dataSize_le, sizeof(dataSize_le), 1, self->theFile);
  } else {
    unsigned int rf64sig = BigEndian32('R','F','6','4');
    unsigned int ds64sig = BigEndian32('d','s','6','4');
    
    /* replace RIFF->RF64 */
    fseek(self->theFile, 0, SEEK_SET);
    fwrite(&rf64sig, sizeof(rf64sig), 1, self->theFile);

    /* riffsize32 */
    fseek(self->theFile, 4, SEEK_SET);
    fwrite(&riffSize_le, sizeof(riffSize_le), 1, self->theFile);

    /* replace JUNK->ds64 */
    fseek(self->theFile, self->junkChunkOffset, SEEK_SET);
    fwrite(&ds64sig, sizeof(ds64sig), 1, self->theFile);

    /* riffSize64, dataSize64, sampleCount64 */
    fseek(self->theFile, self->junkChunkOffset + 8, SEEK_SET);
    fwrite(&riffSize64, sizeof(riffSize64), 1, self->theFile);
    fwrite(&dataSize64, sizeof(dataSize64), 1, self->theFile);
    fwrite(&sampleCount64, sizeof(sampleCount64), 1, self->theFile);

    /* datasize32 */
    fseek(self->theFile, self->dataChunkOffset + 4, SEEK_SET);
    fwrite(&dataSize_le, sizeof(dataSize_le), 1, self->theFile);
  }

  fclose(self->theFile);
  free(self);

  return __TWO_SUCCESS;
}

#ifdef TWO_SUPPORT_BWF
static int CloseWavBWF(
                       WAVEFILEOUT* self,
                       WAVEOUT_LOUDNESSINFO bextData
                       )
{
  int wordData;

  if (!self) return __TWO_ERROR;

  if (self->bextChunkOffset) {
    /* Offset for Loudness Data in bext-chunk: 8: Chunk-Header, 412:prev.Data */
    fseek(self->theFile, self->bextChunkOffset+8+412, SEEK_SET);

    wordData = LittleEndian32(EncodeLoudness(bextData.loudnessVal));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.loudnessRange));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxTruePeakLevel));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxMomentaryLoudnes));
    fwrite(&wordData, 2, 1, self->theFile);

    wordData = LittleEndian32(EncodeLoudness(bextData.maxShortTermLoudness));
    fwrite(&wordData, 2, 1, self->theFile);
  }

  return CloseWav(self);
}
#endif

/*------------- local subs ----------------*/

static __inline unsigned int BigEndian32(char a, char b, char c, char d)
{
#ifdef __TWO_LE
    return
    (unsigned int) d << 24 |
    (unsigned int) c << 16 |
    (unsigned int) b <<  8 |
    (unsigned int) a ;
#else
  return
    (unsigned int) a << 24 |
    (unsigned int) b << 16 |
    (unsigned int) c <<  8 |
    (unsigned int) d ;
#endif
}


static __inline unsigned int LittleEndian32(unsigned int v)
{
#ifdef __TWO_LE
  return v;
#else
  return 
    (v & 0x000000FF) << 24 |
    (v & 0x0000FF00) <<  8 |
    (v & 0x00FF0000) >>  8 |
    (v & 0xFF000000) >> 24 ;
#endif
}


/* signed version of the above */
static __inline unsigned int LittleEndian32s(int v)
{
#ifdef __TWO_LE
  return v;
#else
  return 
    (v & 0x000000FF) << 24 |
    (v & 0x0000FF00) <<  8 |
    (v & 0x00FF0000) >>  8 |
    (v & 0xFF000000) >> 24 ;
#endif
}

 
static __inline short LittleEndian16(short v)
{
#ifdef __TWO_LE
  return v;
#else
  return ((v << 8) & 0xFF00) | ((v >> 8) & 0x00FF);
#endif
}

static __inline TWO_UINT64 LittleEndian64(TWO_UINT64 v)
{
#ifdef __TWO_LE
  return v;
#else
  return 
    (v & 0x00000000000000FF) << 56 |
    (v & 0x000000000000FF00) << 40 |
    (v & 0x0000000000FF0000) << 24 |
    (v & 0x00000000FF000000) <<  8 |
    (v & 0x000000FF00000000) >>  8 |
    (v & 0x0000FF0000000000) >> 24 |
    (v & 0x00FF000000000000) >> 40 |
    (v & 0xFF00000000000000) >> 56 ;
#endif
}

#ifdef TWO_SUPPORT_BWF
static unsigned int EncodeLoudness(float x)
{
  int s = (x>0)-(x<0);
  return (int)( x*100.0f + s*0.5f );
}
#endif

static __inline int __dataSizeChk( WAVEFILEOUT* self, int newbytes )
{
  if (!self) return __TWO_ERROR;

  if ( (self->dataSize + ((TWO_UINT64)newbytes) ) > self->dataSizeLimit64 ) {
    return __TWO_ERROR;
  }

  return __TWO_SUCCESS;
}

#endif  /* __TINYWAVEOUT_C_H__ */



