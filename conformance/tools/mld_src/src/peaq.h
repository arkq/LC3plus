/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef PEAQ_H
#define PEAQ_H

#include <stdio.h>

#define NUM_BANDS 40                          // filterbank bands
#define SUBSAMP_FB 16                         // filterbank subsampling
#define SUBSAMP_EP 6                          // exitation pattern subsampling
#define SUBSAMP_TOT (SUBSAMP_FB * SUBSAMP_EP) // pcm -> excitation pattern subsampling

// constants used by filters influcenced by subsampling
#define SUBSAMP_TIME SUBSAMP_FB                // can be 32 for quicker fallof
#define SPREADING_C ((double)SUBSAMP_TIME)     // frequency domain spreading
#define SMEARING_TAPS (12 * 32 / SUBSAMP_TIME) // time domain smearing 1 filter legnth
#define SMEARING_OFF (SMEARING_TAPS / 2 - 1)   // time domain smearing 1
#define SMEARING_C (6.0 * SUBSAMP_TIME)        // time domain smearing 2

#define CHUNK_SIZE_EP 20                            // processing chunk
#define CHUNK_SIZE_FB (CHUNK_SIZE_EP * SUBSAMP_EP)  // filterbank samples
#define CHUNK_SIZE_PCM (CHUNK_SIZE_FB * SUBSAMP_FB) // pcm samples

#define PEAQ_SAMPLERATE 48000 // expected input sapmling rate

typedef struct peaq Peaq;

// allocate peaq struct, for normal use cases, playback level should be 92
Peaq* peaq_init(double playback_level);

// free peaq struct
void peaq_free(Peaq* pq);

// update peaq with pcm data. The is expected to be in ragne [-1,1]
void peaq_update(Peaq* pq, const float* pcm, int len);

// call finish once after last update to flush internal buffers
void peaq_finish(Peaq* pq);

// print maximum loudness difference between to files
// segment size is the number of frames that should be combined
void  peaq_print_mld(Peaq* pq1, Peaq* pq2, int segment_size, FILE* segment_outfile);
float peaq_get_mld(Peaq* pq1, Peaq* pq2, int segment_size);

#endif
