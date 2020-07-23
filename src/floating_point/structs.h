/******************************************************************************
*                        ETSI TS 103 634 V1.2.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef STRUCTS_H
#define STRUCTS_H

#include "defines.h"

typedef struct {
  LC3_FLOAT r; /* real part */
  LC3_FLOAT i; /* imaginary part */
} Complex;

typedef struct {
  LC3_INT length;
  void *handle;
} Fft;

typedef struct {
  LC3_INT length;
  Fft fft;
} Dct2;

typedef struct {
  LC3_INT length;
  Fft fft;
} Dct3;

typedef struct {
  LC3_INT length;
  Complex *twid1;
  Complex *twid2;
  Fft fft;
} Dct4;

typedef struct {
  LC3_INT length;
  LC3_INT leading_zeros;
  LC3_INT mem_length;
  const LC3_FLOAT *window;
  LC3_FLOAT *mem;
  Dct4 dct;
} Mdct;

typedef struct {
  uint32_t ac_low_fl;
  uint32_t ac_range_fl;
  int BER_detect;
} Decoder_State_fl;

typedef struct {
  LC3_INT bp;
  LC3_INT low;
  LC3_INT range;
  LC3_INT cache;
  LC3_INT carry;
  LC3_INT carry_count;
  uint8_t *ptr;
  LC3_INT *bp_side;
  LC3_INT *mask_side;
} Encoder_State_fl;

typedef struct {
    LC3_INT   nbLostCmpt;
    LC3_INT   prevBfi;
    LC3_INT   prevprevBfi;
    LC3_FLOAT q_d[MAX_LEN];
    LC3_FLOAT q_d_prev[MAX_LEN];
} PlcSetup;

typedef struct PlcAdvSetup PlcAdvSetup;

typedef struct {
    LC3_FLOAT cum_alpha;
    LC3_INT   seed;
} PlcNsSetup;


#endif
