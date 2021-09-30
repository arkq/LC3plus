/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef STRUCTS_H
#define STRUCTS_H

#include "defines.h"
#include "fft/iisfft.h"

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

  LC3_INT32 pc_c_bp;
  LC3_INT32 pc_c_bp_side;
  LC3_INT32 pc_bytes;
  LC3_INT32 pc_b_left;
  LC3_INT32 pc_b_right;
  LC3_INT32 pc_enc;
  LC3_INT32 pc_bfi;
  LC3_INT32 pc_bbi;
  LC3_INT32 pc_be_bp_left;
  LC3_INT32 pc_be_bp_right;
  LC3_INT32 pc_return;
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


typedef struct {
    LC3_FLOAT cum_alpha;
    LC3_INT   seed;
} PlcNsSetup;

typedef struct {
    LC3_INT32 seed;
    LC3_INT32 ns_nbLostCmpt_pc;
    LC3_FLOAT *q_old_res;
    LC3_FLOAT prev_gg;
} pcState;

typedef struct {
    LC3_INT    len;
    LC3_INT    sign;
    LC3_FLOAT* table;
} Cfft;

typedef struct T_IIS_FFT {
    IIS_FFT_DIR sign;
    LC3_INT32   len;
    LC3_FLOAT*  buffer;
    LC3_FLOAT*  sine_table;
    Iisfft      iisfft;
    Cfft        cfft;
} IIS_FFT;

typedef struct T_IIS_FFT* HANDLE_IIS_FFT;

typedef struct {
    Fft      PhEcu_Fft;         /*no counterpart in BASOP */
    Fft      PhEcu_Ifft;        /*no counterpart in BASOP */ 

 
    LC3_FLOAT PhECU_f0hzLtpBin;                      /* BASOP Word16  PhECU_f0hzLtpBinQ7 */
    LC3_FLOAT PhECU_norm_corr;                       /* BASOP Word16  norm_corrQ15 */

    LC3_FLOAT  *PhECU_oold_grp_shape;         /* BASOP Word16  PhECU_oold_grp_shape_fx[MAX_LGW]; */
    LC3_FLOAT  *PhECU_old_grp_shape;          /* BASOP Word16  PhECU_old_grp_shape_fx[MAX_LGW] ; */

    LC3_FLOAT PhECU_L_oold_xfp_w_E;                  /* BASOP  Word32  PhECU_L_oold_xfp_w_E_fx;*/
    LC3_FLOAT PhECU_L_old_xfp_w_E;                   /* BASOP Word32  PhECU_L_old_xfp_w_E_fx; */
    
    LC3_INT32 PhECU_Lprot ;                             /* BASOP Word16  PhECU_Lprot_fx;*/

    LC3_FLOAT *PhECU_xfp;          
    Complex   *PhECU_X_sav_m;      
    LC3_INT32  *PhECU_plocs;    /* BASOP Word16 *PhECU_plocs;    */   /* MAX_PLOCS */ 
    LC3_FLOAT  *PhECU_f0est;         /*BASOP Word32 *PhECU_f0est;*/ 

   LC3_FLOAT *PhECU_mag_chg_1st;     /*  BASOP Word16  PhECU_mag_chg_1st[MAX_LGW];*/
   LC3_FLOAT *PhECU_Xavg;             /*  BASOP Word16  PhECU_Xavg[MAX_LGW] ; */

   LC3_FLOAT PhECU_beta_mute;  /*  BASOP Word16 PhECU_beta_mute*/

   LC3_INT16  PhECU_seed;   /*  BASOP Word16  PhECU_seed_fx;*/

   LC3_INT32  PhECU_LDWIN_OLAP;  /* BASOP Word16  PhECU_LDWIN_OLAP; */
   LC3_INT32  PhECU_t_adv;       /*  BASOP Word16   t_adv; */ 
   
   LC3_INT32  PhECU_short_flag_prev; 
   LC3_INT32  PhECU_time_offs;  
   LC3_INT32  PhECU_num_plocs;
   HANDLE_IIS_FFT handle_fft_phaseecu;
   HANDLE_IIS_FFT handle_ifft_phaseecu;
 
} PlcPhEcuSetup;

typedef struct {
    LC3_INT16 seed;
    LC3_FLOAT gain_c;
    LC3_INT32   lpcorder;
    LC3_FLOAT A[M+1];
    LC3_INT32   fract;
    LC3_INT32   lagw_bw;
    LC3_FLOAT preemphFac;
    LC3_FLOAT *harmonicBuf;
    LC3_FLOAT synthHist[M];
} PlcTdcSetup;

typedef struct {
    LC3_FLOAT *pcmbufHist;
    LC3_INT32   max_len_pcm_plc;
    PlcTdcSetup PlcTdcSetup;
    LC3_FLOAT stabFac;
    LC3_FLOAT cum_fading_slow;
    LC3_FLOAT cum_fading_fast;
    LC3_FLOAT cum_fflcAtten;
    LC3_FLOAT scf_q_old[M];
    LC3_FLOAT scf_q_old_old[M];
    PlcPhEcuSetup PlcPhEcuSetup;   
} PlcAdvSetup;


#endif
