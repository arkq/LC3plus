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

#define BETA_MUTE_FAC 0.5                            /*  % attenuation factor per additional bad frame, FX uses 0.5 (shift right with 1 bit) */
#define BETA_MUTE_FAC_INI 0.5
#define OFF_FRAMES_LIMIT 30                          /*   300 ms for LC3 10 ms   */
#define LGW32k 7
#define LGW16k 5

/*  Tables for attentuation of mag_chg, copied from FX */
/*  Tables are in Q15 */
/*   0.3 dB attenuation per frame for 16 frames, then 6 dB attenuation per frame */
const LC3_INT32 POW_ATT_TABLE1[OFF_FRAMES_LIMIT + 1] = { 32767, 31656, 30581, 29543, 28540, 27571, 26635, 25731, 24857, 24013,
                                                    23198, 22410, 21650, 20915, 20205, 19519,  9759,  4880,  2440,  1220,
                                                      610,   305,   152,    76,    38,    19,    10,     5,     2,     1,
                                                        0 };
/* % 0.4 dB attenuation per frame for 16 frames, then 6 dB attenuation per frame */
const LC3_INT32 POW_ATT_TABLE0[OFF_FRAMES_LIMIT + 1] = { 32767, 31293, 29885, 28540, 27255, 26029, 24857, 23738, 22670, 21650,
                                                    20675, 19745, 18856, 18007, 17197, 16423,  8211,  4106,  2053,  1026,
                                                      513,   257,   128,    64,    32,    16,     8,     4,     2,     1,
                                                        0 };

#ifdef PLC2_FADEOUT_IN_MS 
#if PLC2_FADEOUT_IN_MS == 0
#else 

const LC3_INT32 POW_ATT_TABLE_p3x9_14_7[OFF_FRAMES_LIMIT + 1] = {
 32767,
 31656,  30581,  29543,  28540,  27571,  26635,  25731,  24857,   24013,  /* 9 times .3dB steps , 14  6 dB steps, 7 muted steps */
  12007,   6003,    3002,    1501,    750,    375,     188,     94,     47,      23, 12, 6,      3,      1,
  0,      0,      0,      0,      0,      0,  0 };

const LC3_INT32 POW_ATT_TABLE_p4x9_14_7[OFF_FRAMES_LIMIT + 1] =
{ 32767,
  31293, 29885, 28540, 27255, 26029, 24857, 23738, 22670, 21650,  /* 9 times .4dB steps  , 14  6 dB steps, 7 muted steps */
  10825, 5413,  2706,  1353,  677,   338,  169,  85,  42, 21, 11, 5, 3, 1,
   0, 0,0,0,0,0,0 };

const LC3_INT32 POW_ATT_TABLE_p3x8_6[] = {
 32767,  31656,  30581,  29543,  28540,  27571,  26635,  25731,  12865,   6433,
  3216,   1608,    804,    402,    201,    101,     50,     25,     13,      6,
     3,      2,      1,      0,      0,      0,      0,      0,      0,      0,  0 };
const LC3_INT32  POW_ATT_TABLE_p4x8_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  31293,  29885,  28540,  27255,  26029,  24857,  23738,  11869,   5935,
  2967,   1484,    742,    371,    185,     93,     46,     23,     12,      6,
     3,      1,      1,      0,      0,      0,      0,      0,      0,      0,  0 };

const LC3_INT32  POW_ATT_TABLE_p3x4_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  31656,  30581,  29543,  14772,   7386,   3693,   1847,    923,    462,
   231,    115,     58,     29,     14,      7,      4,      2,      1,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,  0 };
const LC3_INT32  POW_ATT_TABLE_p4x4_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  31293,  29885,  28540,  14270,   7135,   3568,   1784,    892,    446,
   223,    111,     56,     28,     14,      7,      3,      2,      1,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     0 };

const LC3_INT32  POW_ATT_TABLE_p3x2_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  31656,  15828,   7914,   3957,   1979,    989,    495,    247,    124,
    62,     31,     15,      8,      4,      2,      1,      0,      0,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     0 };
const LC3_INT32  POW_ATT_TABLE_p4x2_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  31293,  15647,   7823,   3912,   1956,    978,    489,    244,    122,
    61,     31,     15,      8,      4,      2,      1,      0,      0,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     0 };

const LC3_INT32  POW_ATT_TABLE_p3x1_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  16384,   8192,   4096,   2048,   1024,    512,    256,    128,     64,
    32,     16,      8,      4,      2,      1,      1,      0,      0,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     0 };
const LC3_INT32  POW_ATT_TABLE_p4x1_6[OFF_FRAMES_LIMIT + 1] = {
 32767,  16384,   8192,   4096,   2048,   1024,    512,    256,    128,     64,
    32,     16,      8,      4,      2,      1,      1,      0,      0,      0,
     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,     0 };

const LC3_INT32 *const POW_ATT_TABLES[1 + 12] =
{ NULL,
/*0.3dB col      ,      0.4dB col */
/* 1*/POW_ATT_TABLE_p3x1_6, POW_ATT_TABLE_p4x1_6,       /* 0 0.3dB,  16 6dB, 14mute   */ /*  0.4dB version */  /* old short mute tabs  */
/* 3*/POW_ATT_TABLE_p3x2_6, POW_ATT_TABLE_p4x2_6,       /* 1 0.3dB,  15 6dB ,14mute   */ /* 0.4dB version  */
/* 5*/POW_ATT_TABLE_p3x4_6, POW_ATT_TABLE_p4x4_6,       /* 3 0.3dB,  15 6dB , 12 mute */ /* 0.4dB version  */
/* 7*/POW_ATT_TABLE_p3x8_6, POW_ATT_TABLE_p4x8_6,       /* 7 0.3dB,  15 6dB , 8 mute  */ /* 0.4dB version  */
/* 9*/POW_ATT_TABLE_p3x9_14_7, POW_ATT_TABLE_p4x9_14_7, /* 9 0.3dB,  14 6dB , 7 mute  */ /* 0.4dB version  */   /* opt 120 ms */
/*11*/POW_ATT_TABLE1,       POW_ATT_TABLE0,             /* 15 0.3dB, 14 6dB , 1 mute */  /* 0.4dB version  */   /*  original   curves */
};

#endif
#endif

void plc_phEcu_tba_trans_dect_gains(LC3_INT32 burst_len, LC3_INT32 n_grp, LC3_FLOAT *grp_pow_change,  
                          LC3_FLOAT *stPhECU_beta_mute, LC3_FLOAT *stPhECU_mag_chg_1st,
                          LC3_FLOAT *alpha, LC3_FLOAT *beta, LC3_FLOAT *mag_chg, LC3_FLOAT *ph_dith, LC3_INT32 *tr_dec,
                          LC3_FLOAT *att_val, LC3_INT32 *attDegreeFrames_dbg, LC3_FLOAT *thresh_dbg                         
	                      , LC3_UINT8 plc_fadeout_type
                          ) 
{

    LC3_INT32 i;
    LC3_FLOAT  thresh_tr_dB, max_increase_grp_pow;
    LC3_FLOAT max_increase_grp_pow_lin;
    LC3_FLOAT grp_pow_change_lin[MAX_LGW];
    LC3_FLOAT XavgFadeinFactor;

    LC3_INT32 burst_att_thresh;
    LC3_INT32 att_per_frame_idx;
    LC3_INT32 att_always,   attDegreeFrames;
    const LC3_INT32 *TABLEQ15;
    LC3_INT32 beta_mute_thr;                             /* time threshold in 10 ms frames  to start beta - noise attenuation */
    UNUSED(attDegreeFrames_dbg);

    /* constants setup */
    att_always = 0;  

    XavgFadeinFactor = -1.0;

    /* 10ms constants */
    thresh_tr_dB = 10.0;            /*    dB threshold kept same as for 20ms, even though transient analysis frame size was shortened */
    max_increase_grp_pow = 0;      /*     maximum amplification(dB)  in case of onset transients, offset always deacy */

    max_increase_grp_pow_lin = (LC3_FLOAT)1.0*LC3_POW((LC3_FLOAT)10.0, max_increase_grp_pow / (LC3_FLOAT)10.0)*(LC3_FLOAT)(32767.0 / 32768.0);

    if (plc_fadeout_type != 0)
    {
        i = (PLC2_FADEOUT_LONG_IN_MS - PLC2_FADEOUT_IN_MS_MIN) / PLC2_FADEOUT_RES;   /*a long fading table entry in fade_scheme_tab  */
    } else {
        i = (PLC2_FADEOUT_IN_MS - PLC2_FADEOUT_IN_MS_MIN) / PLC2_FADEOUT_RES;   /* a shorter  fading entry in  fade_scheme_tab  */
    }
    assert(i >= 0 && i <= ((PLC2_FADEOUT_IN_MS_MAX - PLC2_FADEOUT_IN_MS_MIN) / PLC2_FADEOUT_RES) && "fade_scheme_tab index error");

    att_per_frame_idx = fade_scheme_tab[i][0];
    burst_att_thresh  = fade_scheme_tab[i][1];           /* number of 1.0 frames before  muting/mixing phase */
    /* band gain muting may can take place earlier due to a band transient */
    beta_mute_thr = fade_scheme_tab[i][2];          /*   muting of Xavg contribution  start when slow fadeout is over */
 
    attDegreeFrames = 0;
    if (burst_len > burst_att_thresh)
    {
        att_always = 1;

        /*   Added to be able to able to use tables to be aligned with FX */
        /*   Limit attDegreeFrames to OFF_FRAMES_LIMIT */
        attDegreeFrames = burst_len - burst_att_thresh;

        if (attDegreeFrames > OFF_FRAMES_LIMIT)
        {
            attDegreeFrames = OFF_FRAMES_LIMIT;
        }
    }


    set_vec(1.0 * (32767.0/32768.0), mag_chg, n_grp);
    set_vec(0.0, ph_dith, n_grp);

    set_vec(1.0 * (32767.0/32768.0), alpha, n_grp);
    set_vec(0.0, beta, n_grp);
    set_vec_int(0, tr_dec, n_grp);

    set_vec(1.0 * (32767.0/32768.0), att_val, n_grp);

 

    /* transient detection per band  */
    for (i = 0;i < n_grp; i++) {
        if(burst_len == 1)
        {
            /* first bad frame  */
            grp_pow_change_lin[i] = LC3_POW((LC3_FLOAT)10.0, grp_pow_change[i]/(LC3_FLOAT)10.0);

            *stPhECU_beta_mute = BETA_MUTE_FAC_INI;
            *stPhECU_beta_mute = *stPhECU_beta_mute / (LC3_FLOAT)2.0;
             
            /*  transient processing */
            /*   transients may be both rise and decay transients !! */

            if(LC3_FABS(grp_pow_change[i]) >= thresh_tr_dB)
            {
 
                tr_dec[i] = 1;
            }


            /*   magnitude modification */
            att_val[i] = 0.0f;
            if(tr_dec[i] || att_always) {

                att_val[i] = MIN(max_increase_grp_pow_lin, grp_pow_change_lin[i]);  /*  % linear values !! */
                att_val[i] = LC3_SQRT(att_val[i]);
                mag_chg[i] = att_val[i];
                stPhECU_mag_chg_1st[i] = att_val[i];
            }
            else
            {
                mag_chg[i] = 1.0 * (LC3_FLOAT)(32767.0/32768.0);
                stPhECU_mag_chg_1st[i] = (LC3_FLOAT)1.0;
            }
        }
        else
        {
            /*   burst handling based on states */

            assert(burst_len >= 2);      /*  states used here */
            tr_dec[i] = 0;

            {
                assert(att_per_frame_idx >= 1 && att_per_frame_idx <= (10+2));
                TABLEQ15 = POW_ATT_TABLES[att_per_frame_idx]; 
                att_val[i] = (LC3_FLOAT)1.0 * ( (LC3_FLOAT) TABLEQ15[MIN(OFF_FRAMES_LIMIT, attDegreeFrames )]  / (LC3_FLOAT)32768.0);  /* Table idx 0...N-1 therefore no + 1 */
                att_val[i] = att_val[i];
            }
                        
            if ( (att_val[i] != 0) && (att_val[i] * (LC3_FLOAT)32768.0 < (LC3_FLOAT)0.5) )
            {
                att_val[i] = 0.0;            /*   for SNR measurments match in  float lowest possible level to BASOP representation  */
            }

            /*  Apply attenuation */
            mag_chg[i] = stPhECU_mag_chg_1st[i];

            mag_chg[i] = mag_chg[i] * att_val[i]; /*   add additional attenuation from burst attenation  logic */

            if ((mag_chg[i] != 0) && (mag_chg[i] * (LC3_FLOAT)32768.0 < (LC3_FLOAT)0.5))
            {
                mag_chg[i] = 0;             /*   for SNR measurments match in  float lowest possible level to BASOP representation */
            }

            /* note beta_mute decreased once per frame,  not once per band  */
            if (i == 0 && burst_len > beta_mute_thr)
            {
                    *stPhECU_beta_mute = *stPhECU_beta_mute * (LC3_FLOAT)BETA_MUTE_FAC;
            }

            alpha[i] = mag_chg[i];

            if (alpha[i] >= (LC3_FLOAT)(32766.0 / 32768.0))
            {
                beta[i] = 0;  /*  align to BASOP more efficent  use of beta */
            }
            else
            {
                beta[i] = LC3_SQRT((LC3_FLOAT)1.0 - alpha[i]* alpha[i]) * *stPhECU_beta_mute;
            }

            if ( i >= LGW32k-1) {
                beta[i] = beta[i] * (LC3_FLOAT)0.1;
            }
            else if( i >= LGW16k-1)
            {
                beta[i] = beta[i] * (LC3_FLOAT)0.5;
            }


            /*  limit Xavg noise contribution further in case of offset / tr_decay */

            if ((burst_len <= burst_att_thresh) && (stPhECU_mag_chg_1st[i] < (LC3_FLOAT)(32767.0 / 32768.0)))
            {
               XavgFadeinFactor = (LC3_FLOAT)(burst_len - (LC3_FLOAT)1.0) / burst_att_thresh;
        
               XavgFadeinFactor = MIN((LC3_FLOAT)1.0, XavgFadeinFactor);

               beta[i] = beta[i] * XavgFadeinFactor;
  
            }
        }
    }

    if (thresh_dbg != NULL)
    {
        *thresh_dbg = XavgFadeinFactor;
    }

     return;
}

