/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void processPlcUpdate_fl(PlcAdvSetup *PlcAdvSetup, LC3_INT32 frame_length, LC3_FLOAT *syntM, LC3_FLOAT *scf_q,
             LC3_INT32 *nbLostCmpt, LC3_FLOAT *cum_alpha, LC3_INT32 bfi, LC3_INT32 *prevBfi, LC3_INT32 *prevprevBfi)
{
    LC3_FLOAT tmp[MAX_LEN_PCM_PLC];
    
    move_float(tmp, &PlcAdvSetup->pcmbufHist[frame_length], PlcAdvSetup->max_len_pcm_plc - frame_length);
    move_float(&PlcAdvSetup->pcmbufHist[0], tmp, PlcAdvSetup->max_len_pcm_plc - frame_length);
    move_float(&PlcAdvSetup->pcmbufHist[PlcAdvSetup->max_len_pcm_plc - frame_length], syntM, frame_length);
    
    if (bfi != 1)
    {   
        *nbLostCmpt = 0;
        *cum_alpha = 1;
        
        if (PlcAdvSetup)
        {
            move_float(PlcAdvSetup->scf_q_old_old, PlcAdvSetup->scf_q_old, M);
            move_float(PlcAdvSetup->scf_q_old, scf_q, M);
             /* PLC fullband transient detector setting for non-bfi frames */
             PlcAdvSetup->PlcPhEcuSetup.PhECU_short_flag_prev  = 0;    /* fullband transient not active   */
        }
    }
    
    *prevprevBfi = *prevBfi;
    *prevBfi = bfi;
}

void plc_phEcu_processPLCspec2shape(LC3_INT16 prev_bfi, LC3_INT16 bfi, LC3_FLOAT q_d[], LC3_INT32 yLen,
   LC3_FLOAT *stPhECU_oold_grp_shape, LC3_FLOAT *stPhECU_old_grp_shape)
{
   LC3_INT32 i, j, N_grp;
   LC3_INT32 local_prev_bfi;
   LC3_INT32 fs_idx;
   LC3_FLOAT E_tot = 0.0;
   LC3_INT32 l_grp;
   LC3_FLOAT *pX;

   if (bfi != 1)  /* compute only for  bfi== 0 or 2 */
   {
      fs_idx = (LC3_INT32)floor(yLen / 100);
      assert(fs_idx < 5);
      N_grp = xavg_N_grp[fs_idx];
     
      local_prev_bfi = prev_bfi;
      if (local_prev_bfi == 2) {
         local_prev_bfi = 0;
      }


      /* Copy  old to oold grp shape */
      for (i = 0; i < MAX_LGW; i++)
      {
         stPhECU_oold_grp_shape[i] = stPhECU_old_grp_shape[i];
      }

     /* Accumulate DC-coupled  bins to total */
      E_tot = 0;
      pX = q_d;          /*  ptr setup */
      for (i = 0; i < mdct_grp_bins[0]; i++) 
      {
          E_tot +=  sqrf( *pX ); 
          pX++;          
      }

      /* Accumulate middle grps and add to total */
      for (i = 0; i < (N_grp - 1); i++)
      {
         l_grp = mdct_grp_bins[i + 1] - mdct_grp_bins[i]; ;
         stPhECU_old_grp_shape[i] = 0.0;
         for (j = 0; j < l_grp; j++) {
            stPhECU_old_grp_shape[i] +=  sqrf( *pX ); 
            pX++;
         }
         E_tot += stPhECU_old_grp_shape[i];
      }

      /* Accumulate last  subbband and add to total */
      stPhECU_old_grp_shape[(N_grp - 1)] = 0.0;
      l_grp = mdct_grp_bins[N_grp] - mdct_grp_bins[N_grp - 1] - mdct_grp_bins[0];
      assert( (mdct_grp_bins[N_grp] - mdct_grp_bins[0]) <= yLen);
      for (j = 0; j < l_grp; j++)
      {
         stPhECU_old_grp_shape[(N_grp - 1)] += sqrf( *pX );
         pX++;
      }
      E_tot += stPhECU_old_grp_shape[(N_grp - 1)];


      /* Normalize shape */
      for (i = 0; i < (N_grp); i++) {
         if (E_tot > 0.0) {
            stPhECU_old_grp_shape[i] /= E_tot;
         }
         else 
         {
            stPhECU_old_grp_shape[i] = 0.0;
         }
      }
      if (local_prev_bfi == 1) {
         for (i = 0; i < MAX_LGW; i++) {
            stPhECU_oold_grp_shape[i] = stPhECU_old_grp_shape[i];
         }
      }
   }/*bfi*/
   return;
}

void processPlcUpdateSpec_fl(LC3_FLOAT *q_d_prev, LC3_FLOAT *q_d_fl_c, LC3_INT32 yLen)
{
    move_float(q_d_prev, q_d_fl_c, yLen);
}

