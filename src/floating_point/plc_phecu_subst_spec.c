/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#include "defines.h"
#include "functions.h"
#include "constants.h"


static LC3_INT32 own_rand(LC3_INT32 seed);
static Complex valley_magnitude_adj(Complex X_i_in, LC3_INT32 uni_seed, LC3_FLOAT cos_F);
static LC3_INT32 rand_phase(LC3_INT32 seed_in, LC3_FLOAT* cos_F);

void plc_phEcu_subst_spec(LC3_INT32* plocs, LC3_INT32 n_plocs, LC3_FLOAT* f0est, LC3_INT32 time_offs, Complex* X, LC3_INT32 X_len,
                          LC3_FLOAT* mag_chg_gr, LC3_INT32 *seed, LC3_FLOAT* alpha, LC3_FLOAT* beta, LC3_FLOAT* Xavg,
                          LC3_INT32 t_adv_in, LC3_INT32 Lprot, LC3_INT32 delta_corr, LC3_FLOAT *corr_phase_dbg,
                          LC3_FLOAT *X_i_new_re_dbg, LC3_FLOAT *X_i_new_im_dbg) {
    
    LC3_INT32 i, i2, lprotBy2Minus1, one_peak_flag_mask, noise_mag_scale;
    LC3_INT32 t_adv;
   LC3_FLOAT corr_phase[MAX_PLC_NPLOCS] = {0};
   LC3_FLOAT cos_F, mag_chg_local, alpha_local, beta_local, tmp;
    Complex X_i, X_i_new;
    LC3_INT32 segmentLen, e;
    LC3_FLOAT Xph;
    LC3_FLOAT seed_local;
    
    UNUSED(corr_phase_dbg);
    UNUSED(X_i_new_re_dbg);
    UNUSED(X_i_new_im_dbg);

    seed_local = (LC3_FLOAT) *seed;
   
    
   lprotBy2Minus1 = imin(320, Lprot/2 - 1); /* limit to 20 KHz */
    
    
    t_adv = t_adv_in + time_offs;
    
    for (i = 0; i < n_plocs; i++) {
        corr_phase[i] = (LC3_FLOAT)2.0 * (LC3_FLOAT)M_PI * (f0est[i]/Lprot)*(LC3_FLOAT)t_adv;
    }


    // EVOLVE PHASE -----------------
    LC3_INT32 binCounter = 1;
    LC3_INT32 subInd = 0;
    
    one_peak_flag_mask = -1;
    if (n_plocs < 3 && n_plocs > 0) {
        one_peak_flag_mask = 0;
    }
    
    noise_mag_scale = 0;
    if (n_plocs == 0 || time_offs != 0) {
        noise_mag_scale = 1;
    }
    
    if (n_plocs == 0) {
        X[0] = realtoc(0);
        X[X_len-1] = realtoc(0);
    }


   if (n_plocs != 0) {
      for (i = 0; i < n_plocs; i++) {
         LC3_INT32 delta_corr_dn = delta_corr;
         LC3_INT32 delta_corr_up = delta_corr;

         if (i > 0) {
            delta_corr_dn = imin( ((plocs[i] - plocs[i - 1] - 1) / 2), delta_corr_dn);
         }

         if (i < n_plocs - 1) {
            delta_corr_up = imin( ((plocs[i + 1] - plocs[i] - 1) / 2), delta_corr_up);
         }

         segmentLen = (plocs[i] - delta_corr_dn) - binCounter;

         for (i2 = 0; i2 < segmentLen; i2++) {
            seed_local = (LC3_FLOAT)rand_phase((LC3_INT32)seed_local, &cos_F);
           
            X_i = X[binCounter];
            X_i_new = cmul(X_i, cexpi((LC3_FLOAT)M_PI*seed_local / (LC3_FLOAT)32768.0));


            seed_local = (LC3_FLOAT)own_rand((LC3_INT32)seed_local);

            if (noise_mag_scale != 0) {
               X_i = valley_magnitude_adj(X_i_new,(LC3_INT32) seed_local, cos_F);
               X_i_new = X_i;
            }

            mag_chg_local = mag_chg_gr[subInd];
            alpha_local = alpha[subInd];

            if (beta[subInd] != 0) {
               tmp = beta[subInd] * Xavg[subInd];
               if (one_peak_flag_mask == 0) {
                  tmp = 0;
                  X_i_new = realtoc(0);
               }
               X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI*seed_local / (LC3_FLOAT)32768.0)));
            }
            else {
               if (one_peak_flag_mask == 0) {
                  X_i_new = realtoc(0);
               }

               X[binCounter] = cmul(realtoc(mag_chg_local), X_i_new);
            }

            binCounter++;

            if (binCounter >= gwlpr[subInd + 1]) {
               subInd++;
            }
         }
            
            e = plocs[i] + delta_corr_up;
            if (e > lprotBy2Minus1) {
                e = lprotBy2Minus1;
            }
            
            Xph = corr_phase[i];
             
         segmentLen = e - (binCounter - 1);

         for (i2 = 0; i2 < segmentLen; i2++) 
         {
            seed_local = (LC3_FLOAT)own_rand((LC3_INT32)seed_local);
            X_i = X[binCounter];

            {
               LC3_INT32 nrep =(LC3_INT32) LC3_FLOOR(Xph / (2.0f*(LC3_FLOAT)M_PI));

               X_i_new = cmul(X_i, cexpi(Xph - (2.0f*(LC3_FLOAT)M_PI*(LC3_FLOAT)nrep)));
            }


            seed_local = (LC3_FLOAT)own_rand((LC3_INT32)seed_local);

            mag_chg_local = mag_chg_gr[subInd];
            alpha_local = alpha[subInd];
            beta_local = beta[subInd];
            if (beta_local != 0) {

               assert(alpha_local == mag_chg_local);
               tmp = beta_local * Xavg[subInd];

               X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI*seed_local / (LC3_FLOAT)32768.0)));
            }
            else
            {
               X[binCounter] = cmul(realtoc(mag_chg_local), X_i_new);
            }

            binCounter++;

            if (binCounter >= gwlpr[subInd + 1]) {
               subInd++;
            }
         }
      }
   }

   segmentLen = lprotBy2Minus1 - (binCounter - 1);

   for (i = 0; i < segmentLen; i++) {
      seed_local = (LC3_FLOAT)rand_phase((LC3_INT32)seed_local, &cos_F);
        
        X_i = X[binCounter];
        X_i_new = cmul(X_i, cexpi((LC3_FLOAT)M_PI*seed_local/(LC3_FLOAT)32768.0));

        seed_local = (LC3_FLOAT)own_rand((LC3_INT32)seed_local);
        
        if (noise_mag_scale != 0) {
            X_i = valley_magnitude_adj(X_i_new, (LC3_INT32)seed_local, cos_F);
            X_i_new = X_i;
        }
        
        if (one_peak_flag_mask == 0) {
            X_i_new = realtoc(0);
        }
        
        alpha_local = alpha[subInd];
        mag_chg_local = mag_chg_gr[subInd];
        
        if (beta[subInd] != 0) {
            assert(alpha_local == mag_chg_local);
            tmp = beta[subInd]*Xavg[subInd];
            
            if (one_peak_flag_mask == 0) {
                tmp = 0;
            }
            
            X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI*seed_local/(LC3_FLOAT)32768.0)));
        } 
      else 
      {
            X[binCounter] = cmul(realtoc(mag_chg_local), X_i_new);
        }
        
        binCounter++;
        
        if (binCounter >= gwlpr[subInd + 1]) {
            subInd++;
        }
    }

    
   *seed = (LC3_INT32)seed_local;
}

static LC3_INT32 own_rand(LC3_INT32 seed) {
     assert(seed <= 32767 && seed >= -32768);
    LC3_INT32 retSeed = (13849 + (seed + 32768) * 31821) & 65535;
    retSeed -= 32768;
   assert(retSeed <= 32767 && retSeed >= -32768);
    return retSeed;
}

static Complex valley_magnitude_adj(Complex X_i_in, LC3_INT32 uni_seed, LC3_FLOAT cos_F) {
    LC3_FLOAT scale = ((LC3_FLOAT)0.5*(LC3_FLOAT)uni_seed/(LC3_FLOAT)32768.0) + (LC3_FLOAT)0.5*cos_F;
    scale = (LC3_FLOAT)1.0 + (LC3_FLOAT)0.25*scale;
    
    assert(scale <= (LC3_FLOAT)1.25);
    assert(scale >= (LC3_FLOAT)0.75);
    
    return cmul(X_i_in, realtoc(scale));
}

static LC3_INT32 rand_phase(LC3_INT32 seed_in, LC3_FLOAT* cos_F) {
    LC3_FLOAT seed = (LC3_FLOAT)own_rand(seed_in);
    *cos_F = LC3_COS((LC3_FLOAT)M_PI*seed/(LC3_FLOAT)32768.0);
    return (LC3_INT32) seed;
}

