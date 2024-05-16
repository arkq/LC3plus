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
#include "constants.h"

static LC3_INT32 own_rand(LC3_INT32 seed);
static Complex valley_magnitude_adj(Complex X_i_in, LC3_INT32 uni_seed, LC3_FLOAT cos_F);
static LC3_INT32 rand_phase(LC3_INT32 seed_in, LC3_FLOAT* cos_F);

#define ONE_SIDED_SINE_WIDTH    (4)       /* expected pure sine main lobe   maximum width (4+1+4) bins *62.5 hz/bin  => approx 560 Hz total width    */  

static LC3_INT16  plc_phEcu_nonpure_tone_ana(const LC3_INT32* plocs, const LC3_INT32 n_plocs, const  Complex* X,    const  LC3_FLOAT* Xavg, const  LC3_INT32 Lprot);

void plc_phEcu_subst_spec(LC3_INT32* plocs, LC3_INT32 n_plocs, LC3_FLOAT* f0est, LC3_INT32 time_offs, Complex* X, LC3_INT32 X_len,
                          LC3_FLOAT* mag_chg_gr, LC3_INT32 *seed, LC3_FLOAT* alpha, LC3_FLOAT* beta, LC3_FLOAT* Xavg,
                          LC3_INT32 t_adv_in, LC3_INT32 Lprot, LC3_INT32 delta_corr,
                           LC3_INT16    fadeout,            /* need for DC muting */
                           LC3_INT16* nonpure_tone_flag_ptr, 
                          LC3_FLOAT *corr_phase_dbg,
                          LC3_FLOAT *X_i_new_re_dbg, LC3_FLOAT *X_i_new_im_dbg) {
    
    LC3_INT32 i, i2, lprotBy2Minus1, one_peak_flag_mask, noise_mag_scale;
    LC3_INT32 t_adv;
    LC3_FLOAT corr_phase[MAX_PLC_NPLOCS] = {0};
    LC3_FLOAT cos_F, mag_chg_local, alpha_local, beta_local, tmp;
    Complex X_i, X_i_new;
    LC3_INT32 segmentLen, e;
    LC3_FLOAT Xph;
    LC3_FLOAT seed_local;
    LC3_INT32 binCounter = 1, subInd = 0;
	LC3_INT16 fs_idx;

    UNUSED(corr_phase_dbg);
    UNUSED(X_i_new_re_dbg);
    UNUSED(X_i_new_im_dbg);

    seed_local = (LC3_FLOAT) *seed;
   
    lprotBy2Minus1 = imin(320, Lprot/2 - 1); /* limit to 20 KHz */
    
    t_adv = t_adv_in + time_offs;
    
    for (i = 0; i < n_plocs; i++) {
        corr_phase[i] = (LC3_FLOAT)2.0 * (LC3_FLOAT)M_PI_LC3PLUS * (f0est[i]/Lprot)*(LC3_FLOAT)t_adv;
    }

    // EVOLVE PHASE -----------------
    
    one_peak_flag_mask = -1;	
    fs_idx = (LC3_INT16)LC3_FLOOR((LC3_FLOAT)Lprot / 160.0);  /* aquire, fs_idx for 10 ms frame sizes */
	if (n_plocs < 3 && n_plocs > 0) 
	{
		one_peak_flag_mask = 0;  /* initial crude single tone detection, only using n_plocs as a  result from peak_locator() dynamics as input  */

		if ( (*nonpure_tone_flag_ptr < 0 )
				 &&  ( (fs_idx == 2) /*SemiSWB 24 kHz */ ||  (fs_idx >= 4) /* FB 48 kHz */ )  
			)
		{
 	        /* in the first lost frame  analyze spectra to possibly reverse initial pure sine assumption  */
			*nonpure_tone_flag_ptr  = plc_phEcu_nonpure_tone_ana(plocs,  n_plocs,  X,   Xavg,   Lprot );
		} 
	    
		if (  *nonpure_tone_flag_ptr  >  0 ) {
			one_peak_flag_mask = -1; /* actually  revert single pure tone detection  */  /* 0->  mute all surrounding valley bins in evolution ,    0xff  -> generate noise in all valleys */
		}

	}  
    
    noise_mag_scale = 0;
    if (n_plocs == 0 || time_offs != 0) {
        noise_mag_scale = 1;
    }
    
    if (n_plocs == 0) {
        X[0] = realtoc(0);
        X[X_len-1] = realtoc(0);
    }

    /*  binary selection of fadeout scheme  */   
    assert(PLC2_FADEOUT_LONG_IN_MS >= PLC2_FADEOUT_IN_MS_MIN && PLC2_FADEOUT_IN_MS >= PLC2_FADEOUT_IN_MS_MIN);
    assert(PLC2_FADEOUT_LONG_IN_MS <= PLC2_FADEOUT_IN_MS_MAX && PLC2_FADEOUT_IN_MS <= PLC2_FADEOUT_IN_MS_MAX); 
    i = (PLC2_FADEOUT_IN_MS - PLC2_FADEOUT_IN_MS_MIN) / PLC2_FADEOUT_RES;

    if (fadeout != 0)
    {
        i = (PLC2_FADEOUT_LONG_IN_MS - PLC2_FADEOUT_IN_MS_MIN) / PLC2_FADEOUT_RES;
    }

    /* calculate  local burst_len for  securing DC  and fs/2 muting */
    i2 = (time_offs / ((Lprot * 100) / 160)) + 1; /* burst_len */

    if (i2 > (fade_scheme_tab[i][1] + 1))
    {
        /*     start DC  scaling attenuation  */
        X[0].r = alpha[0] * X[0].r;

        /*  start fs/by2   attenuation  */
        X[X_len - 1].r = alpha[(xavg_N_grp[fs_idx] - 1)] * X[X_len - 1].r;
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
            X_i_new = cmul(X_i, cexpi((LC3_FLOAT)M_PI_LC3PLUS*seed_local / (LC3_FLOAT)32768.0));


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
               X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI_LC3PLUS*seed_local / (LC3_FLOAT)32768.0)));
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
               LC3_INT32 nrep =(LC3_INT32) LC3_FLOOR(Xph / (2.0f*(LC3_FLOAT)M_PI_LC3PLUS));

               X_i_new = cmul(X_i, cexpi(Xph - (2.0f*(LC3_FLOAT)M_PI_LC3PLUS*(LC3_FLOAT)nrep)));
            }


            seed_local = (LC3_FLOAT)own_rand((LC3_INT32)seed_local);

            mag_chg_local = mag_chg_gr[subInd];
            alpha_local = alpha[subInd];
            beta_local = beta[subInd];
            if (beta_local != 0) {

               assert(alpha_local == mag_chg_local);
               tmp = beta_local * Xavg[subInd];

               X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI_LC3PLUS*seed_local / (LC3_FLOAT)32768.0)));
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
        X_i_new = cmul(X_i, cexpi((LC3_FLOAT)M_PI_LC3PLUS*seed_local/(LC3_FLOAT)32768.0));

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
            
            X[binCounter] = cadd(cmul(realtoc(alpha_local), X_i_new), cmul(realtoc(tmp), cexpi((LC3_FLOAT)M_PI_LC3PLUS*seed_local/(LC3_FLOAT)32768.0)));
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
    LC3_INT32 retSeed;
    assert(seed <= 32767 && seed >= -32768);
    retSeed = (13849 + (seed + 32768) * 31821) & 65535;
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
    *cos_F = LC3_COS((LC3_FLOAT)M_PI_LC3PLUS*seed/(LC3_FLOAT)32768.0);
    return (LC3_INT32) seed;
}

static LC3_INT16  plc_phEcu_nonpure_tone_ana(const LC3_INT32* plocs, const LC3_INT32 n_plocs, const  Complex* X,   const  LC3_FLOAT* Xavg, const  LC3_INT32 Lprot)
{
 
    LC3_INT16 nonpure_tone_detect;
	LC3_INT16  n_ind, tone_ind, low_ind, high_ind;
	LC3_FLOAT  peak_amp, peak_amp2, valley_amp, x_abs[(1 + 2 * ONE_SIDED_SINE_WIDTH + 2 * 1)];
	LC3_INT16  sineband_ind_low, sineband_ind_high;
	LC3_INT16 i, fs_idx, N_grp;
	LC3_FLOAT tmp, tmp_dB, tot_inc_HF, tot_inc_LF;
	 


	/* use compressed hearing sensitivity curve to allow more deviation in highest and lowest bands */
	/* ROM table 	LC3_FLOAT scATHFx[MAX_LGW - 1]  */

	/* init */
	nonpure_tone_detect = 0;
	tot_inc_HF = 0.0;
	tot_inc_LF = 0.0;

	/*  limit single sine  optimization to when 2 peaks are close enough to represent a single sinusoid */
	if (n_plocs == 2 && (plocs[1] - plocs[0]) >= ONE_SIDED_SINE_WIDTH) /* NB, plocs is an ordered vector */
	{
		nonpure_tone_detect |= 0x1;
	}

	/* local bin wise dynamics analysis, if 2 peaks, we do the analysis based on the location of the largest peak */
	{
		tone_ind = 0;
		plc_phEcu_fft_spec2_sqrt_approx(&(X[plocs[0]]), 1, &peak_amp);  /* get 1st peak amplitude = approx_sqrt(Re^2+Im^2) */


		if ((n_plocs - 2) == 0)
		{
			plc_phEcu_fft_spec2_sqrt_approx(&(X[plocs[1]]), 1, &peak_amp2);  /* get 2nd peak amplitude */
			if (peak_amp2 > peak_amp)
			{
				tone_ind = 1;
				peak_amp = peak_amp2;
			}
		}

		low_ind = MAX(1, plocs[tone_ind] - (ONE_SIDED_SINE_WIDTH + 1));                /* DC is not allowed as valley   */
		high_ind = MIN((Lprot >> 1) - 2, plocs[tone_ind] + (ONE_SIDED_SINE_WIDTH + 1)); /* Fs/2 is not allowed as valley  */

		n_ind = high_ind - low_ind + 1;

		/*  find lowest amplitudes around the  assumed  main  lobe center location */
		plc_phEcu_fft_spec2_sqrt_approx(&(X[low_ind]), n_ind, x_abs);
		valley_amp = peak_amp;
		for (i = 0; i < n_ind; i++) {
			valley_amp = MIN(x_abs[i], valley_amp);
		}

		/* at least  a localized amplitude ratio of 16 (24 dB)  required to declare a pure sinusoid  */
		if (peak_amp < 16 * valley_amp)  /* 1/16 easily implemented in BASOP */
		{

			nonpure_tone_detect |= 0x2;/* not a pure tone due to too low local SNR */

		}
	}
	
	/* analyze LF/ HF bands  energy dynamics vs the  assumed single tone band ( one or two  peaks found) */
	{
		fs_idx = (LC3_INT16)floor(Lprot / 160); /* fs_idx */
		assert(fs_idx < 5);

		/* Xavg , is a vector of rather rough MDCT based band energy estimates in perceptually motivated bands. from approx the last 26 ms of synthesis */

		/*  eval amplitude relations  for assumed tonal  band vs  lower and higher bands */
		N_grp = xavg_N_grp[fs_idx];  /* { 4 NB , 5 WB , 6 SSWB , 7 SWB, 8  FB };  */

	  /* establish band(s) with assumed sinusoid tone */
	  /* if tone freq location is below first  MDCT-band definition,  use first band as location anyway */
		i = 0;                                     /* band                 0   ,    1     ,     2      , 3      , ...*/
		while (plocs[tone_ind] >= gwlpr[i + 1]) {  /*  gwplr=       [      1, 12(750Hz), 20(1250Hz) , 36     , .. */
		                                          /*  dct-inds           "0"...11, 12...19,   20...35,     36 ...     */
			i++;
		}
		sineband_ind_low  = i;
		sineband_ind_high = i; /* typically  in the same band as low  */

		/* a single tone may end up on a band border
		  , handle case when assumed tone is more or less right in between two perceptual bands  +/-4  62.5  Hz */
		if ((sineband_ind_high > 0) &&
			(plocs[tone_ind] - ONE_SIDED_SINE_WIDTH) >= gwlpr[sineband_ind_high + 1]
			) {
			sineband_ind_low = sineband_ind_high - 1;
		}

		if ( (sineband_ind_low < (N_grp - 1)) &&
			 (plocs[tone_ind] + ONE_SIDED_SINE_WIDTH) >= gwlpr[sineband_ind_low + 1]
		    ) {
			sineband_ind_high = sineband_ind_low + 1;
		}
	}



	/* intraframe(26 ms),    weighted LB and HB envelope dynamics/variation  analysis  */
	 /* envelope analysis ,
		require at least two HF or two LF  bands in the envelope taper/roll-off  analysis ,  otherwise skip this  condition  */


	if (nonpure_tone_detect == 0 &&
		(((sineband_ind_high + 2) < N_grp) ||
		((sineband_ind_low - 2) >= 1)
			)
		)
	{
		/* delta taper-off analysis solution, less sensitive to input bandwidth limitation and levels   */

		   /* verify that an assumed  clean sine does not have any odd HF content indications by thresholding the accumulated delta rise in  LF/HF side lobes   */
		for (i = (sineband_ind_high + 1); i < (N_grp - 1); i++) {
			tmp = (Xavg[i + 1] + LC3_EPS) / (Xavg[i] + LC3_EPS);
			tmp_dB = 20.0*LC3_LOGTEN(tmp);
			if ((Xavg[i] + LC3_EPS) > (Xavg[i + 1] + LC3_EPS)) {
				tmp_dB = 0;
			}
			tot_inc_HF += scATHFx[i] * tmp_dB; /*  i is ATH factor between band i, i+1  based on Hearing sensitivity */
		}

		/* verify that an assumed  clean sine does not have any odd LF  content by thresholding the accumulated LF reverse up tilt  */
		for (i = MAX(0, (sineband_ind_low - 1)); i > 0; i--) {
			tmp = (Xavg[i - 1] + LC3_EPS) / (Xavg[i] + LC3_EPS);
			tmp_dB = 20.0*LC3_LOGTEN(tmp);  /* switch to log2() to simplify BASOP  */

			if ((Xavg[i - 1] + LC3_EPS) < (Xavg[i] + LC3_EPS)) {
				tmp_dB = 0;
			}
			tot_inc_LF += scATHFx[i - 1] * tmp_dB;  /* "psycho" scale using  i-1  is ATH factor between band i-1, i , based on Hearing sensitivity */
		}

		if (tot_inc_HF > 4.5){  /* 4.5 dB in log2 is   0.7474   */ 
			nonpure_tone_detect |= 0x10; /* still not a pure tone,  too  great  HF  side increase  */
		}

		if (tot_inc_LF > 4.5)	{ /* 4.5  dB limit in 4.5 = 20log10(x)  corresponds to  limit value 0.7474  in log2(x)  */
			nonpure_tone_detect |= 0x20; /* still not a pure tone,  too  great accumulated   LF  side increase  */
		}

		/* verify that an assumed  clean sine does not have any odd LF+HF  content by thresholding the accumulated LF+HF unexpected  tilt  */
		if ((tot_inc_LF + tot_inc_HF) > 6.0)	{  /* 6  dB limit in 20log10(x)  corresponds to  limit value 1.0 in log2(x)  */
			nonpure_tone_detect |= 0x40;  /* still not a pure tone,  to  great  LF+HF  side  variation/increase  */
		}
	} /* bands available*/

	return nonpure_tone_detect;
}
