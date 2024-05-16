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

#define PEAK_LOCATOR_RES_FX    1  /*  fixed point resolution minimum value  */

static LC3_INT16  plc_phEcu_find_ind_fx(                        /* o  : output maximum  indx 0.. len-1    */
    const LC3_INT16 *inp,      /* i  : vector     */
    const LC3_INT16  len,      /* i  : length     */
    const LC3_INT16  val   /* i  : value to find     */
);

static void plc_phEcu_peak_locator_fxlike(const LC3_INT16 *inp, /* i: vector with values >=0   ,Qx      */
    const LC3_INT16  inp_len,              /* i: length of inp                                 */
    LC3_INT16 *      int_plocs,            /* o:  array of filtered integer plocs           Q0 */
    LC3_INT16 *      n_fsc,                /* o:  total_ number of filtered located highs   Q0 */
    const LC3_INT16  sens,                 /* i  sensitivity,   Qx */
    const LC3_INT16  inp_high,             /* i  global high ,  Qx */
    const LC3_INT16  inp_low               /* i:  global low,  Qx */
);


void plc_phEcu_spec_ana(LC3_FLOAT* xfp, 
   LC3_INT32 xfp_len, 
   const LC3_FLOAT* whr, 
   LC3_FLOAT *pfind_sensPtr, 
   LC3_INT32* plocs, 
   LC3_INT32* n_plocs, 
   LC3_FLOAT* f0est, 
   Complex* x, 
   LC3_INT32* x_len, 
   LC3_FLOAT * f0hzLtpBinPtr, 
   LC3_FLOAT * f0gainLtpPtr, 
   LC3_INT32 bw_idx,
   Fft* PhEcu_fft
) 
{


   LC3_INT32 i, peak_range_1, curr;
   LC3_FLOAT xfp_w[MAX_PLC_LPROT];
 
   LC3_FLOAT Xabs[MAX_LEN] = {0};
   LC3_FLOAT inp_high, inp_low, sens;
   LC3_FLOAT interPos;
   Complex   Xana_p[3];
   LC3_INT32   P_in_plocs;
   LC3_INT32   nSubs;
   LC3_INT32   n_plocs_in;
   LC3_FLOAT phEcu_c_jacob[1];

    LC3_FLOAT fx_fft_scale;
    LC3_FLOAT fft_fs_scale;

    LC3_FLOAT max_xfp_abs;
    LC3_FLOAT PLC2_Q_flt;
    LC3_FLOAT Q_scale_flt;

    LC3_INT16    Xabs_fx[MAX_LEN];
    LC3_INT16    plocs_fx[MAX_LEN];

    LC3_INT16    sens_fx;
    LC3_INT16    inp_high_fx;
    LC3_INT16    inp_low_fx;
    LC3_INT16    n_plocs_fx;

    LC3_FLOAT pfind_sens  ;
    LC3_FLOAT f0hzLtpBin  ; 
    LC3_FLOAT f0gainLtp   ; 
 
   pfind_sens = *pfind_sensPtr;
   f0hzLtpBin = *f0hzLtpBinPtr;
   f0gainLtp = *f0gainLtpPtr;

    for (i = 0; i < xfp_len; i++)
    {
       xfp_w[i] = xfp[i] * whr[i]; /* whr windowing may be split into three segments ,  two loops,  and possibly inplace */
    }
    real_fft_apply(PhEcu_fft, xfp_w,  (LC3_FLOAT *)x); 

    x[xfp_len/2].r = x[0].i;  /* move  the real Fs/2 value to end  */
    x[xfp_len/2].i = 0;       /* safety clear imaginary   Fs/2 value at  end */
    x[0].i = 0.0;             /* safety, make DC  value only real  */

    
    *x_len = xfp_len/2 + 1;
    
   i =(LC3_INT32) LC3_FLOOR(20000.0/PHECU_FRES)+1;
   zero_cmplx( &(x[i]),  *x_len - i);
    
    peak_range_1 = (LC3_INT32) MIN(*x_len, (40000.0 / 100 * 1.6) / 2 + 1);
    
   plc_phEcu_fft_spec2_sqrt_approx(x, peak_range_1, Xabs);

   zero_float(&(Xabs[peak_range_1]), *x_len - peak_range_1);
    
    inp_high = Xabs[0];
    inp_low = Xabs[0];
    
    for (i = 1; i < peak_range_1; i++) {
        inp_high = MAX(inp_high, Xabs[i]);
        inp_low = MIN(inp_low, Xabs[i]);
    }
    
    sens = (inp_high-inp_low)*(1-pfind_sens);

   if (inp_high >  ((LC3_FLOAT) PEAK_LOCATOR_RES_FX)/2.0)
    { 
       {
          /* from ROM constants.c */
          LC3_FLOAT fx_fft_scales[5] = { 6, 7, 7, 8, 8 };           /*NB,WB, sSWB, SWB, FB*/
          fx_fft_scale = LC3_POW(2.0, fx_fft_scales[bw_idx]);  /*% scaling due to up / dn pre shifts in fx FFT */
       }
       {   /* from  ROM constants.c */
          LC3_FLOAT fx_fs_scales[5] = { 1.0, 1.0, 1.5, 1.0, 1.5 };           /*NB,WB, sSWB, SWB, FB*/
          fft_fs_scale = fx_fs_scales[bw_idx];
       }


        max_xfp_abs = (LC3_FLOAT) LC3_FABS(xfp[0]);
        for (i = 1; i < xfp_len; i++) {
            max_xfp_abs = (LC3_FLOAT) MAX(max_xfp_abs, LC3_FABS(xfp[i]));
        } 
        
        if (max_xfp_abs >= 0.5)
        {  
           PLC2_Q_flt = (LC3_FLOAT)LC3_FLOOR(LC3_LOGTWO(32768 / 2 / 2 / max_xfp_abs));
           Q_scale_flt = LC3_POW(2.0, PLC2_Q_flt) / fx_fft_scale / fft_fs_scale; /* basop way using xfp scale */

           /*   C-Float additional safety limit/verification of  the integer xfp based scaling using the available  C-float Xabs max value inp_high as well */
           {
              LC3_FLOAT tmp_scale;
              tmp_scale = LC3_POW(2.0, LC3_FLOOR(LC3_LOGTWO(32768 / 2 / 2 / inp_high)));
              if (Q_scale_flt > tmp_scale) {
                 Q_scale_flt = tmp_scale;              
              }
           }
           /* Round sens, inp_high, inp_low according to BASOP fix-point scaling */

           for (i = 0; i < peak_range_1; i++) {
              Xabs_fx[i] = (LC3_INT16)   LC3_ROUND(Xabs[i]  * Q_scale_flt) ;
           }
           sens_fx       = (LC3_INT16)   LC3_ROUND(sens      * Q_scale_flt) ;    
           inp_high_fx   = (LC3_INT16)   LC3_ROUND(inp_high  * Q_scale_flt) ;
           inp_low_fx    = (LC3_INT16)   LC3_ROUND(inp_low   * Q_scale_flt) ;
           plc_phEcu_peak_locator_fxlike(Xabs_fx, peak_range_1, plocs_fx, &n_plocs_fx, sens_fx, inp_high_fx, inp_low_fx);
           
           *n_plocs = (LC3_INT32)n_plocs_fx;
           for (i = 0; i < *n_plocs; i++) {
              plocs[i] = (LC3_INT32)plocs_fx[i]; /* short Word16  values now stored/saved  as Word32 */
           }
        }      
        else 
        {
                *n_plocs = 0;   /* time domain xfp level near zero */
        }
    }
    else
    {
        *n_plocs = 0; /* Freq domain  Xabs  max level near zero */
    }
    
    for (i = 0; i < *n_plocs; i++) {
      curr = plocs[i];
        if (curr == 0) {
            interPos = plc_phEcu_interp_max(Xabs, 3);  /* returns 0.0 ... 2.0  */
            if (interPos == 2) {
            /* integer peak was at DC, restrict to one of coeffs at [DC or   DC+1] */ 
                interPos = plc_phEcu_interp_max(Xabs, 2); /* returns 0.0 or 1.0 */
            }
            interPos += plocs[i];
        } else if (curr == 1) {
            interPos = plc_phEcu_interp_max(Xabs, 3);
            interPos += plocs[i] - 1;
        } else if (curr == *x_len - 2) {
            interPos = plc_phEcu_interp_max(&Xabs[*x_len - 3], 3);
            interPos += plocs[i] - 1;
        } else if (curr == *x_len - 1) {
         /* integer curr at Fs/2, a real coeff  */
            interPos = plc_phEcu_interp_max(&Xabs[*x_len - 3], 3); /* returns 0.0 ... 2.0  */
         interPos += plocs[i] - 2;    /*  valid for range   ]... 1.0 ... 2.0]  , where 1 is fs/2-1 and 2.0 is Fs/2    */
            if (interPos == 0) {
            /* restrict to one of coeffs at [fs/2-1,  fs/2 ] */
            interPos = plc_phEcu_interp_max(&Xabs[*x_len - 2], 2);    /* returns 0.0 or 1.0 */
                interPos += plocs[i] - 1;
            } 

            if (interPos > (*x_len - 1) ) { /* interPos only defined  up to  Fs/2  */
                interPos = (LC3_FLOAT)(*x_len - 1);
            }
      } else {
            Xana_p[0] = x[plocs[i]-1];
            Xana_p[1] = x[plocs[i]];
            Xana_p[2] = x[plocs[i]+1];
         phEcu_c_jacob[0] = (LC3_FLOAT)PHECU_C_JACOB;
            interPos = plc_phEcu_imax2_jacobsen_mag(Xana_p, phEcu_c_jacob );
            interPos += (LC3_FLOAT) plocs[i];
        }    
        f0est[i] = interPos;
    }
    
    if (*n_plocs >= 2 &&   plocs[0] == 0 &&
      f0est[0] > f0est[1] && plocs[1] <= 2 && Xabs[0] < Xabs[plocs[1]+1]) 
   {
        f0est[0] = f0est[1];
    }
    
    P_in_plocs = plc_phEcu_pitch_in_plocs(plocs, *n_plocs);
    
    if (f0hzLtpBin > 0.0 && P_in_plocs > 0) {
        nSubs = 2;
        n_plocs_in = *n_plocs;
        plc_phEcu_LF_peak_analysis(plocs, n_plocs, f0est, Xabs, &f0hzLtpBin, &f0gainLtp, nSubs);
        
        if (n_plocs_in == *n_plocs) {
            nSubs = 3;
            plc_phEcu_F0_refine_first(plocs, *n_plocs, f0est, *x_len, &f0hzLtpBin, &f0gainLtp, nSubs);
        }
    }
    
    if (f0gainLtp > 0.0 && f0gainLtp < 0.5 && *n_plocs > 14) {
        if (P_in_plocs > 0) {
            *n_plocs = 0;
        }
    }

   return;
}


#define sub(a,b) (a - b)
#define add(a,b) (a + b)
#define s_xor(a,b) (a ^ b)

/* in case a  value (e.g max or min)  is already known , find the first corresponding array  index */
static LC3_INT16 plc_phEcu_find_ind_fx(                        /* o  : output maximum  indx 0.. len-1    */
    const LC3_INT16 *inp,      /* i  : vector     */
    const LC3_INT16  len,      /* i  : length     */
    const LC3_INT16 val   /* i  : value to find     */
)
{
    LC3_INT16  val_ind;
    LC3_INT16 pos;

    val_ind = -1;  

    for(pos = 0; pos < len; pos++)
    {
        if (sub(inp[pos], val) == 0)
        {
            val_ind = pos;
        }
    }

    return   val_ind;
}



/* BASOP  function adapted to compile in float/integer  environment */
/*-----------------------------------------------------------------------------
 * plc_phEcu_peak_locator_fxlike()
 *----------------------------------------------------------------------------*/
static void plc_phEcu_peak_locator_fxlike(const LC3_INT16 *inp, /* i: vector with values >=0   ,Qx      */
    const LC3_INT16  inp_len,              /* i: length of inp                                 */
    LC3_INT16 *      int_plocs,            /* o:  array of filtered integer plocs           Q0 */
    LC3_INT16 *      n_fsc,                /* o:  total_ number of filtered located highs   Q0 */
    const LC3_INT16  sens,                 /* i  sensitivity,   Qx */
    const LC3_INT16  inp_high,             /* i  global high ,  Qx */
    const LC3_INT16  inp_low               /* i:  global low,  Qx */
)
{

   LC3_INT16        j, k, n, idx_high, idx_low;
   LC3_INT16       inp_len_minus1;
   LC3_INT16       pairs_start, pairs_end;
   LC3_INT16       *p_tmp;
   LC3_INT16       prev_delta, curr_delta;
   LC3_INT16       delta_predc, delta_fin;
   LC3_INT16       add_dc_flag, add_fin_flag;
   LC3_INT16       low_val_cand_pairs, val_range;
   LC3_INT16       num_pairs, n_tail_values;
   LC3_INT16       cand_phase_start, cand_idx, prev_low_plus_sens, tmp;
   LC3_INT16       cand_high, prev_low;
   LC3_INT16       *cand_pairs;  /*  actually  [DC ] + pairs + [FS/2]   */

   LC3_INT16       sc_idx[1 + 368 + 1];
   LC3_INT16       cand_pairs_buf[1 + 1 + 368 + 1];
   LC3_INT16       fsc_idx[1 + 368 / 2 + 1];

   
    inp_len_minus1 = sub(inp_len, 1);  /* size of delta=derivative array ,and last index in inp */

    cand_pairs = &cand_pairs_buf[1];  /* ptr init , make space for storing a lowest  amplitude value in location  -1    */
    pairs_start = 1;                 /* adjusted to zero or 1 or 2 when/if,  DC is injected  as sc_idx[0], or initial plateau skipped */

    p_tmp = &(sc_idx[pairs_start]); /* ptr init */


    /*  xor high/low pairs of delta_inp and save sign changes */
    prev_delta = sub(inp[1], inp[0]);  /*  precompute very first delta */

    for(n = 1;  n < inp_len_minus1; n++)
    {   /* sign change analysis */
        curr_delta = sub(inp[n + 1], inp[n]);    /*  n+1 ,n ,   are loop ptrs   */
        if (s_xor(prev_delta, curr_delta) < 0)   /* a "0" delta  treated as a  positive sign */
        {
            *p_tmp++ = n;               /* store sign change bin locations , location n in the inp[] signal */
        }
        prev_delta = curr_delta; 
    }

    k = (LC3_INT16)(p_tmp - &(sc_idx[pairs_start]));

    /* copy sign change location values to a pairs array */
    /* leave one initial sc_idx location open for a potential initial DC value */

    for(j = 0; j < k; j++){
        cand_pairs[j + pairs_start] = inp[sc_idx[j + pairs_start]];     
    }

    /* filter away a potential  single initial/trailing  plateau
      to enable correct analysis for adding DC or  fs/2 bins */

    
    if((sub(k, 2) >= 0) &&
        (sub(cand_pairs[pairs_start], cand_pairs[pairs_start + 1]) == 0)){
        pairs_start = add(pairs_start, 1);
        k = sub(k, 1);
    }

    /* filter away potential single trailing plateu */
    pairs_end = sub(add(pairs_start, k), 1);  /* point  to last established  sign change element  */
    
    if ((sub(k, 2) >= 0) &&
        (sub(cand_pairs[sub(pairs_end, 1)], cand_pairs[pairs_end]) == 0)){
        k = sub(k, 1);
    }
    pairs_end = sub(add(pairs_start, k), 1);  /*  recalc  ptr to last element  */


    /* conditionally add high/lows  on both sides of input (pre_dc or fin) as  candidates  */
    add_dc_flag = 0; 
    add_fin_flag = 0; 


    if(sub(k, 1) == 0) /*  one single sign change found special case */
    {
        if (sub(inp[0], cand_pairs[pairs_start]) != 0)
        {
            add_dc_flag = 1;    /* not plateau    */
        }

        if (sub(cand_pairs[pairs_end], inp[inp_len_minus1]) != 0)
        {
            add_fin_flag = 1;     /* not plateau    */
        }
    }

    if(sub(k, 2) >= 0)
    {
        delta_predc = sub(cand_pairs[pairs_start + 1], cand_pairs[pairs_start]);
        delta_fin = sub(cand_pairs[pairs_end], cand_pairs[pairs_end - 1]);

        /* plateaus are allowed to be detected by xor sign change,
           but still not allowed at the start nor  at the end */

        add_dc_flag = 1;   
        if (sub(inp[0], cand_pairs[pairs_start]) == 0)
        {
            add_dc_flag = 0;      /* plateau down or , plateaus up., --> do not add DC  */
        }

        
        if ((sub(inp[0], cand_pairs[pairs_start]) < 0) && (delta_predc > 0))
        {
            add_dc_flag = -1;    /*UP - up    ... replace */
        }
        
        if ((sub(inp[0], cand_pairs[pairs_start]) > 0) && (delta_predc < 0))
        {
            add_dc_flag = -1;      /* DOWN - down ... % replace */
        }

        add_fin_flag = 1;   
        if (sub(cand_pairs[pairs_end], inp[inp_len_minus1]) == 0)
        {
            add_fin_flag = 0;      /* up - plateau ... */
        }
        
        if ((delta_fin > 0) && (sub(cand_pairs[pairs_end], inp[inp_len_minus1]) < 0))
        {
            add_fin_flag = -1;      /* up - UP ...    % replace , hard to hit  */
        }
        
        if ((delta_fin < 0) && (sub(cand_pairs[pairs_end], inp[inp_len_minus1]) > 0))
        {
            add_fin_flag = -1;    /*down - DOWN ... % replace */
        }

    }

    if(add_dc_flag > 0)
    {  /* add DC */
        pairs_start = sub(pairs_start, 1);
        cand_pairs[pairs_start] = inp[0]; 
        sc_idx[pairs_start] = 0; 
        k = add(k, 1);
    }
    if(add_dc_flag < 0)
    { /*   -1 -->  replace with DC*/
        cand_pairs[pairs_start] = inp[0]; 
        sc_idx[pairs_start] = 0; 
    }

    if(add_fin_flag > 0)
    {  /* add FS/2  */
        pairs_end = add(pairs_end, 1);
        cand_pairs[pairs_end] = inp[inp_len_minus1]; 
        sc_idx[pairs_end] = inp_len_minus1; 
        k = add(k, 1);
    }
    if(add_fin_flag < 0)
    {  /*    -1, replace tail with FS/2*/
        cand_pairs[pairs_end] = inp[inp_len_minus1]; 
        sc_idx[pairs_end]     = inp_len_minus1; 
    }
    /* preliminary cand_pairs now only have  highs , lows , no initial/trailing plateaus */


    /*  we allow the  DC/FsBy2 lows to be used as the candidatelLow  */
    low_val_cand_pairs = inp_low;   
    val_range = sub(inp_high, low_val_cand_pairs); /* used to determine if search is useful at all */

    
    if ((sub(val_range, PEAK_LOCATOR_RES_FX) < 0) ||
        (sub(inp_high, sens) < 0))
    {
        k = 0;   
    }

    
    if ((k == 0) && (sub(val_range, sens) >= 0))
    {
        k = 1;  
    }


    if(sub(k, 2) > 0)
    {
        /*  low, high, low, ... or
            high, low, high, ...*/

        cand_phase_start = pairs_start;       /*assume first candidate   is a high */
        if (sub(cand_pairs[pairs_start], cand_pairs[pairs_start + 1]) < 0)
        {
            cand_phase_start = add(pairs_start, 1);     /* first is a low, --> skip to next higher cand  */
        }

        /*  high, low, high, ... */
        tmp = k;   
        if (sub(cand_phase_start, pairs_start) != 0)
        {
            tmp = sub(tmp, 1);
        }
        num_pairs = tmp / 2; // shr(tmp, 1);
        n_tail_values = sub(tmp, num_pairs * 2); // shl(num_pairs, 1));

        /* filter  preliminary  sign changes into sensitivity filtered sign changes */

        *n_fsc = 0;                        /*   counter of  filtered fsc_idx */
        cand_high = low_val_cand_pairs;      
        cand_idx = -1;                       /*  sentinel location for no high cand found yet. */
        cand_pairs[-1] = low_val_cand_pairs; 

        prev_low = low_val_cand_pairs;    
        prev_low_plus_sens = add(prev_low, sens);

        /* filter loop for   high - low sign change pairs */
        /* idx_high, idx_low are raw pointers into the  cand_pairs and sc_idx arrays */

        for(idx_high = cand_phase_start;  idx_high < (cand_phase_start + 2 * num_pairs); idx_high += 2)
        {
            idx_low = idx_high + 1;  /* loop ptr increase */

            /* new high candidate  larger than previous candidate  and   */
            /* sensitivity still larger  than the the previous low */
            tmp = MAX(cand_high, prev_low_plus_sens);
            if (sub(cand_pairs[idx_high], tmp) > 0)
            {
                cand_idx = idx_high;                 /*   enable or shift candidate position fwd */
            }
            cand_high = cand_pairs[cand_idx];    /* NB, cand_pairs[-1] , has the low_val_cand_pairs value  stored */

            /* now check the fwd  idx_low  of the current  {high,low} pair  */
            prev_low = MIN(cand_pairs[idx_low], prev_low);

            tmp = sub(cand_high, sens);
            if(sub(tmp, cand_pairs[idx_low]) > 0)
            {
                /*  this low  point is now low enough to fix a previous high candidate */

                fsc_idx[*n_fsc] = cand_idx;    /*%  add cand high idx  -> output idx list*/
                *n_fsc = add(*n_fsc, 1);

                prev_low = cand_pairs[idx_low];     /*  use this value  as new low estimate */
                cand_idx = -1;                     /*   no  candidate until next pair or tail  bin, and pt to lowVal */
                cand_high = low_val_cand_pairs;    /*  enable next candidate to be selected immediately  */
            }
            prev_low_plus_sens = add(prev_low, sens);
        } /* { high, low} for loop */

        
        if((n_tail_values == 0) && (cand_idx >= 0))
        {
            /*  no tail  low or high value  to analyze
                still may need to lock a non-locked but qualified candidate */
            fsc_idx[*n_fsc] = cand_idx;  
            *n_fsc = add(*n_fsc, 1);
        }


        /* cand_pairs vector may have a last orphan value */
        if(n_tail_values > 0)
        {
            /*   cand_pairs vector may have a last orphan tail value */
            /*
             logic boils down to   if (nTailValues > 0) && (cand_pairs(n_end) > tmp)
              there is a last  one  trailing high to process

             a) the last high, may be a new high Peak if we have not yet
                locked  the current candidate
             b) if we have locked the last candidate, the last high may also be
                a highpeak if it is high enough from the(newly set previous) valley floor.

               tmp=a||b
            */

            tmp = MAX(cand_high, prev_low_plus_sens);
            tmp = sub(cand_pairs[pairs_end], tmp);
            if(tmp > 0)
            {
                fsc_idx[*n_fsc] = pairs_end;            
                *n_fsc = add(*n_fsc, 1);
            }
            else
            {
               if(cand_idx >= 0)
               { /* we have a previously established high candidate */
                  fsc_idx[*n_fsc] = cand_idx;   
                  *n_fsc = add(*n_fsc, 1);
               }

            }
        }

        /* move high locations info from  fsc_idx , to output  */
        for(j = 0; j < *n_fsc; j++)
        { 
           int_plocs[j] = sc_idx[fsc_idx[j]];    
             
        }

    } /* end of  pairs + [tail] section filtering  */
    else
    {
        /* constant/single  rise or constant decay or very low overall values,   cases */
        *n_fsc = 0;  
        
        tmp = sub(inp_high, sens);
        if((k != 0) && (sub(tmp, low_val_cand_pairs) > 0))
        {
            /*      low,high */
            /*      high,low */
            tmp = plc_phEcu_find_ind_fx(inp, inp_len, inp_high);  
            int_plocs[0] = tmp;    /*  simply  locate the   high peak*/
            *n_fsc = 1;     
            if (tmp < 0)
            {  /*safety in case max value index was not found */
               *n_fsc = 0; 
            }
         }
    }

    return;
}

