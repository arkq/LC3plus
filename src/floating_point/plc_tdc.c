/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/***************************************************************************\
 *   contents/description: Main function for Time domain concealment
\***************************************************************************/

#include <string.h>
#include "functions.h"
#include "constants.h"

static LC3_INT16 TDC_random_short(LC3_INT16 *seed);
static LC3_FLOAT TDC_get_gainp(const LC3_FLOAT x[], const LC3_FLOAT y[], LC3_INT32 n);
static LC3_FLOAT TDC_get_gainc(const LC3_FLOAT x[], const LC3_FLOAT y[], const LC3_FLOAT *gain_p, const LC3_INT32 n, const LC3_INT32 frame_dms);
static void TDC_LPC_synthesis(const LC3_FLOAT a[], LC3_FLOAT x[], LC3_FLOAT y[], LC3_INT32 l, const LC3_FLOAT mem[], LC3_INT32 lpcorder, LC3_FLOAT *buf);
static void TDC_LPC_residu(const LC3_FLOAT *a, LC3_FLOAT *x, LC3_FLOAT *y, LC3_INT32 l, LC3_INT32 lpcorder);
static void TDC_highPassFiltering(const LC3_INT32 L_buffer, LC3_FLOAT exc2[], const LC3_FLOAT hp_filt[], const LC3_INT32 l_fir_fer);
static void TDC_f_preemph(LC3_FLOAT *signal, const LC3_FLOAT *mu, LC3_INT32 L, LC3_FLOAT *mem);
static void TDC_deemph(LC3_FLOAT *signal, const LC3_FLOAT *mu, LC3_INT32 L, const LC3_FLOAT *mem);
const LC3_FLOAT TDC_high_16[TDC_L_FIR_HP]   = { 0.f,     -0.0205f, -0.0651f, -0.1256f, -0.1792f,  0.8028f, -0.1792f, -0.1256f, -0.0651f, -0.0205f,  0.f    };
const LC3_FLOAT TDC_high_32[TDC_L_FIR_HP]   = {-0.0517f, -0.0587f, -0.0820f, -0.1024f, -0.1164f,  0.8786f, -0.1164f, -0.1024f, -0.0820f, -0.0587f, -0.0517f};
const LC3_FLOAT TDC_high_16_harm[TDC_L_FIR_HP]  = { 0.0053f,  0.0000f, -0.0440f,  0.0000f,  0.2637f,  0.5500f,  0.2637f,  0.0000f, -0.0440f,  0.0000f,  0.0053f};
const LC3_FLOAT TDC_high_32_harm[TDC_L_FIR_HP]  = {-0.0053f, -0.0037f, -0.0140f,  0.0180f,  0.2668f,  0.4991f,  0.2668f,  0.0180f, -0.0140f, -0.0037f, -0.0053f};
static void TDC_levinson(LC3_FLOAT *acf, LC3_INT32 len, LC3_FLOAT *out);
static void TDC_copyFLOAT(const LC3_FLOAT * X, LC3_FLOAT * Z, LC3_INT32 n);
static LC3_FLOAT TDC_dotFLOAT(const LC3_FLOAT * X, const LC3_FLOAT * Y, LC3_INT32 n);
static LC3_FLOAT type_2_alpha_long(LC3_INT32 nbLostFramesInRow, LC3_INT32 frame_dms);
const LC3_INT32 beforeNextIncArray[4][4] = {{0,0,0,1},
                                          {0,1,0,1},
                                          {0,1,1,1},
                                          {1,1,1,1}};
const LC3_INT32 nextIncArray[4][4] = {{1,0,0,0},
                                    {1,0,1,0},
                                    {1,0,1,1},
                                    {1,1,1,1}};

void processTdcApply_fl(const LC3_INT32    pitch_int,
                        const LC3_FLOAT  *preemphFac,
                        const LC3_FLOAT* A,
                        const LC3_INT32    lpc_order,
                        const LC3_FLOAT* pcmbufHist,
                        const LC3_INT32    max_len_pcm_plc,
                        const LC3_INT32    N,
                        const LC3_INT32    frame_dms,
                        const LC3_INT32    SampRate,
                        const LC3_INT32    nbLostFramesInRow,
                        const LC3_INT32    overlap,
                        const LC3_FLOAT  *stabFac,
                        LC3_FLOAT        harmonicBuf[MAX_PITCH],
                        LC3_FLOAT        synthHist[M],
                        LC3_INT32*         fract,
                        LC3_INT16*       seed,
                        LC3_FLOAT*       gain_c,
                        LC3_FLOAT*       alpha,
                        LC3_FLOAT*       synth
                        , LC3_UINT8       plc_fadeout_type
                        , LC3_FLOAT*      alpha_type_2_table
                        )
{    
          LC3_FLOAT step, step_n;
          LC3_INT32   i, len, Tc, nbLostCmpt_loc, nextInc, beforeNextInc;
          LC3_FLOAT gain_h, tmp, gain_p;
          LC3_FLOAT *exc2, *exc_buf, *exc, *x_pre, *buf, *pt_exc, *pt1_exc, *synthMemPtr;
          LC3_FLOAT *harmonicBufPtr;
          LC3_FLOAT synth_mem[M];
          const LC3_FLOAT *hp_filt, *high_harm;
          LC3_FLOAT gainInov;
          LC3_FLOAT hpBlendFac;
          char *scratchSpace1st, *scratchSpaceTmp;
          char scratchSpace[(MAX_LEN_PCM_PLC  + MDCT_MEM_LEN_MAX + MAX_LEN_PCM_PLC + 1 + M) * sizeof(LC3_FLOAT)];
          LC3_FLOAT alphaPrev;
          LC3_FLOAT throttle;
          LC3_INT32   frame_dms_idx, nbLostFramesInRow_mod;
          LC3_INT32 plc_fadeout_len = 0;
        
  memset(synth_mem, 0, M * sizeof(LC3_FLOAT));
  memset(scratchSpace, 0, (MAX_LEN_PCM_PLC  + MDCT_MEM_LEN_MAX + MAX_LEN_PCM_PLC + 1 + M) * sizeof(LC3_FLOAT));

  /* len of synthesized signal */
  len = N + overlap;

  nbLostCmpt_loc        = floor(frame_dms/100.0 * (nbLostFramesInRow - 1) + 1);
  frame_dms_idx         = frame_dms / 25 - 1; /* 0,1,2,3 */
  nbLostFramesInRow_mod = (nbLostFramesInRow - 1) % 4;

  beforeNextInc         = beforeNextIncArray[frame_dms_idx][nbLostFramesInRow_mod];
  nextInc               = nextIncArray      [frame_dms_idx][nbLostFramesInRow_mod]; 

  if (plc_fadeout_type >= 1){
    plc_fadeout_len = PLC_FADEOUT_TYPE_1_IN_MS;
  }
  else{
    plc_fadeout_len = PLC_FADEOUT_IN_MS;
  }

  if (nbLostCmpt_loc > plc_fadeout_len/10)
  {
      gain_p = 0;
      *gain_c = 0;
      *alpha = 0;
      memset(synth, 0, len * sizeof(LC3_FLOAT));
      return;
  }

  Tc = pitch_int;
  if (*fract > 0) {
    Tc++;
  }

  /*----------------------------------------------------------------
   * Buffer Initialization for timeDomainConcealment_Apply
   *
   *            1st
   * |--exc_buf--|--x_pre--|
   * |           |--exc2--|
   * |           |--buf (LPC Syn)--|
   *
   *---------------------------------------------------------------*/
    
  scratchSpace1st = scratchSpace;
  exc_buf = (LC3_FLOAT*)scratchSpace1st; scratchSpace1st += (LC3_INT32)sizeof(LC3_FLOAT) * (Tc + N/2 + len);
  exc = exc_buf + (Tc + N/2);

  scratchSpaceTmp = scratchSpace1st;
  x_pre = (LC3_FLOAT*)scratchSpaceTmp; scratchSpaceTmp += (LC3_INT32)sizeof(LC3_FLOAT) * (lpc_order + Tc + N/2 + 1);

  /*---------------------------------------------------------------*
   * LPC Residual                                                  *
   *---------------------------------------------------------------*/
  if (nbLostFramesInRow == 1)
  {
    /* copy buffer to pre-emphasis buffer */
    TDC_copyFLOAT(&(pcmbufHist[max_len_pcm_plc-(lpc_order+Tc+N/2+1)]), &(x_pre[0]), lpc_order+Tc+N/2+1);

    /* apply pre-emphasis to the signal */
    TDC_f_preemph(&(x_pre[1]), preemphFac, lpc_order+Tc+N/2, &x_pre[0]);

    /* copy memory for LPC synth */
    TDC_copyFLOAT(&(x_pre[Tc+N/2+1]), synth_mem, lpc_order);

    /* LPC Residual */
    TDC_LPC_residu(A, &(x_pre[lpc_order+1]), &(exc[-(Tc+N/2)]), Tc+N/2, lpc_order);
  }

  /*---------------------------------------------------------------*
   * Calculate gains                                               *
   *---------------------------------------------------------------*/
  if (nbLostFramesInRow == 1)
  {
    if (pitch_int == Tc)
    {
      gain_p = TDC_get_gainp( &(x_pre[lpc_order+Tc+1]), &(x_pre[lpc_order+1]), N/2 );
    }
    else
    {
      tmp    = TDC_get_gainp( &(x_pre[lpc_order+Tc+1]), &(x_pre[lpc_order+2]), N/2 );
      gain_p = TDC_get_gainp( &(x_pre[lpc_order+Tc+1]), &(x_pre[lpc_order+1]), N/2 );

      if (tmp > gain_p) {
        Tc = pitch_int;
        gain_p = tmp;
        *fract = 0;
      }
    }

    if(gain_p < 0.0f)
    {
      gain_p = 0.0f;
    }

    if(gain_p > 1.0f)
    {
      gain_p = 1.0f;
    }

    *gain_c = 0.0f;

    if (pitch_int == Tc)
    {
      *gain_c = TDC_get_gainc( &(exc[-1]), &(exc[-1-Tc]), &gain_p, N/2, frame_dms );
    }
    else
    {
      tmp     = TDC_get_gainc( &(exc[-1]), &(exc[-1-pitch_int]), &gain_p, N/2, frame_dms );
      *gain_c = TDC_get_gainc( &(exc[-1]), &(exc[-1-Tc])       , &gain_p, N/2, frame_dms );
      *gain_c = MIN(*gain_c, tmp);
    }
  }
  else
  {
      gain_p = *alpha;
  }

  /*---------------------------------------------------------------*
   * Damping factor                                                *
   *---------------------------------------------------------------*/

  alphaPrev = 1;
  if (nbLostFramesInRow > 1)
  {
      alphaPrev = *alpha;
  }

  if (plc_fadeout_type == 2){
    *alpha = alpha_type_2_table[nbLostFramesInRow];
  }
  else{

  if (nextInc != 0)
  {
    switch (nbLostCmpt_loc)
    {
    case 1:
      *alpha = (LC3_FLOAT)sqrt(gain_p);
      if ( *alpha > 0.98f )
      {
        *alpha = 0.98f;
      }
      else if ( *alpha < 0.925f )
      {
        *alpha = 0.925f;
      }
      break;
    case 2:
      *alpha = (0.63f + 0.35f * (*stabFac)) * gain_p;
      if ( *alpha < 0.919f )
      {
        *alpha = 0.919f;
      }
      break;
    default:
      *alpha = (0.652f + 0.328f * (*stabFac)) * gain_p;
    }
  }
  
  if (nbLostCmpt_loc > 3)
  {
      switch (frame_dms)
      {
      case  25: *alpha *= PLC34_ATTEN_FAC_025; break;
      case  50: *alpha *= PLC34_ATTEN_FAC_025; break;
      case  75: *alpha *= PLC34_ATTEN_FAC_075; break;
      case 100: *alpha *= PLC34_ATTEN_FAC_100; break;
      }
  }
  
  if (nbLostCmpt_loc > 5)
  {
      gain_p = *alpha;
  }
  }
  /*---------------------------------------------------------------*
   * Construct the harmonic part                                   *
   *  Last pitch cycle of the previous frame is repeatedly copied. *
   *---------------------------------------------------------------*/

  pt_exc = harmonicBuf;
  pt1_exc = exc - Tc;

  if( nbLostFramesInRow == 1 ) 
  {
    if (*stabFac >= 1)
    {
        TDC_copyFLOAT(pt1_exc, pt_exc, Tc);
    }
    else 
    {
      /* These values are necessary for the last five filtered samples */
      TDC_copyFLOAT(&exc[-Tc], exc, (TDC_L_FIR_HP-1)/2);

      high_harm = TDC_high_32_harm;
      if (SampRate <= 16000)
      {
          high_harm = TDC_high_16_harm;
      }
      
      for( i = 0; i < Tc; i++ )
      {
          pt_exc[i] = TDC_dotFLOAT(&pt1_exc[i-(TDC_L_FIR_HP-1)/2], high_harm, TDC_L_FIR_HP);
      }
    }
  }

  /*---------------------------------------------------------------*
   * Construct the random part of excitation                       *
   *---------------------------------------------------------------*/
  scratchSpaceTmp = scratchSpace1st;
  exc2 = (LC3_FLOAT*)scratchSpaceTmp; scratchSpaceTmp += (LC3_INT32)sizeof(LC3_FLOAT) * (len + TDC_L_FIR_HP - 1);

  for (i = 0; i < len + TDC_L_FIR_HP - 1; i++) {
    exc2[i] = (LC3_FLOAT)TDC_random_short(seed);
  }

  /* high pass noise */
  if (SampRate <= 16000 )
  {
    hp_filt = TDC_high_16;
  } else {
    hp_filt = TDC_high_32;
  }

  if ( nbLostFramesInRow == 1 )
  {
    TDC_highPassFiltering(len, exc2, hp_filt, TDC_L_FIR_HP);
  }
  else
  {
    /* moves from 0 to 1, speed is defined by PLC3_HPBLENDTHROTTLE */
    throttle = (LC3_FLOAT)nbLostCmpt_loc / (nbLostCmpt_loc + PLC3_HPBLENDTHROTTLE);
    hpBlendFac = (1 - *alpha) * throttle;

    for (i = 0; i < len; i++)
    {
      exc2[i] =  hpBlendFac * exc2[i+TDC_L_FIR_HP/2] + (1 - hpBlendFac) * TDC_dotFLOAT(&exc2[i], hp_filt, TDC_L_FIR_HP );
    }
  }

  /* normalize energy */
  gainInov = 1.0f / (LC3_FLOAT)sqrt(TDC_dotFLOAT( exc2, exc2, N ) / (LC3_FLOAT)N + 0.01f );
  gainInov *= (1.1f - 0.75* gain_p);

  /* gains */
  gain_h = alphaPrev;
  tmp = *gain_c * *alpha / alphaPrev;

  /* update steps */
  step = (1.0f/(LC3_FLOAT)N) * (gain_h - *alpha);
  step_n = (1.0f/(LC3_FLOAT)N) * (*gain_c - tmp);

  /*---------------------------------------------------------------*
   * Construct the total excitation                                *
   *---------------------------------------------------------------*/
  harmonicBufPtr = harmonicBuf + ((nbLostFramesInRow - 1) * N) % Tc;

  for ( i = 0; i < len; i++ ) {
    /* harmonic */
    if (harmonicBufPtr - harmonicBuf >= Tc) {
        harmonicBufPtr = harmonicBuf;
    }
    exc[i] = *harmonicBufPtr++;
    exc[i] *= gain_h;

    /* random */
    exc2[i] *= *gain_c * gainInov;

    /* total */
    exc[i] = exc[i] + exc2[i];

    /* update */
    gain_h -= step;
    gain_h = MAX(gain_h, 0);
    *gain_c -= step_n;
    *gain_c = MAX(*gain_c, 0);
  }

  *gain_c = tmp;

  /*----------------------------------------------------------*
   * Compute the synthesis speech                             *
   *----------------------------------------------------------*/
  buf = (LC3_FLOAT*)scratchSpace1st; scratchSpace1st += (LC3_INT32)sizeof(LC3_FLOAT) * (len + lpc_order);
  synthMemPtr = synth_mem;
  if (nbLostFramesInRow != 1)
  {
      synthMemPtr = synthHist;
  }
  
  TDC_LPC_synthesis(A,
                    &exc[0],
                    synth,
                    len,
                    synthMemPtr,
                    lpc_order,
                    buf);

  TDC_copyFLOAT(&synth[N-lpc_order], synthHist, lpc_order);

  /*----------------------------------------------------------*
   * Deemphasis                                               *
   *----------------------------------------------------------*/
  TDC_deemph( synth, preemphFac, len, &pcmbufHist[max_len_pcm_plc-1] );

  /*----------------------------------------------------------*
   * Fade to zero                                             *
   *----------------------------------------------------------*/
  if (beforeNextInc != 0)
  {
    if (nbLostCmpt_loc == plc_fadeout_len/10)
    {
      gain_h = 1;
      step = 1.0f/(LC3_FLOAT)N;
      for ( i = 0; i < N; i++ ) {
        synth[i] *= gain_h;
        gain_h -= step;
      }
      memset(&synth[N], 0, overlap * sizeof(LC3_FLOAT));
    }
  }
}

/* Take only real part */
void processTdcInverseOdft_fl(LC3_FLOAT *in, LC3_INT32 n_bands, LC3_FLOAT *out, LC3_INT32 lpc_order)
{
        LC3_INT32 i, j, k;
        LC3_FLOAT buf[2*MAX_BANDS_NUMBER_PLC];
        Complex sum;
        Complex res;

    /* Buffer for ifft */
    j = 0;
    for (i = 0; i < n_bands - 1; i += 2)
    {
        buf[j] = in[i];
        j++;
    }

    for (i = n_bands - 1; i > 0; i -= 2)
    {
        buf[j] = in[i];
        j++;
    }

    for (i = 0; i < n_bands; i++)
    {
        buf[j] = in[i];
        j++;
    }

    /* ifft */
    for (j = 0; j < n_bands; j++)
    {
        sum.r = 0, sum.i = 0;
        res.r = 0, res.i = 0;
        for (k = 0; k < n_bands; k++)
        {
            res = cexpi((2 * M_PI_LC3PLUS  * (LC3_FLOAT) (j * k)) / (LC3_FLOAT) n_bands);
            res.r = res.r * buf[k];
            res.i = res.i * buf[k];
            sum = cadd(sum, res);
        }

        res = cexpi((LC3_FLOAT) j * M_PI_LC3PLUS / (2.0 * (LC3_FLOAT) n_bands));
        out[j] = (sum.r * res.r - sum.i * res.i);
    }

    out[0] = out[0] * 1.0001;
    if (out[0] == 0)
    {
        out[0] = 1;
        zero_float(&out[1], lpc_order);
    }
}

void processTdcPreemphasis_fl(LC3_FLOAT *in, LC3_FLOAT *pre_emph_factor, LC3_INT32 n_bands)
{
        LC3_INT32 i;

    for (i = 0; i < n_bands; i++)
    {
        in[i] = in[i] * (1.0 - 2.0 * (*pre_emph_factor) * LC3_COS(2.0 * M_PI_LC3PLUS * (0.5 + (LC3_FLOAT) i) / (2.0 * (LC3_FLOAT) n_bands)) + (*pre_emph_factor) * (*pre_emph_factor));
    }
}

void processTdcLpcEstimation_fl(LC3_FLOAT *r, LC3_INT32 fs_idx, LC3_INT32 len, LC3_FLOAT *A, LC3_INT32 frame_dms)
{
        LC3_INT32 i;
        const LC3_FLOAT *lpc_array;
    
    lpc_array = plc_tdc_lpc_all[fs_idx];
    
    if (fs_idx == 0 && frame_dms == 25)
    {
        lpc_array = plc_tdc_lpc_8_25ms;
    }
    
    /* r[0] = r[0] * 1 */
    for (i = 1; i < len; i++)
    {
        r[i] = r[i] * lpc_array[i];
    }

    TDC_levinson(r, len - 1, A);
}

/** random
 *
 * Parameters:
 *    seed        I/O: seed for random number
 *
 * Function:
 *    Signed 16 bits random generator.
 *
 * Returns:
 *    random number
 */
static LC3_INT16 TDC_random_short(LC3_INT16 *seed)
{
   *seed = (LC3_INT16) (*seed * 12821L + 16831L);
   return(*seed);
}

static LC3_FLOAT TDC_get_gainp( /* output: gain of pitch                        */
  const LC3_FLOAT x[],      /* input : input signal                         */
  const LC3_FLOAT y[],      /* input : shifted input signal                 */
  LC3_INT32 n                 /* input : vector length                        */
)
{
        LC3_FLOAT corr, ener;
        LC3_INT16 i;
    
    corr = 0; ener = 1e-6f;

    for (i = 0; i < n; i++)
    {
        corr += x[i]*y[i];
        ener += y[i]*y[i];
    }

    return(corr/ener);
}

static LC3_FLOAT TDC_get_gainc( /* output: gain of code                         */
  const LC3_FLOAT x[],      /* input : input signal                         */
  const LC3_FLOAT y[],      /* input : shifted input signal                 */
  const LC3_FLOAT *gain_p,   /* input : gain of pitch                        */
  const LC3_INT32   n,        /* input : vector length                        */
  const LC3_INT32   frame_dms /* input : frame length in dms                  */
)
{
        LC3_FLOAT gain_c;
        LC3_FLOAT gain_c_max;
        LC3_INT16 i;
    
    gain_c = 0; gain_c_max = 0;

    for (i = 0; i < n; i++)
    {
      gain_c += ( x[-i] - *gain_p * y[-i] ) * ( x[-i] - *gain_p * y[-i] );
    }

    if (frame_dms < 100)
    {
        for (i = 0; i < n; i++)
        {
          gain_c_max += (x[-i] * x[-i]);
        }
        gain_c = MIN(gain_c, gain_c_max);
    }

    gain_c = (LC3_FLOAT)sqrt(gain_c / n );

    return gain_c;
}

static void TDC_highPassFiltering(const LC3_INT32   L_buffer,      /* i:   buffer length                                      */
                       LC3_FLOAT       exc2[],        /* i/o: unvoiced excitation before the high pass filtering */
                       const LC3_FLOAT hp_filt[],     /* i:   high pass filter coefficients                      */
                       const LC3_INT32   l_fir_fer)     /* i:   high pass filter length                            */
{
       LC3_INT32   i;
    
  for( i=0 ; i< L_buffer; i++ ) {
    exc2[i] = TDC_dotFLOAT(&exc2[i], hp_filt, l_fir_fer);
  }
}

static void TDC_LPC_synthesis(
                       const LC3_FLOAT a[],
                       LC3_FLOAT       x[],
                       LC3_FLOAT       y[],
                       LC3_INT32         l,
                       const LC3_FLOAT mem[],
                       LC3_INT32         lpcorder,
                       LC3_FLOAT      *buf
                       )
{
       LC3_FLOAT s, *yy;
       LC3_INT32   i, j;

   /* copy initial filter states into synthesis buffer */
   for (i=0; i < lpcorder; i++)
   {
      buf[i] = mem[i];
   }
   yy = &buf[i];

   for (i = 0; i < l; i++)
   {
      s =  x[i];
      for (j = 1; j <= lpcorder; j++)
      {
          s -= a[j] * yy[i- j];
      }
      y[i] = s;
      yy[i] = y[i];
   }

   return;
}


/** TDC_LPC_residu
 *
 * Parameters:
 *    a           I: LP filter coefficients (Q12)
 *    x           I: input signal (usually speech)
 *    y           O: output signal (usually residual)
 *    l           I: size of filtering
 *    lpcorder    I: Order of LP filter
 *
 * Function:
 *    Compute the LP residual by filtering the input speech through A(z).
 *
 * Returns:
 *    void
 */
static void TDC_LPC_residu(const LC3_FLOAT *a, LC3_FLOAT *x, LC3_FLOAT *y, LC3_INT32 l, LC3_INT32 lpcorder)
{
       LC3_FLOAT s;
       LC3_INT32   i, j;

   for (i = 0; i < l; i++)
   {
      s = x[i];
      for (j = 1; j <= lpcorder; j++)
      {
         s += a[j] * x[i - j];
      }
      y[i] = s;
   }

   return;
}


/** TDC_f_preemph
 *
 * Parameters:
 *    signal       I/O: signal
 *    mu             I: preemphasis factor
 *    L              I: vector size
 *    mem            I: memory (x[-1])
 *
 * Function:
 *    Filtering through 1 - mu z^-1
 *
 *
 * Returns:
 *    void
 */

static void TDC_f_preemph(LC3_FLOAT *signal, const LC3_FLOAT *mu, LC3_INT32 L, LC3_FLOAT *mem)
{
      LC3_INT32 i;

   for (i = L - 1; i > 0; i--)
   {
      signal[i] = signal[i] - *mu * signal[i - 1];
   }

   signal[0] -= *mu * (*mem);

   return;
}

/*
 * TDC_deemph
 *
 * Parameters:
 *    signal       I/O: signal
 *    mu             I: deemphasis factor
 *    L              I: vector size
 *    mem            I: memory (signal[-1])
 *
 * Function:
 *    Filtering through 1/(1-mu z^-1)
 *    Signal is divided by 2.
 *
 * Returns:
 *    void
 */
static void TDC_deemph(LC3_FLOAT *signal, const LC3_FLOAT *mu, LC3_INT32 L, const LC3_FLOAT *mem)
{
      LC3_INT32 i;

   signal[0] = signal[0] + *mu * (*mem);

   for (i = 1; i < L; i++)
   {
      signal[i] = signal[i] + *mu * signal[i - 1];
   }

   return;
}

static void TDC_copyFLOAT(const LC3_FLOAT * X, LC3_FLOAT * Z, LC3_INT32 n)
{
  /* no values to copy */
  if ( (n < 1) || (X == Z) ){
    return;
  }
  /* If overlapping */
  if ( ( (Z > X) && (Z < X+n) ) || ( (Z < X) && (X < Z+n) ) ) {
    memmove(Z, X, sizeof(LC3_FLOAT)*n);
  }
  else{
    memcpy(Z, X, sizeof(LC3_FLOAT)*n);
  }
}

static LC3_FLOAT TDC_dotFLOAT(const LC3_FLOAT * X, const LC3_FLOAT * Y, LC3_INT32 n)
{
      LC3_FLOAT acc;
      LC3_INT32 i;

  acc = 0;
  if (n) {
    acc = X[0]*Y[0];
  }

  for (i=1; i<n; i++) acc += X[i]*Y[i];
    
  return acc;
}

static void TDC_levinson(LC3_FLOAT *acf, LC3_INT32 len, LC3_FLOAT *out)
{ 
        LC3_FLOAT g, v, sum, buf[M], buf2[M];
        LC3_INT32 i, j, k, t;

    g = acf[1] / acf[0];

    out[0] = g;
    v = (1.0 - g * g) * acf[0];

    for (t = 1; t < len; t++)
    {
        sum = 0; j = 0;
        for (i = 1; i <= t; i++)
        {
            sum += out[j] * acf[i];
            j++;
        }

        g = (acf[t + 1] - sum) / v;

        move_float(buf, out, len);
        move_float(buf2, out, len);
        out[0] = g;

        j = 1, k = 0;
        for (i = t - 1; i >= 0; i--)
        {
            out[j] = buf[k] - g * buf2[i];
            j++; k++;
        }

        v = v * (1.0 - g * g);
    }

    move_float(buf, out, len);
    out[0] = 1;

    j = 1;
    for (i = len - 1; i >= 0; i--)
    {
        out[j] = -buf[i];
        j++;
    }
}

static LC3_FLOAT type_2_alpha_long(LC3_INT32 nbLostFramesInRow, LC3_INT32 frame_dms)
{   
    if (nbLostFramesInRow <= 3*100.0/frame_dms){
       return LC3_POW(0.95,(nbLostFramesInRow + (100.0/frame_dms) - 1) * frame_dms/100.0);
    }
    else {
      LC3_INT32 n_shift = (nbLostFramesInRow - 3*100.0/frame_dms) * 50/frame_dms;
      return LC3_POW(0.7,(n_shift + 100.0/frame_dms - 1) * frame_dms/100.0);
    }
}

LC3_FLOAT type_2_fadeout(LC3_INT32 nbLostFramesInRow, LC3_INT32 frame_dms)
{   
    LC3_FLOAT selector = PLC_FADEOUT_TYPE_2_SELECTOR * 2 * 100/frame_dms;
    if (selector >= nbLostFramesInRow){
      return type_2_alpha_long(nbLostFramesInRow, frame_dms);
    }
    else {
      return LC3_POW(0.5,(nbLostFramesInRow + (100.0/frame_dms) - 1) * frame_dms/100.0);
    } 
}
