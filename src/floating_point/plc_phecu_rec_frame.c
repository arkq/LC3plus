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

void plc_phEcu_rec_frame(Complex *X_in,
   LC3_INT32 L,
   LC3_INT32 Lecu,
   const LC3_FLOAT *whr,
   const LC3_FLOAT *winMDCT,
   LC3_INT32 Lprot,
   LC3_FLOAT *xfp,
   LC3_INT32 time_offs,
   LC3_FLOAT *x_out,
   Complex *full_spec_dbg,
   LC3_FLOAT* ifft_out_dbg,
   LC3_FLOAT* xsubst_dbg,
   LC3_INT32 LA_ZEROS,
   LC3_INT32 LA,
   Fft* PhEcu_ifft
)
{

   LC3_INT32 i;

   LC3_FLOAT xrec[2*MAX_LEN];
   LC3_FLOAT xsubst[2*MAX_LEN];
   LC3_FLOAT xsubst_LL[2*MAX_LEN];
   LC3_FLOAT *pXsubst_LL;

   LC3_INT32 fs_idx;

   LC3_FLOAT *pXfp, *pOlaXsubst, *pXOut;
   LC3_INT32 work_part, copy_part, ola_part;

   const LC3_FLOAT *hannOla;
   const LC3_FLOAT *pHannOla;

   UNUSED(time_offs);
   UNUSED(full_spec_dbg);
   UNUSED(ifft_out_dbg);
   UNUSED(xsubst_dbg);
   UNUSED(xsubst_LL);

   fs_idx = FRAME2FS_IDX_10MS(L);
   hannOla = hannOla_wins[fs_idx];

   X_in[0].i = X_in[Lprot / 2].r; /* move fs/2 real to imag part of X_in[0]*/

   real_fft_apply(PhEcu_ifft, (LC3_FLOAT*)X_in, xrec);

   move_float(xsubst, xrec, Lprot);




   {
      for (i = 0; i < Lprot; i++) {

         if (whr[i] != 0) {
            xsubst[i] = xsubst[i] / whr[i];     /*  inverse stored in BASOP */
         }

      }

      assert(xsubst_LL != NULL);
      zero_float(xsubst_LL,  (Lecu-Lprot)/2);  /* initial 2ms */
      zero_float(&(xsubst_LL[ Lecu- (Lecu-Lprot)/2]),  (Lecu-Lprot)/2);  /* tail 2ms */
      {
         /*    position reconstruction  properly  */
         /*          pXsubst_LL = &xsubst_LL[Lecu - Lprot - (Lecu - Lprot) / 2]; */
         pXsubst_LL = &xsubst_LL[(Lecu - Lprot) / 2];
         for (i = 0; i <  Lprot ; i++) {
            *pXsubst_LL++ = xsubst[i];  /*    copy required 14.25 ms into  center  */
         }
      }

   }



   work_part = LA_ZEROS + LA;
   copy_part = (Lecu - Lprot) / 2;
   ola_part = work_part - copy_part;

   pXfp = &xfp[Lprot - work_part];
   for (i = 0; i < copy_part; i++) {
      xsubst_LL[i] = *pXfp++;
   }

   assert(xsubst_LL != NULL);
   pOlaXsubst = &(xsubst_LL[copy_part]);
   pHannOla = hannOla;
   for (i = 0; i < ola_part; i++) {
      *pOlaXsubst = *pOlaXsubst * *pHannOla++;
      pOlaXsubst++;
   }

   pOlaXsubst = &(xsubst_LL[copy_part]);
      for (i = 0; i < ola_part; i++) {
         *pOlaXsubst = *pOlaXsubst + *pXfp++ * *pHannOla--;
         pOlaXsubst++;
      }


   /* clear x_out to start with */
   assert(x_out != NULL);
   zero_float(x_out, L);


   for (i = 0; i < (Lecu - LA_ZEROS); i++) {

      xsubst_LL[i] = xsubst_LL[i] * winMDCT[i];  /*  xsubstLL windowing   up to  16.25 ms  i.e not last 3.75 ms */

   }
     zero_float(&(xsubst_LL[Lecu - LA_ZEROS]),  LA_ZEROS); /* tail 3.75ms always zero */

   /* perform tda */

    /* first half */
   pXsubst_LL = &xsubst_LL[3 * Lecu / 4];
   pXfp = &xsubst_LL[(3 * Lecu / 4) - 1];

   pXOut = x_out;
   for (i = 0; i < Lecu / 4; i++) {
      *pXOut++ = -*pXsubst_LL++ - *pXfp--; /* 3.75 ms mults with 0 . may be skipped, see BASOP */
   }

   /* second half */
   /* */

   pXsubst_LL = &(xsubst_LL[0]);
   pXfp = &xsubst_LL[(Lecu / 2) - 1];
   for (i = 0; i < Lecu / 4; i++) {
      *pXOut++ = *pXsubst_LL++ - *pXfp--;
   }
}

