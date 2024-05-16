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

LC3_FLOAT plc_phEcu_imax2_jacobsen_mag(const Complex *y, LC3_FLOAT *c_jacobPtr) {

    LC3_FLOAT posi;
    const Complex *pY;
    Complex y_m1, y_0, y_p1;
    Complex N;
    Complex D;
    LC3_FLOAT numer, denom;

    /* Jacobsen estimates peak offset relative y_0 using
    *                 X_m1 - X_p1
    *  d = REAL ( ------------------- ) * c_jacob
    *              2*X_0 - X_m1 -Xp1
    *
    *  Where c_jacob is a window  dependent constant
    */

    /* Get the parameters into variables */
    pY = y;
    y_m1 = *pY++;
    y_0 = *pY++;
    y_p1 = *pY++;

    /* prepare numerator real and imaginary parts*/
    N = csub(y_m1, y_p1);

    /* prepare denominator real and imaginary parts */
    D = cmul(cmplx(2.0, 0.0), y_0);
    D = csub(D, y_m1);
    D = csub(D, y_p1);

    /* REAL part of complex division  */
    numer = N.r*D.r + N.i*D.i;
    denom = D.r*D.r + D.i*D.i;

    if (numer != 0 && denom != 0) {
        posi = numer / denom * (*c_jacobPtr);
    } else {
        posi = 0.0; /* flat top,  division is not possible choose center freq */
    }


    posi = fclampf(-1.0, posi, 1.0);
    return posi;
}

/*-------------------------------------------------------------------*
 * imax()
 *
 * Get interpolated maximum position
 *-------------------------------------------------------------------*/

LC3_FLOAT plc_phEcu_interp_max(const LC3_FLOAT *y, LC3_INT32 y_len) {
    LC3_FLOAT posi, y1, y2, y3, y3_y1, y2i;
    LC3_FLOAT ftmp_den1, ftmp_den2;

    /* Seek the extrema of the parabola P(x) defined by 3 consecutive points so that P([-1 0 1]) = [y1 y2 y3] */
    y1 = y[0];
    y2 = y[1];
    
    /* If interp between two values only */
    if (y_len == 2) {
        if (y1 < y2) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    
    y3 = y[2];
    y3_y1 = y3-y1;
    ftmp_den1 =   (y1+y3-2*y2);
    ftmp_den2 =   (4*y2 - 2*y1 - 2*y3);
    
    if(ftmp_den2 == 0.0 || ftmp_den1 == 0.0) {
        return 0.0;  /* early exit with left-most value */
    }
    
    y2i = ((LC3_FLOAT)-0.125) * sqrf(y3_y1) /(ftmp_den1) + y2;
    /* their corresponding normalized locations */
    posi = y3_y1/(ftmp_den2);
    /* Interpolated maxima if locations are not within [-1,1], calculated extrema are ignored */
    if (posi >= (LC3_FLOAT)1.0  || posi <= (LC3_FLOAT)-1.0) {
        posi = y3 > y1 ? (LC3_FLOAT)1.0 : (LC3_FLOAT)-1.0;
    } else {
        if (y1 >= y2i) {
            posi = (y1 > y3) ? (LC3_FLOAT)-1.0 :(LC3_FLOAT) 1.0;
        } else if (y3 >= y2i) {
            posi = (LC3_FLOAT)1.0;
        }
    }
    
    return posi + (LC3_FLOAT)1.0;
}

/*-----------------------------------------------------------------------------
 * fft_spec2_sqrt_approx_ ()
 *
 * Approximation of sqrt(Square magnitude) of fft spectrum
 * if min_abs <= 0.4142135*max_abs
 *     abs = 0.99 max_abs + 0.197*min_abs
 * else
 *     abs = 0.84 max_abs + 0.561*min_abs
 * end
 *
 
 *----------------------------------------------------------------------------*/

void plc_phEcu_fft_spec2_sqrt_approx(const Complex* x, LC3_INT32 x_len, LC3_FLOAT* x_abs) {
    LC3_INT32 i;
    LC3_FLOAT max_abs, min_abs, re, im;
    
    for (i = 0; i < x_len; i++) {
        re = LC3_FABS(x[i].r);
        im = LC3_FABS(x[i].i);
        max_abs = MAX(re, im);
        min_abs = MIN(re, im);
        
        if (min_abs <= (LC3_FLOAT)0.4142135 * max_abs) {
            x_abs[i] = (LC3_FLOAT)0.99*max_abs + (LC3_FLOAT)0.197*min_abs;
        } else {
            x_abs[i] = (LC3_FLOAT)0.84*max_abs + (LC3_FLOAT)0.561*min_abs;
        }
    }

      return;
}

LC3_INT32 plc_phEcu_pitch_in_plocs(LC3_INT32* plocs, LC3_INT32 n_plocs) {
    
    LC3_INT32 i;
    LC3_INT32 p_in_plocs;

    p_in_plocs = 0;
    
    for (i = 0; i < n_plocs; i++) {
        if (plocs[i] > 0 && plocs[i] < 7) {
            p_in_plocs++;
        }
    }

    return p_in_plocs;
}

