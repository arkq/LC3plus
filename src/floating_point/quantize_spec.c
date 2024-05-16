/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

LC3_INT32 find_last_nz_pair(LC3_INT32 x[], LC3_INT32 length);
LC3_INT32 find_last_nz_pair(LC3_INT32 x[], LC3_INT32 length)
{
    LC3_INT32 last_nz, lobs[4]; LC3_INT32 stage, i;

    lobs[0] = 4;

    lobs[1] = (length >> 1); /* length/2 */

    lobs[2] = (lobs[1]+ (length >> 2));

    lobs[3] = (lobs[2]+ (length >> 3));


    last_nz = 0;

    i = length;

    for (stage = 3; stage >= 0; --stage)
    {
        /* unmapped kernel */
        for (; i >= lobs[stage]; i -= 2)
        {
            if (x[i - 2] != 0)
            {
                last_nz = MAX(last_nz, i);
            }
            if (x[i - 1] != 0)
            {
                last_nz = MAX(last_nz, i);
            }
        }
        if (last_nz > 0)
        {
            break;
        }
    }

    return MAX(last_nz, 2);
}


void processQuantizeSpec_fl(LC3_FLOAT x[], LC3_FLOAT gain, LC3_INT xq[], LC3_INT nt, LC3_INT totalBits, LC3_INT* nbits, LC3_INT* nbits2, LC3_INT fs,
                            LC3_INT* lastnzout, LC3_INT* codingdata, LC3_INT* lsbMode, LC3_INT mode, LC3_INT target, LC3_INT hrmode)
{

    LC3_INT rateFlag, i, lastnz2, m, maxlev, k;
    LC3_INT nbits_lsb;
    LC3_INT c;
    LC3_INT a, b, lev1, sym, t, pki;
    LC3_INT a1_msb, b1_msb;
    LC3_INT lastnz = 1, nt_half;
    LC3_FLOAT offset = 0.375;
	LC3_INT32 bits, bits2;
    LC3_FLOAT inv_gain;
        
    assert(target >= 0);

    nbits_lsb = 0;

    nt_half = nt >> 1;
    rateFlag = 0; c = 0;

    if (hrmode)
    {
        offset = 0.5;
    }

    /* Quantization */
    inv_gain = 1.0 / gain;
    
    for (i = 0; i < nt; i++) {
        if (x[i] > 0)
        {
            xq[i] =   (LC3_INT32) ( x[i] * inv_gain + offset);
        }
        else
        {
            xq[i] = -((LC3_INT32) (-x[i] * inv_gain + offset));
        }
        if (hrmode == 0) {
            assert(xq[i] <= 32767 && xq[i] >= -32768);
        }
    }

    /* Rate flag */
    if (fs != 96000 && (totalBits > (160 + FS2FS_IDX(fs) * 160)))
    {
        rateFlag = 512;
    }

    /* Init */
    if (fs != 96000 && (mode == 0 && (totalBits >= (480 + FS2FS_IDX(fs) * 160))))
    {
        mode = 1;
    }

    /* Last non-zero 2-tuple */
    for (i = nt - 2; i >= 2; i = i - 2) {
        if (xq[i + 1] != 0 || xq[i] != 0) {
            lastnz = i + 1;
            break;
        }
    }

    if (mode < 0)
    {
        lastnz2 = lastnz + 1;
    }
    else
    {
        lastnz2 = 2;
    }

    bits = bits2 = 0;

    /* Calculate number of estimated bits */

    for (k = 0; k < lastnz; k = k + 2)
    {
        t = c + rateFlag;
        if (k > nt_half)
        {
            t += 256;
        }

        codingdata[0] = t;

        a = abs(xq[k]);
        b = abs(xq[k + 1]);
        m = MAX(a, b);

        if (m == 0)
        {
            maxlev = -1;
        }
        else
        {
            maxlev = 29 - (clz_func(MAX(m, 3)) - 1);
        }

        codingdata[1] = maxlev;

        if (mode <= 0) {
            bits = bits + (MIN(a, 1) << 11);
            bits = bits + (MIN(b, 1) << 11);
        }

        lev1 = 0;

        while (MAX(a, b) >= 4)
        {
            pki  = ari_spec_lookup_fl[t + lev1 * 1024];
            bits = bits + ari_spec_bits_fl[pki][16];

            if (lev1 == 0 && mode > 0)
            {
                nbits_lsb += 2;
            }
            else
            {
                bits = bits + 2 * 2048;
            }

            a    = a >> 1;
            b    = b >> 1;
            lev1 = MIN(lev1 + 1, 3);
        }

        pki           = ari_spec_lookup_fl[t + lev1 * 1024];
        sym           = a + 4 * b;
        codingdata[2] = sym;
        codingdata   += 3;
        bits          = bits + ari_spec_bits_fl[pki][sym];

        if (mode > 0)
        {
            a1_msb = abs(xq[k]);
            b1_msb = abs(xq[k + 1]);

            if (lev1 > 0)
            {
                a1_msb = a1_msb >> 1;
                b1_msb = b1_msb >> 1;

                if (a1_msb == 0 && xq[k] != 0)
                {
                    nbits_lsb++;
                }

                if (b1_msb == 0 && xq[k + 1] != 0)
                {
                    nbits_lsb++;
                }
            }

            bits = bits + (MIN(a1_msb, 1) << 11);
            bits = bits + (MIN(b1_msb, 1) << 11);
        }

        if (mode >= 0 && (abs(xq[k]) != 0 || abs(xq[k + 1]) != 0) && bits <= target * 2048)
        {
            lastnz2 = k + 2;
            bits2   = bits;
        }

        lev1 = lev1 - 1;
        
        if (lev1 <= 0)
        {
            t = 1 + (a + b) * (lev1 + 2);
        }
        else
        {
            t = 13 + lev1;
        }

        c = (c & 15) * 16 + t;
    }

    *nbits = (bits + 2047) >> 11; // Exactly same as ceil((LC3_FLOAT)*nbits / 2048.0);

    if (mode >= 0)
    {
        *nbits2 = (bits2 + 2047) >> 11;  //ceil((LC3_FLOAT)*nbits2 / 2048.0);
    }
    else
    {
        *nbits2 = *nbits;
    }

    if (mode > 0)
    {
        *nbits  += nbits_lsb;
        *nbits2 += nbits_lsb;
    }

    /* Truncation of high-frequency coefficients */
    for (i = lastnz2; i <= lastnz; i++)
    {
        xq[i] = 0;
    }

    /* Truncation of LSBs */
    if (mode > 0 && *nbits > target)
    {
        *lsbMode = 1;
    }
    else
    {
        *lsbMode = 0;
    }

    *lastnzout = lastnz2;
}
