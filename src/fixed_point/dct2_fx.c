/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

void idct16_fx(const Word16 *in, Word16 *out)
{
    Dyn_Mem_Deluxe_In(
        Word16 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
        Word16 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;
    );

    a8  = add(mult_r(in[1], 1136), mult_r(in[15], -11529)); /* Sπ/32/√8   -S15π/32/√8 */
    a9  = add(mult_r(in[9], 8956), mult_r(in[7], -7350));   /* S9π/32/√8  -S7π/32/√8  */
    a10 = add(mult_r(in[5], 5461), mult_r(in[11], -10217)); /* S5π/32/√8  -S11π/32/√8 */
    a11 = add(mult_r(in[13], 11086), mult_r(in[3], -3363)); /* S13π/32/√8 -S3π/32/√8  */
    a12 = add(mult_r(in[3], 11086), mult_r(in[13], 3363));  /* C3π/32/√8   C13π/32/√8 */
    a13 = add(mult_r(in[11], 5461), mult_r(in[5], 10217));  /* C11π/32/√8  C5π/32/√8  */
    a14 = add(mult_r(in[7], 8956), mult_r(in[9], 7350));    /* C7π/32/√8   C9π/32/√8  */
    a15 = add(mult_r(in[15], 1136), mult_r(in[1], 11529));  /* C15π/32/√8  Cπ/32/√8   */

    b4  = add(mult_r(in[2], 2260), mult_r(in[14], -11363)); /* Sπ/16/√8  -S7π/16/√8 */
    b5  = add(mult_r(in[10], 9633), mult_r(in[6], -6436));  /* S5π/16/√8 -S3π/16/√8 */
    b6  = add(mult_r(in[6], 9633), mult_r(in[10], 6436));   /* C3π/16/√8  C5π/16/√8 */
    b7  = add(mult_r(in[14], 2260), mult_r(in[2], 11363));  /* C7π/16/√8  Cπ/16/√8  */
    b8  = add(a9, a8);
    b9  = sub(a8, a9);
    b10 = sub(a11, a10);
    b11 = add(a10, a11);
    b12 = add(a13, a12);
    b13 = sub(a12, a13);
    b14 = sub(a15, a14);
    b15 = add(a14, a15);

    a0  = add(mult_r(in[0], 8192), mult_r(in[8], 8192));    /*  Cπ/4/√8  Cπ/4/√8  */
    a1  = add(mult_r(in[8], -8192), mult_r(in[0], 8192));   /* -Cπ/4/√8  Cπ/4/√8  */
    a2  = add(mult_r(in[4], 4433), mult_r(in[12], -10703)); /*  Sπ/8/√8 -S3π/8/√8 */
    a3  = add(mult_r(in[12], 4433), mult_r(in[4], 10703));  /*  C3π/8/√8 Cπ/8/√8  */
    a4  = add(b5, b4);
    a5  = sub(b4, b5);
    a6  = sub(b7, b6);
    a7  = add(b6, b7);
    a8  = b8;                                            move16();
    a9  = add(mult_r(b9, -30274), mult_r(b14, 12540));   /* -Cπ/8  C3π/8 */
    a10 = add(mult_r(b10, -12540), mult_r(b13, -30274)); /* -Sπ/8 -S3π/8 */
    a11 = b11;                                           move16();
    a12 = b12;                                           move16();
    a13 = add(mult_r(b13, 12540), mult_r(b10, -30274));  /* C3π/8 -Cπ/8 */
    a14 = add(mult_r(b14, 30274), mult_r(b9, 12540));    /* S3π/8  Sπ/8 */
    a15 = b15;                                           move16();

    b0  = add(a3, a0);
    b1  = add(a2, a1);
    b2  = sub(a1, a2);
    b3  = sub(a0, a3);
    b4  = a4;                                         move16();
    b5  = add(mult_r(a5, -23170), mult_r(a6, 23170)); /* -Cπ/4 Cπ/4 */
    b6  = add(mult_r(a6, 23170), mult_r(a5, 23170));  /*  Cπ/4 Cπ/4 */
    b7  = a7;                                         move16();
    b8  = add(a11, a8);
    b9  = add(a10, a9);
    b10 = sub(a9, a10);
    b11 = sub(a8, a11);
    b12 = sub(a15, a12);
    b13 = sub(a14, a13);
    b14 = add(a13, a14);
    b15 = add(a12, a15);

    a0  = add(b7, b0);
    a1  = add(b6, b1);
    a2  = add(b5, b2);
    a3  = add(b4, b3);
    a4  = sub(b3, b4);
    a5  = sub(b2, b5);
    a6  = sub(b1, b6);
    a7  = sub(b0, b7);
    a10 = add(mult_r(b10, -23170), mult_r(b13, 23170)); /* -Cπ/4 Cπ/4 */
    a11 = add(mult_r(b11, -23170), mult_r(b12, 23170)); /* -Cπ/4 Cπ/4 */
    a12 = add(mult_r(b12, 23170), mult_r(b11, 23170));  /*  Cπ/4 Cπ/4 */
    a13 = add(mult_r(b13, 23170), mult_r(b10, 23170));  /*  Cπ/4 Cπ/4 */

    out[0]  = add(b15, a0); move16();
    out[1]  = add(b14, a1); move16();
    out[2]  = add(a13, a2); move16();
    out[3]  = add(a12, a3); move16();
    out[4]  = add(a11, a4); move16();
    out[5]  = add(a10, a5); move16();
    out[6]  = add(b9, a6);  move16();
    out[7]  = add(b8, a7);  move16();
    out[8]  = sub(a7, b8);  move16();
    out[9]  = sub(a6, b9);  move16();
    out[10] = sub(a5, a10); move16();
    out[11] = sub(a4, a11); move16();
    out[12] = sub(a3, a12); move16();
    out[13] = sub(a2, a13); move16();
    out[14] = sub(a1, b14); move16();
    out[15] = sub(a0, b15); move16();

    Dyn_Mem_Deluxe_Out();
}

void dct32_fx(const Word32 *in, Word32 *out)
{
    Dyn_Mem_Deluxe_In(Word32 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
                      Word32 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;);

    a0  = L_add(in[15], in[0]);
    a1  = L_add(in[14], in[1]);
    a2  = L_add(in[13], in[2]);
    a3  = L_add(in[12], in[3]);
    a4  = L_add(in[11], in[4]);
    a5  = L_add(in[10], in[5]);
    a6  = L_add(in[9], in[6]);
    a7  = L_add(in[8], in[7]);
    a10 = L_sub(in[5], in[10]);
    a11 = L_sub(in[4], in[11]);
    a12 = L_sub(in[3], in[12]);
    a13 = L_sub(in[2], in[13]);

    b0  = L_add(a7, a0);
    b1  = L_add(a6, a1);
    b2  = L_add(a5, a2);
    b3  = L_add(a4, a3);
    b4  = L_sub(a3, a4);
    b5 = L_sub_sat(a2, a5);
    b6 = L_sub_sat(a1, a6);
    b7  = L_sub(a0, a7);
    b8  = L_sub(in[7], in[8]);
    b9  = L_sub(in[6], in[9]);
    b10 = L_add(Mpy_32_16(a10, -23170), Mpy_32_16(a13, 23170)); /* -Cπ/4 Cπ/4 */
    b11 = L_add(Mpy_32_16(a11, -23170), Mpy_32_16(a12, 23170)); /* -Cπ/4 Cπ/4 */
    b12 = L_add(Mpy_32_16(a12, 23170), Mpy_32_16(a11, 23170));  /*  Cπ/4 Cπ/4 */
    b13 = L_add(Mpy_32_16(a13, 23170), Mpy_32_16(a10, 23170));  /*  Cπ/4 Cπ/4 */
    b14 = L_sub(in[1], in[14]);
    b15 = L_sub(in[0], in[15]);

    a0 = L_add(b3, b0);
    a1 = L_add(b2, b1);
    a2 = L_sub(b1, b2);
    a3 = L_sub_sat(b0, b3);
    a4 = b4;
    move16();
    a5 = L_add(Mpy_32_16(b5, -23170), Mpy_32_16(b6, 23170)); /* -Cπ/4 Cπ/4 */
    a6 = L_add(Mpy_32_16(b6, 23170), Mpy_32_16(b5, 23170));  /*  Cπ/4 Cπ/4 */
    a7 = b7;
    move16();
    a8  = L_add(b11, b8);
    a9  = L_add(b10, b9);
    a10 = L_sub(b9, b10);
    a11 = L_sub(b8, b11);
    a12 = L_sub(b15, b12);
    a13 = L_sub(b14, b13);
    a14 = L_add(b13, b14);
    a15 = L_add(b12, b15);

    out[0] = L_add(Mpy_32_16(a0, 8192), Mpy_32_16(a1, 8192));
    move16(); /*  Cπ/4/√8   Cπ/4/√8  */
    out[8] = L_add(Mpy_32_16(a1, -8192), Mpy_32_16(a0, 8192));
    move16(); /* -Cπ/4/√8   Cπ/4/√8  */
    out[4] = L_add(Mpy_32_16(a2, 4433), Mpy_32_16(a3, 10703));
    move16(); /*  Sπ/8/√8   Cπ/8/√8  */
    out[12] = L_add(Mpy_32_16(a3, 4433), Mpy_32_16(a2, -10703));
    move16(); /*  C3π/8/√8 -S3π/8/√8 */
    b4 = L_add(a5, a4);
    b5 = L_sub(a4, a5);
    b6 = L_sub_sat(a7, a6);
    b7 = L_add(a6, a7);
    b8 = a8;
    move16();
    b9  = L_add(Mpy_32_16(a9, -30274), Mpy_32_16(a14, 12540));   /* -Cπ/8  Sπ/8 */
    b10 = L_add(Mpy_32_16(a10, -12540), Mpy_32_16(a13, -30274)); /* -Sπ/8 -Cπ/8 */
    b11 = a11;
    move16();
    b12 = a12;
    move16();
    b13 = L_add(Mpy_32_16(a13, 12540), Mpy_32_16(a10, -30274)); /* C3π/8 -S3π/8 */
    b14 = L_add(Mpy_32_16(a14, 30274), Mpy_32_16(a9, 12540));   /* S3π/8  C3π/8 */
    b15 = a15;
    move16();

    out[2] = L_add(Mpy_32_16(b4, 2260), Mpy_32_16(b7, 11363));
    move16(); /* Sπ/16/√8   Cπ/16/√8  */
    out[10] = L_add(Mpy_32_16(b5, 9633), Mpy_32_16(b6, 6436));
    move16(); /* S5π/16/√8  C5π/16/√8 */
    out[6] = L_add(Mpy_32_16(b6, 9633), Mpy_32_16(b5, -6436));
    move16(); /* C3π/16/√8 -S3π/16/√8 */
    out[14] = L_add(Mpy_32_16(b7, 2260), Mpy_32_16(b4, -11363));
    move16(); /* C7π/16/√8 -S7π/16/√8 */

    a8  = L_add_sat(b9, b8);
    a9  = L_sub_sat(b8, b9);
    a10 = L_sub_sat(b11, b10);
    a11 = L_add_sat(b10, b11);
    a12 = L_add_sat(b13, b12);
    a13 = L_sub_sat(b12, b13);
    a14 = L_sub_sat(b15, b14);
    a15 = L_add_sat(b14, b15);

    out[1] = L_add(Mpy_32_16(a8, 1136), Mpy_32_16(a15, 11529));
    move16(); /* Sπ/32/√8    Cπ/32/√8   */
    out[9] = L_add(Mpy_32_16(a9, 8956), Mpy_32_16(a14, 7350));
    move16(); /* S9π/32/√8   C9π/32/√8  */
    out[5] = L_add(Mpy_32_16(a10, 5461), Mpy_32_16(a13, 10217));
    move16(); /* S5π/32/√8   C5π/32/√8  */
    out[13] = L_add(Mpy_32_16(a11, 11086), Mpy_32_16(a12, 3363));
    move16(); /* S13π/32/√8  C13π/32/√8 */
    out[3] = L_add(Mpy_32_16(a12, 11086), Mpy_32_16(a11, -3363));
    move16(); /* C3π/32/√8  -S3π/32/√8  */
    out[11] = L_add(Mpy_32_16(a13, 5461), Mpy_32_16(a10, -10217));
    move16(); /* C11π/32/√8 -S11π/32/√8 */
    out[7] = L_add(Mpy_32_16(a14, 8956), Mpy_32_16(a9, -7350));
    move16(); /* C7π/32/√8  -S7π/32/√8  */
    out[15] = L_add(Mpy_32_16(a15, 1136), Mpy_32_16(a8, -11529));
    move16(); /* C15π/32/√8 -S15/32/√8  */

    Dyn_Mem_Deluxe_Out();
}

void idct32_fx(const Word32 *in, Word32 *out)
{
    Dyn_Mem_Deluxe_In(Word32 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
                      Word32 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;);

    a8  = L_add(Mpy_32_16(in[1], 1136), Mpy_32_16(in[15], -11529)); /* Sπ/32/√8   -S15π/32/√8 */
    a9  = L_add(Mpy_32_16(in[9], 8956), Mpy_32_16(in[7], -7350));   /* S9π/32/√8  -S7π/32/√8  */
    a10 = L_add(Mpy_32_16(in[5], 5461), Mpy_32_16(in[11], -10217)); /* S5π/32/√8  -S11π/32/√8 */
    a11 = L_add(Mpy_32_16(in[13], 11086), Mpy_32_16(in[3], -3363)); /* S13π/32/√8 -S3π/32/√8  */
    a12 = L_add(Mpy_32_16(in[3], 11086), Mpy_32_16(in[13], 3363));  /* C3π/32/√8   C13π/32/√8 */
    a13 = L_add(Mpy_32_16(in[11], 5461), Mpy_32_16(in[5], 10217));  /* C11π/32/√8  C5π/32/√8  */
    a14 = L_add(Mpy_32_16(in[7], 8956), Mpy_32_16(in[9], 7350));    /* C7π/32/√8   C9π/32/√8  */
    a15 = L_add(Mpy_32_16(in[15], 1136), Mpy_32_16(in[1], 11529));  /* C15π/32/√8  Cπ/32/√8   */

    b4  = L_add(Mpy_32_16(in[2], 2260), Mpy_32_16(in[14], -11363)); /* Sπ/16/√8  -S7π/16/√8 */
    b5  = L_add(Mpy_32_16(in[10], 9633), Mpy_32_16(in[6], -6436));  /* S5π/16/√8 -S3π/16/√8 */
    b6  = L_add(Mpy_32_16(in[6], 9633), Mpy_32_16(in[10], 6436));   /* C3π/16/√8  C5π/16/√8 */
    b7  = L_add(Mpy_32_16(in[14], 2260), Mpy_32_16(in[2], 11363));  /* C7π/16/√8  Cπ/16/√8  */
    b8  = L_add(a9, a8);
    b9  = L_sub(a8, a9);
    b10 = L_sub(a11, a10);
    b11 = L_add(a10, a11);
    b12 = L_add(a13, a12);
    b13 = L_sub(a12, a13);
    b14 = L_sub(a15, a14);
    b15 = L_add(a14, a15);

    a0 = L_add(Mpy_32_16(in[0], 8192), Mpy_32_16(in[8], 8192)); /*  Cπ/4/√8  Cπ/4/√8  */
    a1 = L_add(Mpy_32_16(in[8], -8192), Mpy_32_16(in[0], 8192)); /* -Cπ/4/√8  Cπ/4/√8  */
    a2 = L_add(Mpy_32_16(in[4], 4433), Mpy_32_16(in[12], -10703)); /*  Sπ/8/√8 -S3π/8/√8 */
    a3 = L_add(Mpy_32_16(in[12], 4433), Mpy_32_16(in[4], 10703));  /*  C3π/8/√8 Cπ/8/√8  */
    a4 = L_add(b5, b4);
    a5 = L_sub(b4, b5);
    a6 = L_sub(b7, b6);
    a7 = L_add(b6, b7);
    a8 = b8;
    move32();
    a9  = L_add(Mpy_32_16(b9, -30274), Mpy_32_16(b14, 12540)); /* -Cπ/8  C3π/8 */
    a10 = L_add(Mpy_32_16(b10, -12540), Mpy_32_16(b13, -30274)); /* -Sπ/8 -S3π/8 */
    a11 = b11;
    move32();
    a12 = b12;
    move32();
    a13 = L_add(Mpy_32_16(b13, 12540), Mpy_32_16(b10, -30274)); /* C3π/8 -Cπ/8 */
    a14 = L_add(Mpy_32_16(b14, 30274), Mpy_32_16(b9, 12540));   /* S3π/8  Sπ/8 */
    a15 = b15;
    move32();

    b0 = L_add(a3, a0);
    b1 = L_add(a2, a1);
    b2 = L_sub(a1, a2);
    b3 = L_sub(a0, a3);
    b4 = a4;
    move32();
    b5 = L_add(Mpy_32_16(a5, -23170), Mpy_32_16(a6, 23170)); /* -Cπ/4 Cπ/4 */
    b6 = L_add(Mpy_32_16(a6, 23170), Mpy_32_16(a5, 23170));  /*  Cπ/4 Cπ/4 */
    b7 = a7;
    move32();
    b8  = L_add(a11, a8);
    b9  = L_add(a10, a9);
    b10 = L_sub(a9, a10);
    b11 = L_sub(a8, a11);
    b12 = L_sub(a15, a12);
    b13 = L_sub(a14, a13);
    b14 = L_add(a13, a14);
    b15 = L_add(a12, a15);

    a0  = L_add(b7, b0);
    a1  = L_add(b6, b1);
    a2  = L_add(b5, b2);
    a3  = L_add(b4, b3);
    a4  = L_sub(b3, b4);
    a5  = L_sub(b2, b5);
    a6  = L_sub(b1, b6);
    a7  = L_sub(b0, b7);
    a10 = L_add(Mpy_32_16(b10, -23170), Mpy_32_16(b13, 23170)); /* -Cπ/4 Cπ/4 */
    a11 = L_add(Mpy_32_16(b11, -23170), Mpy_32_16(b12, 23170)); /* -Cπ/4 Cπ/4 */
    a12 = L_add(Mpy_32_16(b12, 23170), Mpy_32_16(b11, 23170));  /*  Cπ/4 Cπ/4 */
    a13 = L_add(Mpy_32_16(b13, 23170), Mpy_32_16(b10, 23170));  /*  Cπ/4 Cπ/4 */

    out[0] = L_add(b15, a0);
    move32();
    out[1] = L_add(b14, a1);
    move32();
    out[2] = L_add(a13, a2);
    move32();
    out[3] = L_add(a12, a3);
    move32();
    out[4] = L_add(a11, a4);
    move32();
    out[5] = L_add(a10, a5);
    move32();
    out[6] = L_add(b9, a6);
    move32();
    out[7] = L_add(b8, a7);
    move32();
    out[8] = L_sub(a7, b8);
    move32();
    out[9] = L_sub(a6, b9);
    move32();
    out[10] = L_sub(a5, a10);
    move32();
    out[11] = L_sub(a4, a11);
    move32();
    out[12] = L_sub(a3, a12);
    move32();
    out[13] = L_sub(a2, a13);
    move32();
    out[14] = L_sub(a1, b14);
    move32();
    out[15] = L_sub(a0, b15);
    move32();

    Dyn_Mem_Deluxe_Out();
}

void idct32_32_fx(const Word32 *in, Word32 *out)
{
    Dyn_Mem_Deluxe_In(Word32 a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15;
                      Word32 b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15;);

    /*round(sin(pi*(1:16)/32)/sqrt(8)*2^31) = 
    74419526	148122351	220398677	290552444	357908031	421816769	481663180	536870912
    586908283	631293407	669598830	701455651	726557070	744661347	755594128	759250125
    */
    /*ound(cos(pi*(1:16)/32)/sqrt(8)*2^31) = 
    755594128	744661347	726557070	701455651	669598830	631293407	586908283	536870912
    481663180	421816769	357908031	290552444	220398677	148122351	74419526	0
    */

    a8  = L_add(Mpy_32_32(in[1], 74419526), Mpy_32_32(in[15], -755594128)); /* Sπ/32/√8   -S15π/32/√8 */
    a9  = L_add(Mpy_32_32(in[9], 586908283), Mpy_32_32(in[7], -481663180)); /* S9π/32/√8  -S7π/32/√8  */
    a10 = L_add(Mpy_32_32(in[5], 357908031), Mpy_32_32(in[11], -669598830)); /* S5π/32/√8  -S11π/32/√8 */
    a11 = L_add(Mpy_32_32(in[13], 726557070), Mpy_32_32(in[3], -220398677)); /* S13π/32/√8 -S3π/32/√8  */

    a12 = L_add(Mpy_32_32(in[3], 726557070), Mpy_32_32(in[13], 220398677)); /* C3π/32/√8   C13π/32/√8 */
    a13 = L_add(Mpy_32_32(in[11], 357908031), Mpy_32_32(in[5], 669598830)); /* C11π/32/√8  C5π/32/√8  */
    a14 = L_add(Mpy_32_32(in[7], 586908283), Mpy_32_32(in[9], 481663180));  /* C7π/32/√8   C9π/32/√8  */
    a15 = L_add(Mpy_32_32(in[15], 74419526), Mpy_32_32(in[1], 755594128));  /* C15π/32/√8  Cπ/32/√8   */

    b4  = L_add(Mpy_32_32(in[2], 148122351), Mpy_32_32(in[14], -744661347)); /* Sπ/16/√8  -S7π/16/√8 */
    b5  = L_add(Mpy_32_32(in[10], 631293407), Mpy_32_32(in[6], -421816769)); /* S5π/16/√8 -S3π/16/√8 */
    b6  = L_add(Mpy_32_32(in[6], 631293407), Mpy_32_32(in[10], 421816769));  /* C3π/16/√8  C5π/16/√8 */
    b7  = L_add(Mpy_32_32(in[14], 148122351), Mpy_32_32(in[2], 744661347));  /* C7π/16/√8  Cπ/16/√8  */
    b8  = L_add(a9, a8);
    b9  = L_sub(a8, a9);
    b10 = L_sub(a11, a10);
    b11 = L_add(a10, a11);
    b12 = L_add(a13, a12);
    b13 = L_sub(a12, a13);
    b14 = L_sub(a15, a14);
    b15 = L_add(a14, a15);

    a0 = L_add(Mpy_32_32(in[0], 536870912), Mpy_32_32(in[8], 536870912)); /*  Cπ/4/√8  Cπ/4/√8  */
    a1 = L_add(Mpy_32_32(in[8], -536870912), Mpy_32_32(in[0], 536870912)); /* -Cπ/4/√8  Cπ/4/√8  */
    a2 = L_add(Mpy_32_32(in[4], 290552444), Mpy_32_32(in[12], -701455651)); /*  Sπ/8/√8 -S3π/8/√8 */
    a3 = L_add(Mpy_32_32(in[12], 290552444), Mpy_32_32(in[4], 701455651));  /*  C3π/8/√8 Cπ/8/√8  */
    a4 = L_add(b5, b4);
    a5 = L_sub(b4, b5);
    a6 = L_sub(b7, b6);
    a7 = L_add(b6, b7);
    a8 = b8;
    move32();
    a9  = L_add(Mpy_32_32(b9, -1984016189), Mpy_32_32(b14, 821806413)); /* -Cπ/8  C3π/8 */
    a10 = L_add(Mpy_32_32(b10, -821806413), Mpy_32_32(b13, -1984016189)); /* -Sπ/8 -S3π/8 */
    a11 = b11;
    move32();
    a12 = b12;
    move32();
    a13 = L_add(Mpy_32_32(b13, 821806413), Mpy_32_32(b10, -1984016189)); /* C3π/8 -Cπ/8 */
    a14 = L_add(Mpy_32_32(b14, 1984016189), Mpy_32_32(b9, 821806413));   /* S3π/8  Sπ/8 */
    a15 = b15;
    move32();

    b0 = L_add(a3, a0);
    b1 = L_add(a2, a1);
    b2 = L_sub(a1, a2);
    b3 = L_sub(a0, a3);
    b4 = a4;
    move32();
    b5 = L_add(Mpy_32_32(a5, -1518500250), Mpy_32_32(a6, 1518500250)); /* -Cπ/4 Cπ/4 */
    b6 = L_add(Mpy_32_32(a6, 1518500250), Mpy_32_32(a5, 1518500250));  /*  Cπ/4 Cπ/4 */
    b7 = a7;
    move32();
    b8  = L_add(a11, a8);
    b9  = L_add(a10, a9);
    b10 = L_sub(a9, a10);
    b11 = L_sub(a8, a11);
    b12 = L_sub(a15, a12);
    b13 = L_sub(a14, a13);
    b14 = L_add(a13, a14);
    b15 = L_add(a12, a15);

    a0  = L_add(b7, b0);
    a1  = L_add(b6, b1);
    a2  = L_add(b5, b2);
    a3  = L_add(b4, b3);
    a4  = L_sub(b3, b4);
    a5  = L_sub(b2, b5);
    a6  = L_sub(b1, b6);
    a7  = L_sub(b0, b7);
    a10 = L_add(Mpy_32_32(b10, -1518500250), Mpy_32_32(b13, 1518500250)); /* -Cπ/4 Cπ/4 */
    a11 = L_add(Mpy_32_32(b11, -1518500250), Mpy_32_32(b12, 1518500250)); /* -Cπ/4 Cπ/4 */
    a12 = L_add(Mpy_32_32(b12, 1518500250), Mpy_32_32(b11, 1518500250));  /*  Cπ/4 Cπ/4 */
    a13 = L_add(Mpy_32_32(b13, 1518500250), Mpy_32_32(b10, 1518500250));  /*  Cπ/4 Cπ/4 */

    out[0] = L_add(b15, a0);
    move32();
    out[1] = L_add(b14, a1);
    move32();
    out[2] = L_add(a13, a2);
    move32();
    out[3] = L_add(a12, a3);
    move32();
    out[4] = L_add(a11, a4);
    move32();
    out[5] = L_add(a10, a5);
    move32();
    out[6] = L_add(b9, a6);
    move32();
    out[7] = L_add(b8, a7);
    move32();
    out[8] = L_sub(a7, b8);
    move32();
    out[9] = L_sub(a6, b9);
    move32();
    out[10] = L_sub(a5, a10);
    move32();
    out[11] = L_sub(a4, a11);
    move32();
    out[12] = L_sub(a3, a12);
    move32();
    out[13] = L_sub(a2, a13);
    move32();
    out[14] = L_sub(a1, b14);
    move32();
    out[15] = L_sub(a0, b15);
    move32();

    Dyn_Mem_Deluxe_Out();
}
