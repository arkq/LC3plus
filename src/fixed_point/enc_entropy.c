/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static Word32 ac_enc_mux_st2VQ_cws(                    /*  o:  max 25 bits total codeword */
                                   const Word32 L_szA, /*  i:  max 22 bits */
                                   const Word32 L_szB, /*  i:  max 4 bits  */
                                   const Word32 L_cwA, const Word32 L_cwB);

void processEncoderEntropy(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits, Word16 targetBytes,
                           Word16 L_spec, Word16 BW_cutoff_bits, Word16 tns_numfilters,
                           Word16 lsbMode, Word16 lastnz, Word16 *tns_order, Word16 fac_ns_idx, Word16 gg_idx,
                           Word16 BW_cutoff_idx, Word16 *ltpf_idx, Word32 *L_scf_idx, Word16 bfi_ext, Word16 fs_idx
#ifdef CR9_C_ADD_1p25MS
                           , LC3PLUS_FrameDuration frame_dms, Word16* Tx_ltpf
#endif
                           )
{
    Word16  tmp;
    Word32  L_tmp;
    Word16  submode_LSB, submode_MSB, gain_MSBs;
    Word32  L_gain_LSB;
    Counter n;
    UWord8 *ptr;

    Word16  lastnzTrigger[5] = {63, 127, 127, 255, 255};

#ifdef DYNMEM_COUNT
    struct _dynmem
    {
        Word16  tmp;
        Word32  L_tmp;
        Word16  submode_LSB, submode_MSB, gain_MSBs;
        Word32  L_gain_LSB;
        Counter n;
        UWord8 *ptr;
        Word16  lastnzTrigger[5];
    };
    Dyn_Mem_In("processEncoderEntropy", sizeof(struct _dynmem));
#endif

    UNUSED(L_gain_LSB); UNUSED(gain_MSBs); UNUSED(submode_MSB); UNUSED(submode_LSB); UNUSED(L_tmp);

    /* Init */
    *bp_side   = shr_pos(sub(nbbits, 1), 3);
    *mask_side = shl(1, sub(8, sub(nbbits, shl_pos(*bp_side, 3))));
    ptr        = bytes;

    basop_memset(bytes, 0, targetBytes * sizeof(*bytes));

    /* Cutoff-detection */
    IF (BW_cutoff_bits > 0)
    {
        write_indice_backward(ptr, bp_side, mask_side, BW_cutoff_idx, BW_cutoff_bits);
    }

    /* Encode last non-zero tuple */
    tmp = getLastNzBits_fx(L_spec);

    IF (sub(bfi_ext, 1) == 0)
    {
        write_indice_backward(ptr, bp_side, mask_side, lastnzTrigger[fs_idx], tmp);
    }
    ELSE
    {
        write_indice_backward(ptr, bp_side, mask_side, sub(shr_pos(lastnz, 1), 1), tmp);
    }

    /* Mode bit */
    write_bit_backward(ptr, bp_side, mask_side, lsbMode);

    /* Encode global-gain */
    write_indice_backward(ptr, bp_side, mask_side, gg_idx, 8);

    /* TNS on/off flag */
    FOR (n = 0; n < tns_numfilters; n++)
    {
        write_bit_backward(ptr, bp_side, mask_side, s_min(tns_order[n], 1));
    }

    /* LTPF on/off*/
#ifdef CR9_C_ADD_1p25MS
#ifdef NEW_SIGNALLING_SCHEME_1p25
    IF( sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0 )
    {
        writeLtpData_fx(ptr, bp_side, mask_side, ltpf_idx, Tx_ltpf); /* LTP-flag and LTPFflag  and interleaved lag  */
    }
    ELSE
    {
        write_bit_backward(ptr, bp_side, mask_side, ltpf_idx[0]); /* LTP-flag only */
    }
#endif
#else
    write_indice_backward(ptr, bp_side, mask_side, ltpf_idx[0], 1);
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
    /* Encode SCF VQ parameters */
    writeSNSData_fx(bytes, bp_side, mask_side, frame_dms, L_scf_idx);
#else
    /* Encode SCF VQ parameters - 1st stage (10 bits) */
    write_indice_backward(ptr, bp_side, mask_side, extract_l(L_scf_idx[0]), 5); /* stage1 LF   5 bits */
    write_indice_backward(ptr, bp_side, mask_side, extract_l(L_scf_idx[1]), 5); /* stage1 HF   5  bits  */

    /* Encode SCF VQ parameters - 2nd stage side-info (3-4 bits) */
    submode_MSB = shr_pos(extract_l(L_scf_idx[2]), 1);        /*  explicit tx */
    write_bit_backward(ptr, bp_side, mask_side, submode_MSB); /* submode MSB  1 explicit bit */
    submode_LSB = s_and(extract_l(L_scf_idx[2]), 0x1);        /* for joint coding with shapeCw */
    gain_MSBs   = extract_l(L_scf_idx[3]);                    /* all gain bits */
    L_gain_LSB  = L_and(L_scf_idx[3], 0x1L);
    gain_MSBs   = shr(gain_MSBs, sns_gainLSBbits[L_scf_idx[2]]);

    ASSERT(gain_MSBs >= 0 && gain_MSBs < (1 << sns_gainMSBbits[L_scf_idx[2]])); /* ASSERT  max 2 MSB(s) in gain bits */

    write_indice_backward(ptr, bp_side, mask_side, gain_MSBs,
                          sns_gainMSBbits[L_scf_idx[2]]);                 /* adjgain or MSBs of adjGains   1-2 bits  */
    write_bit_backward(ptr, bp_side, mask_side, extract_l(L_scf_idx[4])); /*  shape  LS 1 bit */

    /* Encode SCF VQ parameters - 2nd stage data (24-25 bits) */
    IF (submode_MSB == 0)
    { /* regular,regular_lf*/
        ASSERT(submode_MSB == 0);

        L_tmp = L_add(L_gain_LSB, 0); /* gain-LSB 0,1 for regular_lf,  offset is 0 */
        if (submode_LSB == 0)
        {
            L_tmp = L_add(L_scf_idx[6],
                          sns_MPVQ_Sz[1][1]); /* shape B pos offset is 2 , upshifted two positions , 0..11 -> 2..13 */
        }
        /* regular mode A,B indexes multiplexed, total 24.x bits  MPVQ codeword section A + codeword for section B */
        L_tmp = ac_enc_mux_st2VQ_cws(sns_MPVQ_Sz[0][0],                               /*  max 21.3  bits*/
                                     UL_addNsD(sns_MPVQ_Sz[0][1], sns_MPVQ_Sz[1][1]), /*  max log2(14)  bits */
                                     L_scf_idx[5] /* shapeA */, L_tmp /*   shapeB joint with adjGainLSB */);
        /* regular mode  mode shape  index   total  1+23.9999 bits    MPVQ codeword  */
        ASSERT(L_tmp < (1L << 25));
        write_indice_backward(ptr, bp_side, mask_side, extract_l(L_tmp), 13);                /*  multiplex 13  bits  */
        write_indice_backward(ptr, bp_side, mask_side, extract_l(L_shr_pos(L_tmp, 13)), 12); /* multiplex 12 bits  */
    }
    ELSE
    { /* outlier near, outlier far */
        ASSERT(submode_MSB == 1);
        L_tmp = L_scf_idx[5];  move32(); /* outlier near  section assumed */
        if (submode_LSB != 0)
        {                                                   /* outl_far */
            L_tmp = L_add(L_shl_pos(L_tmp, 1), L_gain_LSB); /*  add lsb bit of Gain */
            L_tmp = L_add(L_tmp, sns_MPVQ_Sz[2][0]);        /*  outlier far section offset added */
        }

        ASSERT(L_tmp < (1L << 24));
        /* outlier mode shape  index   total  23.8536 ( +  ~.14 ) bits as   MPVQ codeword  */
        write_indice_backward(ptr, bp_side, mask_side, extract_l(L_tmp), 12);            /*  multiplex 12  bits  LSB*/
        write_indice_backward(ptr, bp_side, mask_side, extract_l(L_shr(L_tmp, 12)), 12); /* multiplex 12 bits  MSBs */
    }
#endif

    /* LTPF data */
#ifdef NEW_SIGNALLING_SCHEME_1p25
    test(); test();
    IF( (sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) != 0)  && (ltpf_idx[0] != 0) )
#else
    test();
    IF( ltpf_idx[0] != 0 )
#endif
    {
          write_indice_backward(ptr, bp_side, mask_side, ltpf_idx[1], 1);
          write_indice_backward(ptr, bp_side, mask_side, ltpf_idx[2], 9);
    }


    /* Encoder noise-fac */
    write_indice_backward(ptr, bp_side, mask_side, fac_ns_idx, 3);

#ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#endif
}

static __forceinline Word32
ac_enc_mux_st2VQ_cws(                    /*  o:  max 25 bits total codeword */
                     const Word32 L_szA, /*  i:  max 22 bits */
                     const Word32 L_szB, /*  i:  max 4 bits  0..13  */
                     const Word32 L_cwA,
                     const Word32 L_cwB) /*  [0..13} corresponding to gains{0,1}, shapeB{0..11}  or */
{

    Word32 L_cwTx;
    /* L_cw_tx =   L_cwB(21.z bits) * L_szA(3.y bits)   + L_cwA(21.x bits)); */
    L_cwTx = (Word32)UL_Mpy_32_32(
        (UWord32)L_cwB, (UWord32)L_szA); /* non-fractional 16x32 -> 32  may possibly also be used if available */
    L_cwTx = L_add(L_cwTx, L_cwA);

    ASSERT((L_szA * L_szB) <= 1 << 25);          /* multiplexing only allowed up to 25 bits  (+ leading sign)  */
    ASSERT(L_cwTx >= 0 && L_cwTx <= 0x01ffFFff); /*  max 25 bits allowed */
    UNUSED(L_szB);

    return L_cwTx;
}

#ifdef NEW_SIGNALLING_SCHEME_1p25
void writeLtpData_fx(
    UWord8 *ptr,
    Word16 *bp_side,
    Word16 *mask_side,
    Word16* ltpf_idx,
    Word16* Tx_ltpf
)
{
    Dyn_Mem_Deluxe_In(
    Word16 tmp;
    Word16 bitsTx;);

    tmp = s_min(*Tx_ltpf, 1);           /*phaseA==0, phaseB==1*/
    test();
    IF(ltpf_idx[0] == 0)
    {
        write_indice_backward(ptr, bp_side, mask_side, 0, 2);  /* "00" */

        *Tx_ltpf = 0;
        /*  *Tx_ltpf A/B state forced to zero or  kept at zero */
        /*  decoder will discard any sofar in phase received phaseA MSB bits */
    }
    ELSE IF(ltpf_idx[1] == 0)
    {
        /* no current LTPF activation,
           lag transmitted for PLC,  or for next frame LTPF activation */
        ASSERT(ltpf_idx[0] != 0);

        /* A "010"  3 bits Hdr transmitted */
        /* B "011"  3 bits Hdr transmitted*/

        write_indice_backward(ptr, bp_side, mask_side, 1, 2); /* "01"*/
        write_indice_backward(ptr, bp_side, mask_side, tmp, 1);  /* phase A or phaseB  */

        test();
        IF(*Tx_ltpf == 0)
        {  /* phase A transmission */
            ASSERT(tmp == 0);
            ASSERT((ltpf_idx[2] & ~(0x01ff)) == 0); /* only 9 bits info allowed  within  ltpf_idx[2] */
            tmp = shr_pos(ltpf_idx[2], 5);  /* shift_out LSBS, send 4 MSBs  */
            *Tx_ltpf = s_or(0x200, ltpf_idx[2]);  move16(); /*  remember full lag,  as phaseB sentinel in bit 10,   */
        }
        ELSE
        {   /* phase B */
            ASSERT(tmp == 1);
            ASSERT(*Tx_ltpf > 511);  /*sentinel in b10 should have been set in previous phaseA frame */
            tmp = shr_pos(s_and(*Tx_ltpf, 0x001f), 1);  /* B send 4* LSBs,  1 bit truncated */
            *Tx_ltpf = 0;         move16();             /*  clear sentinel in b10  and the old  remebered lag value */
        }
        write_indice_backward(ptr, bp_side, mask_side, tmp, 4); /* 4 bits  lag info Tx, when LTPF is deactivated  */
    }
    ELSE
    {   /* LTPF activated */
        /* A "10"  2 bits Hdr */
        /* B "11"  2 bits Hdr */
        ASSERT(ltpf_idx[0] != 0 && ltpf_idx[1] != 0);

        tmp = s_or(0x02, tmp);
        write_indice_backward(ptr, bp_side, mask_side, tmp, 2);
        test();
        IF(*Tx_ltpf == 0)
        {
            bitsTx = 4; move16();
            ASSERT((ltpf_idx[2] & ~(0x01ff)) == 0); /* only 9 bits info allowed  within  ltpf_idx[2] */
            tmp = shr_pos(ltpf_idx[2], 5);          /* shift away 5 LBS,  LTPF active, send phaseA 4 MSBs  */

            *Tx_ltpf = s_or(0x0200, ltpf_idx[2]);  /*  remember full lag in state *Tx_ltpf,  add phaseB Tx state sentinel in bit 10   */
        }
        ELSE
        {
            bitsTx = 5;          move16();
            tmp = s_and(*Tx_ltpf, 0x001f);        /* LTPF active   B send 5 LSBs , full regular resolution  */
            *Tx_ltpf = 0;         move16();         /* clear sentinel in b10  and also the old  remebered lag value */
        }

        write_indice_backward(ptr, bp_side, mask_side, tmp, bitsTx);  /* ltp==1 , ltpf==1,  4(A,MSBs)  or 5(b, LSBs) bits */
    }
    Dyn_Mem_Deluxe_Out();
}


void writeSNSData_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, LC3PLUS_FrameDuration frame_dms, Word32* L_scf_idx)
{
    Dyn_Mem_Deluxe_In(
        Word16 submode_LSB, submode_MSB, gain_MSBs;
    Word32 L_gain_LSB;
    Word16 write_legacy_sns_vq_bits;
    Word16 tmp;
    Word32 L_tmp;
    Word16 aux_idx;         /* aux value: location the LS, or the s0 sign bit     */
    Word16 shape_idx;        /* st2 shape 0 .. 5  , where [2.3.4.5] are fixed FESS shapes */
    Word16 fixenv_shape_idx;
    Word16 gain_idx;         /* idx of 2-3 bits valued  gains */
    Word16 n5k;              /* number of unit pulses k for the width N=5 LF region of splitLF */
    );
 
    write_legacy_sns_vq_bits = 1;
    if (sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0) /*  9,10,29/30 */
    {
        write_legacy_sns_vq_bits = 0;
    }
 
    test();
    IF(write_legacy_sns_vq_bits != 0)  /* 10 bits (5+5) for stage1  + a total of 28 bits for stage2 */
    {
        /* Encode SCF VQ parameters - 1st stage (10 bits) */
        write_indice_backward(bytes, bp_side, mask_side, extract_l(L_scf_idx[0]), 5); /* stage1 LF   5 bits */
        write_indice_backward(bytes, bp_side, mask_side, extract_l(L_scf_idx[1]), 5); /* stage1 HF   5  bits  */

        /* Encode SCF VQ parameters - 2nd stage side-info (3-4 bits) */
        submode_MSB = shr_pos(extract_l(L_scf_idx[2]), 1); /*  explicit tx */
        write_bit_backward(bytes, bp_side, mask_side, submode_MSB);    /* submode MSB  1 explicit bit */
        submode_LSB = s_and(extract_l(L_scf_idx[2]), 0x1); /* for joint coding with shapeCw */
        gain_MSBs = extract_l(L_scf_idx[3]);                 /* all gain bits */
        L_gain_LSB = L_and(L_scf_idx[3], 0x1L);
        gain_MSBs = shr(gain_MSBs, sns_gainLSBbits[L_scf_idx[2]]);

        ASSERT(gain_MSBs >= 0 && gain_MSBs < (1 << sns_gainMSBbits[L_scf_idx[2]])); /* ASSERT  max 2 MSB(s) in gain bits */

        write_indice_backward(bytes, bp_side, mask_side, gain_MSBs,
            sns_gainMSBbits[L_scf_idx[2]]);             /* adjgain or MSBs of adjGains   1-2 bits  */
        write_bit_backward(bytes, bp_side, mask_side, extract_l(L_scf_idx[4])); /*  shape  LS 1 bit */

        /* Encode SCF VQ parameters - 2nd stage data (24-25 bits) */
        test();
        IF(submode_MSB == 0)
        { /* regular,regular_lf*/
            ASSERT(submode_MSB == 0);

            L_tmp = L_add(L_gain_LSB, 0); /* gain-LSB 0,1 for regular_lf,  offset is 0 */
            if (submode_LSB == 0)
            {
                L_tmp = L_add(L_scf_idx[6],
                    sns_MPVQ_Sz[1][1]); /* shape B pos offset is 2 , upshifted two positions , 0..11 -> 2..13 */
            }
            /* regular mode A,B indexes multiplexed, total 24.x bits  MPVQ codeword section A + codeword for section B */
            L_tmp = ac_enc_mux_st2VQ_cws(sns_MPVQ_Sz[0][0],                                 /*  max 21.3  bits*/
                UL_addNsD(sns_MPVQ_Sz[0][1], sns_MPVQ_Sz[1][1]), /*  max log2(14)  bits */
                L_scf_idx[5] /* shapeA */, L_tmp /*   shapeB joint with adjGainLSB */);
            /* regular mode  mode shape  index   total  1+23.9999 bits    MPVQ codeword  */
            ASSERT(L_tmp < (1L << 25));
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_tmp), 13);                /*  multiplex 13  bits  */
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr_pos(L_tmp, 13)), 12); /* multiplex 12 bits  */
        }
        ELSE
        { /* outlier near, outlier far */
            ASSERT(submode_MSB == 1);
            L_tmp = L_scf_idx[5];
            move32(); /* outlier near  section assumed */
            if (submode_LSB != 0)
            {                                                       /* outl_far */
                L_tmp = L_add(L_shl_pos(L_tmp, 1), L_gain_LSB); /*  add lsb bit of Gain */
                L_tmp = L_add(L_tmp, sns_MPVQ_Sz[2][0]);          /*  outlier far section offset added */
            }

            ASSERT(L_tmp < (1L << 24));
            /* outlier mode shape  index   total  23.8536 ( +  ~.14 ) bits as   MPVQ codeword  */
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_tmp), 12);              /*  multiplex 12  bits  LSB*/
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr(L_tmp, 12)), 12); /* multiplex 12 bits  MSBs */
        }
    }

    test();
    IF(write_legacy_sns_vq_bits == 0) /* 9 or 10 bits for stage1,    19-20 bits for stage2  */
    {
        /* SNS-VQ 1st stage is jointly multiplexed into 9 bits or 10 bits */
        /* input  scf_idx[0],  has the  stage1 index  for one of BCA   , if B+stage 1 only */
        tmp = extract_l(L_scf_idx[0]);
        aux_idx = extract_l(L_scf_idx[1]);      /* aux value: we have the LS, or the s0 sign bit     */
        shape_idx = extract_l(L_scf_idx[2]);        /* st2 shape 0 .. 5  , where [2.3.4.5] are fixed FESS shapes */
        gain_idx = extract_l(L_scf_idx[3]);         /* idx of 2-3 bits valued  gains */

#ifdef DEBUG
        if (shape_idx == 0 && L_scf_idx[0] >= 0)
        { /*aux==2 -->  split mode  st1B (aux==2), +  a 2 bit gain +  P(5,6)10.96b + P(8,2)7b   */
            ASSERT(shape_idx == 0 && gain_idx < 4);
        }
        if ((shape_idx == 1) && L_scf_idx[0] >= 0)
        { /*  regular mode st1B  , +  a 3 bit gain +   P(15,5) */
            ASSERT(shape_idx == 1 && gain_idx >= 0 && gain_idx < 8);
        }

        if ((shape_idx >= 2) && L_scf_idx[0] >= 0)
        { /*  regular mode st1B  , +  a 3 bit gain +   fixenv */
            ASSERT(shape_idx >= 2 && shape_idx <= 5 && gain_idx >= 0 && gain_idx < 8);
        }
#endif
        /*        b0-b8     b9
        segm ,  idx9b , stop bit,  comment use
        -----+--------+---------
         A   | 510,511| n/a,   2 entries,  9 bit total
       ------+--------+--------
         B   | 0--169 |   1  , 170 entries,  10 bit total
       ------+--------+--------
         C   | 170-339|   1  , 170 entries,   10 bit total
       ------+--------+--------+------------*/

        test(); test(); test(); test();
        IF(L_sub(L_scf_idx[0], 510) >= 0)
        {   /* stage1A */
            ASSERT(L_scf_idx[0] < (1 << 9));
            write_indice_backward(bytes, bp_side, mask_side, tmp, 9);             /*  writes values  510 or 511 */
            shape_idx = -9;   move16();                               /* use only the coarse stage 1A , no more bits to send */
        }
        ELSE IF((sub(tmp, 2 * 170) < 0) && (shape_idx < 0))
        {
            write_indice_backward(bytes, bp_side, mask_side, tmp, 9); /* 1B  [ 0..169]  or   1C   [ 170 ... 339 ] */

            /* write the stop bit sentinel value "1",  so that the demux dec_entropy_fx)  can stop already at (9+1)=10 bits */
            write_indice_backward(bytes, bp_side, mask_side, 1, 1); /*   dec_entropy_fx()  will read this a 1 "stop" bit   */
            shape_idx = -10; move16();
        }
        ELSE IF(shape_idx == 0) {
            /* aux info as a  part of  (9+1) 10 initial bits */
            /*        b0-b8     b9
             segm ,  idx9b , stop bit,  comment use
             -----+--------+---------
              B*  | 0--169 |   0  ,  --> aux=0, 170,  2b+17b for stage2 'LR_SplitLF',  29 bit total
            ------+--------+--------+-------
              B*  | 170-339|   0  ,  --> aux=1, 170,  2b+17b for stage2 'LR_SplitLF',  29 bit total
            ------+--------+--------+-------
           */
           /* 29 bit total LR_splitLF */
#ifdef DEBUG
            ASSERT(L_scf_idx[0] >= 0 && L_scf_idx[0] < 170); /* st1B range */
            ASSERT(aux_idx >= 0 && aux_idx <= 1);            /* aux_bit  can only be 1 or 0  */
#endif
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_mac0(L_scf_idx[0], aux_idx, 170)), 9);  /*  aux_bit  defined by region 0.169 or 170..339 in decoder */
            write_indice_backward(bytes, bp_side, mask_side, 0, 1); /*  "stop" bit. always zero for  'LR_splitLF' path */

            write_indice_backward(bytes, bp_side, mask_side, gain_idx, 2);   /* always 2bits == 4 gain levels for the splitLF mode */

            n5k = 6; move16();
            test();
            if (L_scf_idx[5] < 0)
            {
                ASSERT(L_scf_idx[5] == -8);
                n5k = 8;  move16();        /* 8 LF pulses ,  no HF pulses */
            }

            test();
            IF(sub(n5k, 6) == 0)
            {  /*  multiplex remaining  10 + 1+6 bit    */
                write_indice_backward(bytes, bp_side, mask_side, extract_l(L_scf_idx[4]), 10);         /*      P(5,6)=10.94  in LS+10 bits,  */
                write_indice_backward(bytes, bp_side, mask_side, extract_l(L_and(L_scf_idx[5], 0x1)), 1);    /* LS for P(8,2)=7  */
                write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr_pos(L_scf_idx[5], 1)), 6);    /* mPVQ(8,2) in 6   */
            }
            ELSE
            {   /* LF is PVQ(N=5,K=8), HF is all zero , multiplexed as top  section in stage2 10b  ,  + 7b  */
                /* scf_idx[4] is in the range  [0 ... (SNSLR_NPVQ_L5K8 >>1)[,   [0 .. 2945[  ,  11.52 bits */
                /* 985 - 1024 = 39 entries where we now use  32 (5 bits) */
                /* index_to_send = (SNSLR_NPVQ_L5K6 >> 1) + (scf_idx[4] & 0x001f) , range [985 .. 1017[ , i.e.   7 values are not sent:[1018 ... 1023] */
                ASSERT(L_scf_idx[4] < (SNSLR_NPVQ_L5K8 >> 1));
                write_indice_backward(bytes, bp_side, mask_side, extract_l(L_add((SNSLR_NPVQ_L5K6 >> 1), L_and(L_scf_idx[4], 0x001f))), 10); /* 5 lsb's as top in the 10b  block*/
                ASSERT((L_scf_idx[4] >> 5) < (1 << 7));
                write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr_pos(L_scf_idx[4], 5)), 7);                             /* 7 msb's of mpvq(5,8)*/
                /* unused region of 7*2^7 = 1792= 10.8 bits  may be used for BER  detection */
            }
        }
        ELSE IF(shape_idx == 1)
        {
            /*       b0-b8     b9
            segm ,  idx9b , stop bit,  comment/use
            ------+--------+--------+------------
              B*  | 340-509|   1     --> aux=1, 170, 3b+17b for stage2 'LR_full', 30 bit total
            ------+--------+--------+-------
              B*  | 340-509|   0     --> aux=0, 170, 3b+17b for stage2 'LR_full', 30 bit total
            ------+--------+--------+-------
            */
            ASSERT(L_scf_idx[0] >= 0 && L_scf_idx[0] < 170); /* validate  st1B* range */

            write_indice_backward(bytes, bp_side, mask_side, add(tmp, 2 * 170), 9); /* stage1B*   signal */
            write_indice_backward(bytes, bp_side, mask_side, aux_idx, 1);          /*  auxbit value transmitted  in the stop bit location */

            write_indice_backward(bytes, bp_side, mask_side, gain_idx, 3); /* 30b  always 8  gain levels in 3 bits for the full mode  */
            /* the  next 17 bit index is used to decode submode   1==full  or  one of submode 2,3,4,5 == fix */
            ASSERT(L_scf_idx[4] >= 0 && L_scf_idx[4] < ((SNSLR_NPVQ_L15K5 >> 1)));

            L_tmp = L_scf_idx[4];  move32(); /* full 17 bit index */
            /* typically P(15,5) without the leading sign,   16.6593 bits, written   */
            /* subset of 7 lsb's written first to match Cflt writing of all 17 bits at once */
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_tmp), 17 - 10);
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr_pos(L_tmp, 7)), 10);
        }
        ELSE
        {  /* Fixed  shapes 2(fixenv0) and 3(fixenv1) , 4(fixenv2),   5(fixenv3) */
            ASSERT(shape_idx >= 2 && shape_idx <= 5);
        /*  fixenv 0,1,2,  :  s0(in aux)  + shift(2bits)+ 11 signs  , section size 8192
            fixenv 3       :  s0(in aux)  + shift(2bits)+ 9 signs   , section size 2048
        */
        fixenv_shape_idx = extract_l(L_scf_idx[4]);
        ASSERT(fixenv_shape_idx == (shape_idx - 2)); /* scf_idx[4] has the fixed env index*/
        ASSERT(L_scf_idx[0] >= 0 && L_scf_idx[0] < 170); /* validate st1B* range */

        write_indice_backward(bytes, bp_side, mask_side, add(tmp,  2 * 170), 9); /* stage1B*   signal */
        write_indice_backward(bytes, bp_side, mask_side, aux_idx, 1);   /*  auxbit transmitted  in the stop bit location */

        write_indice_backward(bytes, bp_side, mask_side, gain_idx, 3); /* for this 30bit path it is always 3 gain bits    */

        IF(sub(shape_idx, 5) < 0) /* 2,3,4 */
        {
            ASSERT(L_scf_idx[5] >= 0 && L_scf_idx[5] < (1 << 13));
            ASSERT((SNSLR_NPVQ_L15K5 >> 1) + fixenv_shape_idx * (1 << 13) + L_scf_idx[5] < (1 << 17));
            /* offset is PVQ(15,5)  without  leading sign   */
            /* int32_t tmp_idx = (SNSLR_NPVQ_L15K5 >> 1) + env_shape_idx * (1 << 13) + scf_idx[5]; */
            L_tmp = L_add(L_mac0((SNSLR_NPVQ_L15K5 >> 1), fixenv_shape_idx, 1 << 13), L_scf_idx[5]); /* full 17 bits*/
        }
        ELSE
        {
            ASSERT(shape_idx == 5); /* 5*/
            ASSERT(fixenv_shape_idx == 3);
            ASSERT(L_scf_idx[5] < (1 << 11));
            ASSERT(((SNSLR_NPVQ_L15K5 >> 1) + 3 * (1 << 13) + L_scf_idx[5]) < (1 << 17));
            /* offset is PVQ(15,5)  without  leading sign   */
            L_tmp = L_add((SNSLR_NPVQ_L15K5 >> 1) + 3 * (1 << 13), L_scf_idx[5]); /* full 17 bits*/
        }
            /* subset of 7 lsb's written first to match Cflt writing of all 17 bits at once */
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_tmp), 17 - 10);
            write_indice_backward(bytes, bp_side, mask_side, extract_l(L_shr_pos(L_tmp, 7)), 10);
        }
    }
    Dyn_Mem_Deluxe_Out();
}
#endif /* CR9_C_ADD_1p25MS_LRSNS  */



#ifdef __cplusplus
} /* NAMESPACE_VERSION */
#endif
