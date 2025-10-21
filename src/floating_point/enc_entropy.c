/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

#ifdef CR9_C_ADD_1p25MS

#  ifdef NEW_SIGNALLING_SCHEME_1p25
void writeLtpData_fl( LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT32* ltpf_idx, LC3_INT16* Tx_ltpf );
# endif

void writeSNSData_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3PLUS_FrameDuration frame_dms, LC3_INT32* scf_idx);
#else /* CR9_C_ADD_1p25MS */
void writeSNSData_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT32* scf_idx);
#endif /* CR9_C_ADD_1p25MS */


static const LC3_INT32 gainMSBbits[4] = {1, 1, 2, 2};
static const LC3_INT32 gainLSBbits[4] = {0, 1, 0, 1};


#ifdef CR9_C_ADD_1p25MS
#  ifdef NEW_SIGNALLING_SCHEME_1p25
void writeLtpData_fl(
    LC3_UINT8* ptr,
    LC3_INT* bp_side,
    LC3_INT* mask_side,
    LC3_INT32* ltpf_idx,
    LC3_INT16* Tx_ltpf
)
{
    LC3_INT32 tmp;
    LC3_INT32 bitsTx;
    /*        Hdr, information               , bits used
              00 , no lag info , no phase info  sum=2
              010, PhaseA,LTPF=0, lagAbits=4 ,  sum=7    : pitch-PLC may be activated with old lag , 4 MSbs
              011, PhaseB,LTPF=0, lagBbits=4*,  sum=7*   : pitch-PLC may be activated with fresh lag, 4* =  reduced lag resolution in Quantized ltpf_idx domain for PLC
              10 , PhaseA,LTPF=1, lagAbits=4 ,  sum=6    : LTPF may be activated with old lag
              11 , PhaseB,LTPF=1, lagBbits=5 ,  sum=7    : LTPF activated with fresh lag
*/
 

    tmp = MIN(*Tx_ltpf, 1);           /*phaseA==0, phaseB==1*/

    if (ltpf_idx[0] == 0)
    {
        write_bit_backward_fl(ptr, bp_side, mask_side, 0); /* "00" */
        write_bit_backward_fl(ptr, bp_side, mask_side, 0);

        *Tx_ltpf = 0;
        /*  *Tx_ltpf A/B state forced to zero or  kept at zero */
        /*  decoder will discard any sofar in phase received phaseA MSB bits */
    }
    else if (ltpf_idx[1] == 0)
    {
        /* no current LTPF activation,
           lag transmitted for PLC,  or for next frame LTPF activation */
        assert(ltpf_idx[0] != 0);

        /* A "010"  3 bits Hdr transmitted */
        /* B "011"  3 bits Hdr transmitted*/

        write_uint_backward_fl(ptr, bp_side, mask_side, 1, 2); /* "01"*/
        write_bit_backward_fl(ptr, bp_side, mask_side, tmp);  /* phase A or phaseB  */

        if (*Tx_ltpf == 0)
        {  /* phase A transmission */
            assert(tmp == 0);
            assert((ltpf_idx[2] & ~(0x01ff)) == 0); /* only 9 bits info allowed  within  ltpf_idx[2] */
            tmp = (ltpf_idx[2] >> 5 );  /* shift_out LSBS, send 4 MSBs  */
            *Tx_ltpf = (0x200 | ltpf_idx[2]);   /*  remember full lag,  as phaseB sentinel in bit 10,   */
        }
        else
        {   /* phase B */
            assert(tmp == 1);
            assert(*Tx_ltpf > 511);  /*sentinel in b10 should have been set in previous phaseA frame */
            tmp = ((*Tx_ltpf & 0x001f) >> 1);  /* B send 4* LSBs,  1 bit truncated */
            *Tx_ltpf = 0;                   /*  clear sentinel in b10  and the old  remebered lag value */
        }
        write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 4); /* 4 bits  lag info Tx, when LTPF is deactivated  */
    }
    else
    {   /* LTPF activated */
        /* A "10"  2 bits Hdr */
        /* B "11"  2 bits Hdr */
        assert(ltpf_idx[0] != 0 && ltpf_idx[1] != 0);

        tmp |= 0x02;
        write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 2);

        if (*Tx_ltpf == 0)
        {
            bitsTx = 4;
            assert((ltpf_idx[2] & ~(0x01ff)) == 0); /* only 9 bits info allowed  within  ltpf_idx[2] */
            tmp = (ltpf_idx[2] >> 5);          /* shift away 5 LBS,  LTPF active, send phaseA 4 MSBs  */

            *Tx_ltpf = 0x0200 | ltpf_idx[2];  /*  remember full lag in state *Tx_ltpf,  add phaseB Tx state sentinel in bit 10   */
        }
        else
        {
            bitsTx = 5;
            tmp = (*Tx_ltpf & 0x001f);        /* LTPF active   B send 5 LSBs , full regular resolution  */
            *Tx_ltpf = 0;                       /* clear sentinel in b10  and also the old  remebered lag value */
        }

        write_uint_backward_fl(ptr, bp_side, mask_side, tmp, bitsTx);  /* ltp==1 , ltpf==1,  4(A,MSBs)  or 5(b, LSBs) bits */
    }
}
#endif /*  NEW_SIGNALLING_SCHEME_1p25 */
void writeSNSData_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3PLUS_FrameDuration frame_dms, LC3_INT32* scf_idx)
#else
void writeSNSData_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT32* scf_idx)
#endif
{
#  ifdef CR9_C_ADD_1p25MS_LRSNS

    LC3_INT32 submodeMSB, submodeLSB, tmp, gainMSB, gainLSB;
    LC3_INT16 write_legacy_sns_vq_bits;
    LC3_INT32  aux_idx;
    LC3_INT32 shape_idx;
    LC3_INT32 env_shape_idx;
    LC3_INT32 gain_idx;
    LC3_INT16 n5k;

#  else
        LC3_INT32 submodeMSB, submodeLSB, tmp, gainMSB, gainLSB;
#  endif

#if defined(CR9_C_ADD_1p25MS) && !defined( CR9_C_ADD_1p25MS_LRSNS)
UNUSED(frame_dms);
#endif

#ifdef CR9_C_ADD_1p25MS_LRSNS
    write_legacy_sns_vq_bits = 1;
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) /*  9,10,29/30 */
    {
        write_legacy_sns_vq_bits = 0;
    }
    if (write_legacy_sns_vq_bits != 0)
    {
#endif
        /* SNS-VQ 1st stage */
        write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0], 5);
        write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[1], 5);

        /* SNS-VQ 2nd stage side-info (3-4 bits) */
        submodeMSB = scf_idx[2] / 2;
        submodeLSB = scf_idx[2] & 1;
        write_bit_backward_fl(ptr, bp_side, mask_side, submodeMSB);
        gainMSB = scf_idx[3] >> (gainLSBbits[scf_idx[2]]);
        gainLSB = scf_idx[3] & 1;
        write_uint_backward_fl(ptr, bp_side, mask_side, gainMSB, gainMSBbits[scf_idx[2]]);
        write_bit_backward_fl(ptr, bp_side, mask_side, scf_idx[4]);

        /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */
        if (submodeMSB == 0)
        {
            if (submodeLSB == 0)
            {
                tmp = scf_idx[6] + 2;
            }
            else
            {
                tmp = gainLSB;
            }

            tmp = tmp * 2390004 + scf_idx[5];
            write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 25);
        }
        else
        {
            tmp = scf_idx[5];

            if (submodeLSB != 0)
            {
                tmp = 2 * tmp + gainLSB + 15158272;
            }
            write_uint_backward_fl(ptr, bp_side, mask_side, tmp, 24);
        }
#ifdef CR9_C_ADD_1p25MS_LRSNS
    }
    if (write_legacy_sns_vq_bits == 0)
    {
        /* SNS-VQ 1st stage is jointly multiplexed into 9 bits or 10 bits */
        /* input  scf_idx[0],  has the  stage1 index  for one of BCA  or B*  */
        aux_idx = scf_idx[1];      /* aux value: we have the LS, or the s0 sign bit     */
        shape_idx = scf_idx[2];        /* st2 shape 0 .. 5  , where [2.3.4.5] are fixed FESS shapes */
        gain_idx  = scf_idx[3];         /* idx of 2-3 bits valued  gains */

        if (shape_idx == 0 && scf_idx[0] >= 0)
        { /*aux==2 -->  split mode  st1B (aux==2), +  a 2 bit gain +  P(5,6)10.96b + P(8,2)7b */
            assert(shape_idx == 0 && gain_idx < 4);
        }
        if ((shape_idx == 1) && scf_idx[0] >= 0)
        { /*  regular mode st1B  , +  a 3 bit gain +   P(15,5) */
            assert(shape_idx == 1 && gain_idx >= 0 && gain_idx < 8);
        }
        if ((shape_idx >= 2) && scf_idx[0] >= 0)
        { /*  regular mode st1B  , +  a 3 bit gain +   fixenv */
            assert(shape_idx >= 2 && shape_idx <= 5 && gain_idx >= 0 && gain_idx < 8);
        }

        /*        b0-b8     b9
        segm ,  idx9b , stop bit,  comment use
        -----+--------+---------
         A   | 510,511| n/a,   2 entries,  9 bit total
       ------+--------+--------
         B   | 0--169 |   1  , 170 entries,  10 bit total
       ------+--------+--------
         C   | 170-339|   1  , 170 entries,   10 bit total
       ------+--------+--------+------------
        */

        if (scf_idx[0] >= 510)
        {   /* stage1A */
            assert(scf_idx[0] < (1 << 9));
            write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0], 9);   /*  currently writes    510 or 511 */
            shape_idx = -9;                                           /* only crude stage 1 A , no more bits to send */
        }
        else if ( (scf_idx[0] < 2 * 170) && (shape_idx < 0) )
        {
            if (scf_idx[0] < 170)
            {
                write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0] + 0 * 170, 9); /* 1B */
            }
            else if (scf_idx[0] < 2 * 170)
            {
                write_uint_backward_fl(ptr, bp_side, mask_side, (scf_idx[0] - 170) + 1 * 170, 9);  /* 1C  */
            }
            /* write the stop bit value  sentinel value "1",  so that the demux dec_entropy  can stop already at (9+1)=10 bits */

            write_uint_backward_fl(ptr, bp_side, mask_side, 1, 1); /*   dec_entropy will read a 1 "stop" bit   */

            /* pitch_rx and ltpf_rx handled separately */
        }
        else if (shape_idx == 0)
        {
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
            assert(scf_idx[0] >= 0 && scf_idx[0] < 170); /* st1B range */
            assert(aux_idx >= 0 && aux_idx <= 1);        /* aux_bit  can only be 1 or 0  */
            write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0] + aux_idx * 170, 9);  /* aux_bit  defined by region 0..339 region in decoder */
            write_uint_backward_fl(ptr, bp_side, mask_side, 0, 1); /*  "stop" bit. always zero for  'LR_splitLF' */

            write_uint_backward_fl(ptr, bp_side, mask_side, gain_idx, 2);             /* always 2bits == 4 gain levels for the splitLF mode */

            n5k = 6;
            if (scf_idx[5] < 0) {
                assert(scf_idx[5] == -8 );
                n5k = 8;
            }

            if (n5k == 6)
            {  /*  multiplex remaining  10 + 1+6 bit    */
                write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[4], 10);        /* 29b   P(5,6)=10.94  in LS+10 bits,  */
                write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[5] & 0x1, 1);   /* LS for P(8,2)=7  */
                write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[5] >> 1, 6);    /* mPVQ(8,2) in 6   */
            }
            else
            {   /* LF is PVQ(N=5,K=8), HF is all zero , multiplexed as top  section in stage2 10b  ,  + 7b  */
                /* scf_idx[4] is in the range  [0 ... (SNSLR_NPVQ_L5K8 >>1)[,   [0 .. 2945[  ,  11.52 bits */
                /* 985 - 1024 = 39 entries where we now use only 32 (5 bits) */
                /* index_to_send = (SNSLR_NPVQ_L5K6 >> 1) + (scf_idx[4] & 0x001f) , range [985 .. 1017[ , i.e.   7 values are not sent:[1018 ... 1023] */
                write_uint_backward_fl(ptr, bp_side, mask_side, (SNSLR_NPVQ_L5K6 >> 1) + (scf_idx[4] & 0x001f), 10); /* 5 lsb's as top in the 10b */
                assert((scf_idx[4] >> 5) < (1 << 7));
                write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[4] >> 5, 7);                             /* 7 msb's of mpvq(5,8)*/
            }
        }
        else if (shape_idx == 1)
        {
            /*       b0-b8     b9
            segm ,  idx9b , stop bit,  comment use
            ------+--------+--------+------------
              B*  | 340-509|   1     --> aux=1, 170, 3b+17b for stage2 'LR_full', 30 bit total
            ------+--------+--------+-------
              B*  | 340-509|   0     --> aux=0, 170, 3b+17b for stage2 'LR_full', 30 bit total
            ------+--------+--------+-------
            */
            assert(scf_idx[0] >= 0 && scf_idx[0] < 170); /* st1B* range */

            write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0] + 2 * 170, 9); /* stage1B*   signal */
            write_uint_backward_fl(ptr, bp_side, mask_side, aux_idx, 1);              /*  auxbit transmitted  in the stop bit location */

            write_uint_backward_fl(ptr, bp_side, mask_side, gain_idx, 3); /*30b  always 8  gain levels in 3 bits for the full mode  */
            /* the  next 17 bit index is used to decode submode   1==full  or  one of submode 2,3,4,5 == fix */
            assert(scf_idx[4] >= 0 && scf_idx[4] < ((SNSLR_NPVQ_L15K5 >> 1)));
            write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[4], 17); /* P(15,5) without the leading sign,   16.6593 bits, written   */
        }
        else
        {  /* Fixed  shapes 2(ones, env0) and 3(env1) , 4(env2), and  later 5(env3) */
            assert(shape_idx >= 2 && shape_idx <= 5);
            /*  env 0,1,2,  :  s0(in aux)  + shift(2bits)+ 11 signs  , section size 8192
                env 3       :  s0(in aux)  + shift(2bits)+ 9 signs   , section size 2048
            */
            env_shape_idx = scf_idx[4];
            assert(env_shape_idx == (shape_idx - 2)); /*scf_idx[4] has the fixed env index*/

            assert(scf_idx[0] >= 0 && scf_idx[0] < 170); /* st1B* range */

            write_uint_backward_fl(ptr, bp_side, mask_side, scf_idx[0] + 2 * 170, 9); /* stage1B*   signal */
            write_uint_backward_fl(ptr, bp_side, mask_side, aux_idx, 1);              /*  auxbit transmitted  in the stop bit location */

            write_uint_backward_fl(ptr, bp_side, mask_side, gain_idx, 3); /* for 30bit   its is always 3 gain bits    */

            if (shape_idx < 5) /*2,3,4*/
            {
                assert(scf_idx[5] >= 0 && scf_idx[5] < (1 << 13));
                assert((SNSLR_NPVQ_L15K5 >> 1) + env_shape_idx * (1 << 13) + scf_idx[5] < (1 << 17));
                /* offset is PVQ(15,5)  without  leading sign   */
                /* int32_t tmp_idx = (SNSLR_NPVQ_L15K5 >> 1) + env_shape_idx * (1 << 13) + scf_idx[5]; */
                write_uint_backward_fl(ptr, bp_side, mask_side, (SNSLR_NPVQ_L15K5 >> 1) + env_shape_idx * (1 << 13) + scf_idx[5], 17);
            }
            else
            {
                assert(shape_idx == 5);
                assert(env_shape_idx == 3);
                assert(scf_idx[5] < (1 << 11));
                assert(((SNSLR_NPVQ_L15K5 >> 1) + 3 * (1 << 13) + scf_idx[5]) < (1 << 17));
                /* offset is (15,5)  without  leading sign   */
                write_uint_backward_fl(ptr, bp_side, mask_side, (SNSLR_NPVQ_L15K5 >> 1) + 3 * (1 << 13) + scf_idx[5], 17);
            }
        }
    }
#endif
}

void processEncoderEntropy_fl(LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT numbytes, LC3_INT bw_cutoff_bits,
                              LC3_INT bw_cutoff_idx, LC3_INT lastnz, LC3_INT N, LC3_INT lsbMode, LC3_INT gg_idx, LC3_INT num_tns_filters,
                              LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* scf_idx, LC3_INT fac_ns_idx,
                              LC3_INT bfi_ext, LC3_INT fs_idx
#ifdef CR9_C_ADD_1p25MS
                              , LC3PLUS_FrameDuration frame_dms, LC3_INT16* Tx_ltpf
#endif
                              )
{
    LC3_UINT8* ptr;
    LC3_INT    i;

    LC3_INT16 lastnzTrigger[5] = {63, 127, 127, 255, 255};

#if defined(CR9_C_ADD_1p25MS) && !defined(CR9_C_ADD_1p25MS_LRSNS)
    UNUSED(Tx_ltpf);
#endif
    *bp_side   = numbytes - 1;
    *mask_side = 1;
    ptr        = bytes;

    /* Bandwidth */
    if (bw_cutoff_bits > 0) {
        write_uint_backward_fl(ptr, bp_side, mask_side, bw_cutoff_idx, bw_cutoff_bits);
    }

    /* Last non zero touple */
    if (bfi_ext == 1) {
        write_uint_backward_fl(ptr, bp_side, mask_side, lastnzTrigger[fs_idx], getLastNzBits (N));
    }
    else
    {
        write_uint_backward_fl(ptr, bp_side, mask_side, lastnz / 2 - 1, getLastNzBits (N));
    }

    /* LSB mode bit */
    write_bit_backward_fl(ptr, bp_side, mask_side, lsbMode);

    /* Global gain */
    write_uint_backward_fl(ptr, bp_side, mask_side, gg_idx, 8);

    /* TNS activation flag */
    for (i = 0; i < num_tns_filters; i++) {
        write_bit_backward_fl(ptr, bp_side, mask_side, MIN(1, tns_order[i]));
    }

    /* LTPF activation flag */
#ifdef CR9_C_ADD_1p25MS
#   ifdef NEW_SIGNALLING_SCHEME_1p25
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        writeLtpData_fl(ptr, bp_side, mask_side, ltpf_idx, Tx_ltpf);     /* LTP and LTPF-active and interleaved lag-idx */
    }
    else
    {
        write_bit_backward_fl(ptr, bp_side, mask_side, ltpf_idx[0]);    /* ltp tx bit */
    }
#   endif
#else
    write_bit_backward_fl(ptr, bp_side, mask_side, ltpf_idx[0]);
#endif /* CR9_C_ADD_1p25MS  */

    /* SNS data*/
#ifdef CR9_C_ADD_1p25MS
    writeSNSData_fl(ptr, bp_side, mask_side, frame_dms, scf_idx);
#else
    writeSNSData_fl(ptr, bp_side, mask_side, scf_idx);
#endif

    /* LTPF data */
#ifdef CR9_C_ADD_1p25MS
#   ifdef NEW_SIGNALLING_SCHEME_1p25
    if ((frame_dms != LC3PLUS_FRAME_DURATION_1p25MS ) && (ltpf_idx[0] == 1) )
        {
        write_uint_backward_fl( ptr, bp_side, mask_side, ltpf_idx[1], 1 );
        write_uint_backward_fl( ptr, bp_side, mask_side, ltpf_idx[2], 9 );
    }
#   endif
#else /* CR9_C_ADD_1p25MS */
    if (ltpf_idx[0] == 1) {
        write_uint_backward_fl(ptr, bp_side, mask_side, ltpf_idx[1], 1);
        write_uint_backward_fl(ptr, bp_side, mask_side, ltpf_idx[2], 9);
    }
#endif /* CR9_C_ADD_1p25MS */

    /* Noise factor */
    write_uint_backward_fl(ptr, bp_side, mask_side, fac_ns_idx, 3);
}

void write_uint_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT val, LC3_INT numbits)
{
    LC3_INT k, bit;

    for (k = 0; k < numbits; k++) {
        bit = val & 1;
        write_bit_backward_fl(ptr, bp_side, mask_side, bit);
        val = val >> 1;
    }
}

void write_bit_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT bit)
{
    if (bit != 0) {
        ptr[*bp_side] = ptr[*bp_side] | *mask_side;
	}

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side << 1;
    }
}
