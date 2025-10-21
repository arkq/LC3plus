/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

 

#ifdef CR9_C_ADD_1p25MS_LRSNS

static void readSNSData_fx(UWord8* ptr,
    Word16* bfiPtr,
    Word16* mask_side,
    Word16* bp_side,
    Word16* ltpf_idx_2_lrsns,
    Word32* L_scf_idx,
    LC3PLUS_FrameDuration frame_dms);

#endif 

static Word16 read_indice(UWord8 *ptr, Word16 *bp, Word16 *mask, Word16 numbits)
{
    Dyn_Mem_Deluxe_In(
        Word16  indice, bit;
        Counter i;
    );

    indice = read_bit(ptr, bp, mask);

    FOR (i = 1; i < numbits; i++)
    {
        bit    = read_bit(ptr, bp, mask);
        indice = add(indice, lshl_pos(bit, i));
    }

    Dyn_Mem_Deluxe_Out();
    return indice;
}

static Word16 ac_dec_split_st2VQ_CW(                     /* local BER flag */
                                    const Word32 L_cwRx, /* max 25 bits */
                                    const Word32 L_szA, const Word32 L_szB, Word32 *L_cwA, Word32 *L_cwB,
                                    Word16 *submodeLSB);

void processDecoderEntropy_fx(UWord8 *bytes, Word16 *bp_side, Word16 *mask_side, Word16 nbbits,
                              Word16 L_spec, Word16 fs_idx, Word16 BW_cutoff_bits, Word16 *tns_numfilters,
                              Word16 *lsbMode, Word16 *lastnz, Word16 *bfi, Word16 *tns_order, Word16 *fac_ns_idx,
                              Word16 *gg_idx, Word16 *BW_cutoff_idx, Word16 *ltpf_idx, Word32 *L_scf_idx,
                              LC3PLUS_FrameDuration frame_dms
#ifdef CR9_C_ADD_1p25MS
                              , Word16 rx_status[2], Word16* mem_continuation
#  ifdef NEW_SIGNALLING_SCHEME_1p25
                               ,Word16 *ltpfinfo_frame_cntr_fx  /* set here , but  also increased outside  by  bfi for the channel */
#  endif
#endif

                              )
{
#ifdef CR9_C_ADD_1p25MS_LRSNS
    Dyn_Mem_Deluxe_In(
        Word16 L;
    Word32 tmp32;
    Word16 gain_e, gain;
    Counter n;
    UWord8 * ptr;
    Word16  ltpf_idx_2_lrsns[3];
    Word16 bfiSNS;
    );
#else
    Dyn_Mem_Deluxe_In(
        Word16 L, submodeLSB;
        Word32 tmp32, tmp32lim;
        Word16 gain_e, gain, submodeMSB, BER_detect;
        Counter n;
        UWord8 * ptr; );
#endif

#if !defined(LRSNS_PC_SIGNAL_FIX) && defined(CR9_C_ADD_1p25MS_LRSNS)
    UNUSED(bfiSNS);
#endif

    ptr        = bytes;
    *bp_side   = shr_pos(sub(nbbits, 1), 3);
    *mask_side = shl(1, sub(8, sub(nbbits, shl_pos(*bp_side, 3))));

    /* Cutoff-detection */
    IF (BW_cutoff_bits > 0)
    {
        *BW_cutoff_idx = read_indice(ptr, bp_side, mask_side, BW_cutoff_bits);
        /* check for bitflips */
        IF (sub(fs_idx, *BW_cutoff_idx) < 0)
        {
            *BW_cutoff_idx = fs_idx;
            *bfi           = 1;  move16();
            Dyn_Mem_Deluxe_Out();
            return;
        }
    }
    ELSE
    {
        *BW_cutoff_idx = fs_idx;
    }

    /* Number of TNS filters */
# ifdef CR9_C_ADD_1p25MS
    IF (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
        *tns_numfilters = 0;
    } ELSE
#ifdef CR9_C_ADD_1p25MS_LRSNS
    {
#endif
# endif
        IF (sub(*BW_cutoff_idx, 3) >= 0 && frame_dms >= LC3PLUS_FRAME_DURATION_5MS)
        {
            *tns_numfilters = 2;  move16();
        }
        ELSE
        {
            *tns_numfilters = 1;  move16();
        }
#ifdef CR9_C_ADD_1p25MS_LRSNS
    }
#endif

    /* Decode number of ntuples */
    L       = getLastNzBits_fx(L_spec);
    n       = read_indice(ptr, bp_side, mask_side, L);
    n       = add(n, 1);
    *lastnz = shl_pos(n, 1);
    IF (sub(*lastnz, L_spec) > 0)
    {
        *bfi = 1;  move16();
        Dyn_Mem_Deluxe_Out();
        return;
    }

    /* Mode bit */
    *lsbMode = read_bit(ptr, bp_side, mask_side);

    /* Decode global-gain */
    *gg_idx = read_indice(ptr, bp_side, mask_side, 8);  move16();
    tmp32  = L_shl_pos(L_mult0(*gg_idx, 0x797D), 7);  /* 6Q25; 0x797D -> log2(10)/28 (Q18) */
    gain_e = add(extract_l(L_shr_pos(tmp32, 25)), 1); /* get exponent */
    gain   = round_fx(BASOP_Util_InvLog2(L_or(tmp32, 0xFE000000)));
    assert(gain >= 0); /* JSv, check if shr_pos(gain,1)  is more appropriate) */
    gain   = shr_r(gain, 1);
    gain_e = add(gain_e, 1);

    /* Decode TNS on/off flag */
    tns_order[1] = 0; move16(); /* fix problem with uninitialized memory */
    FOR (n = 0; n < *tns_numfilters; n++)
    {
        tns_order[n] = read_bit(ptr, bp_side, mask_side);  move16();
    }

    /* LTPF on/off */
#ifdef NEW_SIGNALLING_SCHEME_1p25
    IF(sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) != 0)
    {
        ltpf_idx[0] = read_bit(ptr, bp_side, mask_side); move16();
    }
    ELSE
    {    /* read one of {2, 6, 7} bits into ltp/ltpf/lag  variable  ltpf_idx[ 0 ...  2]  */
          readLtpData_fx(ptr, bfi, mask_side, bp_side, ltpf_idx, rx_status, ltpfinfo_frame_cntr_fx, mem_continuation);
    }
#else
    ltpf_idx[0] = read_indice(ptr, bp_side, mask_side, 1);  move16();
#endif


/* Decode SNS VQ parameters - 1st stage (10 bits) */
#ifdef CR9_C_ADD_1p25MS_LRSNS
    ltpf_idx_2_lrsns[0] = ltpf_idx[0]; move16();
    ltpf_idx_2_lrsns[1] = ltpf_idx[1]; move16();
#ifdef LRSNS_PC_SIGNAL_FIX  
    bfiSNS = 0;  move16(); /* Local BFI flag for Errors SNS bit area */

    readSNSData_fx(ptr, &bfiSNS, mask_side, bp_side, ltpf_idx_2_lrsns, L_scf_idx, frame_dms);
    IF ( bfiSNS != 0 )
    {   /* corrupt SNSbits triggers PLC through global PLC flag *bfi==1.
           *bfi==2 and bfiSNS == 0 maintains bfi==2 for PC(Partial Concealmnet)
        */
        *bfi = 1;
        Dyn_Mem_Deluxe_Out();
        return;
    }
    IF ( sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS ) == 0 )
    {   /* for 1.25ms  and previously detected bit errors -->  handle frame as a completely corrupt bad frame  */
        IF( sub(*bfi, 2) == 0)
        {
            *bfi = 1;
            Dyn_Mem_Deluxe_Out();
            return;
        }
    }

# else 
    readSNSData_fx(ptr, bfi, mask_side, bp_side, ltpf_idx_2_lrsns, L_scf_idx, frame_dms);

    IF(*bfi != 0)
    {
        *bfi = 1;
        Dyn_Mem_Deluxe_Out();
        return;
    }
#endif
#else

#ifdef ENABLE_HR_MODE
    L = read_indice(ptr, bp_side, mask_side, 5 + 5);
    L_scf_idx[0] = L_deposit_l(s_and(L, 0x1F)); /* stage1 LF  5  bits */
    L_scf_idx[1] = L_deposit_l(shr_pos(L, 5)); /* stage1 HF  5 bits  */
#else
    L_scf_idx[0] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 5)); /* stage1 LF  5  bits */
    L_scf_idx[1] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 5)); /* stage1 HF  5 bits  */
#endif

    /* Decode SNS VQ parameters - 2nd stage side-info (3-4 bits) */
    submodeMSB   = read_bit(ptr, bp_side, mask_side); /* submodeMSB 1 bit */
    L_scf_idx[2] = L_deposit_l(shl_pos(submodeMSB, 1));
    ASSERT(sns_gainMSBbits[L_scf_idx[2]] > 0);
    L_scf_idx[3] = L_deposit_l(
        read_indice(ptr, bp_side, mask_side, sns_gainMSBbits[L_scf_idx[2]])); /* gains or gain MSBs  1-2 bits  */
    L_scf_idx[4] = read_bit(ptr, bp_side, mask_side);                         /*  shape LS 1 bit */

    /* Decode SNS VQ parameters - 2nd stage data (24-25 bits) */
    IF (submodeMSB == 0)
    { /* shape_j = 0, or 1  */
        /* regular mode A,B indexes integer multiplexed, total 24.x bits  MPVQ codeword section A and  codeword for
         * section B */
        /* regular mode  mode shape  index   total  24.9999 bits    MPVQ codeword  */
        tmp32 = L_deposit_l(read_indice(ptr, bp_side, mask_side, 13));
        tmp32 = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 12), 13));  move16(); /*for ber state   */
        BER_detect =
            ac_dec_split_st2VQ_CW(       /* local BER flag */
                                  tmp32, /* L_cwRx  max 25 bits */
                                  sns_MPVQ_Sz[0][0], UL_addNsD(sns_MPVQ_Sz[0][1], sns_MPVQ_Sz[1][1]), /* 12+2 = 14 */
                                  (&L_scf_idx[5]),                                                    /* shape A */
                                  (&L_scf_idx[6]), /* shape B or  gain LSB */
                                  &submodeLSB      /* total submode update below  */
            );
        IF (submodeLSB != 0)
        { /* add gainLSB bit */
            L_scf_idx[3] = L_add(L_shl_pos(L_scf_idx[3], 1), L_scf_idx[6]);
            L_scf_idx[6] = -2L;
        }
    }
    ELSE
    { /* shape_j = 2 or 3  */
        ASSERT(submodeMSB == 1);
        /* outlier mode shape  index   total  23.8536 +  19.5637 (19.5637 < (log2(2.^24 -2.^23.8537))    bits    MPVQ
         * codeword  */

        tmp32 = L_deposit_l( read_indice( ptr, bp_side, mask_side, 12 ) );
        tmp32 = L_or( tmp32, L_shl_pos( read_indice( ptr, bp_side, mask_side, 12 ), 12 ) );

        L_scf_idx[5] = tmp32;
        move32(); /*shape outl_near or outl_far */
        submodeLSB = 0;
        move16();
        BER_detect = 0;
        move16();
        tmp32lim = L_add( sns_MPVQ_Sz[2][0], L_shl_pos( sns_MPVQ_Sz[3][0], 1 ) );
        IF( L_sub( tmp32, tmp32lim ) >= 0 )
        {
            BER_detect = 1;
            move16();
        }
        ELSE
        {
            tmp32 = L_sub( tmp32, sns_MPVQ_Sz[2][0] ); /*  a potential high index is computed */
            IF( tmp32 >= 0 )
            {
                submodeLSB = 1;
                move16();
                ASSERT( tmp32 >= 0 && tmp32 < (Word32) ( 2 * sns_MPVQ_Sz[3][0] ) );
                L_scf_idx[3] = L_add( L_shl_pos( L_scf_idx[3], 1 ), L_and( tmp32, 0x1 ) ); /* add LSB_gain bit to gain MSBs */
                L_scf_idx[5] = L_shr_pos( tmp32, 1 );                                      /* MPVQ index with offset and gainLSB removed */
                L_scf_idx[6] = -2L;
                move32();
            }
            ELSE
            {
                L_scf_idx[6] = -1L;
                move32();
            }
        }
    }
    L_scf_idx[2] =
        L_add( L_scf_idx[2], L_deposit_l( submodeLSB ) ); /* decoder internal signal shape_j = submode 0..3 to VQ */

    IF( BER_detect > 0 )
    {
        *bfi = 1;        move16();
        Dyn_Mem_Deluxe_Out();
        return;
    }
#endif

#ifdef CR9_C_ADD_1p25MS

#ifdef NEW_SIGNALLING_SCHEME_1p25
    IF( sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) != 0 )
    {
         ltpf_idx[1] = 0; move16();
         ltpf_idx[2] = 0; move16();

         test();
         IF( ltpf_idx[0] != 0 )
         {
             L = read_indice(ptr, bp_side, mask_side, 1 + 9);
             ltpf_idx[1] = s_and(L, 1);                     move16();
             ltpf_idx[2] = shr_pos(L, 1);                   move16();
         }
    }
#else
    IF (ltpf_idx[0] == 1) {
        if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
            Word32 rx_current_status = read_bit(ptr, bp_side, mask_side);
            IF (rx_current_status == rx_status[0]) {
                IF (rx_current_status == 0) {
                    rx_status[0] = 1;                                       move16();
                    rx_status[1] = read_indice(ptr, bp_side, mask_side, 5);   move16();
                    if (*mem_continuation == 0)
                    {
                        *mem_continuation = 1;
                    }
                } ELSE {
                    /* rx current status 1 */
                    rx_status[0] = 0;                                       move16();
                    ltpf_idx[1] = shr_pos (rx_status[1],4);
                    ltpf_idx[2] = read_indice(ptr, bp_side, mask_side, 5);
                    ltpf_idx[2] = add(ltpf_idx[2], shl_pos(s_and(rx_status[1],15), 5));
                    rx_status[1] = 0;                                       move16();
                    *mem_continuation = 0;
                }
            } ELSE { /* frame loss condtion */
                rx_status[0] = 0;
#ifdef FIX_PLC_CONFORM_ISSUES
                rx_status[1] = read_indice( ptr, bp_side, mask_side, 5 );
#else
                rx_status[1] = 0;
#endif       
            }
        } ELSE {
            L = read_indice(ptr, bp_side, mask_side, 1+9);            move16();
            ltpf_idx[1] = s_and(L, 1);                              move16();
            ltpf_idx[2] = shr_pos(L, 1);                            move16();
        }
    }    
    ELSE
    {
        ltpf_idx[1] = 0;  move16();
        ltpf_idx[2] = 0;  move16();
        rx_status[0] = 0;
        rx_status[1] = 0;
        *mem_continuation = 0;
    }
#endif
#else /* CR9_C_ADD_1p25MS */
    /* LTPF data */
    IF (ltpf_idx[0] != 0)
    {
#ifdef ENABLE_HR_MODE
        L = read_indice(ptr, bp_side, mask_side, 1+9);          move16();
        ltpf_idx[1] = s_and(L, 1);                              move16();
        ltpf_idx[2] = shr_pos(L, 1);                            move16();
#else
        ltpf_idx[1] = read_indice(ptr, bp_side, mask_side, 1);  move16();
        ltpf_idx[2] = read_indice(ptr, bp_side, mask_side, 9);  move16();

#endif
    }
    ELSE
    {
         ltpf_idx[1] = 0; move16();
         ltpf_idx[2] = 0; move16();
    }
#endif  /* new signalling */

    /* Decode noise-fac */
    *fac_ns_idx = read_indice( ptr, bp_side, mask_side, 3 );  move16();

    Dyn_Mem_Deluxe_Out();
}


#  ifdef ENABLE_PADDING
Word32 paddingDec_fx( UWord8* bytes, Word16 nbbits, Word16 L_spec, Word16 BW_cutoff_bits, Word16 ep_enabled, Word16* total_padding, Word16* np_zero )
{
    Word16 lastnz_threshold;
    Word16 padding_len_bits, padding_len;

    Word16 bp_side;
    Word16 nbbytes = shr( nbbits, 3 );

    Word16 mask_side;
    UWord8* ptr = bytes;

    Word16 lastnz;
    Word16 nbits = getLastNzBits_fx( L_spec );
    if ( sub( nbbits, nbits ) < 0 )
    {
        return 1;
    }
    *np_zero = 0;

    *total_padding = 0;

    bp_side = shr_pos( sub( nbbits, 1 ), 3 );
    mask_side = shl( 1, sub( 8, sub( nbbits, shl_pos( bp_side, 3 ) ) ) );

    test();
    IF( sub( bp_side, 19 ) < 0 || sub( bp_side, LC3PLUS_MAX_BYTES ) >= 0 )
    {
        return 1;
    }

    ptr = bytes;

    IF( BW_cutoff_bits > 0 )
    {
        read_indice( ptr, &bp_side, &mask_side, BW_cutoff_bits );
        move16();
    }

    lastnz = read_indice( ptr, &bp_side, &mask_side, nbits );
    move16();

    lastnz_threshold = sub( shl( 1, nbits ), 2 );

    WHILE( lastnz == lastnz_threshold )
    {
        padding_len_bits = sub( sub( 12, nbits ), BW_cutoff_bits );

        /*Read padding length*/
        padding_len = read_indice( ptr, &bp_side, &mask_side, padding_len_bits );
        move16();

        /* Read 4 reserved bits */
        read_indice( ptr, &bp_side, &mask_side, 4 );
        move16();

        IF( ep_enabled == 0 )
        {
            /* Discard padding length bytes */
            bp_side = sub( bp_side, padding_len );
            *total_padding = add( add( *total_padding, padding_len ), 2 );
            move16();
        }
        ELSE
        {
            *total_padding = add( *total_padding, 2 );
            move16();
            *np_zero = add( *np_zero, padding_len );
            move16();
        }

        /* test if we have less than 20 bytes left; if so frame is broken */
        IF( sub( sub( nbbytes, add( *total_padding, *np_zero ) ), 20 ) < 0 )
        {
            return 1;
        }

        /* Read bandwidth bits */
        IF( BW_cutoff_bits > 0 )
        {
            read_indice( ptr, &bp_side, &mask_side, BW_cutoff_bits );
            move16();
        }

        lastnz = read_indice( ptr, &bp_side, &mask_side, nbits );
        move16();
    }

    IF( ep_enabled != 0 )
    {
        *total_padding = add( *total_padding, *np_zero );
        move16();
    }
    return 0;
}
#  endif

static Word16 ac_dec_split_st2VQ_CW(                      /* local BER flag */
                                            const Word32 L_cwRx, /* max 25 bits */
                                            const Word32 L_szA,
                                            const Word32 L_szB,
                                            Word32* L_cwA,
                                            Word32* L_cwB,
                                            Word16* submodeLSB )
{
    /* demultiplex:  L_cwRx =   L_cwB(21.z bits) * L_szA(3.y bits)   + L_cwA(21.x bits)); */
    Word16 start, fin, ind;
    Word32 L_tmp, L_max_size;
    Counter i;

    L_max_size = (Word32) UL_Mpy_32_32( (UWord32) L_szB, (UWord32) L_szA ); /*  may be tabled  */

    /* section B  ind larger than 13  out of the possible  14 =   0..13  */
    IF( L_sub( L_cwRx, L_max_size ) >= 0 )
    {
        *L_cwA = L_deposit_l( 0 );
        *L_cwB = L_deposit_l( 0 );
        *submodeLSB = 0;
        move16();
        return (Word16) 1; /* set berFlag and exit */
    }

    /*initial binary split of cw,  select top or low half */
    start = 0;
    move16();

    ASSERT( ( L_szB & 0x1L ) == 0 ); /* this middle split only works if  L_szB is even  */
    if ( L_sub( L_cwRx, L_shr_pos( L_max_size, 1 ) ) >= 0 )
    {
        start = L_shr_pos( L_szB, 1 ); /* top half start index */
    }

    /*linear loop over a low  or a  high section */
    ind = start;
    move16();
    L_tmp = L_negate( L_cwRx ); /* search from negative side */

    L_tmp = L_add( L_tmp, (Word32) UL_Mpy_32_32( UL_deposit_l( (UWord16) start ), (UWord32) L_szA ) );
    /* start is 0 or 7 */ /*non-fractional mult is   (int)start * L_szA */

    /* a short linear run  over  ceil(szB/2) =  7   values  */

    fin = add( start, shr_pos( L_szB, 1 ) );
    FOR( i = start; i < fin; i++ )
    {
        ind = add( ind, 1 );
        L_tmp = L_add( L_tmp, L_szA );
        if ( L_tmp > 0 )
        {
            ind = sub( ind, 1 ); /* passed criteria point, keep index    */
        }
    }

    *L_cwB = L_deposit_l( ind );
    *L_cwA = L_sub( L_cwRx, (Word32) UL_Mpy_32_32( UL_deposit_l( (UWord16) ind ),
                                                   (UWord32) L_szA ) ); /* non-fractional mult;   (int)ind * L_szA */

    ASSERT( *L_cwA >= 0 && *L_cwA < L_szA );
    ASSERT( *L_cwB >= 0 && *L_cwB < L_szB );

    *submodeLSB = 0;
    *L_cwB = L_sub( *L_cwB, 2 );
    if ( *L_cwB < 0 )
    {
        *submodeLSB = 1;
        move16();
    }
    *L_cwB = L_mac0( *L_cwB, 2, *submodeLSB ); /* add back gain ind if needed */

    return 0; /* no BER */
}

#ifdef CR9_C_ADD_1p25MS_LRSNS

void readSNSData_fx(UWord8* ptr,
    Word16* bfiPtr,
    Word16* mask_side,
    Word16* bp_side,
    Word16* ltpf_idx_2_lrsns,
    Word32* L_scf_idx,
    LC3PLUS_FrameDuration frame_dms)
{
    Dyn_Mem_Deluxe_In(
        Word32 i, tmp32, tmp32lim;
    Word16  submodeMSB, submodeLSB;
    Word16  L, BER_detect;
    Word16  read_legacy_sns_vq_bits_fx;
    Word16  shape_idx, gain_idx, aux_idx, tmp_shape, stop_bit;
    Word16  plc_trigger_SNS1, plc_trigger_SNS2;

    );

    BER_detect = 0;  move16();
    plc_trigger_SNS1 = 1; move16();
#ifdef    LRSNS_10MS_BFISIGNAL_FIX
    plc_trigger_SNS2 = 1; move16();
#else 
    plc_trigger_SNS2 = 2; move16();
#endif 
    read_legacy_sns_vq_bits_fx = 1;  move16();
    IF(sub(frame_dms, LC3PLUS_FRAME_DURATION_1p25MS) == 0)
    {
        read_legacy_sns_vq_bits_fx = 0;   move16();  /*   decode   9,  10, or   29/30 bits */
    }

    IF(read_legacy_sns_vq_bits_fx != 0)
    {
        /* Decode SNS VQ parameters - 1st stage (10 bits) */
        L = read_indice(ptr, bp_side, mask_side, 5 + 5);
        L_scf_idx[0] = L_deposit_l(s_and(L, 0x1F)); /* stage1 LF  5  bits */
        L_scf_idx[1] = L_deposit_l(shr_pos(L, 5));  /* stage1 HF  5 bits  */


        /* Decode SNS VQ parameters - 2nd stage side-info (3-4 bits) */
        submodeMSB = read_bit(ptr, bp_side, mask_side); /* submodeMSB 1 bit */
        L_scf_idx[2] = L_deposit_l(shl_pos(submodeMSB, 1));
        ASSERT(sns_gainMSBbits[L_scf_idx[2]] > 0);
        L_scf_idx[3] = L_deposit_l(
            read_indice(ptr, bp_side, mask_side, sns_gainMSBbits[L_scf_idx[2]])); /* gains or gain MSBs  1-2 bits  */
        L_scf_idx[4] = read_bit(ptr, bp_side, mask_side);                        /*  shape LS 1 bit */

        /* Decode SNS VQ parameters - 2nd stage data (24-25 bits) */
        IF(submodeMSB == 0)
        {   /* shape_j = 0, or 1  */
            /* regular mode A,B indexes integer multiplexed, total 24.x bits  MPVQ codeword section A and  codeword for
             * section B */
             /* regular mode  mode shape  index   total  24.9999 bits    MPVQ codeword  */

            tmp32 = L_deposit_l(read_indice(ptr, bp_side, mask_side, 13));
            tmp32 = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 12), 13));
            move16(); /*for ber state   */
            BER_detect =
                ac_dec_split_st2VQ_CW(                                                                       /* local BER flag */
                    tmp32,                                                                /* L_cwRx  max 25 bits */
                    sns_MPVQ_Sz[0][0], UL_addNsD(sns_MPVQ_Sz[0][1], sns_MPVQ_Sz[1][1]), /* 12+2 = 14 */
                    (&L_scf_idx[5]),                                                    /* shape A */
                    (&L_scf_idx[6]),                                                    /* shape B or  gain LSB */
                    &submodeLSB                                                           /* total submode update below  */
                );
            IF(submodeLSB != 0)
            { /* add gainLSB bit */
                L_scf_idx[3] = L_add(L_shl_pos(L_scf_idx[3], 1), L_scf_idx[6]);
                L_scf_idx[6] = -2L;
            }
        }
        ELSE
        { /* shape_j = 2 or 3  */
            ASSERT(submodeMSB == 1);
        /* outlier mode shape  index   total  23.8536 +  19.5637 (19.5637 < (log2(2.^24 -2.^23.8537))    bits    MPVQ
         * codeword  */
        tmp32 = L_deposit_l(read_indice(ptr, bp_side, mask_side, 12));
        tmp32 = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 12), 12));

        L_scf_idx[5] = tmp32;   move32(); /*shape outl_near or outl_far */
        submodeLSB = 0;   move16();
        BER_detect = 0;   move16();
        tmp32lim = L_add(sns_MPVQ_Sz[2][0], L_shl_pos(sns_MPVQ_Sz[3][0], 1));

        IF (L_sub(tmp32, tmp32lim) >= 0)
        {
            BER_detect = 1;  move16();
        }
        ELSE
        {
            tmp32 = L_sub(tmp32, sns_MPVQ_Sz[2][0]); /*  a potential high index is computed */
            IF (tmp32 >= 0)
            {
                submodeLSB = 1;  move16();
                ASSERT(tmp32 >= 0 && tmp32 < (Word32)(2 * sns_MPVQ_Sz[3][0]));
                L_scf_idx[3] = L_add(L_shl_pos(L_scf_idx[3], 1), L_and(tmp32, 0x1)); /* add LSB_gain bit to gain MSBs */
                L_scf_idx[5] = L_shr_pos(tmp32, 1); /* MPVQ index with offset and gainLSB removed */
                L_scf_idx[6] = -2L;  move32();
            }
            ELSE
            {
                L_scf_idx[6] = -1L;  move32();
            }
        }
    }
    L_scf_idx[2] =
        L_add(L_scf_idx[2], L_deposit_l(submodeLSB)); /* decoder internal signal shape_j = submode 0..3 to VQ */


#ifdef   LRSNS_10MS_BFISIGNAL_FIX
    IF( BER_detect > 0)
    {
                *bfiPtr = 1;  move16();
                Dyn_Mem_Deluxe_Out();
                return;
    }
#else 
     *bfiPtr = BER_detect; move16();
     IF(*bfiPtr != 0)
     {
                Dyn_Mem_Deluxe_Out();
                return;
     }
#endif 
    }
    ELSE
    {
       ASSERT(read_legacy_sns_vq_bits_fx == 0 );
       /* lrsns   9/10/29/30  */
        /* init auxiliary demuxing variabls */
        shape_idx = -1; move16();
        gain_idx = -1;  move16();
        aux_idx = -1;   move16();
        tmp_shape = -1; move16();
        stop_bit = -1;  move16();

        FOR(i = 0; i < SCF_MAX_PARAM; i++)
        {
            L_scf_idx[i] = L_sub(-32000, i);  move32(); /* init parameters to be fwded to LRSNS VQ demultiplexor */
        }

        /* start actual Q-mode and fractional bits demultiplexing */

        /*  SNS-VQ 1st stage , 3 sections of 7.4 bits is stored in the first 9 bits */
        L_scf_idx[0] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 9));

        IF(L_sub(L_scf_idx[0], 510) >= 0)
        {
            ASSERT(L_scf_idx[0] < 512);
            L_scf_idx[0] -= 510;     move32();   /* sent only  idx 0,1 */
            shape_idx = -9;          move16();
            L_scf_idx[2] = shape_idx;  move32(); /* actual signal to  LRSNS decoder vector reconstruction  */
        }
        ELSE
        {
            /* read stop bit, 10th bit */
            stop_bit = read_indice(ptr, bp_side, mask_side, 1);

            test(); test();
            IF((L_sub(L_scf_idx[0] , (2 * 170)) < 0) && (stop_bit != 0))
            {
                /*keep values  0...339 in sns_vq_idx[0] , so that B vs C can be determined later in DecLR_fx function */
                L_scf_idx[2] = -10L;      move32();
                L_scf_idx[3] = L_deposit_l(ltpf_idx_2_lrsns[0]);  /*LTP active flag */
                L_scf_idx[4] = L_deposit_l(ltpf_idx_2_lrsns[1]);  /*LTPF active flag  */
            }
            ELSE
            { /* stage1B* + stage2*/
                /*0...169 in sns_vq_idx[0]*/
                test(); test();
                IF((L_sub(L_scf_idx[0], (2 * 170)) < 0) && (stop_bit == 0))
                {
                    aux_idx = 0;   move16();   /* typically a leading or first sign is stored in aux_idx */

                    IF(L_sub(L_scf_idx[0], 170) >= 0)
                    {
                        aux_idx = 1; move16();
                        L_scf_idx[0] = L_sub(L_scf_idx[0], 170);  move32();
                    }
                    L_scf_idx[1] = aux_idx;  move32();   /* forward aux bit for , LR_Split_LF, 29 bits  */

                    shape_idx = 0; move16(); /* point to splitLF parsing */
                    L_scf_idx[2] = L_deposit_l(shape_idx);           move32();

                    gain_idx = read_indice(ptr, bp_side, mask_side, 2);
                    L_scf_idx[3] = L_deposit_l(gain_idx); move32();

                    /* stage2 shape demux for LR_splitLF */

                    L_scf_idx[4] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 10));/* 10bits mPVQ(N=5,K=6) */

                    IF(L_sub(L_scf_idx[4] , (SNSLR_NPVQ_L5K6 >> 1) + (1 << 5)) >= 0) /* some limited  bit error detection possible here  */
                    {
                        *bfiPtr = plc_trigger_SNS1;  move16();
                        IF(*bfiPtr != 0)
                        {
#ifdef   LRSNS_10MS_BFISIGNAL_FIX
                            ASSERT(*bfiPtr == 1);
#endif
                            Dyn_Mem_Deluxe_Out();
                            return;
                        }
                    }

                    /* determine section of splitLF   mpvq(5,6)+P(8,2)+P(2,0) or  mpvq(5,8)+P(10,0)  */
                    IF(L_sub(L_scf_idx[4] , (SNSLR_NPVQ_L5K6 >> 1)) < 0)
                    {
                        tmp_shape = read_indice(ptr, bp_side, mask_side, 1);/* LS (8,2) */
                        L_scf_idx[5] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 6)); move32();  /* mPVQ(8,2) */
                        L_scf_idx[5] = L_mac0(L_shl_pos(L_scf_idx[5], 1), 1, tmp_shape);  move32(); /* P(8,2) LS put as lsb */
                    }
                    ELSE
                    {
                        L_scf_idx[4] = L_sub(L_scf_idx[4], (SNSLR_NPVQ_L5K6 >> 1));  move32();     /* 5 lsbs of mpvq (5,8) */

                        L_scf_idx[5] = L_deposit_l(read_indice(ptr, bp_side, mask_side, 7)); move32();/*  7 msbs of mPVQ(5,8) */

                        L_scf_idx[4] = L_or(L_shl_pos(L_scf_idx[5], 5),  L_scf_idx[4]);    move32();   /*  binary merge MSB's and LSBs  */
                        L_scf_idx[5] = -8L;   move32(); /* signal to sns_decoder split_LF subshape  to decode 8 lf pulses, and no hf pulses */

                        IF(L_sub(L_scf_idx[4] , (SNSLR_NPVQ_L5K8 >> 1)) >= 0)
                        {
                            *bfiPtr = plc_trigger_SNS2;  move16();
                            IF( *bfiPtr != 0 )
                            {  
#ifdef   LRSNS_10MS_BFISIGNAL_FIX
                                ASSERT(*bfiPtr == 1);
#endif
                                Dyn_Mem_Deluxe_Out();
                                return;
                            }
                        }
                    }
                }
                ELSE IF(L_sub(L_scf_idx[0], 2 * 170) >= 0)
                {
                    aux_idx = stop_bit; move16();
                    L_scf_idx[0] = L_sub(L_scf_idx[0], 2 * 170); move32();
                    L_scf_idx[1] = L_deposit_l(aux_idx);          move32();


                    shape_idx = 1; move16();/* point to full  parsing */
                    L_scf_idx[2] = shape_idx; move32();/* LR_full , 30 bits  */

                    gain_idx = read_indice(ptr, bp_side, mask_side, 3);
                    L_scf_idx[3] = gain_idx; move32();

                    /* stage2 shape demux for LR_full  */
                    tmp32 = L_deposit_l(read_indice(ptr, bp_side, mask_side, 10));     /* 10 LSBs. total 16.666  bits mPVQ(N=15,K=5) */
                    tmp32 = L_or(tmp32, L_shl_pos(read_indice(ptr, bp_side, mask_side, 7), 10)); /*7 MSBs */
                    L_scf_idx[4] = tmp32; move32();
                    IF(L_sub(L_scf_idx[4],  (SNSLR_NPVQ_L15K5 >> 1)) >= 0)
                    {  /* fixenv shapes demultiplexing */
                        L_scf_idx[5] = L_sub(L_scf_idx[4], (SNSLR_NPVQ_L15K5 >> 1));   move32();

                        IF(L_sub(L_scf_idx[5] , 3 * (1 << 13)) < 0)
                        {   /*fix_env's "0,1,2"  with 2 shiftbits and 11 remaining sign bits s1..s11  */
                           L_scf_idx[4] = 0L;                  move32();
                           WHILE(L_sub(L_scf_idx[5], (1 << 13)) >= 0)
                           {
                                L_scf_idx[5] = L_mac0(L_scf_idx[5], -1, (1 << 13));  move32();
                                L_scf_idx[4] = L_add(L_scf_idx[4], 1); move32();
                            }
                            assert(L_scf_idx[4] >= 0 && L_scf_idx[4] <= 3);
                            assert(L_scf_idx[5] >= 0 && L_scf_idx[5] < (1 << 13));
                        }
                        ELSE IF(L_sub(L_scf_idx[5],  3 * (1 << 13) + (1 << 11)) < 0)
                        {
                            L_scf_idx[4] = 3L;  move32(); /*smaller fix_env "3"  with 2 shiftbits and 9 remaining sign bits s1..s9  */
                            L_scf_idx[5] = L_mac0(L_scf_idx[5] , -1,  3 * (1 << 13));   move32();
                            assert(L_scf_idx[5] >= 0 && L_scf_idx[5] < (1 << 11));
                        }
                        ELSE
                        {
                            /* unused section indicate bit error */
                            *bfiPtr = plc_trigger_SNS2; move16();
                            test();
                            IF( *bfiPtr != 0 )
                            {     
                                Dyn_Mem_Deluxe_Out();
                                return;
                            }
                        }
                        shape_idx = add(extract_l(L_scf_idx[4]), 2); move32();
                        L_scf_idx[2] = L_deposit_l(shape_idx);           move32();
                    } /* fixenv */
                } /*full*/
            } /*stage1B* + stage2 */
        } /*10+ bits*/
    }
#ifdef  LRSNS_PC_SIGNAL_FIX  
    assert(*bfiPtr == 0 || *bfiPtr == 1); /* local SNS BFI-flag  output check */
#endif 
    Dyn_Mem_Deluxe_Out();
} /*read SNS*/
#endif /* LRSNS */



#ifdef NEW_SIGNALLING_SCHEME_1p25

void readLtpData_fx(
    UWord8* ptr,
    Word16* bfiPtr,
    Word16* mask_side,
    Word16* bp_side,
    Word16* ltpf_idx,
    Word16* rx_status,
    Word16*  ltpfinfo_frame_cntr_fx,
    Word16* mem_continuation
)
{
    Word16 rx_current_status = -1;
    Word16 tmp, MSBs, LSBs;

    ltpf_idx[2] = -1;  move16();   /* -1 indicates incomplete lag,  conditionally decoded if phase is B , and consecutive A/B has arrived */

    tmp = read_indice(ptr, bp_side, mask_side, 2);

    test();
    IF(tmp == 0)
    {
        ltpf_idx[0] = 0;   move16(); /* ltp ltpf/lag was not transmitted */
        ltpf_idx[1] = 0;   move16(); /* ltpf  activation bit zeroed    */

        /* *ltp_bits_fx = 2; */  /* note: ltpbits  bitbudget not really used in decoder */

        rx_status[0] = -32768;  move16();  /* set unknown phase A , due to rxLTP==0 */
        rx_status[1] = -1;      move16(); /* set unknown phase A MSBs content       */
        *ltpfinfo_frame_cntr_fx = -32768;  move16();
        ASSERT(ltpf_idx[2] < 0); /* ltpf_idx[2] = -1; , no ready lag available */
#ifdef FIX_LTPF_1p25
        *mem_continuation = 0; move16();  /* also kill lag continuation state */
#endif
    }
    ELSE IF(sub(tmp, 1) == 0)
    {
        ltpf_idx[0] = 1;
        ltpf_idx[1] = 0;  /* LTP=1, LTPF=0,  inactive ltpf */
        rx_current_status = read_bit(ptr, bp_side, mask_side);

        test();
        IF(rx_current_status == 0)
        {
            rx_status[0] = 0; move16(); /* phaseA */
            rx_status[1] = read_indice(ptr, bp_side, mask_side, 4);   /* read four MSBs, and store in rx_status[1] */
#    ifdef FIX_LTPF_1p25
            test();
            if (*mem_continuation == 0)
            {
                *mem_continuation = 1;
            }
#    endif
            *ltpfinfo_frame_cntr_fx = 0;  /* handle longer loss bursts  */
        }
        ELSE
        { /* LSB part of delta coded lag information */
           ASSERT(rx_current_status == 1);
           LSBs = shl(read_indice(ptr, bp_side, mask_side, 4), 1);  /* NB  LSB  is on purpose always zero, truncation on encoder side   */
           IF ( rx_status[1] < 0 )
            {
                *bfiPtr = 1; move16();
                return;
            }
           ltpf_idx[2] = s_or(shl(rx_status[1], 5), LSBs);

           /* check  frame cntr info to not combine oldA with a newB */
            IF(sub(*ltpfinfo_frame_cntr_fx, 1) != 0)
            {
                ltpf_idx[1] = 0; move16();   /*turn of LTPF, even number of bfi frames may have happened */
                ltpf_idx[2] = -1; move16(); /*  indicate bfi burst and corrupt lagLSBs to PLC and ltpf_decoder_fx       */
            }
#    ifdef FIX_LTPF_MEM_CONTINUATION
            else
            {
                *mem_continuation = 0; move16();
            }
 #    endif
            rx_status[0] = -32768; move16();
            *ltpfinfo_frame_cntr_fx = -32678; move16();
        }
    }
    ELSE
    {   /*2 or 3*/
        ltpf_idx[0] = 1;  move16();
        ltpf_idx[1] = 1;  move16();     /* active ltpf */

        IF(sub(tmp, 2) == 0)
        {
            /* phaseA */
            MSBs = read_indice(ptr, bp_side, mask_side, 4);
            rx_status[0] = 0; move16();
            rx_status[1] = MSBs; move16();/* remember the four MSBs */
#    ifdef FIX_LTPF_1p25
            test();
            if (*mem_continuation == 0)
            {
                *mem_continuation = 1; move16();
            }
#    endif
            *ltpfinfo_frame_cntr_fx = 0; move16();
        }
        ELSE
        {
            ASSERT(tmp == 3); /*   phaseB */
            LSBs = read_indice(ptr, bp_side, mask_side, 5);  /*  all 5 LSBs  available*/
            IF ( rx_status[1] < 0 )
            {
                *bfiPtr = 1; move16();
                return;
            }
            ltpf_idx[2] = s_or(shl(rx_status[1], 5), LSBs);

            /* check  frame cntr info to not combine oldA  MSBs with a newB LSBs  */
            IF(sub(*ltpfinfo_frame_cntr_fx, 1) != 0)
            {
                ltpf_idx[1] = 0;    move16();   /*  turn off LTPF activation,  ltpf_idx[2] is not read */
                ltpf_idx[2] = -1;    move16();  /*  indicate bfi burst and corrupt lagLSBs to PLC and ltpf_decoder_fx       */
            }
            *ltpfinfo_frame_cntr_fx = -32678; move16(); /*cntr init in phaseA*/
            rx_status[0] = -32768; move16();        /* phase init in phaseA*/

#    ifdef FIX_LTPF_MEM_CONTINUATION
            *mem_continuation = 0;  move16();
#    endif
        }
    }
}

#endif
