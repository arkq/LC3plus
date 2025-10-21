/******************************************************************************
*                        ETSI TS 103 634 V1.6.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"

static const LC3_INT32 gainMSBbits[4] = {1, 1, 2, 2};

static void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit);
static void read_uint_fl(LC3_INT nbits, LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* val);

#ifdef CR9_C_ADD_1p25MS_LRSNS
void readSNSData_fl(LC3_UINT8* ptr, LC3_INT32* bfi, LC3_INT32* mask_side_local, LC3_INT32* bp_side_local, const LC3_INT32* ltpf_idx_2_lrsns, LC3_INT32* sns_vq_idx, LC3PLUS_FrameDuration frame_dms, LC3_INT32 plc_trigger_SNS1, LC3_INT32 plc_trigger_SNS2);
#else
void readSNSData_fl(LC3_UINT8* ptr, LC3_INT32* bfi, LC3_INT32* mask_side_local, LC3_INT32* bp_side_local, LC3_INT32* sns_vq_idx, LC3PLUS_FrameDuration frame_dms, LC3_INT32* plc_trigger_SNS1, LC3_INT32* plc_trigger_SNS2);
#endif

 

#ifdef NEW_SIGNALLING_SCHEME_1p25
void readLtpData_fl(LC3_UINT8* ptr, LC3_INT32* bfiPtr, LC3_INT32* mask_side, LC3_INT32* bp_side, LC3_INT32* ltpf_idx, LC3_INT32* rx_status, LC3_INT32* ltpfinfo_frame_cntr, LC3_INT16* mem_continuation);
#endif
 
void read_bit_fl(LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* bit)
{
    if (ptr[*bp_side] & *mask_side) {
        *bit = 1;
    } else {
        *bit = 0;
    }

    if (*mask_side == 128) {
        *mask_side = 1;
        *bp_side   = *bp_side - 1;
    } else {
        *mask_side = *mask_side * 2;
    }
}

void read_uint_fl(LC3_INT nbits, LC3_UINT8* ptr, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT* val)
{
    LC3_INT bit, i;

    read_bit_fl(ptr, mask_side, bp_side, val);

    for (i = 1; i < nbits; i++) {
        read_bit_fl(ptr, mask_side, bp_side, &bit);
        *val = *val + (bit << i);
    }
}

#ifdef ENABLE_PADDING
LC3_INT paddingDec_fl(LC3_UINT8* bytes, LC3_INT nbbits, LC3_INT L_spec, LC3_INT bw_cutoff_bits, LC3_INT ep_enabled, LC3_INT* total_padding, LC3_INT *np_zero)
{
    LC3_INT lastnz_threshold;
    LC3_INT val, padding_len_bits, padding_len;
    LC3_INT bp_side;
    LC3_INT    mask_side;
    LC3_UINT8* ptr = bytes;

    LC3_INT nbbytes = nbbits >> 3;
    LC3_INT lastnz;
    LC3_INT bw_cutoff_idx;
    LC3_INT nbits = getLastNzBits (L_spec);

    if (nbits > nbbits)
    {
        return 1;
    }

    *np_zero = 0;

    *total_padding = 0;

    bp_side   = (nbbits - 1) >> 3;
    mask_side = 1 << (8 - (nbbits - (bp_side << 3)));

    if (bp_side < 19 || bp_side >= LC3PLUS_MAX_BYTES) {
        return 1;
    }

    ptr = bytes;

    if (bw_cutoff_bits > 0) {
        read_uint_fl(bw_cutoff_bits, ptr, &mask_side, &bp_side, &bw_cutoff_idx);
    }

    read_uint_fl(nbits, ptr, &mask_side, &bp_side, &lastnz);

    lastnz_threshold = (1 << nbits) - 1 - 1;

    while (lastnz == lastnz_threshold) {
        padding_len_bits = 16 - nbits - bw_cutoff_bits - 4;

        /*Read padding length*/
        read_uint_fl(padding_len_bits, ptr, &mask_side, &bp_side, &padding_len);

        /* Read 4 reserved bits */
        read_uint_fl(4, ptr, &mask_side, &bp_side, &val);

        if (ep_enabled == 0)
        {
            /* Discard padding length bytes */
            bp_side        = bp_side - padding_len;
            *total_padding = *total_padding + padding_len + 2;
        }
        else
        {
            *total_padding = *total_padding + 2;
            *np_zero       = *np_zero + padding_len;
        }

        /* check if minimum payload size is reached */
        if ((nbbytes - (*total_padding + *np_zero)) < 20) {
            return 1;
        }

        /* Read bandwidth bits */
        if (bw_cutoff_bits > 0) {
            read_uint_fl(bw_cutoff_bits, ptr, &mask_side, &bp_side, &bw_cutoff_idx);
        }

        read_uint_fl(nbits, ptr, &mask_side, &bp_side, &lastnz);
    }

    if (ep_enabled != 0)
    {
        *total_padding = *total_padding + *np_zero;
    }

    return 0;
}
#endif


#ifdef NEW_SIGNALLING_SCHEME_1p25
void readLtpData_fl(
    LC3_UINT8* ptr,
    LC3_INT32* bfiPtr,
    LC3_INT32* mask_side,
    LC3_INT32* bp_side,
    LC3_INT32* ltpf_idx,
    LC3_INT32* rx_status,
    LC3_INT32* ltpfinfo_frame_cntr,
    LC3_INT16* mem_continuation
)
{
    LC3_INT32 rx_current_status = -1;
    LC3_INT32 tmp, MSBs, LSBs;

    /*        Hdr, information               , bits used
              00 , no lag info , no phase info   sum=2
              010, PhaseA,LTPF=0, lagAbits=4 ,   sum=7    : PLC may be activated, 4 MSbs
              011, PhaseB,LTPF=0, lagBbits=4*,   sum=7*   : PLC may be activated, 4* =  reduced lag resolution in Q_ltpf_Idx domain for PLC-activation
              10 , PhaseA,LTPF=1, lagAbits=4 ,   sum=6    : LTPF activated
              11 , PhaseB,LTPF=1, lagBbits=5 ,   sum=7    : LTPF activated
    */

    ltpf_idx[2] = -1; /* no ready lag available, conditionally decoded if phase is B, and consecutive A/B has arrived */

    read_uint_fl(2, ptr, mask_side, bp_side, &tmp);
    if (tmp == 0) /* "00" */
    {
        ltpf_idx[0] = 0;    /* ltp ltpf/lag was not transmitted */
        ltpf_idx[1] = 0;    /* ltpf  activation bit zeroed    */

        rx_status[0] = -32768;   /* set unknown phase A , due to rx LTP==0 */
        rx_status[1] = -1;       /* set unknown phase A MSBs content       */
        *ltpfinfo_frame_cntr = -32768;
        assert(ltpf_idx[2] < 0); /* ltpf_idx[2] = -1; , no ready lag available */
#    ifdef FIX_LTPF_1p25
        *mem_continuation = 0;   /* also kill lag continuation state */
#    endif
    }
    else if (tmp == 1) /* "01" */
    {
        ltpf_idx[0] = 1;
        ltpf_idx[1] = 0;  /* LTP=1, LTPF=0,  inactive ltpf */
        read_bit_fl(ptr, mask_side, bp_side, &rx_current_status);

        if (rx_current_status == 0)
        {
            rx_status[0] = 0;   /* phaseA */
            read_uint_fl(4, ptr, mask_side, bp_side, &(rx_status[1])); /* read four MSBs, and store in rx_status[1] */
#    ifdef FIX_LTPF_1p25
            if (*mem_continuation == 0)
            {
             *mem_continuation = 1;
            }
#    endif
            *ltpfinfo_frame_cntr = 0;  /*same as rx_status    [0] */
        }
        else
        { /* LSB part of delta coded lag information */
            assert(rx_current_status == 1);
            read_uint_fl(4, ptr, mask_side, bp_side, &LSBs);
            LSBs = (LSBs << 1);   /* NB  Least Signifcant bit  is on purpose always zero, truncation on encoder side   */
            if ( rx_status[1] < 0 )
            {
                *bfiPtr = 1;
                return;
            }
            ltpf_idx[2] = ((rx_status[1] << 5) | LSBs); /* bitwise OR */

            /* check  frame cntr info to not combine oldA with a newB */
            if (*ltpfinfo_frame_cntr != 1)
            {
                rx_status[0] = -32768;   /*even number of bfi frames may have happened */
                ltpf_idx[1] = 0;
                ltpf_idx[2] = -1;     /*   send signal of non-decoded lag to PLC and LTPF decoder  */
            }
#    ifdef FIX_LTPF_MEM_CONTINUATION
            else
            {
                *mem_continuation = 0;
            }
 #    endif
            rx_status[0] = -32768;
            *ltpfinfo_frame_cntr = -32678;
        }
    }
    else
    {   /* 2 or 3 */
        ltpf_idx[0] = 1;
        ltpf_idx[1] = 1;   /* active ltpf */

        if (tmp == 2) /* 2="10" */
        {
            /* phaseA */
            read_uint_fl(4, ptr, mask_side, bp_side, &MSBs);

            rx_status[0] = 0;
            rx_status[1] = MSBs; /* remember the four MSBs */
#    ifdef FIX_LTPF_1p25
            if (*mem_continuation == 0)
            {
                *mem_continuation = 1;
            }
#    endif
            *ltpfinfo_frame_cntr = 0;
            assert(ltpf_idx[2] < 0); /* ltpf_idx[2] = -1; , no ready lag available */
        }
        else
        {   /*  3="11" */
            assert(tmp == 3); /*   phaseB */
            read_uint_fl(5, ptr, mask_side, bp_side, &LSBs);  /*  all 5 LSBs  available*/
            if ( rx_status[1] < 0 )
            {
                *bfiPtr = 1;
                return;
            }
            ltpf_idx[2] = ((rx_status[1] << 5) | LSBs); /* bitwise OR  */

            /* check  frame cntr info to not combine oldA with a newB */
            if (*ltpfinfo_frame_cntr != 1)
            {
                ltpf_idx[1] = 0;     /*  turn off LTPF activation for now, so that  ltpf_idx[2] is not read */
                ltpf_idx[2] = -1;    /*  send signal to PLC and ltpf_decoder, that phase B could not be decoded     */
            }
            *ltpfinfo_frame_cntr = -32678;  /*cntr init in phaseA*/
            rx_status[0] = -32768;         /* phase init in phaseA*/

#    ifdef FIX_LTPF_MEM_CONTINUATION
            *mem_continuation = 0;
#    endif
        }
    }

}
#endif


#ifdef CR9_C_ADD_1p25MS_LRSNS
void readSNSData_fl(LC3_UINT8* ptr, LC3_INT32* bfi, LC3_INT32* mask_side_local, LC3_INT32* bp_side_local, const LC3_INT32* ltpf_idx_2_lrsns,
    LC3_INT32* sns_vq_idx, LC3PLUS_FrameDuration frame_dms,
    LC3_INT32 plc_trigger_SNS1, LC3_INT32 plc_trigger_SNS2)
#else
void readSNSData_fl(LC3_UINT8* ptr, LC3_INT32* bfi, LC3_INT32* mask_side_local, LC3_INT32* bp_side_local, LC3_INT32* sns_vq_idx,
    LC3PLUS_FrameDuration frame_dms, LC3_INT32* plc_trigger_SNS1, LC3_INT32* plc_trigger_SNS2)
#endif
{
    LC3_INT32 tmp, submodeMSB, idxBorGainLSB, submodeLSB;
#ifndef CR9_C_ADD_1p25MS
    (void) frame_dms;
#endif
#ifdef CR9_C_ADD_1p25MS_LRSNS
    LC3_INT32 read_legacy_sns_vq_bits;
#endif
    

    UNUSED(frame_dms);

#ifdef CR9_C_ADD_1p25MS_LRSNS
    read_legacy_sns_vq_bits = 1;
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {
        read_legacy_sns_vq_bits = 0;  /*   decode   9,  10, or   29/30 bits */
    }

    if (read_legacy_sns_vq_bits != 0)
    {
#  endif
        /* SNS-VQ 1st stage */
        read_uint_fl(5, ptr, mask_side_local, bp_side_local, &sns_vq_idx[SNS_IDX_LF]);
        read_uint_fl(5, ptr, mask_side_local, bp_side_local, &sns_vq_idx[SNS_IDX_HF]);

        /* SNS-VQ 2nd stage side-info (3-4 bits) */
        read_bit_fl(ptr, mask_side_local, bp_side_local, &submodeMSB);

        read_uint_fl(gainMSBbits[submodeMSB * 2], ptr, mask_side_local, bp_side_local, &sns_vq_idx[SNS_IDX_GAIN]);
        read_bit_fl(ptr, mask_side_local, bp_side_local, &sns_vq_idx[SNS_IDX_LS_INDA]);

        /* SNS-VQ 2nd stage VQ decoding (24-25 bits) */
        if (submodeMSB == 0)
        {
            read_uint_fl(25, ptr, mask_side_local, bp_side_local, &tmp);
            if (tmp >= 33460056)
            {
#ifdef CR9_C_ADD_1p25MS_LRSNS
                *bfi = plc_trigger_SNS1;
#else
                *bfi = *plc_trigger_SNS1;
#endif
                if (*bfi)
                {
                    return;
                }
            }
            idxBorGainLSB = floor(tmp / 2390004);
            sns_vq_idx[SNS_IDX_A] = tmp - idxBorGainLSB * 2390004;

            if (idxBorGainLSB < 2)
            {
                submodeLSB = 1;
                sns_vq_idx[SNS_IDX_GAIN] = sns_vq_idx[SNS_IDX_GAIN] * 2 + idxBorGainLSB;
                sns_vq_idx[SNS_IDX_BORGAINLSB] = -2;
            }
            else
            {
                submodeLSB = 0;
                sns_vq_idx[SNS_IDX_BORGAINLSB] = idxBorGainLSB - 2;
            }
        }
        else
        {
            read_uint_fl(24, ptr, mask_side_local, bp_side_local, &tmp);
            if (tmp >= 16708096)
            {
#ifdef CR9_C_ADD_1p25MS_LRSNS
                *bfi = plc_trigger_SNS2;
#else
                *bfi = *plc_trigger_SNS2;
#endif
                if (*bfi)
                {
                    return;
                }
            }

            if (tmp >= 15158272)
            {
                submodeLSB = 1;
                tmp -= 15158272;
                sns_vq_idx[SNS_IDX_GAIN] = sns_vq_idx[SNS_IDX_GAIN] * 2 + (tmp & 1);
                sns_vq_idx[SNS_IDX_A] = floor(tmp / 2);
                sns_vq_idx[SNS_IDX_BORGAINLSB] = -2;
            }
            else
            {
                submodeLSB = 0;
                sns_vq_idx[SNS_IDX_A] = tmp;
                sns_vq_idx[SNS_IDX_BORGAINLSB] = -1;
            }
        }

        sns_vq_idx[SNS_IDX_SHAPEJ] = submodeMSB * 2 + submodeLSB;

#ifdef CR9_C_ADD_1p25MS_LRSNS
    }

    if (read_legacy_sns_vq_bits == 0)
    {
        LC3_INT32 shape_idx = -1;
        LC3_INT32 gain_idx = -1;
        LC3_INT32 aux_idx = -1;
        LC3_INT32 tmp_shape = -1;
        LC3_INT32 stop_bit = -1;
        /* SNS-VQ 1st stage  in 9-10  bits */
        read_uint_fl(9, ptr, mask_side_local, bp_side_local, &sns_vq_idx[SNS_IDX_LF]);
        if (sns_vq_idx[SNS_IDX_LF] >= 510) /* stage 1A */
        {
            assert(sns_vq_idx[SNS_IDX_LF] < 512);
            sns_vq_idx[SNS_IDX_LF] -= 510;    /* send only  idx 0,1 */
            sns_vq_idx[SNS_IDX_HF] = -32768;  /* unused */
            shape_idx = -9;
            sns_vq_idx[2] = shape_idx; /* actual signal to  LR SNS vector reconstruction  */
        }
        else
        {
            /* read stop bit */
            read_uint_fl(1, ptr, mask_side_local, bp_side_local, &stop_bit);
            sns_vq_idx[SNS_IDX_HF] = -32768;  /* unused */

            if (sns_vq_idx[SNS_IDX_LF] < (2 * 170) && stop_bit != 0)
            {
                /*B or C , keep values  0...339 in sns_vq_idx[0] , so that stage1 B vs stage1 C can be determined later in the DecLR function */
                sns_vq_idx[2] = -10;
                sns_vq_idx[3] = ltpf_idx_2_lrsns[0];  /* forward LTP active flag */
                sns_vq_idx[4] = ltpf_idx_2_lrsns[1];  /* forward LTPF active flag  */
            }
            else
            { /* stage1B + stage2 */
                /*0...169 in sns_vq_idx[0]*/
                if (sns_vq_idx[SNS_IDX_LF] < (2 * 170) && stop_bit == 0)
                {
                    aux_idx = 0;
                    if (sns_vq_idx[SNS_IDX_LF] >= (170))
                    {
                        aux_idx = 1;
                        sns_vq_idx[SNS_IDX_LF] -= 170;
                    }
                    sns_vq_idx[SNS_IDX_HF] = aux_idx;   /* aux bit for , LR_Split_LF, 29 bits  */

                    shape_idx = 0; /* point to splitLF parsing */
                    sns_vq_idx[2] = shape_idx;

                    read_uint_fl(2, ptr, mask_side_local, bp_side_local, &gain_idx);
                    sns_vq_idx[3] = gain_idx;

                    /* stage2 shape demux for LR_splitLF */
                    read_uint_fl(10, ptr, mask_side_local, bp_side_local, &sns_vq_idx[4]); /* 10bits mPVQ(N=5,K=6) */

                    if (sns_vq_idx[4] >= (SNSLR_NPVQ_L5K6 >> 1) + (1 << 5)) /* some limited  bit error detection possible here  */
                    {
                        *bfi = plc_trigger_SNS1;
                        if (*bfi != 0)
                        {
                            return;
                        }
                    }

                    /* determine section of splitLF   mpvq(5,6)+P(8,2)+P(2,0) or  mpvq(5,8)+P(10,0)  */
                    if (sns_vq_idx[4] < (SNSLR_NPVQ_L5K6 >> 1))
                    {
                        read_uint_fl(1, ptr, mask_side_local, bp_side_local, &tmp_shape);      /* LS (8,2) */
                        read_uint_fl(6, ptr, mask_side_local, bp_side_local, &sns_vq_idx[5]); /* mPVQ(8,2) */
                        sns_vq_idx[5] = (sns_vq_idx[5] << 1) + tmp_shape; /* P(8,2) LS put as lsb */
                    }
                    else
                    {
                        sns_vq_idx[4] = sns_vq_idx[4] - (SNSLR_NPVQ_L5K6 >> 1);   /* 5 LSBs of mpvq (5,8) */

                        read_uint_fl(7, ptr, mask_side_local, bp_side_local, &sns_vq_idx[5]); /* 7 msbs of mPVQ(5,8) */
                        sns_vq_idx[4] = (sns_vq_idx[5] << 5) | sns_vq_idx[4]; /* merge MSB's and LSBs  */
                        sns_vq_idx[5] = -8; /* signal to sns_decoder split_LF subshape  to decode 8 lf pulses,  and no hf pulses */

                        if (sns_vq_idx[4] >= (SNSLR_NPVQ_L5K8 >> 1)) {
                            *bfi = plc_trigger_SNS1;
                            if (*bfi != 0) {
                                return;
                            }
                        }
                    }
                }
                else if (sns_vq_idx[SNS_IDX_LF] >= (2 * 170))
                {
                    aux_idx = stop_bit;
                    sns_vq_idx[SNS_IDX_LF] -= (2 * 170);
                    sns_vq_idx[1] = aux_idx;

                    shape_idx = 1; /* point to full  parsing */
                    sns_vq_idx[2] = shape_idx; /* LR_full , 30 bits  */

                    read_uint_fl(3, ptr, mask_side_local, bp_side_local, &gain_idx);
                    sns_vq_idx[3] = gain_idx;

                    /* stage2 shape demux for LR_full  */
                    read_uint_fl(17, ptr, mask_side_local, bp_side_local, &sns_vq_idx[4]); /* 16.666  bits mPVQ(N=15,K=5) */

                    if (sns_vq_idx[4] >= (SNSLR_NPVQ_L15K5 >> 1))
                    { /* fixenv shapes demultiplexing */
                        sns_vq_idx[5] = (sns_vq_idx[4] - (SNSLR_NPVQ_L15K5 >> 1));
                        if (sns_vq_idx[5] < (3 * (1 << 13)))
                        {   /*fix_env's "0,1,2"  with 2 shiftbits and 11 remaining sign bits s1..s11  */
                            sns_vq_idx[4] = 0;
                            while (sns_vq_idx[5] >= (1 << 13)) {
                                sns_vq_idx[5] = sns_vq_idx[5] - (1 << 13);
                                sns_vq_idx[4] += 1;
                            }
                            assert(sns_vq_idx[4] >= 0 && sns_vq_idx[4] <= 3);
                            assert(sns_vq_idx[5] >= 0 && sns_vq_idx[5] < (1 << 13));
                        }
                        else if (sns_vq_idx[5] < (3 * (1 << 13) + (1 << 11)))
                        {
                            sns_vq_idx[4] = 3;  /*smaller fix_env "3"  with 2 shiftbits and 9 remaining sign bits s1..s9  */
                            sns_vq_idx[5] = sns_vq_idx[5] - (3 * (1 << 13));
                            assert(sns_vq_idx[5] >= 0 && sns_vq_idx[5] < (1 << 11));
                        }
                        else
                        { /* bit error */
                            *bfi = plc_trigger_SNS2;
                            if (*bfi != 0)
                            {
                                return;
                            }
                        }
                        shape_idx = sns_vq_idx[4] + 2;
                        sns_vq_idx[2] = shape_idx;
                    } /* fixenv */
                }/*full*/
            }/*stage1B* + stage2 */
        } /*10+ bits*/
    }
#  endif
#ifdef  LRSNS_PC_SIGNAL_FIX  
    assert(*bfi == 0 || *bfi == 1 ); /* local SNS BFI-flag  output check */
#endif 
}

void processDecoderEntropy_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT N, LC3_INT fs_idx,
                              LC3_INT bw_cutoff_bits, LC3_INT* bfi, LC3_INT* gg_idx, LC3_INT* sns_vq_idx, LC3_INT* fac_ns_idx,
                              LC3_INT* tns_numfilters, LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* bw_cutoff_idx, LC3_INT* lastnz,
                              LC3_INT* lsbMode, LC3PLUS_FrameDuration frame_dms
#ifdef CR9_C_ADD_1p25MS
                              , LC3_INT32 rx_status[2], LC3_INT16* mem_continuation
#ifdef NEW_SIGNALLING_SCHEME_1p25
                               , LC3_INT32 *ltpfinfo_frame_cntr  /* set here , but  also increased outside  by  bfi for the channel */
#endif
#endif
                              )
{

#  ifdef CR9_C_ADD_1p25MS_LRSNS
        LC3_INT32 plc_trigger_bw, plc_trigger_last_nz, plc_trigger_SNS1, plc_trigger_SNS2, bit,
        i, ltpf_tmp[3], bp_side_local, mask_side_local, rx_current_status, ltpf_idx_2_lrsns[3];
    LC3_UINT8 * ptr;
#ifdef LRSNS_PC_SIGNAL_FIX  
    LC3_INT bfiSNS;
#endif

    //UNUSED(rx_status);
    //UNUSED(mem_continuation);
    UNUSED(rx_current_status);
#  else
        LC3_INT32 plc_trigger_bw, plc_trigger_last_nz, plc_trigger_SNS1, plc_trigger_SNS2, bit,
        i, ltpf_tmp[3], bp_side_local, mask_side_local, rx_current_status;
    LC3_UINT8 * ptr;
    UNUSED(rx_current_status);
#  endif

#ifdef NEW_SIGNALLING_SCHEME_1p25
        UNUSED(rx_current_status);
#endif

    *bp_side = -1;
    bp_side_local   = numbytes - 1; /* Matlab offset by 1 */
    mask_side_local = 1;
    *mask_side = -1;
    ptr        = bytes;
    *lsbMode = -1;
    *lastnz = -1;

    plc_trigger_bw      = 1; /* Bandwidth */
    plc_trigger_last_nz = 1; /* Last non-zero tuple */

    plc_trigger_SNS1 = 1; /* SNS-VQ 2nd stage MPVQ data (24-25 bits) */
#ifdef    LRSNS_10MS_BFISIGNAL_FIX
    plc_trigger_SNS2 = 1; /* SNS-VQ 2nd stage MPVQ data (10-16 bits) */
#else 
    plc_trigger_SNS2 = 2; /* SNS-VQ 2nd stage MPVQ data (10-16 bits) */
#endif 
    /* Bandwidth */
    if (bw_cutoff_bits > 0) {
        read_uint_fl(bw_cutoff_bits, ptr, &mask_side_local, &bp_side_local, bw_cutoff_idx);

        if (fs_idx < *bw_cutoff_idx) {
            *bfi = plc_trigger_bw;

            if (*bfi) {
                return;
            }
        }
    } else {
        *bw_cutoff_idx = fs_idx;
    }

    /* Number of TNS filters */
#ifdef CR9_C_ADD_1p25MS
    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS) {
        *tns_numfilters = 0;
    } else {
#endif
        if (*bw_cutoff_idx < 3 || frame_dms == LC3PLUS_FRAME_DURATION_2p5MS) {
            *tns_numfilters = 1;
        } else {
            *tns_numfilters = 2;
        }
#ifdef CR9_C_ADD_1p25MS
    }
#endif

    /* Last non-zero tuple */
    read_uint_fl(getLastNzBits (N), ptr, &mask_side_local, &bp_side_local, lastnz);
    *lastnz = (*lastnz + 1) * 2;

    if (*lastnz > N) {
        *bfi = plc_trigger_last_nz;
        if (*bfi) {
            return;
        }
    }

    /* LSB mode bit */
    read_bit_fl(ptr, &mask_side_local, &bp_side_local, lsbMode);

    /* Global gain */
    read_uint_fl(8, ptr, &mask_side_local, &bp_side_local, gg_idx);

    /* TNS activation flag */
    for (i = 0; i < *tns_numfilters; i++) {
        read_bit_fl(ptr, &mask_side_local, &bp_side_local, &bit);
        tns_order[i] = bit;
    }

    /* LTPF activation flag */
#ifdef NEW_SIGNALLING_SCHEME_1p25
    ltpf_tmp[1] = 0;  /* ltpf activation idx */
    ltpf_tmp[2] = 0;  /* quantized lag idx */
    if (frame_dms != LC3PLUS_FRAME_DURATION_1p25MS)
    {
        read_bit_fl(ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[0]);
    }
    else
    {    /* read one of {2, 6, 7} bits into ltp/ltpf/lag  variable  ltpf_idx[ 0 ...  2]  */
        readLtpData_fl(ptr, bfi, &mask_side_local, &bp_side_local, ltpf_tmp, rx_status, ltpfinfo_frame_cntr, mem_continuation);
    } /* !  LC3PLUS_FRAME_DURATION_1p25MS */
#else
    read_bit_fl(ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[0]);
#endif

    /* read SNS data */
# ifdef CR9_C_ADD_1p25MS_LRSNS
    ltpf_idx_2_lrsns[0] = ltpf_tmp[0];  /* raw LTP flag   input to LRSNS */
    ltpf_idx_2_lrsns[1] = ltpf_tmp[1];  /* raw LTPF flag  input to LRSNS */
#   ifdef LRSNS_PC_SIGNAL_FIX  
    bfiSNS = 0; /* Local BFI flag for Errors SNS bit area */
    readSNSData_fl(ptr, &bfiSNS, &mask_side_local, &bp_side_local, ltpf_idx_2_lrsns, sns_vq_idx, frame_dms, plc_trigger_SNS1, plc_trigger_SNS2);
    if (bfiSNS != 0 ) 
    {   /* corrupt SNSbits triggers PLC through  global PLC flag.  
          *bfi==2 and  bfiSNS == 0 maintains bfi==2 for PC    
        */ 
        *bfi = 1;
        return;
    }

    if (frame_dms == LC3PLUS_FRAME_DURATION_1p25MS)
    {   /* for 1.25ms  and previously detected bit errors -->  handle frame as a completely corrupt bad frame  */
        if (*bfi == 2) 
        {
            *bfi = 1;
            return;
        }
    }

#   else
    readSNSData_fl(ptr, bfi, &mask_side_local, &bp_side_local, ltpf_idx_2_lrsns, sns_vq_idx, frame_dms, plc_trigger_SNS1, plc_trigger_SNS2);
    if (*bfi != 0)
    {
        *bfi = 1;
        return;
    }
#   endif
# else
    readSNSData_fl(ptr, bfi, &mask_side_local, &bp_side_local, sns_vq_idx, frame_dms, &plc_trigger_SNS1, &plc_trigger_SNS2);
# endif

    /* LTPF data */
#ifdef CR9_C_ADD_1p25MS
#  ifdef NEW_SIGNALLING_SCHEME_1p25
    if ( frame_dms != LC3PLUS_FRAME_DURATION_1p25MS )
    {
        ltpf_tmp[1] = 0;
        ltpf_tmp[2] = 0;
        if (ltpf_tmp[0] == 1)
        {
            read_bit_fl(ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[1]);
            read_uint_fl(9, ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[2]);
        }
    }
#  endif
#endif

#ifndef CR9_C_ADD_1p25MS
    if (ltpf_tmp[0] == 1)
    {
        read_bit_fl(ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[1]);
        read_uint_fl(9, ptr, &mask_side_local, &bp_side_local, &ltpf_tmp[2]);
    }
    else
    {
        ltpf_tmp[1] = 0;
        ltpf_tmp[2] = 0;
    }
#endif /* CR9_C_ADD_1p25MS */

    for (i = 0; i < 3; i++) {
        ltpf_idx[i] = ltpf_tmp[i];
    }

    /* Noise factor */
    read_uint_fl(3, ptr, &mask_side_local, &bp_side_local, fac_ns_idx);

    *bp_side = bp_side_local;
    *mask_side = mask_side_local;
}
