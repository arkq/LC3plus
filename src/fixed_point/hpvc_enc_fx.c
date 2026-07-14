/******************************************************************************
*                        ETSI TS 103 634 V1.7.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "functions.h"
#include "setup_enc_lc3plus.h"

#ifdef LL_INCL_HPVC
#  ifdef CR14_A_ADD_LOSSLESS_MODE

#  ifdef __cplusplus
namespace NAMESPACE_VERSION
{
#  endif

#ifdef DEBUG
    extern int32_t frame;
    extern int32_t dbg_frame;
    extern int32_t frame_dbg;
#endif



    /*decl of helper functions */
void hpvc_get_quasi_uniform_sizes(
        Word16 Np,
        Word16 Ns,
        Word16* nLeafs /* o: vector of length Ns with leaf widthd  sizes */
    );


void hpvc_get_log_sizes(
    Word16 splitRule,
    Word16 Np,
    Word16 Ns,
    Word16* nLeafs   /* o: vector of length Ns with leaf width  sizes acording to splitRule */
);



Word16  hpvc_analyze_start_fx(Word32* truncated_data , Word16 lastnz, Word16* L1_signal, HpvcEncCfg* hpvcEncCfgPtr )
{
    Counter i,j,k;
    Word16  n_blocks, n_blocks_low, n_blocks_high;
    Word16 cnt_low,cnt_high, mode;
    Word16 cnt_zeros_L1_4[2], cnt_zeros_L1_8[2];
    Word16 nFreqBlocks[2];
    Word32 *L_ptr;
    Word32 L_acc, L_acc2;
    Word16  low_encodable, high_encodable;

    mode = -1; move16();
    Word16 coef0, coef1;

    coef0 = hpvc_adjust_startcoefs(hpvcEncCfgPtr->startCoefListNom, lastnz, LL_HPVC_N_SIGNAL,  N_SIGNAL_LOG, hpvcEncCfgPtr->startCoefList);
    coef1 = hpvcEncCfgPtr->startCoefList[1];
    ASSERT( ((coef1 - coef0) % LL_HPVC_N_SIGNAL) == 0); //only check if blocksize is correct. It is OK if coef1 = coef0, means only HF analyzed


    /*    (0 ... startcoef0-1)   |   startcoef0 ... startcoef1-1    |   startcoef1... lastnz-1 */

    IF( sub( lastnz, coef0) > 0 )
    {

        n_blocks = shr_pos(sub(lastnz, coef0), N_SIGNAL_LOG);
        n_blocks_high = shr_pos(sub(lastnz, coef1), N_SIGNAL_LOG);

        n_blocks_low = s_max(0, sub(n_blocks, n_blocks_high));

        nFreqBlocks[0] = n_blocks_low;  move16();
        nFreqBlocks[1] = n_blocks_high; move16();
        cnt_zeros_L1_4[0] = 0; move16();
        cnt_zeros_L1_4[1] = 0; move16();
        cnt_zeros_L1_8[0] = 0; move16();
        cnt_zeros_L1_8[1] = 0; move16();

        L_ptr = &(truncated_data[coef0]); /* ptr init */
        Word16 *tmp_ptr = &(L1_signal[0]);
        FOR(k = 0; k < 2; k++)
        {
            FOR(i = 0; i < nFreqBlocks[k]; i++)
            {
                L_acc = 0; move32();
                FOR (j = 0; j < (LL_HPVC_N_SIGNAL / 2); j++)
                {
                    L_acc = L_add(L_acc, L_abs(*L_ptr++));
                }
                if (L_acc == 0)
                {
                    cnt_zeros_L1_4[k] = add(cnt_zeros_L1_4[k], 1);
                }

                L_acc2 = 0; move32();
                FOR (j = 0; j < (LL_HPVC_N_SIGNAL / 2); j++)
                {
                    L_acc2 = L_add(L_acc2, L_abs(*L_ptr++));
                }

                if (L_acc2 == 0)
                {
                    cnt_zeros_L1_4[k] = add(cnt_zeros_L1_4[k], 1);
                }

                L_acc = L_add(L_acc, L_acc2);
                if (L_acc == 0)
                {
                    cnt_zeros_L1_8[k] = add(cnt_zeros_L1_8[k], 1);
                }
                L_acc = L_min(L_acc, 32767);  /* saturate to Word16 positive  domain */
                *tmp_ptr++ = extract_l(L_acc);

            }
        }

        /* check low section */
        cnt_low = 0; move16();
        FOR (i = 0; i < n_blocks_low; i++)
        {
            if (sub(L1_signal[i] , LL_HPVC_KP_MAX)  < 0  )  /* kpMax for Nsignal */
            {
                cnt_low = add(cnt_low, LL_HPVC_N_SIGNAL);
            }
        }
        low_encodable = cnt_low;

        /* laplacian boost low section */
        if (low_encodable > 0  &&  cnt_zeros_L1_4[0] > 0  /* TBD SET LIM LOW 4 */ )
        {
            cnt_low = add(cnt_low, LL_HPVC_N_SIGNAL / 2);/* TBD SET INCREASE 4  */
        }
        if (low_encodable > 0 && cnt_zeros_L1_8[0] > 0  /* TBD SET LIM LOW 8 */)
        {
            cnt_low = add(cnt_low, LL_HPVC_N_SIGNAL );/* TBD SET INCREASE 8  */
        }


        /* check high section */
        cnt_high = 0;  move16();
        FOR (i = n_blocks_low; i < n_blocks; i++)
        {
            if (sub(L1_signal[i], LL_HPVC_KP_MAX) < 0)   /* kpMax for Nsignal */
            {
                cnt_high = add(cnt_high, LL_HPVC_N_SIGNAL);
            }
        }

        high_encodable = cnt_high;

        /* laplacian boost high  section */
        if ((high_encodable >0) &&  (cnt_zeros_L1_4[1]  > 0)   /* TBD SET LIM2  4 */)
        {
            cnt_high = add(cnt_high, LL_HPVC_N_SIGNAL / 2);/* TBD SET INCREASE  */
        }
        if ((high_encodable > 0) && (cnt_zeros_L1_8[1] > 0)   /* TBD SET LIM2_8 */)
        {
            cnt_high = add(cnt_high, LL_HPVC_N_SIGNAL );/* TBD SET INCREASE  */
        }



        if (n_blocks_low < 2)
        {
            cnt_low = (2 * LL_HPVC_N_SIGNAL + 1);
        }
        if ((low_encodable > 0) && (cnt_low > (2*LL_HPVC_N_SIGNAL) ))  /*  TUNE  LOW  LIMIT  */
        {
            mode = 0;  move16();
        }
        if(n_blocks_high < 2)
        {
            cnt_high = (2 * LL_HPVC_N_SIGNAL + 1);
        }
        logic16();
        if (((low_encodable == 0) || cnt_low <= (2*LL_HPVC_N_SIGNAL)) && (high_encodable > 0) && (cnt_high > (2*LL_HPVC_N_SIGNAL) ))   /*   TUNE  HIGH  LIMIT  */
        {
            mode = 1; move16();
        }

    }

    hpvcEncCfgPtr->mode = mode; move16();

    hpvcEncCfgPtr->startCoef = coef0; move16();
    hpvcEncCfgPtr->startCoefNom = hpvcEncCfgPtr->startCoefListNom[0]; move16();  /* set LF */
    IF (mode > 0)
    {
        hpvcEncCfgPtr->startCoefNom = hpvcEncCfgPtr->startCoefListNom[1]; move16(); /* HF */
        hpvcEncCfgPtr->startCoef = coef1; move16();
    }

    return mode; /*-1, 0, 1  */
}

/*  exact calc of N_mpvq */
Word32 hpvc_leaf_sz_Nmpvq(  /* o:   N_MPVQ(N,K) ) in   Q0  */
    Word16    N,            /* i:  Dimension */
    Word16    K             /* i:  L1-norm  */
)
{
    Word16  u_shift ;
    Word32  L_tmp;
    const   UWord32 *hpMemPtr;
    Word16 Ntab;

    ASSERT(N >= 0 && K >= 0 && (HPVC_tabledKMAX[MIN(N,65)] != -1)  && (K <= HPVC_tabledKMAX[MIN(N,65)]));



    Ntab = s_min((LL_HPVC_NP_MAX / 2) + 1, N); /* clamp N==128 to 65 , saving ROM */
    ASSERT(HPVC_tabledKMAX[Ntab] >= 0);

    IF(K == 0)
    {
        L_tmp = 0;
    }
    ELSE IF(K <= HPVC_tabledKMAX[Ntab])
    {
        hpMemPtr = &(HPVC_MPVQ_offs_ptr[Ntab][K]);  move32(); /* adaptive ptr init */

        ASSERT((hpMemPtr != NULL) && "NOT ALL dimesions (Ns) are supported as leaves  in the realized  HPVC tree ");

        u_shift = 1; move16();  /* not a last A offset item in MPVQ offsets */
        if (sub(K, HPVC_tabledKMAX[Ntab]) == 0)  /* last column */
        {
            u_shift = 0; move16();/*  a last U item in MPVQ offsets */
        }

        L_tmp = L_add(1L, L_shr_pos((Word32)(*hpMemPtr++), 1));
        L_tmp = L_add(L_tmp, L_shr_pos((Word32)(*hpMemPtr), u_shift));
    }
    ELSE
    {
        L_tmp = INT32_MAX;  /* invalid leaf, return maxInt for now */
    }


    return L_tmp;
}

/*  near exact bit calculation of  S_mpvq */
/* cost about  :   16(log2_LC)  cycles  */
Word16 hpvc_leaf_rate_Smpvq(  /* o:  log2( N_MPVQ(N,K) ) in   S6Q9   */
       Word16    N,           /* i:  Dimension */
       Word16    K            /* i:  L1-norm  */
)
{
    Word16 leaf_bitrateQ9;
    Word32 L_tmp;
    Word16 Ntab;



    ASSERT(N != 0);


    Ntab = s_min((LL_HPVC_NP_MAX / 2) + 1, N);
    leaf_bitrateQ9 = INT16_MAX; move16();   /* K might be too high  */

    IF( K != 0  && (K <= HPVC_tabledKMAX[Ntab] )  )
    {
        L_tmp = hpvc_leaf_sz_Nmpvq(N, K);   /* N_mpvq from table offsets */

        L_tmp = L_add((31L << 25) + (1L << 15), BASOP_Util_Log2_LC(L_tmp));   /* input x in Q31,  output   is log2(x), 6Q25 */
        leaf_bitrateQ9 = extract_h(L_tmp);  /* extract in Q9,  rounding (1L<<15) added above  */
        ASSERT(leaf_bitrateQ9 >= 0 && leaf_bitrateQ9 <= ((31 << 9)+1)); /* stay within  30.99999 bits */
    }

    if(K == 0)
    {
        leaf_bitrateQ9 = 0; move16();   /* K might be 0 */
    }


    ASSERT(HPVC_tabledKMAX[Ntab] > 0);
    if( K > HPVC_tabledKMAX[Ntab] )
    {
        leaf_bitrateQ9 = -32768;  move16();  /* sentinel signaling  invalid leaf */
    }



    return   leaf_bitrateQ9;
}

Word16  hpvc_abs_sum(Word16 *in, Word16 len)
{
    Counter i;
    Word32 Lacc;
    Word16 l1norm;



    /*   max L1_norm is ( 128*128 )  */
    Lacc= L_mult0(1, abs_s(in[0])) ;
    FOR(i = 1; i < len; i++)
    {
        Lacc = L_mac0(Lacc, 1, abs_s(in[i]));
    }

    Lacc = L_min(Lacc, 32767);
    l1norm = extract_l(Lacc);


    return l1norm;
}

/***********************************/
/* HPVC encoder tree_parsing logic */
/**********************************/

/* split */

/* HPVC split rule parsing  common between encoder and decoder */
Word16 MPVQ_HPVC_SplitSetup( /* o: Ns */
    Word16 NpIn, Word16 Kp, Word16 *NsHdrPtr)
{
    Word16 NsSafe, NpIdx;
    const Word16* splitVecK;


    /*N_SIGNAL LOG may change   , (64) to 128   */

    NpIdx = sub(14 - N_SIGNAL_LOG, norm_s(NpIn));       /*NSignal=64, -->  (64,(1)128 */


    ASSERT(NpIn == (1 << (NpIdx + N_SIGNAL_LOG)));
    ASSERT(NpIdx >= 0);

    splitVecK = hpvc_kMax[NpIdx];  /* point to the relevant kMaxVec vector for block width NpIn */;

    NsSafe = 1; move16();   /* assume noSplit */
    IF(sub(Kp, splitVecK[1]) > 0)   /* position  “1” contains the limit for MPVQ (single idx)  vs HPVC-tree coding */
    {
        NsSafe = 3; move16();
        WHILE(sub(Kp, splitVecK[NsSafe]) > 0)
        {
            NsSafe = add(NsSafe, 1);
        }
    }
    /* this loop can also be tabled with  5 x kpMax(36) elements ~= 180 values */
    /* Ns now has the number of normal HPVC segments for an input vector of width Np */

    *NsHdrPtr = 0; move16();  /* an LS and a single leaf will be encoded, NO  Hdr  at all   */
    if(NsSafe != 1)
    {  /* 3+  segments */
       *NsHdrPtr = 1; move16();  /* assume a single non-split Hdr */
    }

     /* The Header of the flat tree might need an additional split  */
     /* Ns is in the range 9 to 16, and Np = 128 in this implementation */
     /* Ns belongs to 9,10,12,13,14,15,16   */
     /* examples: log2(Nmpvq(Ns=9..16,K=36)) > 31bits ) */
    ASSERT(NsSafe <= 16);
    if (sub(Kp, hpvc_noHdrSplitK[NsSafe]) > 0)
    {
        ASSERT(NsSafe > 8 && NpIn == 128);
        *NsHdrPtr = 3; move16();   /* one NsHdr split value is always sufficient,   9/3 => {3,3,3}, ...,   16/3  => {5, 5, 6}   */
        /* NsHdrSafe widths belongs to {3,4,5,6} */
    }

    /* *NsHdrPtr now has the number of segments for the header vector of total width Ns */
    return NsSafe ;
}


Word16  splitRuleAdapt_Ns_NsHdr( /*o: updated Ns flat leaves */
    Word16 splitRule,   /*i: rule to adapt */
    Word16 Np, /*i: for asserts only*/
    Word16 Kp, /*i: */
    Word16 NsSafe, /*i: >  1   */
    Word16 *NsHdr)  /*o: updated hdr splits */
{

    Word16 NsRedKeepQ14[4] = { 16384/*uniSafe 100%*/, 14746/*log 90%*/, 14746/* revlog 90%*/, 12288 /* UniAggressive 75%*/ };  /*Q14, 32768=keep 100% , 29491=keep 90%   26241= keep 80 %*/
    Word16 Ns;
    UNUSED(Np);  /*   */

    ASSERT(NsSafe >=3 &&  "Adaptation only possible with NsSafe>=3");


    Ns = NsSafe; move16();

    IF (splitRule > 0)
    { /* log(1), revlog(2), uni_optimistic(3) */
       Ns = extract_l(L_shr_pos(L_mult0(shl_pos(NsSafe, 10), NsRedKeepQ14[splitRule]), 10 + 14)); /*mult + trunc */
       Ns = s_min(sub(NsSafe, 1), Ns);
       Ns = s_max(3, Ns);
       if (splitRule == 3)
       {
           ASSERT(NsSafe > Ns && "Adaptation can not move down to less segments then NsSafe");
       }
    }

    ASSERT(Ns >= 3 && "Adaptation can not move down to less segments than 3 ");


    /* <Np>, <Kp> -->"NsSafe,NsHdrSafe",
       <LS>, <SplitRule>->{newNs, ->newNsHdr},   <Hdr><HdrLeaves x newNsHdr>, <FlatLeaves newNs>    */

    /*  reevaluate Hdr split number is a necessity  ,  as new lower Ns may not need a header split any longer */
    *NsHdr = 1;
    if (sub(Kp, hpvc_noHdrSplitK[Ns]) > 0)
    {
        *NsHdr = 3; move16();
    }

#ifdef DEBUG
    if (sub(Ns,1) == 0)
    {
        *NsHdr = 0; move16();
    }
#endif

    return Ns;
}


/* calculate_bits for a segment width Np in  coeffs L_X , with known {Kp, Ns and NsHdr}  */
Word32   MPVQ_HPVC_BitEst(  /*o:  bits for the long segment in Q9 */
    Word16 Np,
    Word16 Kp,
    Word16 Ns,
    Word16 NsHdr,
    Word32* L_X)
{
    /* temporary variables for widths and L1-norms */
    Counter f, i, j;


    Word32  L_bitsQ9;
    Word16  kLeafs[LL_HPVC_NS_MAX];    /*flat split Ks  NsMax is 16(LL_HPVC_NS_MAX ); */
    Word16  nLeafs[LL_HPVC_NS_MAX];
    Word16  kHdrLeafs[LL_HPVC_NSHDR_MAX];   /* NsHdrMax is 3;*/
    Word16  nHdrLeafs[LL_HPVC_NSHDR_MAX];
    Word16  x_abs[LL_HPVC_NP_MAX];   /* needed for  leaf K calculations,  Np_max is 128 */
    Word16 * ptr;
    Word16 tmpQ9;
    Word16* ptrFlatLeaves;



    logic16();
    IF( Kp == 0 || Np == 0)
    {
        BASOP_sub_sub_end();
        return 0L;  /* no bits for no pulses */
    }

    IF( sub(Ns, 1)  == 0)
    {  /* no segment splitting at all required : only the LS bit + MPVQ - idx will eventually be created */
        ASSERT(NsHdr == 0);
        kLeafs[0] = Kp;
        nLeafs[0] = Np;
    }
    ELSE
    {
        /* create the quasi-uniform  tree   segments */
        hpvc_get_quasi_uniform_sizes(Np, Ns, nLeafs);

        /* sum up the L1-norm for each of the Ns segments.   NB not  limited by blocks of 8 coefs */

        /* absolute from L_X */
        FOR(f = 0; f < Np; f++)
        {
             x_abs[f] = extract_l(L_abs(L_X[f])); /* we assume that the coeffs of segments of L_X has less than HPVC_KMAX values */
             ASSERT(L_X[f] == ((L_X[f] << 16) >> 16));
        }


        ptr = &(x_abs[0]);  /* ptr init */
        FOR(i = 0; i < Ns; i++)
        {
              kLeafs[i] = 0; move16();
              FOR(j = 0; j < nLeafs[i]; j++)
              {
                  kLeafs[i] = add(kLeafs[i], *ptr++);
              }
        }

        IF(NsHdr == 1)
        {
            nHdrLeafs[0] = Ns; move16();
            kHdrLeafs[0] = Kp; move16();
        }
        ELSE
        {    /* a Header split is required */
             ASSERT(NsHdr == 3);
             ASSERT(Ns >= 8 && Ns <= LL_HPVC_NS_MAX);

             /* find the quasi - uniform header widths for a total of NsHdr segments */

             hpvc_get_quasi_uniform_sizes(Ns, NsHdr, nHdrLeafs);

             /*  sum up the L1 - norm for each of the NsHdr segments.  */

             ptrFlatLeaves = &(kLeafs[0]);/*  sequential position in segments kLeafs[0 … Ns-1] e.g [0..15]   */
             FOR(i = 0; i < NsHdr; i++)
             {
                 kHdrLeafs[i] = 0;  move16();
                 FOR(j = 0; j < nHdrLeafs[i]; j++)  /* e.g[5, 5, 6] , for a header width 16 split in three */
                 {
                     kHdrLeafs[i] = add(kHdrLeafs[i], *ptrFlatLeaves++);
                 }
             }
        }/* % end header split was required */
    }  /* % end   initial split coding was required*/



    /*
        % we now have the segments widths and the L1 - norm desired
        %   sum up the cost in bits for leaves
        % Find the total bit rate estimate for the(Np, Kp) tree with segments { Ns, NsHdr }
    */

    L_bitsQ9 = L_add(1L << 9, 0);   /* add cost for the always present initial Leading Sign LS */

    /*   header splitting  was required, one or more header symbols are produced */

    IF( NsHdr > 1) /* the split hdr  top node */
    {
        ASSERT(NsHdr == LL_HPVC_NSHDR_MAX);
        tmpQ9 = hpvc_leaf_rate_Smpvq(NsHdr, Kp);
        L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);       /*accumulate bits to L word*/
    }

    FOR(i = 0; i < NsHdr; i++)
    {
        /* header leaves represented by  NsHdr half - pyramids */
        /* non-split Header  summed up as a single first Hdr leaf  */
        tmpQ9 = hpvc_leaf_rate_Smpvq(nHdrLeafs[i], kHdrLeafs[i]);
        L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);       /*accumulate bits to L word*/
    }


    /* sum  up the normal(nonSplit  Header) flat tree  HPVC leaf segments */
    FOR(i = 0; i < Ns; i++)
    {
        /*  regular leaf segments  represented by  Ns half - pyramids */
        tmpQ9    = hpvc_leaf_rate_Smpvq(nLeafs[i], kLeafs[i]);
        L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);
    }
    return L_bitsQ9;

}



Word32   MPVQ_HPVC_BitEstMulti(  /*o:  bits for the long segment in Q9 */
    Word16 splitRule, /* splitrule, 0 uni,1 log,2 revlog, 3 uniform_optimistic  */
    Word16 Np,
    Word16 Kp,
    Word16 Ns,            /* actual Ns transmitted */
    Word16 NsHdr,          /*actual NsHdr transmitted */

    Word32* L_X)
{
    /* temporary variables for widths and L1-norms */
    Counter f, i, j;


    Word32  L_bitsQ9;
    Word16  kLeafs[LL_HPVC_NS_MAX];    /*flat split Ks  NsMax is 16; */
    Word16  nLeafs[LL_HPVC_NS_MAX];
    Word16  kHdrLeafs[LL_HPVC_NSHDR_MAX];   /* NsHdrMax is 3;*/
    Word16  nHdrLeafs[LL_HPVC_NSHDR_MAX];
    Word16  x_abs[LL_HPVC_NP_MAX];   /* needed for  leaf K calculations,  Np_max is 128 */
    Word16 * ptr;
    Word16 tmpQ9;
    Word16* ptrFlatLeaves;
    Word16 invalidLeafFlag;


    IF(Kp == 0  )
    {
        return 0L;  /* no bits for no pulses */
    }

    invalidLeafFlag = 0; move16();
    IF(sub(Ns, 1) == 0)
    {  /* no segment splitting at all required : only the LS bit + MPVQ - idx will eventually be created */
        ASSERT(NsHdr == 0);
        kLeafs[0] = Kp;
        nLeafs[0] = Np;
        ASSERT(Kp <= HPVC_tabledKMAX[Np]);
    }
    ELSE
    {
        IF(splitRule <= 0 || splitRule==3 )
        {
        /* create the quasi-uniform  tree   segments */
             hpvc_get_quasi_uniform_sizes(Np, Ns, nLeafs);
        }
        ELSE
        {
            hpvc_get_log_sizes(splitRule,Np, Ns, nLeafs);    /* table lookup for now f==1 ..> log, f==2 revlog  */
            ASSERT(splitRule <= LL_HPVC_SPLITRULE_GLOBAL_MAX);
        }



        /* sum up the L1-norm for each of the Ns segments.
          NB not  limited to blocks of Nsignal  coefs */

        /* absolute and to Word16  from Word32  L_X */
        FOR(f = 0; f < Np; f++)
        {
                x_abs[f] = extract_l(L_abs(L_X[f])); /* we assume that the coeffs of segments of L_X has less than HPVC_KMAX values */
                ASSERT(L_X[f] == ((L_X[f] << 16) >> 16));
        }


           ptr = &(x_abs[0]);  /* ptr init to Word16 */
           FOR(i = 0; i < Ns; i++)
           {
                 kLeafs[i] = 0;           move16();
                 FOR(j = 0; j < nLeafs[i]; j++)
                 {
                     kLeafs[i] = add(kLeafs[i], *ptr++);
                 }
                 ASSERT(HPVC_tabledKMAX[nLeafs[i]] != -1);   /*all leaf n's shoud be tabled */

                 if ( kLeafs[i] > HPVC_tabledKMAX[nLeafs[i]])   /* verify that  the  tree   is still possible  for splitRule 1, 2,3   */
                 {
                     invalidLeafFlag = 1; move16();
                 }
           }
           ASSERT(ptr[-1] == x_abs[Np - 1]);

           IF(sub(NsHdr,1) == 0 || invalidLeafFlag!=0 /* quicker exit */)
           {
               nHdrLeafs[0] = Ns; move16();
               kHdrLeafs[0] = Kp; move16();
           }
           ELSE
           {    /* a Header split is required */
                ASSERT(NsHdr == LL_HPVC_NSHDR_MAX);
                ASSERT(Ns >= 8 && Ns <= LL_HPVC_NS_MAX);

                /* Hdr find the quasi  - uniform header widths for a total of  NsHdr segments */
                hpvc_get_quasi_uniform_sizes(Ns, NsHdr, nHdrLeafs);

                /*  sum up the L1 - norm for each of the NsHdr segments.  */

                ptrFlatLeaves = &(kLeafs[0]);/*  sequential position in segments kLeafs[0 … Ns-1] e.g [0..15]   */
                FOR(i = 0; i < NsHdr; i++)
                {
                    kHdrLeafs[i] = 0;  move16();
                    FOR(j = 0; j < nHdrLeafs[i]; j++)  /* e.g[5, 5, 6] , for a header width 16 split in three */
                    {
                        kHdrLeafs[i] = add(kHdrLeafs[i], *ptrFlatLeaves++);
                    }
                }
           }/*   end header split was required */
    }  /*   end   initial split coding was required*/


    L_bitsQ9 = INT32_MIN; move32();  /* negative value is  a  sentinel for an invalid tree */
    IF( invalidLeafFlag == 0)
    {


        /*
            % we now have the segments widths and the L1 - norm desired
            % time to sum up the cost in bits for leaves
            % Find the total bit rate estimate for the(Np, Kp, ) tree with segments { Ns, NsHdr }
        */

        L_bitsQ9 = L_add(1L << 9, 0);   /* add cost for the always present initial Leading Sign LS */

        /*   header splitting  was required, one or more header symbols are produced */

        IF(NsHdr > 1)  /* the split hdr  top node */
        {
            ASSERT(NsHdr == LL_HPVC_NSHDR_MAX );
            tmpQ9 = hpvc_leaf_rate_Smpvq(NsHdr, Kp);
            L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);       /* accumulate bits to L word*/
        }

        FOR(i = 0; i < NsHdr; i++)
        {
            ASSERT((HPVC_tabledKMAX[MIN(65, nHdrLeafs[i])] > 0 && kLeafs[i] <= HPVC_tabledKMAX[MIN(65, nHdrLeafs[i])]) && " Hdr leaf size too large");
            /*   header leaves represented by  NsHdr half - pyramids */
            /* non-split Header  summed up as a single first Hdr leaf  */
            tmpQ9 = hpvc_leaf_rate_Smpvq(nHdrLeafs[i], kHdrLeafs[i]);
            L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);       /*accumulate bits to L word*/
        }

        /* sum  up the normal(nonSplit  Header) flat tree  HPVC leaf segments */
        FOR(i = 0; i < Ns; i++)
        {
            ASSERT((HPVC_tabledKMAX[MIN(65, nLeafs[i])] > 0 &&  kLeafs[i] <= HPVC_tabledKMAX[MIN(65, nLeafs[i])]) && " Flat leaf size too large");

            /*  regular leaf segments  represented by  Ns half - pyramids */
            tmpQ9 = hpvc_leaf_rate_Smpvq(nLeafs[i], kLeafs[i]);
            L_bitsQ9 = L_msu0(L_bitsQ9, -1, tmpQ9);
        }
    }

    return L_bitsQ9;
}



/*  hpvc_segment  accept or reject based on estimated bit rate  */
Word16  hpvc_segment_accept_reject(     /* o: allTCX_flag   and  Updated Tx_dec  */
    Word32* L_bits_shortQ9,             /* i: tcx_bitrates per L1_signal block of coeffs  */
    Word32* L_bits_longQ9,              /* i: selected long segment bit rates has a correct bitrate, others "-1" */
    Word16* L1_signal,                   /* L1norm per L1-signal segment */
    Word16 Nqm_ana,                     /* i: Nq_mana, number of computed blocks in L1_signal */
    Word16* Tx_dec,                     /*i/o:  initial tcx_hpvc transmit decisions ,  size at least Nqm_ana */
    Word32* L_totalBits_fulltcxQ9,       /*o:   full tcx only  sum of bits  no signalling overhead  added     */
    Word32* L_totalBits_mixedQ9    /*o:    mixed coding sum of bits including signalling overhead     */
#ifdef HPVC_APPLY_TREE_LIMIT
    , Word16 nomTreeLimit
#endif
#ifdef HPVC_MAXTREE_LIMIT
    , Word16 maxTreeLimit
#endif
)
{
    Counter f, i;
    Word16 idxNpLen, nBlocksNp,idxNpLin;
    Word32 L_totalShortBitsQ9;
    Word16 f_adv; /* step in for loop*/
    Word32 L_totalBitsQ9;
    Word16 allTCX_flag;
    Word16 Kp;

    Word32 L_longBitsLocalQ9;
    Word32 L_shortBitsLocalQ9;
    Word16 NpIdxBitQ9[1 + LL_HPVC_NB_NP] ;

    /* an initial  a bit too optimistic Np cost  table */






#ifdef HPVC_APPLY_TREE_LIMIT
    Word32 L_treeBits[4 * HPVC_NOMTREE_COUNT_FB];
    Word16 treePosition[4 * HPVC_NOMTREE_COUNT_FB];
    Word16 treeK[4 * HPVC_NOMTREE_COUNT_FB];
    Word16 tree_nb;
    tree_nb = 0; move16();
#endif

    basop_memcpy(NpIdxBitQ9, NpTabBitsQ9, sizeof(Word16) * (1 + LL_HPVC_NB_NP));

    Word16 KpIdxCostQ9 = 2667;   /*512*log2(1+36)*/

    L_totalShortBitsQ9 = 0;


#ifdef  LL_HPVC_GLOBAL_FRAC

    L_totalBitsQ9 = LL_HPVC_GLOBCOST_Q9_MIXED; move16();       /*  mix of  long and  blocks */
#else
    L_totalBitsQ9 = LL_HPVC_GLOBCOST_Q9; move16();       /* incl sum of both short and long blocks */
#endif
    f_adv = 1; move16();
    allTCX_flag = 1; move16();

#ifdef  LL_HPVC_GLOBAL_FRAC
    *L_totalBits_fulltcxQ9 =  L_add(LL_HPVC_GLOBCOST_Q9_TCX,0);        /*  short blocks only  blocks */
#else
    *L_totalBits_fulltcxQ9 = 0; move32();
#endif

    FOR(f = 0; f < Nqm_ana;  /* NB f += f_adv,  i.e. adaptive f increase in tail of loop */)
    {
        ASSERT( L_bits_shortQ9[f] >= 0 );

#ifdef LL_HPVC_FORCE_NP_TCX_ENC
            Tx_dec[f] = 2;
#endif

        /* mixed sum */
        IF( Tx_dec[f] == 2)
        {
            /* a L1_signal block  with too high L1 - norm to be HPVC coded */
            L_totalShortBitsQ9 = L_mac0(L_totalShortBitsQ9, NpIdxBitQ9[0], 1);
            L_totalShortBitsQ9 = L_add(L_totalShortBitsQ9, L_bits_shortQ9[f]);    /* short block  bit-sum among the  mixed blocks */
            f_adv = 1;  move16();/* move to next Nsignal instance, typicaly 4 tuples fwd */

            L_totalBitsQ9 = L_mac0(L_totalBitsQ9, NpIdxBitQ9[0], 1);

            L_totalBitsQ9 = L_add(L_totalBitsQ9, L_bits_shortQ9[f]);  /* the mixed total sum*/

            /* full tcx coding  sum  , without signaling cost added */
           *L_totalBits_fulltcxQ9  = L_add(*L_totalBits_fulltcxQ9, L_bits_shortQ9[f]);
        }
        ELSE
        {
            idxNpLen = Tx_dec[f]; move16();                              /* 2,8,16,32,64,128 , or 2, 16,32,64,128  or 2,  32,64,128*/
            nBlocksNp = shr_pos(idxNpLen, N_SIGNAL_LOG);                /*% 0(TCX), (1, 2, 4, 8, 16 ) -> MPVQ / HPVC */
            idxNpLin = sub(15 - N_SIGNAL_LOG, norm_s(idxNpLen));         /*    0       1, 2, 3, 4, 5  */

            /* recalc Kp */
            Kp = 0; move16();
            FOR(i = 0; i < nBlocksNp; i++)
            {
                Kp = add(Kp, L1_signal[f + i]);
                ASSERT((f + i) < Nqm_ana);
            }
#ifdef  LL_HPVC_KP_FRAC
            KpIdxCostQ9 =  hpvc_KpTabBitsQ9[Kp];
#endif
            IF(nBlocksNp > 0) /* HPVC segment analysis  for Np=8,16,32,64,128 */
            {
#ifdef DEBUG
                if (L_bits_longQ9[f] < 0) {
                    printf("\n ASSERT will fail: f=%d, Tx_dec[f]=%d, nBlocksNp=%d, L_bits_longQ9[f]=%ld\n",
                           (int)f, (int)Tx_dec[f], (int)nBlocksNp, (long)L_bits_longQ9[f]);
                }
#endif
                ASSERT(L_bits_longQ9[f] >= 0);
                ASSERT( 1 << (idxNpLin-1) == nBlocksNp);
                L_longBitsLocalQ9 = L_add(KpIdxCostQ9, NpIdxBitQ9[idxNpLin]);  /*L1-norm sig cost+ width in NP signal cost */
                L_longBitsLocalQ9 = L_add(L_longBitsLocalQ9, L_bits_longQ9[f]);



                /* sum up cost of sub-alternative TCX bits over the blocks over total length Np , idxNp segments */
                L_shortBitsLocalQ9 = 0;  move16();
                FOR(i = 0; i < nBlocksNp; i++) {
                    L_shortBitsLocalQ9 = L_add(L_shortBitsLocalQ9, L_bits_shortQ9[f + i]);
                }
                *L_totalBits_fulltcxQ9 = L_add(*L_totalBits_fulltcxQ9, L_shortBitsLocalQ9); /*no signaling cost*/
                /* store before signalling cost is added to TCX*/
#ifdef LL_HPVC_SEGMENT_FILE_LOG
                if (fp_segDbg) fprintf(fp_segDbg, "%d,%d,%d,%ld,%ld\n", (int)f, (int)Nqm_ana, (int)nBlocksNp, (long)L_shortBitsLocalQ9, (long)L_longBitsLocalQ9);
#endif

                L_shortBitsLocalQ9 = L_mac0 (L_shortBitsLocalQ9, NpIdxBitQ9[0], nBlocksNp); /* add signaling cost */



                f_adv = nBlocksNp;   /* 1,2,4,8,16 */

                /* possibly revert to TCX, short block coding  ,  no additional testing of shorter Np blocks here !! */

                /* compare with signaling cost and a bit_gain offset included */
                Word16  logShiftNpVal = add(N_SIGNAL_LOG - 1, idxNpLin);
                Word32 L_tmpOffset = L_shl_pos(LL_HPVC_REQ_BIT_GAIN_PER_CEOFF_Q9, logShiftNpVal); /* make bit offset per Np number of coeffs */
                Word32 L_tmp_longPlusOffset = L_add(L_longBitsLocalQ9, L_tmpOffset );

                IF (L_sub(L_tmp_longPlusOffset, L_shortBitsLocalQ9) >= 0)
                {
                    /* switch to TCX-tuple coding for this NP section */
                    assert(f_adv == nBlocksNp); /* rejected Np length in blocks  */
                    FOR (i = 0; i < nBlocksNp; i++)
                    {
                        Tx_dec[f + i] = 2;  move16();    /* block(s) of short TCX tuples signaled  */

                    }
                    L_totalBitsQ9 = L_add(L_totalBitsQ9, L_shortBitsLocalQ9);
                }
                ELSE
                {
                    L_totalBitsQ9 = L_add(L_totalBitsQ9, L_longBitsLocalQ9);     /*  assume Np coded as hpvc  */
#ifdef HPVC_APPLY_TREE_LIMIT
                    L_treeBits[tree_nb]   = L_sub(L_shortBitsLocalQ9, L_longBitsLocalQ9);
                    treePosition[tree_nb] = f;
                    treeK[tree_nb] = Kp;
                    tree_nb = add(tree_nb, 1);
#endif
                }
            }

            if ( sub(Tx_dec[f],2) != 0)
            {
                allTCX_flag = 0; move16(); /* no possibility to revert the global signaling */
            }

        }  /* % end of ‘for’ f */
        f = add(f, f_adv);  /*  the adaptive step fwd in the FOR  loop */
    }

#ifdef HPVC_APPLY_TREE_LIMIT
    {
        Word16 n_diff = sub(tree_nb, nomTreeLimit);
#ifdef HPVC_MAXTREE_LIMIT
        Word16 n_diff_max = sub(tree_nb, maxTreeLimit);
        ASSERT(n_diff >= n_diff_max);
#endif
        logic16();
        IF((tree_nb > 0) && (n_diff > 0))
        {
            FOR(i = 0; i < n_diff; i++)
            {
                Counter j;
                Word16 idx = 0;
                Word32 L_lowest_bit_gain = L_add(L_treeBits[0], 0);
                FOR(j = 1; j < tree_nb; j++)
                {
                    if (L_sub(L_treeBits[j], L_lowest_bit_gain) < 0)
                    {
                        idx = j;  move16();
                    }
                    L_lowest_bit_gain = L_add(L_treeBits[idx], 0);
                }

                /* remove lowest bitgain tree */
                f = treePosition[idx];
                nBlocksNp = shr_pos(Tx_dec[f], N_SIGNAL_LOG);
                Word16 kLim[3] = { -1, 16, 20 };

#ifdef HPVC_MAXTREE_LIMIT
                IF(n_diff_max > 0 || sub(treeK[idx], kLim[nBlocksNp]) >= 0)   /* 20 for 128, 16 for 64 */
#else
                IF(sub(treeK[idx], kLim[nBlocksNp]) >= 0)   /* 20 for 128, 16 for 64 */
#endif
                {
#ifdef HPVC_MAXTREE_LIMIT
                    n_diff_max = sub(n_diff_max, 1);
#endif
                    Counter n;
                    FOR(n = 0; n < nBlocksNp; n++)
                    {
                        Tx_dec[f + n] = 2;  move16();    /* block(s) of short TCX tuples signaled */
                    }
                    L_totalBitsQ9 = L_add(L_totalBitsQ9, L_lowest_bit_gain);  /* add loss */
                }
                L_treeBits[idx] = L_add(INT32_MAX, 0);  /* sentinel to be able to select a new low value, next idx */
            }
        }
    }
#endif

    *L_totalBits_mixedQ9 = L_add(L_totalBitsQ9, 0);

    /*finally revert to full==all TCX if added total TCX/HPVC signalling cost added too much signaling overhead */

#ifdef LL_HPVC_REQ_BIT_GAIN_GLOBAL_Q9
    Word32 L_mix_bit_gain = L_sub(*L_totalBits_fulltcxQ9, *L_totalBits_mixedQ9);
    if (L_sub(L_mix_bit_gain, LL_HPVC_REQ_BIT_GAIN_GLOBAL_Q9) < 0)
#else
    if (L_sub(*L_totalBits_mixedQ9, *L_totalBits_fulltcxQ9) >= 0 ) /* 0 bit margin */
#endif
    {
        allTCX_flag = 1;   move16();   /* signal to skip the mixed tcx/hpvc coding outside this rate loop */

    }


    return  allTCX_flag;  /*  [allTCX_Flag, Tx_dec, L_totalBits_mixedQ9 ] */
}



/*
   Low - complexity segmentation: accumulate only fwd
*/
void  hpvc_segment_tcx_or_hpvc(     /* o:  Tx_dec, selected Np from block f and forward */
    Word16* L1_signal,               /* i: L1signal   precomputed L1 norms in blocks of Nsignal==8 */
    Word16 Nqm_ana,                 /* i: Nqmana, number of computed blocks in L1_signal */
    Word16* Tx_dec                   /*o:   initial tcx_hpvc transmit decisions ,  size at least Nqm_ana */
)
{
    Counter f, i;
    Word16 Np, h, endp_p1;


    f = 0;   move16();    /* current block index */
    Np = -1;  move16();   /*  */

    WHILE( sub(f, Nqm_ana) < 0)
    {
        IF( sub( L1_signal[f], KmaxHPVC[0] ) > 0 )
        {
            Tx_dec[f] = 2;  move16(); /* to high dynamics Set decision to TCX coding */
            /* increment for next f,  go to next block - using Np = 8 % the Nsignal stride */
            endp_p1 = add(f, 1);
        }
        ELSE
        {    /*accumulate and test block lengths versus KmaxHPVC[0...4] */
             Np = hpvc_find_valid_Np(f, L1_signal, Nqm_ana); /*  analyze all HPVC blocks alternatives */
                /* 1xNp, 2xNp, 4Np..16xNp */
             Tx_dec[f] = Np;  /* signal to transmit hpvc codeword(s) representing Np coeffs.*/
             h = sub(14 - N_SIGNAL_LOG, norm_s(Np));    /* h = log2(Np) - log2(Nsignal) */

             endp_p1 = add(f, shl_pos(1, h));

             /* safety mark-up the merged signaling decision*/
             FOR(i = add(f,1); i < endp_p1; i++)
             {           /* store a negated Np value as sentinel for an already started Np section*/
                 Tx_dec[i] = negate(Np);
                 ASSERT(i < Nqm_ana);
             }

        }  /*end else */
        f = endp_p1; /* update  f for next starting point*/
    }  /* end while */


}  /* end func */

Word16  hpvc_find_valid_Np( /* o:  Np, selected Np from block f and forward */
    Word16 f,               /* i: f         current signaling position */
    Word16* L1_signal,      /* i: L1signal   precomputed L1 norms in blocks of Nsignal==8 ( or 16) */
    Word16 Nqm_ana          /* i: Nqm_ana, number  of computed blocks in L1_signals */
)
{
    Counter i;
    Word16* L1_8 = L1_signal; /* % for notation convenience as we now use Nsignal == 8 */
    Word16 Np, h;
    Word32 L_acc;  /* we assume  L1_signal can have  values higher than hpvc max L1, but less than 32767 ....  */
    Word16 end_block_plus1_nb;


    /* set base case for  block “f” */
    ASSERT(L1_8[f] <= KmaxHPVC[0]);  /*  upper limit with split */
    Np = LL_HPVC_N_SIGNAL;	  move16();       /*  length in coeffs of a of a single block  */
    h = 0;                    move16();       /*  note h = log2(Np) - 3; */
    L_acc = L_deposit_l(L1_8[f]) ;            /*  accumulate the single block  at f  */

    /* analyze:  2xNp, 4xNp, 8xNp, 16xNp */
    /* only re-accumulate future blocks */



//#define FOR_LOOP_VERSION
#ifdef FOR_LOOP_VERSION
    Word16 next_f ;  /* next block to accumulate */
    Word16 continue_flag;

    next_f = add(f,1);  /* next block to accumulate */
    continue_flag = Nqm_ana; move16();
    ASSERT(continue_flag > 0);

    FOR(h = 1;  h <= TOPH; h++)          /* [ 2^(N_SIGNAL_LOG + "1")==16 ...  2^(N_SIGNAL_LOG + "4")==128 ] */
    {
        end_block_plus1_nb = add(f, shl_pos(1, h));   /*  in block number range [0...(Nmq_ana-1)]  , f+2^h */
        if(sub(end_block_plus1_nb, Nqm_ana) > 0)   /* (Nmq_ana-1) is last valid block number */
        {
            continue_flag = 0; move16();   /* Np hypothesis goes beyond end , stay with previous found Np */
        }

        /*
           L1fullsum  = sum( L1_8, f ,            f + 2^h - 1 ); % full reaccumulation
           L1acc_sum += sum( L1_8, f + 2^(h - 1), f + 2^h - 1 );  % only new blocks
        */
        /* init       h=0  L1_acc = acc[f ...  f]           0=floor((2^h -1)) */
        /* first time h=1  L1_acc + acc[f+1, f+1]    1=(2^h -1) */
        /* 2nd time   h=2  L1_acc + acc[f+2, f+3]    3=(2^h -1) */
        /* 3rd time   h=3  L1_acc + acc[f+4, f+7]    7=(2^h -1),  4= 2^(h-1)*/
        /* 4th time   h=4  L1_acc + acc[f+8, f+15]   15(2^h -1) */

        FOR(i = next_f; i < MIN(end_block_plus1_nb, continue_flag);  i++)
        {
            L_acc = L_mac0(L_acc, L1_8[i], 1);
        }
        next_f = end_block_plus1_nb;  move16(); /* for next loop */

        logic16();
        /*  test for end,  or  upper K-limit with  splitting allowed */
        IF( ( continue_flag == 0)  || (L_sub(L_acc, L_deposit_l(KmaxHPVC[h]) ) > 0)  )
        {
            BREAK;  /* % hypothesis breaks, stay with previous found Np */
        }

        Np = shl_pos(1, add(h, 3));   /* accept this longer accumulation length,   FOR loop may proceed */

        ASSERT(Np ==  LL_HPVC_N_SIGNAL || Np == 16 || Np == 32 || Np == 64 || Np == 128);
        ASSERT(KmaxHPVC[h - 1] <= KmaxHPVC[h]);
    }  /* % FOR */
    /* a valid NP is available */

#else


    Word16 Klim_while_exit;
    Word16 next_f;
    Word16 end_lim;

    Klim_while_exit = 0;

    next_f = add(f,1);
    logic16(); logic16();
    WHILE(Klim_while_exit == 0 && (sub(next_f, Nqm_ana) < 0 ) && sub(h, TOPH)<= 0 )  /*  test vs  [ Klim, end_point, h ]   limits */
    {
        end_block_plus1_nb = add(f, shl_pos(1, h));

        end_lim = s_min(end_block_plus1_nb, Nqm_ana);
        FOR(i = next_f; i < end_lim ;  i++)
        {
            L_acc = L_mac0(L_acc, L1_8[i], 1);
        }
        next_f = end_block_plus1_nb;  move16(); /* for next loop *//*% next  blocks to test */


        logic16();
        if ( (L_sub(L_acc, L_deposit_l(KmaxHPVC[h])) > 0 )  ||  (sub(end_block_plus1_nb, Nqm_ana)>0) )
        {
            Klim_while_exit = 1; move16();  /* break out of the loop */
        }

        Np = shl_pos(1, add(N_SIGNAL_LOG, sub(h, Klim_while_exit)) );  /*  actual “correct ” Np result so far result,  while_exit -> 1 step back */


        ASSERT( Np >= (1<<N_SIGNAL_LOG) && Np <= (1 << (N_SIGNAL_LOG+TOPH)));


        h = add(h,1);            /*   update Np length hypothesis for next iter */

        logic16(); logic16();  /* while condition logic recalc cost in next iter  */
    }
#endif /* LOOP version */


    return Np;   /* f needs to be updated outside the function, according to Np result */
}


Word16  hpvc_adjust_startcoefs(      /* o:  first adjusted startcoef */
    Word16 *startCoefListNom,       /* i: startCoefNom0 and startCoefNom1  */
    Word16 lastnz,
    Word16 Nsignal,
    Word16 log2_Nsignal,
    Word16 *startCoefList)      /* o:  lastnz adjusted startcoef state 2 variables */
{

    /* align the nominal starting coef(s) by blocks to the determined lastnz  value */
    /* shift down by Nsignal/2 coefs or up by Nsignal/2 coefs of  first   possible */
    /* NB startCoef0 and  startCoef1  are constants pre-separated by integer n  ( x Np) blocks , n > 0 */

    Word16 startCoef, NsignalBy2, diff_nom, coef_nom_delta;
    Word16 blocks_nom, blocks_delta;

    startCoef = -1;   /* dbg invalid sentinel , eventually we always set positive values */


    coef_nom_delta = sub(startCoefListNom[1], startCoefListNom[0]);
    blocks_delta = shr_pos(coef_nom_delta, log2_Nsignal);
    ASSERT(coef_nom_delta == (blocks_delta) << log2_Nsignal);

    NsignalBy2 = shr_pos(Nsignal, 1);

    diff_nom = sub(lastnz, startCoefListNom[0]);
    blocks_nom = shr_pos(diff_nom, log2_Nsignal);
    Word16 rem_diff = sub(diff_nom, shl_pos(blocks_nom, log2_Nsignal)); /* remaining coef diff  0...63 */

    logic16();


    IF(blocks_nom > 0 && rem_diff == 0)
    {
        /* nominal coef0 value is an exact hit  */
        startCoef = startCoefListNom[0]; move16();
    }
    ELSE   /* Adjustment needed to synch with teh tail  lastnz  */ /* subtract or add 32 coeffs */
    {
        IF(sub(rem_diff, NsignalBy2) > 0)  /* assume increase of LF-region */
        {
            startCoef = sub(startCoefListNom[0], sub(Nsignal,rem_diff));
        }
        ELSE      /* decrease i.e. move LF limit up  */
        {
            ASSERT(rem_diff <= NsignalBy2);
            startCoef = add(startCoefListNom[0], s_max(0, rem_diff));
        }




    }

    IF( sub(add(startCoef, Nsignal), lastnz) > 0)
    {   /* not  a single block   fits , within shifts  */
        /* i.e. shift of Nsignal/2 was not enough, */
        startCoef = add(lastnz, Nsignal); move16();/* set start   above lastnz   */

        startCoefList[0] = startCoef;  move16();
        startCoefList[1] = add(startCoef, coef_nom_delta); /* initially keep the original delta  */
    }
    ELSE
    {

        /* Set up final list and separation */
        startCoefList[0] = startCoef;
        startCoefList[1] = add(startCoef, coef_nom_delta); /* initially keep the original delta  */


        /*  finally adjust the  coef0 to coef1   delta distance
           if high band is not encodable due to lastnz being too low
               and
           if the lastnz value does allow splitting the full band  into at least two segments */
        Word16 tmp_hf_tail = sub(lastnz, Nsignal);
        IF(sub(startCoefList[1], tmp_hf_tail) > 0)
        {    /* no full block fits in the so far set up  high band */

             if (coef_nom_delta >= shl_pos(Nsignal,1)) /* if 2 bands or more */
             {
                startCoefList[1] = tmp_hf_tail; /* allow a chance  for using the high region only  */

            }
        }
    }

    ASSERT((startCoef & 0x1) == 0); /* must be even */
    return startCoef;
}


 /* function for getting the estimated bitrate of HPVC "long" segments */
Word16  hpvc_collect_hpvc_bitrates( /* o:  number of analyzed  HPVC trees  */
    Word16 *L1_signal,          /* i: f       vector of  L1-norms per L1signal block of length 8   */
    Word32 *Xqm,           /* i: signal to be encoded , need to get individual leaves K values  */
    Word16* Tx_dec,        /* i: initial decision  [2,8,16,32,64,128]  */
    Word16* Tx_splitRule,  /* o:   can also reduce  final Ns and final  NsHdr                 */
    Word16  Nqm_ana,       /* i: Nqm_ana,  number  of computed blocks in L1_signals  Tx_dec*/
    Word32 *L_bitsLongQ9   /* o: vector of individual tree bit rates */
)
{
    Counter f, f_adv,i;
    Word16 Np, Kp, nBlocksNp, hpvcTreeCount;
    Word16 NsSafe, NsHdrSafe;
    Word16 Ns, NsHdr;
    Word32 *L_Xptr;

    UNUSED(NsSafe);
    UNUSED(NsHdrSafe);


    hpvcTreeCount = 0;  move16();
    f_adv = 1;    move16();
    for (f = 0; f < Nqm_ana; f += f_adv)
    {
        L_bitsLongQ9[f] = -1; move32();
        Tx_splitRule[f] = -1;   move16(); /* assume unencoded uniform split*/

        f_adv = 1;   move16();

        Tx_splitRule[f] = -2; /* inactive split  for tcx */
        ASSERT(Tx_dec[f] > 0);

        IF(Tx_dec[f] != 2)
        {
            hpvcTreeCount = add(hpvcTreeCount, 1);

            {
                Word16 tmp = 0;
                for (i = NP_LOG_MIN; i <= NP_LOG_MAX; i++)
                {
                    tmp |= (Tx_dec[f] == (1<<i));
                }
                ASSERT(tmp != 0);
            }


            /* an assumed  long block coding segment with MPVQ / HPVC */
            Np = Tx_dec[f];
            nBlocksNp = shr_pos(Np, NP_LOG_MIN);
            Kp = 0; move16();

            FOR(i = 0; i < nBlocksNp; i++)
            {
                Kp = add(Kp, L1_signal[f + i]);
                ASSERT((f + i) < Nqm_ana);
            }

            ASSERT(Kp <= LL_HPVC_KP_MAX);

            L_Xptr = &(Xqm[f << N_SIGNAL_LOG ]);  /* Ptr init   in C to vector Xqm[f * 8 …(f * 8 + Np - 1)] */
            Tx_splitRule[f] = -1; move16();

            NsSafe = MPVQ_HPVC_SplitSetup(Np, Kp, &NsHdrSafe);   /* the common, Ns, NsHdr setup aligned with decoder */
            Ns=NsSafe;
            NsHdr = NsHdrSafe;


            L_bitsLongQ9[f] = MPVQ_HPVC_BitEst(Np, Kp, NsSafe, NsHdrSafe, L_Xptr); /* Tree cost for splitRule 0 (uni safe),  rules 1...3 later  */

            Word32 L_bitsQ9[1+LL_HPVC_SPLITRULE_GLOBAL_MAX];

            Word16 splitRule;
            Word32 L_BitCostsQ9[1+LL_HPVC_SPLITRULE_GLOBAL_MAX];

            splitRule = -1;

            IF( NsSafe > 1)
            {

                Word32 L_tmp;

                Word16 tmp = sub(15 - 4, norm_s(Np));  /* linear  idx  from Np 8., 16. ...128 */
                Word16 splitRuleMax = splitRuleMaxPerNp[tmp];
                Word16 splitRuleNsMin = splitRuleNsMinPerNp[tmp];

                if (NsSafe == LL_HPVC_SPLITRULE_GLOBAL_NSMIN)
                {
                    splitRuleMax =  s_min(2, splitRuleMax);  /* Aggressive uniform split not making sense, only {uniSafe, log, revlog}  allowed  */
                }
                ASSERT(splitRuleMax >= -1 && splitRuleMax <= LL_HPVC_SPLITRULE_GLOBAL_MAX);

                /* set up actual Np, KP dependent splitRule cost table */
                L_BitCostsQ9[0] = 0;    /*  move32  */
                L_tmp = 0;
                if (splitRuleMax == 3 )
                {
                    L_tmp = L_add(L_tmp, 2L << 9); /* 2 bits i Q9*/
                }
                if (splitRuleMax == 2)
                {
                    L_tmp = L_add(L_tmp, 2L << 9); /* 2 initial bits  Q9*/
                }
                if (splitRuleMax == 1)
                {
                    L_tmp = L_add(L_tmp, 1L << 9); /* 1 bit  Q9*/
                }
                FOR(i = 0; i <= splitRuleMax;i++)
                {
                    L_BitCostsQ9[i] = L_tmp;
                }
                if (splitRuleMax == 2)
                {
                    L_BitCostsQ9[0] = L_sub(L_BitCostsQ9[0], 1 << 9);  /* adjust   bits for safe uniform  only 1 bit   Q9*/
                }


                L_bitsQ9[0] = L_add(L_BitCostsQ9[0], L_bitsLongQ9[f]);   /*  add splitRule Tx   cost */
                L_bitsLongQ9[f] = L_bitsQ9[0];

                Ns = NsSafe;
                NsHdr = NsHdrSafe;

                Word32 L_minBits = L_bitsQ9[0];
                ASSERT(L_minBits >= 0);

                IF(splitRuleMax > 0 && NsSafe >= splitRuleNsMin)
                {
                    splitRule =  0;    /* safe uniform split is always available  */

                    Word16 sR;
                    FOR(sR = 1; sR <= splitRuleMax; sR++)    /* 0=UniSafe,   aggressive {1=log, 2=revLog, 3=Uni }   */
                    {

                        Ns = splitRuleAdapt_Ns_NsHdr(sR, Np, Kp, NsSafe, &NsHdr); /* more aggressive splits than safe is allowed  */
                        /* NB!!  NsHdr may be updated when NS is reduced */
                        L_bitsQ9[sR] = MPVQ_HPVC_BitEstMulti(sR, Np, Kp, Ns, NsHdr, L_Xptr);

                        IF(L_bitsQ9[sR] > 0) /* negative bits  means a  pyramid is too large to be indexed */
                        {
                            Word32 L_tmp = L_add(L_BitCostsQ9[sR], L_bitsQ9[sR]);       /* 1 or 2, or X   bits for signalling */
                            L_bitsQ9[sR] = L_tmp;

                            if (L_sub(L_tmp, L_minBits) < 0 )
                            {
                                splitRule = sR;     move16();
                            }
                            L_minBits = L_min(L_bitsQ9[splitRule], L_minBits);
                            ASSERT(L_minBits>=0);
                        }

                    }

                    L_bitsLongQ9[f] = L_bitsQ9[splitRule]; move32();
                    Tx_splitRule[f] = splitRule;

                }/* IF RULE  */
            }/*NSSafe>1*/






            f_adv = shr_pos(Np, N_SIGNAL_LOG);
        }
    } /* for f */


    return hpvcTreeCount;
}




void hpvc_get_quasi_uniform_sizes(
    Word16 Np,
    Word16 Ns,
    Word16* nLeafs /* o: vector of length Ns with leaf widthd  sizes */
)
{
    Counter i;
    Word16 nMin;
    Word16 *ptr;
    Word32 L_NpRem;

    /* create the quasi-uniform  tree   segments */

    Word16 factorQ16 = (Word16)((2.0*32768.0) / ((double)Ns));
    /* implement as slighly inexact fractional division using tabled factors  1/3 ... 1/16 , and a post sanity check */
    nMin = L_shr_pos(L_mult0(add(Np, add(Np, 1)), factorQ16), 16 + 1);  /*+1 for rounding  */
#ifdef DEBUG
    if ( ((nMin+1)*Ns)<= Np )
    {
        nMin =  nMin+1;
        ASSERT(((nMin - 1)*Ns < Np) && ((nMin)*Ns <= Np) && "inexact division sanity check failed");
        ASSERT(0 && "inexact division Np/Ns occurs ");  /* dbg do we ever end up here ? */
    }
#endif


    FOR(i = 0; i < Ns; i++)
    {
        nLeafs[i] = nMin;
    }
    L_NpRem = L_msu0(L_deposit_l(Np), Ns, nMin);


    ptr = &(nLeafs[Ns - 1]);

    FOR(i = 0; i < L_NpRem; i++)
    {   /* fill up dimension from the top  */
        *ptr = add(*ptr, 1);
        ptr--;
    }
    /* nLeafs updated */
}

 Word16  hpvc_find_first_nz(Word16* ptr)    /* vectors should not be all zero */
{
    Counter i;

    Word16  f_nz = -1;

    FOR(i = 0; i < LL_HPVC_NP_MAX ; i++)
    {
        IF(*ptr++ != 0)
        {
            f_nz = i;  move16();
            BREAK;
        }

    }
    ASSERT(f_nz != -1 );

    return f_nz;
}



static void  createHdrL1part(
    Word16 *vec,
    Word16 Ns,
    Word16* nSub,
    Word16 *hdrVec)
{
    Counter i, j;
    Word16* ptr;
    Word16 abs_sum;

    /* pre-process Hdr leafs */
    /* move hdrLeaf L1-norms  a  Hdr  */
    ptr = vec;
    FOR(i = 0; i < Ns; i++)
    {
        abs_sum = 0;
        FOR(j = 0; j < nSub[i]; j++)
        {
            abs_sum = add(abs_sum, abs_s(*ptr++));
        }
        hdrVec[i] = abs_sum; move16();
    }
}

static void setHdrLSmodifyLeaf2HP(Word16 NsHdr, Word16* nSub, Word16* hdrVals, Word16* xm)
{
    /* set signs in header and  modify Leafs to a half pyramid */
    Counter i;
    Word16 f_nz, acc_pos, tmp_pos;

    acc_pos = 0; move16();
    FOR(i = 0; i < NsHdr; i++)
    {
        IF( hdrVals[i] != 0 )
        {
            f_nz = hpvc_find_first_nz(&(xm[acc_pos]));
            ASSERT(f_nz < nSub[i]);
            tmp_pos = add(acc_pos, f_nz);
            if (xm[tmp_pos] < 0)
            {
                hdrVals[i] = negate(hdrVals[i]);  /* move leading sign to header */
            }
            xm[tmp_pos] = abs_s(xm[tmp_pos]); /* make leaf  into a half pyramid */
        }
        acc_pos = add(acc_pos, nSub[i]);
    }
}




void hpvc_get_log_sizes(
    Word16 splitRule,  /* 1,2 allowed */
    Word16 Np,
    Word16 Ns,
    Word16* nLeafs   /* o: vector of length Ns with leaf width  sizes acording to splitRule */
)
{
    Counter i;
    const Word16 *(*ptr1) = NULL;
    const Word16 *ptr;

    /* create deterministic log tree   segments */
    SWITCH(Np)
    {
       case 128:
           ptr1 = logSzNp128;    BREAK;
       case 64: ptr1 = logSzNp64;       BREAK;
       case 32: ptr1 = logSzNp32;       BREAK;
       case 16: ptr1 = logSzNp16;       BREAK;
       case 8: ptr1 =  logSzNp8;        BREAK;
       default:        ASSERT(0);    BREAK;
    }
    ASSERT( ptr1[Ns] != NULL);

    ptr = ptr1[Ns];  /* move ptr to target sz vector */
    ASSERT( ptr[0]>0  &&  ptr[Ns - 1] > 0);

    IF( sub(splitRule, 2) == 0)
    {
        FOR( i = 0; i < Ns; i++ )
        {
            nLeafs[i] = ptr[(Ns-1) - i];
        }
    }
    ELSE
    {
        ASSERT(splitRule == 1);
        basop_memcpy(nLeafs, ptr, sizeof(Word16)*Ns);
    }
}



Word16  hpvc_enumerate_trees(     /* o: n_Trees */
    Word32* Xm_ptr,               /* i: signal to be encoded , need to obtain individual split leaves K values for Np>8  */
    Word16* Tx_dec,               /* i: initial decision  [2,8,16,32,64,128]  */
    Word16* Tx_splitRule,         /* i:  decision from rate selection step */
    Word16* L1_signal,            /* i:  L1-norms */
    Word16 Nqm_ana,               /* i: Nqm_ana,  number  of computed blocks in L1_signals  Tx_dec*/
    HpvcTreeEnumCfg* hpvc_tree)   /* o: vector of tree bit cfgs  and idxes  [Nqmana][1+4+16] Word32 idxs or series ?   */
{
    Counter f, i;
    Word16 nTrees, nCoeffs;
    Word16 f_adv;
    Word16   nBlocksNp;
    HpvcTreeEnumCfg *tree;
    Word16 xm[LL_HPVC_NP_MAX];
    Word16 hdrVals[LL_HPVC_NS_MAX];
    Word16 nFlatLeafs[LL_HPVC_NS_MAX];    /*  sizes of leafs */
    Word16 nHdrLeafs[LL_HPVC_NSHDR_MAX];      /*  sizes of split hdr leafs */
    Word16 splitHdrVals[LL_HPVC_NSHDR_MAX];
    UWord32 h_mem[1 + LL_HPVC_KP_MAX + 1];  /* MPVQ  A and U offsets created during encoding  of a leaf*/

    UWord32 UL_idx, UL_LS, UL_sz;
    Word16  tmpK;
    Word16 flatK;
    Word16 acc_pos;





    nTrees = 0;
    nCoeffs = 0;
    f_adv = 1;

    tree = &(hpvc_tree[0]);

    for (f = 0; f < Nqm_ana; f += f_adv)
    {

        f_adv = 1; move16();

        IF(Tx_dec[f] != 2)
        {
            tree->Np = Tx_dec[f];  move16();
            nBlocksNp = shr_pos(tree->Np, N_SIGNAL_LOG);
            f_adv = nBlocksNp; move16();

            ASSERT((Tx_dec[f] > 0) && (Tx_dec[f] == (nBlocksNp << N_SIGNAL_LOG )));

            tree->Kp = 0; move16();
            FOR(i = 0; i < nBlocksNp; i++)
            {
                tree->Kp = add(tree->Kp, L1_signal[f + i]);
                ASSERT((f + i) < Nqm_ana);
            }
            tree->NsSafe = MPVQ_HPVC_SplitSetup(tree->Np, tree->Kp, &(tree->NsHdrSafe));


            tree->Xqm_ptr = &(Xm_ptr[f << N_SIGNAL_LOG ]);/* Ptr init   in C to vector Xqm[f * 8 …(f * 8 + Np - 1)] */
            tree->start_coeff_nb = f << N_SIGNAL_LOG ;



            /* move the vector to a Word16 temporary  vector "xm"  */
            FOR(i = 0; i < tree->Np; i++)
            {
                xm[i] = extract_l(tree->Xqm_ptr[i]);
                tree->xDbg[i] = xm[i];
                tree->xDeEnumDbg[i] = -1;
            }

            tree->splitRule = Tx_splitRule[f];  move16();  /* rule found in a loop by estimated bitrate comparison */



            tree->Ns = tree->NsSafe;
            tree->NsHdr = tree->NsHdrSafe;  // Ns==1 --> 0  , Ns=3...16 , -->(1 or 3  based on Kp)
            IF(tree->splitRule > 0)
            {   /* splitRule bit(s) can change the applied Ns, applied NsHdr */

                tree->Ns = splitRuleAdapt_Ns_NsHdr(tree->splitRule, tree->Np, tree->Kp, tree->NsSafe, &(tree->NsHdr));

            }

            /*below this line use Ns and NsHdr */

            IF( tree->Kp != 0)
            { /* one or more idx's will actually be transmitted */


                IF(tree->NsHdr == 0)
                {   /* enumerate the single leaf  */
                    ASSERT(tree->Ns == 1);
                    //printf("\n tree with single leaf Np=%d,  Kp==%d ", tree->Np, tree->Kp);
                    vec2mind_fx(tree->Np, tree->Kp, xm, &UL_LS, &UL_idx, &UL_sz, h_mem);  /* actual enumeration */
                    ASSERT(UL_idx <= INT32_MAX && UL_sz <= INT32_MAX && UL_idx < UL_sz);

                    tree->flatLeafIdx[0] = L_add((Word32)UL_idx, 0);
                    tree->flatLeafSz[0] = L_add((Word32)UL_sz, 0);
                    tree->LS = extract_l((Word32)UL_LS);

                }
                ELSE IF(tree->NsHdr > 0)
                {   /* enumerate hdr + flat split */
                    ASSERT(tree->NsHdr == 1 || tree->NsHdr == 3);

                    logic16();
                    IF(tree->splitRule <= 0 || sub(tree->splitRule,3) == 0)
                    {
                        hpvc_get_quasi_uniform_sizes(tree->Np, tree->Ns, nFlatLeafs);
                    }
                    ELSE
                    {
                        hpvc_get_log_sizes(tree->splitRule, tree->Np, tree->Ns, nFlatLeafs);
                    }

                    createHdrL1part(xm, tree->Ns, nFlatLeafs, hdrVals);

                    setHdrLSmodifyLeaf2HP(tree->Ns, nFlatLeafs, hdrVals, xm); /* sign moved up  to hdrVals,  xm sections made into proper HalfPyramids */

                    IF(tree->NsHdr == 1)
                    {
                        /* hdr fits in a single leaf */
                        vec2mind_fx(tree->Ns, tree->Kp, hdrVals, &UL_LS, &UL_idx, &UL_sz, h_mem);  /* actual single idx hdr enumeration */
                        ASSERT(UL_idx <= INT32_MAX && UL_sz <= INT32_MAX && UL_idx < UL_sz);

                        tree->hdrIdx = L_add((Word32)UL_idx, 0);
                        tree->hdrSz = L_add((Word32)UL_sz, 0);
                        tree->LS = extract_l((Word32)UL_LS);
                    }
                    ELSE
                    {
                        ASSERT(tree->NsHdr == LL_HPVC_NSHDR_MAX );
                        /* printf("\n tree with  split header "); */
                        /* split hdr  +   flat split */
                        /* hdr half pyramid  to be split into 3 uniform parts  */

                        hpvc_get_quasi_uniform_sizes(tree->Ns, tree->NsHdr, nHdrLeafs);

                        createHdrL1part(hdrVals, tree->NsHdr, nHdrLeafs, splitHdrVals);

                        setHdrLSmodifyLeaf2HP(tree->NsHdr, nHdrLeafs, splitHdrVals, hdrVals); /* sign moved to splitHdrVals, hdrVals sections made into proper HalfPyramids */

                        ASSERT(tree->NsHdr == LL_HPVC_NSHDR_MAX);
                        ASSERT(tree->Kp <= LL_HPVC_KP_MAX);
                        vec2mind_fx(tree->NsHdr, tree->Kp, splitHdrVals, &UL_LS, &UL_idx, &UL_sz, h_mem);  /* actual splitHdr enumeration */

                        tree->hdrIdx = L_add((Word32)UL_idx, 0); /* To RANGE  */

                        tree->hdrSz = L_add((Word32)UL_sz, 0);   /* To RANGE  */
                        tree->LS = extract_l((Word32)UL_LS);      /* To RANGE , as a bit  */

                        acc_pos = 0; move16();
                        FOR(i = 0; i < tree->NsHdr; i++)
                        {
                            tmpK = abs_s(splitHdrVals[i]);

                            tree->splitHdrLeafSz[i] = L_sub(0,1L);  /* required sentinel for NoTx (tmpK==0) in  AriCodec */

                            IF( tmpK != 0 )
                            {
                                vec2mind_fx(nHdrLeafs[i], tmpK, &(hdrVals[acc_pos]), &UL_LS, &UL_idx, &UL_sz, h_mem);  /* actual splitHdr leaf enumeration */
                                ASSERT(UL_idx <= INT32_MAX && UL_sz <= INT32_MAX && UL_idx < UL_sz);
                                tree->splitHdrLeafIdx[i] = L_add((Word32)UL_idx, 0);
                                tree->splitHdrLeafSz[i] = L_add((Word32)UL_sz, 0);
                            }

                            acc_pos = add(acc_pos, nHdrLeafs[i]);
                        }
                    }

                    /*  enum the flat tree  part */
                    Word16 acc_pos = 0; move16();
                    FOR(i = 0; i < tree->Ns; i++)
                    {
                        flatK = abs_s(hdrVals[i]);
                        tree->flatLeafSz[i] = L_sub(0, 1L);  /* "-1" is a sentinel to ari_codec ,  for teh case flatK==0 */
                        IF(flatK != 0)
                        {
                            ASSERT((HPVC_tabledKMAX[MIN(65, nFlatLeafs[i])] > 0 && flatK <= HPVC_tabledKMAX[MIN(65, nFlatLeafs[i])]) && "Flat leaf size too large");
                            vec2mind_fx(nFlatLeafs[i], flatK, &(xm[acc_pos]), &UL_LS, &UL_idx, &UL_sz, h_mem);  /* actual splitHdr leaf enumeration */
                            ASSERT(UL_idx <= INT32_MAX && UL_sz <= INT32_MAX && UL_idx < UL_sz);
                            tree->flatLeafIdx[i] = L_add((Word32)UL_idx, 0);
                            tree->flatLeafSz[i] = L_add((Word32)UL_sz, 0);
                            ASSERT(UL_LS == 0);     /* first nz has to be positive at this stage */
                        }

                        acc_pos = add(acc_pos, nFlatLeafs[i]);
                    }

                } /* hdr or splitHdr */
            } /* Kp!=0 */


            nTrees = add(nTrees, 1);
            nCoeffs = add(nCoeffs, tree->Np);
            tree++;

        }
    }

    return nTrees;
}


/* increased tables  for HPVC vs. MPVQ */
static PvqEntry_fx get_size_hpvc_calc_offset_fx(   /* o : size, dim, k_val        */
    Word16 dim_in,   /* i : dimension                */
    Word16 k_val_in, /* i : nb unit pulses           */
    UWord32* h_mem   /* o : offsets                  */
)
{

    PvqEntry_fx entry;
    Word16 kp1;
    Word16 kp2;
    Word16 dim_in_tmp;

#  ifdef DYNMEM_COUNT
    Dyn_Mem_In("get_size_hpvc_calc_offset_fx", sizeof(struct {
        Counter i;
        PvqEntry_fx entry;
        Word16 kp1;
        Word16 kp2;
        Word16 dim_in_tmp;
    }));
#  endif

    entry.dim = dim_in;  move16();
    entry.k_val = k_val_in;  move16();

    entry.index = L_deposit_l(0);
    entry.lead_sign_ind = 0; move16();

    ASSERT(dim_in <= LL_HPVC_NP_MAX);

    dim_in_tmp = s_min((LL_HPVC_NP_MAX/2)+1, dim_in);   /* a table reduction warping */

    ASSERT(HPVC_tabledKMAX[dim_in_tmp] != 0);

    /* tabled values for worst case Ks */   /* made into table lookup for N=8,16,32,64,128 */
    kp1 = add(k_val_in, 1);

    kp2 = add(k_val_in, 2);
    basop_memcpy(/*dst*/ h_mem, HPVC_MPVQ_offs_ptr[dim_in_tmp], kp2 * sizeof(UWord32));

    /* MPVQ special handling of last  U offset  in k+1 column  */
    if ( sub(k_val_in, HPVC_tabledKMAX[dim_in_tmp]) != 0 )
    {
        h_mem[kp1] = UL_lshr(h_mem[kp1], 1); /* (A+1)/2 , convert from  A(K+1) to  U(K+1)  domain */
    }
    entry.size = UL_addNsD(1U, UL_addNsD(h_mem[kp1], UL_lshr(h_mem[k_val_in], 1)));  /* MPVQ size calc. 1 + U(K+1) + (A(K)>>1) */

#  ifdef DYNMEM_COUNT
    Dyn_Mem_Out();
#  endif


    return entry;
}


Word16 hpvc_leaf_dec_deidx_fx(  /* out BER detected 1 , ok==0 */
    Word16* y,                  /* o:   decoded vector (non-scaled int)  */
    const Word16 k_val,         /* i:   number of allocated pulses       */
    const Word16 dim,           /* i:   Length of vector                 */
    const Word16 LS_ind,        /* i; lS index              1 bit        */
    const UWord32 UL_MPVQ_ind   /* i; MPVQ  index                        */
)
{
    Dyn_Mem_Deluxe_In(
        Word16 BER_flag;
    UWord32 h_mem[1 + LL_HPVC_KP_MAX + 1];   /* dynamic buffer offsets may change from starting dimension in decoding  */
    PvqEntry_fx entry; );


    BER_flag = 0;    move16();

    entry = get_size_hpvc_calc_offset_fx(dim, k_val, h_mem);

    entry.lead_sign_ind = LS_ind;  move16();
    entry.index = L_deposit_l(0); /* only  in case dim == 1 */

    ASSERT(dim >= 3);
    entry.index = UL_MPVQ_ind;    /* a move to the struct, can be optimized away */

    /* safety check in case of bit errors */
    if( L_sub(entry.index, entry.size) >= 0 )
    {
            BER_flag = 1;   move16();
    }

    mpvq_deindex_fx(&entry, h_mem, y);  /* actual deindexing,  h_mem is a buffer that will change for every N decrease    */


    Dyn_Mem_Deluxe_Out();

    return BER_flag;
}

#  ifdef __cplusplus
} /* NAMESPACE_VERSION */
#  endif
#  endif

#endif
