/******************************************************************************
*                        ETSI TS 103 634 V1.5.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include "stdint.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"


/* channel coder specific constants and macros */
#define RS16_CW_LEN_MAX 15

#define FEC_N_MODES 4
#define FEC_N_SYNDROMES_MAX 6
#define FEC_N_ERR_POS_MAX 3
#define FEC_N_ELP_COEFF_MAX 4
#define FEC_N_ERR_SYMB_MAX 3
#define FEC_N_MODE_DETECTION_CW 6

#define SYNDROME_IDX(mode_index, cw_index) (((mode_index)*FEC_N_MODE_DETECTION_CW + (cw_index)) * FEC_N_SYNDROMES_MAX)
#define ELP_IDX(mode_index, cw_index) (((mode_index)*FEC_N_MODE_DETECTION_CW + (cw_index)) * FEC_N_ELP_COEFF_MAX)
#define ERR_POS_IDX(mode_index, cw_index) (((mode_index)*FEC_N_MODE_DETECTION_CW + (cw_index)) * FEC_N_ERR_POS_MAX)
#define ERR_SYMB_IDX(mode_index, cw_index) (((mode_index)*FEC_N_MODE_DETECTION_CW + (cw_index)) * FEC_N_ERR_SYMB_MAX)
#define DEG_ELP_IDX(mode_index, cw_index) ((mode_index)*FEC_N_MODE_DETECTION_CW + (cw_index))

#define FEC_TOTAL_SYNDROME_SIZE (FEC_N_SYNDROMES_MAX * FEC_N_MODES * FEC_N_MODE_DETECTION_CW)
#define FEC_TOTAL_ELP_SIZE (FEC_N_ELP_COEFF_MAX * FEC_N_MODES * FEC_N_MODE_DETECTION_CW)
#define FEC_TOTAL_ERR_POS_SIZE (FEC_N_ERR_POS_MAX * FEC_N_MODES * FEC_N_MODE_DETECTION_CW)
#define FEC_TOTAL_ERROR_SIZE (FEC_N_ERR_SYMB_MAX * FEC_N_MODES * FEC_N_MODE_DETECTION_CW)
#define FEC_TOTAL_DEG_ELP_SIZE (FEC_N_MODES * FEC_N_MODE_DETECTION_CW)

#define ERROR_REPORT_BEC_MASK   ((0x0FFF)>>1)
#define ERROR_REPORT_EP1_OK     ((0x1000)>>1)
#define ERROR_REPORT_EP2_OK     ((0x2000)>>1)
#define ERROR_REPORT_EP3_OK     ((0x4000)>>1)
#define ERROR_REPORT_EP4_OK     ((0x8000)>>1)                    
#define ERROR_REPORT_ALL_OK (ERROR_REPORT_EP1_OK | ERROR_REPORT_EP2_OK | ERROR_REPORT_EP3_OK | ERROR_REPORT_EP4_OK)

/* debugging switches */

/* constants concerning mode detection */
#define EP_RISK_THRESH_NS_M 21990
#define EP_RISK_THRESH_NS_E -23
#define EP_RISK_THRESH_OS_M 25166
#define EP_RISK_THRESH_OS_E -10

#define SIMPLE_FLOAT_1_MANTISSA 16384

#define FEC_STATIC static

/* DISCLAIMER: Strict instrumentation of GF16 arithmetic would have to take into account
 * the initial conversion of the arguments from LC3_UINT8 to LC3_INT16 (one move16() per argument).
 * Behind this is the assumption that one would store GF16 elements in LC3_INT16 for strict BASOP
 * implementation.
 */
#define GF16_MUL(a, b) gf16_mult_table[(a) | (b << 4)]
#define GF16_MUL0(a, b) gf16_mult_table[(a) | (b)]
#define GF16_ADD(a, b) ((a) ^ (b))

/* tables for finite field arithmetic */
/* tables for arithmetic in GF(16)  *
 * generator polynomial: 19
 * unit group generator (g): 2
 */

static const LC3_UINT8 gf16_mult_table[256] = {
    /* gf16_mult_table[a | (b << 4)] contains the product of a and b in GF(16) */
    0,  0,  0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
    13, 14, 15, 0,  2,  4,  6,  8, 10, 12, 14, 3,  1,  7,  5,  11, 9,  15, 13, 0,  3,  6,  5,  12, 15, 10, 9,  11, 8,
    13, 14, 7,  4,  1,  2,  0,  4, 8,  12, 3,  7,  11, 15, 6,  2,  14, 10, 5,  1,  13, 9,  0,  5,  10, 15, 7,  2,  13,
    8,  14, 11, 4,  1,  9,  12, 3, 6,  0,  6,  12, 10, 11, 13, 7,  1,  5,  3,  9,  15, 14, 8,  2,  4,  0,  7,  14, 9,
    15, 8,  1,  6,  13, 10, 3,  4, 2,  5,  12, 11, 0,  8,  3,  11, 6,  14, 5,  13, 12, 4,  15, 7,  10, 2,  9,  1,  0,
    9,  1,  8,  2,  11, 3,  10, 4, 13, 5,  12, 6,  15, 7,  14, 0,  10, 7,  13, 14, 4,  9,  3,  15, 5,  8,  2,  1,  11,
    6,  12, 0,  11, 5,  14, 10, 1, 15, 4,  7,  12, 2,  9,  13, 6,  8,  3,  0,  12, 11, 7,  5,  9,  14, 2,  10, 6,  1,
    13, 15, 3,  4,  8,  0,  13, 9, 4,  1,  12, 8,  5,  2,  15, 11, 6,  3,  14, 10, 7,  0,  14, 15, 1,  13, 3,  2,  12,
    9,  7,  6,  8,  4,  10, 11, 5, 0,  15, 13, 2,  9,  6,  4,  11, 1,  14, 12, 3,  8,  7,  5,  10,
};

static const LC3_UINT8 rs16_elp_deg2_table[256] = {
    /* If the polynomial x^2 + ax + b has distinct non-zero roots z1 and z2 in GF(16),  *
     * then table entry a + 16*b contains log_g(z1) | log_g(z2) << 4, and otherwise it  *
     * contains 0.                                                                      */
    0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   165, 0,   0,   0,   0,
    105, 195, 0,   210, 0,   225, 0,   180, 120, 0,   0,   121, 0,   16,  0,   211, 0,   0,   181, 0,   0,   106,
    196, 226, 0,   0,   0,   214, 64,  0,   199, 0,   0,   0,   0,   0,   49,  184, 0,   154, 0,   229, 0,   227,
    182, 0,   0,   32,  0,   0,   0,   197, 0,   0,   122, 0,   212, 152, 0,   203, 0,   158, 128, 0,   0,   0,
    98,  113, 218, 0,   0,   0,   53,  0,   0,   65,  0,   0,   185, 110, 215, 80,  0,   0,   200, 0,   50,  0,
    0,   0,   0,   130, 205, 115, 0,   0,   160, 190, 145, 0,   0,   0,   0,   0,   0,   100, 0,   0,   168, 198,
    0,   183, 33,  0,   0,   48,  228, 213, 0,   0,   0,   0,   0,   0,   0,   0,   164, 0,   179, 0,   224, 104,
    0,   194, 149, 0,   0,   209, 0,   0,   0,   189, 99,  84,  0,   129, 0,   0,   0,   144, 0,   0,   234, 114,
    0,   0,   82,  0,   0,   0,   0,   217, 202, 0,   112, 52,  232, 0,   97,  0,   0,   0,   126, 0,   81,  201,
    0,   36,  216, 186, 0,   0,   0,   96,  0,   0,   0,   0,   0,   88,  0,   0,   0,   103, 0,   148, 178, 0,
    208, 193, 0,   58,  0,   0,   0,   0,   0,   161, 206, 0,   116, 0,   101, 0,   0,   56,  146, 176, 0,   0,
    147, 162, 222, 0,   132, 0,   0,   0,   0,   0,   177, 117, 192, 0,
};

static const LC3_UINT16 rs16_elp_deg3_table[256] = {
    /* If the polynomial x^3 + ax + b has distinct roots z1, z2 and z3 in GF(16),                       *
     * then table entry a + 16*b contains z1) | z2 << 4 | z3 << 8, and otherwise it                     *
     * contains 0.                                                                                      */
    0, 0, 0,    0,    0,    0,    0,    0,    0,    0,    0, 0,    0, 0,    0,    0,   1889, 0, 0,    0,    0,    0,
    0, 0, 0,    0,    0,    0,    0,    0,    0,    0,    0, 0,    0, 2977, 0,    0,   0,    0, 0,    3990, 1859, 0,
    0, 0, 0,    0,    0,    0,    3521, 0,    0,    0,    0, 0,    0, 0,    0,    0,   1874, 0, 3718, 0,    0,    0,
    0, 0, 0,    2433, 0,    0,    1619, 0,    0,    0,    0, 3495, 0, 0,    0,    0,   0,    0, 4065, 0,    0,    0,
    0, 0, 0,    3255, 0,    0,    0,    1602, 0,    3735, 0, 0,    0, 0,    3238, 801, 0,    0, 0,    0,    0,    0,
    0, 0, 0,    3510, 0,    0,    0,    0,    1345, 3975, 0, 0,    0, 0,    0,    0,   0,    0, 3778, 0,    0,    0,
    0, 0, 0,    0,    0,    0,    0,    0,    0,    0,    0, 0,    0, 0,    2947, 0,   0,    0, 0,    0,    0,    0,
    0, 0, 3476, 0,    4005, 0,    3461, 0,    0,    0,    0, 0,    0, 0,    0,    0,   0,    0, 0,    0,    0,    0,
    0, 0, 0,    0,    0,    3748, 0,    0,    2962, 0,    0, 0,    0, 4035, 0,    0,   4020, 0, 0,    0,    0,    0,
    0, 0, 0,    0,    0,    0,    0,    0,    0,    0,    0, 0,    0, 0,    3221, 0,   0,    0, 0,    0,    0,    2690,
    0, 0, 0,    3795, 0,    0,    0,    4050, 0,    0,    0, 0,    0, 3204, 3765, 0,   0,    0, 0,    0,    2707, 0,
    0, 0, 0,    0,    0,    0,    0,    0,    0,    0,    0, 0,    0, 0,
};

static const LC3_UINT8 gf16_g_pow[16] = {1, 2, 4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13, 9, 1};
/* g_pow[i] contains g^i*/

static const LC3_UINT8 gf16_log_g[16] = {255, 0, 1, 4, 2, 8, 5, 10, 3, 14, 9, 7, 6, 13, 11, 12};
/* log_g[n] contains contains the value i such that g^i = n for n=1, 2, ..., 15, log_g[0] is set to 255 */

static const LC3_UINT8 gf16_inv_table[16] = {255, 1, 9, 14, 13, 11, 7, 6, 15, 2, 12, 5, 10, 4, 3, 8};
/* gf16_inv_table[n] contains the multiplicative inverse of n in GF(16) (1/0 is set to 255)*/

/* RS16 generating polynomials (from lowest to highest coefficient without leading 1)*/

static const LC3_UINT8 rs16_gp_d3[] = {8, 6};
static const LC3_UINT8 rs16_gp_d5[] = {7, 8, 12, 13};
static const LC3_UINT8 rs16_gp_d7[] = {12, 10, 12, 3, 9, 7};

/* FEC mode signaling polynomials */

#define EP_SIG_POLY_DEG 12

static const LC3_UINT8 sig_polys[4][15] = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                                        {7, 15, 5, 6, 14, 9, 1, 3, 12, 10, 13, 3, 2, 0, 0},
                                        {7, 11, 14, 1, 2, 3, 12, 11, 6, 15, 7, 6, 12, 0, 0},
                                        {6, 15, 12, 2, 9, 15, 2, 8, 12, 3, 10, 5, 4, 0, 0}};

static const LC3_UINT8 sig_poly_syndr[4][6] = {
    {0, 0, 0, 0, 0, 0}, {0, 4, 5, 11, 5, 8}, {0, 5, 9, 0, 1, 7}, {0, 12, 5, 12, 9, 8}};

/* bit count table for error report (0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111) */

static const LC3_UINT8 rs16_bit_count_table[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

/* List of RS16 generators by Hamming distance */

static const LC3_UINT8 *const rs16_gp_by_hd[8] = {NULL, NULL, NULL, rs16_gp_d3, NULL, rs16_gp_d5, NULL, rs16_gp_d7};

/* fec config data */

static const LC3_UINT8 hamming_distance_by_mode0[] = {1, 3, 3, 5, 7};
static const LC3_UINT8 hamming_distance_by_mode1[] = {1, 1, 3, 5, 7};

static const LC3_UINT8 crc1_bytes_by_mode0[] = {0, 3, 2, 2, 2};
static const LC3_UINT8 crc1_bytes_by_mode1[] = {0, 3, 3, 3, 3};
static const LC3_UINT8 crc2_bytes_by_mode[]  = {0, 0, 2, 2, 2};

/* fec mode risk table */
typedef struct
{
    LC3_UINT32 mantissa;
    LC3_INT16  exponent;
} simple_float;

static const simple_float risk_table_f[4][4] = {{{16384, 0}, {16384, 0}, {16384, 0}, {16384, 0}},
                                                {{16384, -8}, {26880, -1}, {16384, 0}, {16384, 0}},
                                                {{16384, -16}, {26880, -9}, {20475, -2}, {16384, 0}},
                                                {{16384, -24}, {26880, -17}, {20475, -10}, {19195, -4}}};
/* bit error limits for slot size 40 */
static LC3_INT16 const low_br_max_bit_errors_by_mode[] = {0, 0, 3, 9, 18};

/*
corresponding float values:
    {1.f, 1.f, 1.f, 1.f},
    {0.00390625f, 0.820312f, 1.f, 1.f},
    {1.52588e-05f, 0.00320435f, 0.312424f, 1.f},
    {5.96046e-08f, 1.2517e-05f, 0.00122041f, 0.0732243f}
*/

/* internal encoder routines */

FEC_STATIC void fec_interleave_pack(LC3_UINT8 *out, LC3_UINT8 *in, LC3_INT16 n_nibbles, LC3_INT16 n_codewords);

FEC_STATIC void rs16_enc(LC3_UINT8 *iobuf, LC3_INT16 codeword_length, LC3_INT16 hamming_distance, LC3_INT16 fec_mode,
                         LC3_INT16 signal_mode);

/* internal decoder routines */

FEC_STATIC void fec_deinterleave_unpack(LC3_UINT8 *out, LC3_UINT8 *in, LC3_INT16 n_nibbles, LC3_INT16 n_codewords);

FEC_STATIC LC3_INT16 fec_data_preproc(LC3_INT16 mode, LC3_INT16 epmr, LC3_UINT8 *iobuf, LC3_UINT8 *cw_buf, LC3_INT16 data_bytes,
                                   LC3_INT16 slot_bytes, LC3_INT16 pc_split);

FEC_STATIC void fec_data_postproc(LC3_INT16 mode, LC3PLUS_EpModeRequest *epmr, LC3_UINT8 *iobuf, LC3_INT16 data_bytes, LC3_UINT8 *cw_buf,
                                  LC3_INT16 slot_bytes, LC3_INT16 pc_split, LC3_INT32 *bfi);

FEC_STATIC LC3_INT32 rs16_detect_and_correct(LC3_UINT8 *iobuf, LC3_INT32 n_symb, LC3_INT32 n_codewords, LC3PLUS_EpModeRequest *epmr, LC3_INT16 *error_report,
                                       LC3_INT32 *bfi, LC3_UINT8 *array_of_trust, LC3_INT32 ccc_flag_flag, LC3_INT16 *n_pccw);

FEC_STATIC void rs16_calculate_six_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg);

FEC_STATIC void rs16_calculate_four_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg);

FEC_STATIC void rs16_calculate_two_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg);

FEC_STATIC LC3_INT8 rs16_calculate_elp(LC3_UINT8 *elp, LC3_UINT8 *syndromes, LC3_INT16 hamming_distance);

FEC_STATIC LC3_INT16 rs16_factorize_elp(LC3_UINT8 *error_locations, LC3_UINT8 *elp, LC3_INT16 deg_elp, LC3_INT16 max_pos);

FEC_STATIC void rs16_calculate_errors(LC3_UINT8 *errors, LC3_UINT8 *err_pos, LC3_UINT8 *syndromes, LC3_INT8 deg_elp, LC3_INT8 t);

/* auxiliary routines */

FEC_STATIC LC3_INT16 crc1(LC3_UINT8 *data, LC3_INT16 data_size, LC3_INT16 epmr, LC3_UINT8 *hash_val, LC3_INT16 hash_size, LC3_INT16 check);

FEC_STATIC LC3PLUS_EpModeRequest fec_estimate_epmr_from_cw0(LC3_UINT8 *cw0, LC3_INT8 *t, LC3_UINT8 *syndromes, LC3_UINT8 *elp, LC3_INT8 *deg_elp,
                                            LC3_UINT8 *err_pos, LC3_UINT8 *err_symb, LC3_INT16 n_codewords, LC3_INT16 n_symb);

FEC_STATIC void dw0_bitswap(LC3_UINT8 *dw0, LC3_INT16 mode, LC3_INT16 slot_bytes);

FEC_STATIC LC3PLUS_EpModeRequest cw0_get_epmr(LC3_UINT8 *cw0, LC3_INT16 epmr_position);

FEC_STATIC LC3PLUS_EpModeRequest dw0_get_epmr(LC3_UINT8 *dw0, LC3_INT16 mode, LC3_INT16 slot_size);

FEC_STATIC LC3_INT16 crc2(LC3_UINT8 *data, LC3_INT16 data_size, LC3_UINT8 *hash_val, LC3_INT16 hash_size, LC3_INT16 check);

FEC_STATIC simple_float simple_float_mul(simple_float op1, simple_float op2);

FEC_STATIC LC3_INT16 simple_float_cmp(simple_float op1, simple_float op2);

FEC_STATIC LC3_INT16 get_total_crc_size(LC3_INT16 slot_bytes, LC3_INT16 fec_mode, LC3_INT16 pc_split);

FEC_STATIC LC3_INT16 get_n_codewords(LC3_INT16 slot_bytes);

FEC_STATIC LC3_INT16 get_codeword_length(LC3_INT16 n_codewords, LC3_INT16 slot_nibbles, LC3_INT16 codeword_index);



LC3_INT16 fec_get_n_pccw(LC3_INT16 slot_bytes, LC3_INT16 fec_mode, LC3_INT16 ccc_flag)
{
    LC3_INT16 n_pccw;

    if (fec_mode ==  3)
    {
        n_pccw = (LC3_INT16) (0.080447761194030 * slot_bytes - 1.791044776119394 + 0.5);
    }
    else if (fec_mode == 4)
    {
        n_pccw = (LC3_INT16) (0.066492537313433 * slot_bytes - 1.970149253731338 + 0.5);
    }
    else
    {
        n_pccw = 0;
    }

    if (ccc_flag == 1 || slot_bytes < 80)
    {
        n_pccw = 0;
    }

    return n_pccw;
}

FEC_STATIC LC3_INT16 get_total_crc_size(LC3_INT16 slot_bytes, LC3_INT16 fec_mode, LC3_INT16 pc_split)
{
        LC3_INT16 n_crc;

    n_crc = crc1_bytes_by_mode1[fec_mode];
    if (slot_bytes == 40)
    {
        n_crc = crc1_bytes_by_mode0[fec_mode];
    }

    if (pc_split > 0)
    {
        n_crc = n_crc + crc2_bytes_by_mode[fec_mode];
    }

     
    
    return n_crc;
}

FEC_STATIC LC3_INT16 get_n_codewords(LC3_INT16 slot_bytes)
{
    return (2*slot_bytes + 14)/15;
}

FEC_STATIC LC3_INT16 get_codeword_length(LC3_INT16 n_codewords, LC3_INT16 slot_nibbles, LC3_INT16 codeword_index)
{
    return (slot_nibbles - codeword_index - 1) / n_codewords + 1;
}

/* Encoder */

LC3_INT16 fec_get_data_size(LC3_INT16 fec_mode, LC3_INT16 ccc_flag, LC3_INT16 slot_bytes)
/* not time critical */
{
        LC3_INT16 n_codewords, payload_size;

    n_codewords = get_n_codewords(slot_bytes);

    assert(n_codewords == (2 * slot_bytes + RS16_CW_LEN_MAX - 1) / RS16_CW_LEN_MAX);
    payload_size = slot_bytes;

    if (fec_mode > 0)
    {
        if (fec_mode == 1)
        {
            payload_size --;
        }
        else
        {
            payload_size -= (fec_mode - 1) * n_codewords;
        }
        if (slot_bytes == 40)
        {
            payload_size -= crc1_bytes_by_mode0[fec_mode];
        }
        else
        {
            payload_size -= crc1_bytes_by_mode1[fec_mode];
        }

        if (ccc_flag == 0 && fec_mode > 2 && slot_bytes >= 80)
        {
            payload_size -= crc2_bytes_by_mode[fec_mode];
        }
    }
    
     

    return payload_size;
}

LC3_INT16 fec_get_n_pc(LC3_INT16 fec_mode, LC3_INT16 n_pccw, LC3_INT16 slot_bytes)
/* not time critical */
{
        LC3_INT16 n_codewords, pc_split;
        LC3_INT32    i;

    n_codewords = get_n_codewords(slot_bytes);

    assert(n_codewords == (2 * slot_bytes + RS16_CW_LEN_MAX - 1) / RS16_CW_LEN_MAX);

    pc_split = - 2*n_pccw*(fec_mode - 1);

    if (fec_mode == 1 || slot_bytes < 80)
    {
        pc_split = 0;
    }
    else
    {
        for (i = 0; i < n_pccw; i++)
        {
            pc_split += (2 * slot_bytes + i) / n_codewords;
        }
    }
    
     

    return pc_split;
}

/* functions for EPMR handling */
FEC_STATIC void dw0_bitswap(LC3_UINT8 *dw0, LC3_INT16 mode, LC3_INT16 slot_bytes)
/* swap epmr bits with bits that will be positioned at 30 and 32 in code word 0 */
{
        LC3_UINT8 tmp;
        LC3_INT32    ind0, ind1, position;

    position = get_codeword_length(get_n_codewords(slot_bytes), 2*slot_bytes, 0) - 1;

    if (slot_bytes == 40)
    {
        ind0 = 2*crc1_bytes_by_mode0[mode] - 1;
    }
    else
    {
        ind0 = 2*crc1_bytes_by_mode1[mode] - 1;
    }

    ind1 = position - hamming_distance_by_mode0[mode] + 1;

    /* swap bits 2 and 3 of dw0[ind0] with bits 0 and 1 of dw0[ind1] */
    tmp = (dw0[ind0] >> 2) & 3;
    dw0[ind0] = dw0[ind0] & 3;
    dw0[ind0] = dw0[ind0] | ((dw0[ind1] & 3) << 2);
    dw0[ind1] = dw0[ind1] & 12;
    dw0[ind1] = dw0[ind1] | tmp;

     
}

FEC_STATIC LC3PLUS_EpModeRequest cw0_get_epmr(LC3_UINT8 *cw0, LC3_INT16 position)
{
    return (LC3PLUS_EpModeRequest)(cw0[position] & 3);
}

FEC_STATIC LC3PLUS_EpModeRequest dw0_get_epmr(LC3_UINT8 *dw0, LC3_INT16 mode, LC3_INT16 slot_size)
{
        LC3_INT32    ncrc1;
        LC3PLUS_EpModeRequest epmr;

    ncrc1 = crc1_bytes_by_mode1[mode];

    if (slot_size == 40)
    {
        ncrc1 = crc1_bytes_by_mode0[mode];
    }

    epmr = (LC3PLUS_EpModeRequest)(dw0[2 * ncrc1 - 1] >> 2);
    
     

    return epmr;
}


FEC_STATIC LC3_INT16 fec_data_preproc(LC3_INT16 mode, LC3_INT16 epmr, LC3_UINT8 *iobuf, LC3_UINT8 *cw_buf, LC3_INT16 data_bytes,
                                   LC3_INT16 slot_bytes, LC3_INT16 pc_split)
{
        LC3_INT16 data_offset, n_crc1, n_crc2;
        LC3_INT32    i, j;

    data_offset = 2*(slot_bytes - data_bytes);

    /* extract and reverse data*/
    j = 2*slot_bytes - 1;
    for (i = 0; i < data_bytes; i++)
    {
        cw_buf[j--] = iobuf[i] & 15;
        cw_buf[j--] = iobuf[i] >> 4;
    }

    /* add crc hashes */
    if (slot_bytes == 40)
    {
        n_crc1 = crc1_bytes_by_mode0[mode];
    }
    else
    {
        n_crc1 = crc1_bytes_by_mode1[mode];
    }

    if (pc_split > 0 && mode > 1)
    {
        n_crc2 = crc2_bytes_by_mode[mode];
    }
    else
    {
        n_crc2 = 0;
    }

    if (n_crc2)
    {
        crc2(cw_buf + data_offset + 2 * data_bytes - pc_split, pc_split, cw_buf + data_offset - 2 * n_crc2, n_crc2, 0);
    }
    if (n_crc1)
    {
        crc1(cw_buf + data_offset, 2 * data_bytes - pc_split, epmr, cw_buf + data_offset - 2 * (n_crc1 + n_crc2), n_crc1,
             0);
    }

    data_offset -= 2* (n_crc1 + n_crc2);

    dw0_bitswap(cw_buf + data_offset, mode, slot_bytes);
    
     

    return data_offset;
}

void fec_encoder(LC3_INT16 mode, LC3_INT16 epmr, LC3_UINT8 *iobuf, LC3_INT16 data_bytes, LC3_INT16 slot_bytes, LC3_INT16 n_pccw)
{
    LC3_INT16  n_codewords, codeword_length, hd, redundancy_nibbles, cw_offset, dw_offset, pc_split;
    LC3_INT32     i, j;
    LC3_UINT8 cw_buf[2 * FEC_SLOT_BYTES_MAX];

    cw_offset = 0;
    dw_offset = 0;
    pc_split  = 0;

    n_codewords = get_n_codewords(slot_bytes);

    /* some sanity checks */
    {
        LC3_INT32 tmp = slot_bytes;
        
        assert((slot_bytes >= FEC_SLOT_BYTES_MIN && slot_bytes <= FEC_SLOT_BYTES_MAX) &&
               "fec_encoder: slot_bytes out of range");
        tmp -= mode == 1 ? 1 : n_codewords * (mode - 1);                                 // reed solomon redundancy
        tmp -= slot_bytes == 40 ? crc1_bytes_by_mode0[mode] : crc1_bytes_by_mode1[mode]; // crc1
        tmp -= (n_pccw > 0) && (mode > 1) ? crc2_bytes_by_mode[mode] : 0;                // crc2
        assert(data_bytes == tmp && "fec_encoder: inconsistent payload size");
        assert(n_codewords - n_pccw >= 6);
    }

    /* data preproc: re-ordering and hash extension */
    pc_split = fec_get_n_pc(mode, n_pccw, slot_bytes);

    dw_offset = fec_data_preproc(mode, epmr, iobuf, cw_buf, data_bytes, slot_bytes, pc_split);

    /* encoding of first data word*/
    hd                 = hamming_distance_by_mode0[mode];
    redundancy_nibbles = hd - 1;
    codeword_length    = get_codeword_length(n_codewords, 2 * slot_bytes, 0);

    assert(codeword_length == (2 * slot_bytes - 1) / n_codewords + 1);

    for (j = redundancy_nibbles; j < codeword_length; (j++, dw_offset++))
    {
        cw_buf[j] = cw_buf[dw_offset];
    }

    rs16_enc(cw_buf, codeword_length, hd, mode, 1);

    cw_offset += codeword_length;

    /* encoding of remaining data words */
    hd                 = hamming_distance_by_mode1[mode];
    redundancy_nibbles = hd - 1;

    for (i = 1; i < n_codewords; i++)
    {
        codeword_length = get_codeword_length(n_codewords, 2*slot_bytes, i);

        for (j = redundancy_nibbles; j < codeword_length; (j++, dw_offset++))
        {
            cw_buf[cw_offset + j] = cw_buf[dw_offset];
        }

        rs16_enc(cw_buf + cw_offset, codeword_length, hd, mode, i < 6);

        cw_offset += codeword_length;
    }

    assert(cw_offset == 2 * slot_bytes && dw_offset == 2 * slot_bytes);

    fec_interleave_pack(iobuf, cw_buf, 2 * slot_bytes, n_codewords);

     
}

FEC_STATIC void rs16_enc(LC3_UINT8 *iobuf, LC3_INT16 codeword_length, LC3_INT16 hamming_distance, LC3_INT16 fec_mode,
                         LC3_INT16 signal_mode)
/* expects (data polynomial) * x^(hamming_distance - 1) in iobuf */
{
        LC3_UINT8 const *gp;
        LC3_UINT8        shift_buffer[RS16_CW_LEN_MAX + 1], lc;
        LC3_INT32          i, j, deg_gp;

    memset(shift_buffer, 0, sizeof(shift_buffer));
    gp     = rs16_gp_by_hd[hamming_distance];
    deg_gp = hamming_distance - 1;

    if (hamming_distance > 1)
    {
        assert(codeword_length > deg_gp);

        /* initialize redundancy part to zero */
        memset(iobuf, 0, deg_gp);

        /* initialize shift_buffer */
        memmove(shift_buffer + 1, iobuf + codeword_length - deg_gp, deg_gp);

        /* calculate remainder */
        for (i = codeword_length - deg_gp - 1; i >= 0; i--)
        {
            shift_buffer[0] = iobuf[i];
            lc              = shift_buffer[deg_gp] << 4;

            for (j = deg_gp - 1; j >= 0; j--)
            {
                shift_buffer[j + 1] = GF16_ADD(shift_buffer[j], GF16_MUL0(gp[j], lc));
            }
        }

        /* add remainder to shifted data polynomial */
        for (i = 0; i < deg_gp; i++)
        {
            iobuf[i] = shift_buffer[i + 1];
        }

        /* add signaling polynomial */
        if (signal_mode)
        {
            assert(codeword_length > EP_SIG_POLY_DEG);
            for (i = 0; i <= EP_SIG_POLY_DEG; i++)
            {
                iobuf[i] = GF16_ADD(iobuf[i], sig_polys[fec_mode - 1][i]);
            }
        }
    }

     
}

FEC_STATIC void fec_interleave_pack(LC3_UINT8 *out, LC3_UINT8 *in, LC3_INT16 n_nibbles, LC3_INT16 n_codewords)
{
        LC3_INT16 out_offset, cw_offset, codeword_length;
        LC3_INT32   i, j;

    out_offset = 0;
    cw_offset  = 0;

    /* initialize output buffer to zero */
    memset(out, 0, n_nibbles >> 1);

    /* interleave and pack codewords */
    for (i = 0; i < n_codewords; i++)
    {
        codeword_length = get_codeword_length(n_codewords, n_nibbles, i);

        for (j = 0; j < codeword_length; j++)
        {
                out_offset = n_nibbles - 1 - j*n_codewords - i;
            out[out_offset >> 1] |= in[cw_offset] << ((out_offset & 1) << 2);
            cw_offset = cw_offset + 1;
        }
    }

     
    assert(cw_offset == n_nibbles);
}

/* Decoder */
FEC_STATIC void fec_data_postproc(LC3_INT16 mode, LC3PLUS_EpModeRequest *epmr, LC3_UINT8 *obuf, LC3_INT16 data_bytes, LC3_UINT8 *cw_buf,
                                  LC3_INT16 slot_bytes, LC3_INT16 pc_split, LC3_INT32 *bfi)
{
        LC3_INT16 i;
        LC3_INT16 n_crc1, n_crc2;
        LC3_INT16 cw_buf_len;
        LC3PLUS_EpModeRequest tmp_epmr;

    n_crc1 = crc1_bytes_by_mode1[mode];
    if (slot_bytes == 40)
    {
        n_crc1 = crc1_bytes_by_mode0[mode];
    }

    n_crc2 = 0;
    if (pc_split > 0)
    {
        n_crc2 = crc2_bytes_by_mode[mode];
    }

    assert(n_crc1 == (slot_bytes == 40 ? crc1_bytes_by_mode0[mode] : crc1_bytes_by_mode1[mode]));
    assert(n_crc2 == ((pc_split > 0) && (mode > 1) ? crc2_bytes_by_mode[mode] : 0));

    cw_buf_len = 2 * (data_bytes + n_crc1 + n_crc2);

    if ((mode - 1))
    {
        /* reverse bit-swap */
        dw0_bitswap(cw_buf, mode, slot_bytes);
        tmp_epmr = dw0_get_epmr(cw_buf, mode, slot_bytes);

        if (crc1(cw_buf + ((n_crc1 + n_crc2) << 1), ((data_bytes << 1) - pc_split), tmp_epmr, cw_buf, n_crc1, 1))
        {
            *bfi = 1;
             
            return;
        }
        else
        {
            *epmr = tmp_epmr;
        }
    }

    if (pc_split > 0 && *bfi != 2)
    {
        if (crc2(cw_buf + (((data_bytes + (n_crc1 + n_crc2)) << 1) - pc_split), pc_split,
                 cw_buf + (n_crc1 << 1), n_crc2, 1))
        {
            *bfi = 2;
        }
    }

    for (i = 0; i < data_bytes; i++)
    {
        obuf[i] = (LC3_UINT8)(cw_buf[cw_buf_len - 2 * i - 1] | (cw_buf[cw_buf_len - 2 * i - 2] << 4));
    }

     
}

LC3_INT32 fec_decoder(LC3_UINT8 *iobuf, LC3_INT16 slot_bytes, LC3_INT32 *data_bytes, LC3PLUS_EpModeRequest *epmr, LC3_INT16 ccc_flag, LC3_INT16 *n_pccw,
                LC3_INT32 *bfi, LC3_INT16 *be_bp_left, LC3_INT16 *be_bp_right, LC3_INT16 *n_pc, LC3_INT16 *m_fec)
{
        LC3_UINT8 cw_buf[2 * FEC_SLOT_BYTES_MAX];
        LC3_UINT8 array_of_trust[MAX_LEN];
        LC3_INT16  i, j;
        LC3_INT16  cw_offset, dw_offset;
        LC3_INT16  n_codewords, redundancy_nibbles, codeword_length;
        LC3_INT16  mode, error_report;
        LC3_INT16  n_crc;
        LC3_INT16  first_bad_cw;
        LC3_INT16  pc_split;
    
    UNUSED(n_crc);
    

    if (*bfi == 1)
    {
        return ERROR_REPORT_BEC_MASK;
    }

    if (slot_bytes < FEC_SLOT_BYTES_MIN || slot_bytes > FEC_SLOT_BYTES_MAX)
    {
        *bfi = 1;
         
        return ERROR_REPORT_BEC_MASK;
    }

    if (ccc_flag == 0)
    {
        *be_bp_left = -1;
        *be_bp_right = -1;
    }

    n_codewords = get_n_codewords(slot_bytes);

    /* extract and de-interleave nibbles */
    fec_deinterleave_unpack(cw_buf, iobuf, 2 * slot_bytes, n_codewords);

    /* mode detection and error correction */
    mode = rs16_detect_and_correct(cw_buf, 2 * slot_bytes, n_codewords, epmr, &error_report, bfi, array_of_trust,
                                   ccc_flag, n_pccw);

    /* for normal slots the maximal number of bit errors is limited */
#ifndef APPLY_MAX_ERRORS
    if (slot_bytes == 40 && mode > 0)
    {
    if ((error_report & ERROR_REPORT_BEC_MASK) > low_br_max_bit_errors_by_mode[mode])
    {
        error_report &= ERROR_REPORT_BEC_MASK;
            mode = -1;
            *bfi = 1;
        }
        else
        {
            if ((error_report & ERROR_REPORT_BEC_MASK) > low_br_max_bit_errors_by_mode[2])
            {
                error_report &= ~ERROR_REPORT_EP2_OK;
            }
            if ((error_report & ERROR_REPORT_BEC_MASK) > low_br_max_bit_errors_by_mode[3])
            {
                error_report &= ~ERROR_REPORT_EP3_OK;
            }
        }
    }
#endif
    
    if (*bfi == 1)
    {
        *data_bytes = 0;

         
        return error_report;
    }

    /* initialization for decoding */
    *data_bytes = fec_get_data_size(mode, ccc_flag, slot_bytes);
    pc_split    = fec_get_n_pc(mode, *n_pccw, slot_bytes);
    n_crc       = get_total_crc_size(slot_bytes, mode, pc_split);

    /* decoding of first code word */
    redundancy_nibbles = hamming_distance_by_mode0[mode] - 1;
    codeword_length    = get_codeword_length(n_codewords, slot_bytes + slot_bytes, 0);

    dw_offset = 0;
    cw_offset = 0;

    for (j = redundancy_nibbles; j < codeword_length; j++)
    {
        cw_buf[dw_offset++] = cw_buf[j];
    }
    cw_offset = cw_offset + codeword_length;

    /* decoding of remaining code words */
    redundancy_nibbles = hamming_distance_by_mode1[mode] - 1;

    for (i = 1; i < n_codewords; i++)
    {
        codeword_length = get_codeword_length(n_codewords, slot_bytes + slot_bytes, i);

        for (j = redundancy_nibbles; j < codeword_length; j++)
        {
            cw_buf[dw_offset++] = cw_buf[j + cw_offset];
        }

        cw_offset = cw_offset + codeword_length;
    }

    /* data postproc: hash validation and re-ordering */

    fec_data_postproc(mode, epmr, iobuf, *data_bytes, cw_buf, slot_bytes, pc_split, bfi);

    if (*bfi == 1)
    {
        *data_bytes = 0;

        error_report &= ERROR_REPORT_BEC_MASK;

         
        return error_report;
    }

    if (*bfi == 2)
    {
        first_bad_cw            = 0;
        array_of_trust[*n_pccw] = 0;
        while (array_of_trust[first_bad_cw] != 0)
        {
            first_bad_cw = first_bad_cw + 1;
        }
        if (first_bad_cw == *n_pccw)
        {
            /* this is the case when CRC failed */
            *be_bp_left = 0;
        }
        else
        {
            *be_bp_left = 4*fec_get_n_pc(mode, first_bad_cw, slot_bytes);
        }

        for (i = *n_pccw - 1; i >= 0; i--)
        {
            if (!array_of_trust[i])
            {
                break;
            }
        }
        if (i < 0)
        {
            i = *n_pccw - 1;
        }
        *be_bp_right = 4*fec_get_n_pc(mode, i + 1, slot_bytes) - 1;
    }

    if (ccc_flag == 0)
    {
        *n_pc  = pc_split;
        *m_fec = mode;
    }

     
    return error_report;
}

FEC_STATIC void fec_deinterleave_unpack(LC3_UINT8 *out, LC3_UINT8 *in, LC3_INT16 n_nibbles, LC3_INT16 n_codewords)
{
        LC3_INT16 in_offset, out_offset, codeword_length;
        LC3_INT32   i, j;

    in_offset  = 0;
    out_offset = 0;

    /* unpack nibbles in input buffer and deinterleave codewords */
    for (i = 0; i < n_codewords; i++)
    {
        codeword_length = get_codeword_length(n_codewords, n_nibbles, i);
        for (j = 0; j < codeword_length; (j++, out_offset++))
        {
            in_offset = n_nibbles - 1 - j*n_codewords - i;
            out[out_offset] = (in[in_offset >> 1] >> ((in_offset & 1) << 2)) & 15;
        }
    }

     
    assert(out_offset == n_nibbles);

}

FEC_STATIC LC3PLUS_EpModeRequest fec_estimate_epmr_from_cw0(LC3_UINT8 *cw0, LC3_INT8 *t, LC3_UINT8 *syndromes, LC3_UINT8 *elp, LC3_INT8 *deg_elp,
                                            LC3_UINT8 *err_pos, LC3_UINT8 *err_symb, LC3_INT16 n_codewords, LC3_INT16 n_symb)
{
        LC3_INT32    epmr_lowest_risk_exp;
        LC3_INT32    start, inc, i, n_candidates;
        LC3_INT32    first_codeword_length;
        LC3_INT32    mode_counter;
        LC3PLUS_EpModeRequest  epmr;

    epmr_lowest_risk_exp   = 0;
    first_codeword_length = get_codeword_length(n_codewords, n_symb, 0);
    start                 = 2;
    inc                   = 1;
    n_candidates          = 0;

    /* test if first code word decodes in mode 0 or 1 without error correction */
    if ((syndromes[SYNDROME_IDX(0, 0)] | syndromes[SYNDROME_IDX(0, 0) + 1]) == 0 ||
        (syndromes[SYNDROME_IDX(1, 0)] | syndromes[SYNDROME_IDX(1, 0) + 1]) == 0)
    {
        epmr_lowest_risk_exp = risk_table_f[1][0].exponent;
    }
    /* test if first code word decodes in mode 2 or 3 with lower risk */
    if (deg_elp[DEG_ELP_IDX(2, 0)] <= t[2])
    {
        if (risk_table_f[2][deg_elp[DEG_ELP_IDX(2, 0)]].exponent <= -8)
        {
            n_candidates++;
            start = 2;
        }
    }

    if (deg_elp[DEG_ELP_IDX(3, 0)] <= t[3])
    {
        if (risk_table_f[3][deg_elp[DEG_ELP_IDX(3, 0)]].exponent <= -8)
        {
            n_candidates++;
            start = 3;
        }
    }

    if (n_candidates > 1)
    {
        /* decide on order if mode 2 and 3 are considered */
        if (simple_float_cmp(risk_table_f[2][deg_elp[DEG_ELP_IDX(2, 0)]], risk_table_f[3][deg_elp[DEG_ELP_IDX(3, 0)]]) <
            0)
        {
            start = 2;
            inc   = 1;
        }
        else
        {
            start = 3;
            inc   = -1;
        }
    }

    for (mode_counter = start, i = 0; i < n_candidates; mode_counter += inc, i++)
    {
        if (risk_table_f[mode_counter][deg_elp[DEG_ELP_IDX(mode_counter, 0)]].exponent < epmr_lowest_risk_exp)
        {
            if (!rs16_factorize_elp(err_pos + ERR_POS_IDX(mode_counter, 0), elp + ELP_IDX(mode_counter, 0),
                                    deg_elp[DEG_ELP_IDX(mode_counter, 0)], first_codeword_length - 1))
            {
                /* code word is decodable with error correction */
                epmr_lowest_risk_exp = risk_table_f[mode_counter][deg_elp[DEG_ELP_IDX(mode_counter, 0)]].exponent;

                rs16_calculate_errors(err_symb + ERR_SYMB_IDX(mode_counter, 0), err_pos + ERR_POS_IDX(mode_counter, 0),
                                      syndromes + SYNDROME_IDX(mode_counter, 0), deg_elp[DEG_ELP_IDX(mode_counter, 0)],
                                      t[mode_counter]);

                for (i = 0; i < deg_elp[DEG_ELP_IDX(mode_counter, 0)]; i++)
                {
                    cw0[err_pos[ERR_POS_IDX(mode_counter, 0) + i]] = GF16_ADD(
                        cw0[err_pos[ERR_POS_IDX(mode_counter, 0) + i]], err_symb[ERR_SYMB_IDX(mode_counter, 0) + i]);
                }
                break;
            }
        }
    }

    epmr = cw0_get_epmr(cw0, first_codeword_length - 1);

    if (epmr_lowest_risk_exp > -16)
    {
        epmr += 4;
    }
    if (epmr_lowest_risk_exp > -8)
    {
        epmr += 4;
    }

     
    return epmr;
}

FEC_STATIC LC3_INT32 rs16_detect_and_correct(LC3_UINT8 *iobuf, LC3_INT32 n_symb, LC3_INT32 n_codewords, LC3PLUS_EpModeRequest *epmr, LC3_INT16 *error_report,
                                       LC3_INT32 *bfi, LC3_UINT8 *array_of_trust, LC3_INT32 ccc_flag, LC3_INT16 *n_pccw)
{

        LC3_INT16        mode_broken[4];
        LC3_INT16        error_report_ep_ok[4];
        LC3_INT16       i, cw_counter, mode_counter, cw_offset;
        LC3_INT16       codeword_length;
        LC3_INT16       mode;
        LC3_INT16       mode_candidates[4];
        LC3_INT16       n_mode_candidates;
        LC3_INT16       broken_cw, n_broken_cw;
        LC3_INT16       j, idx_min;
        LC3_INT16       n_pccw0;
        simple_float val_min_f;
        LC3_INT16       tmp;
        LC3_INT16       epmr_position;
        simple_float dec_risk_f[FEC_N_MODES];
        simple_float risk_min_f;
        simple_float ep_risk_thresh;
        LC3_INT32 epmr_dec_fail_increment;
        LC3_UINT8 const *hamming_distance;
        LC3_UINT8        syndromes[FEC_TOTAL_SYNDROME_SIZE];
        LC3_UINT8        elp[FEC_TOTAL_ELP_SIZE];
        LC3_UINT8        err_pos[FEC_TOTAL_ERR_POS_SIZE];
        LC3_UINT8        err_symb[FEC_TOTAL_ERROR_SIZE];
        LC3_INT8         t[FEC_N_MODES];
        LC3_INT8         deg_elp[FEC_TOTAL_DEG_ELP_SIZE];
        LC3_UINT8        blacklist[FEC_N_MODES];
        LC3_INT32          rop;

    void (*syndr_calc[3])(LC3_UINT8 *, LC3_UINT8 *, LC3_INT32);
    rop = 0;

    /* initialization */
    blacklist[0]        = 0;
    blacklist[1]        = 0;
    blacklist[2]        = 0;
    blacklist[3]        = 0;
    mode_broken[0]        = 0;
    mode_broken[1]        = 0;
    mode_broken[2]        = 0;
    mode_broken[3]        = 0;
    error_report_ep_ok[0] = ERROR_REPORT_EP1_OK;
    error_report_ep_ok[1] = ERROR_REPORT_EP2_OK;
    error_report_ep_ok[2] = ERROR_REPORT_EP3_OK;
    error_report_ep_ok[3] = ERROR_REPORT_EP4_OK;
    hamming_distance    = &hamming_distance_by_mode0[1];
    mode                = -1;
    n_mode_candidates   = 0;
    risk_min_f.mantissa = SIMPLE_FLOAT_1_MANTISSA;
    risk_min_f.exponent = 0;
    
    if (n_symb <= 80)
    {
    ep_risk_thresh.mantissa = EP_RISK_THRESH_NS_M;
    ep_risk_thresh.exponent = EP_RISK_THRESH_NS_E;
    }
    else
    {
    ep_risk_thresh.mantissa = EP_RISK_THRESH_OS_M;
    ep_risk_thresh.exponent = EP_RISK_THRESH_OS_E;
    }
    
    syndr_calc[0] = &rs16_calculate_two_syndromes;
    syndr_calc[1] = &rs16_calculate_four_syndromes;
    syndr_calc[2] = &rs16_calculate_six_syndromes;
    
    for (i = 0; i < FEC_N_MODES; i++)
    {
    t[i] = (hamming_distance[i] -1)/2;
    }
    
    *error_report = 0;
    *bfi          = 0;
    
    /* mode detection (stage 1) */
    codeword_length = get_codeword_length(n_codewords, n_symb, 0);
    
    epmr_position = codeword_length - 1;
    
    rs16_calculate_two_syndromes(syndromes + SYNDROME_IDX(0, 0), iobuf, codeword_length - 1);
    
    if ((syndromes[0 + SYNDROME_IDX(0, 0)] | syndromes[1 + SYNDROME_IDX(0, 0)]) == 0)
    {
    
    /* data validation for fec mode 1 */
    *epmr = cw0_get_epmr(iobuf, epmr_position);
    
    dw0_bitswap(iobuf + 2, 1, n_symb / 2);
    
    if (!crc1(iobuf + 8, n_symb - 8, *epmr, iobuf + 2, 3, 1))
    {
            *error_report |= ERROR_REPORT_ALL_OK;
        mode = 0;
        
         
        rop = mode + 1;
        goto CLEANUP;
    }
    else
    {
        /* reverse bit swap */
        dw0_bitswap(iobuf + 2, 1, n_symb / 2);
        
        *epmr += 4;
    }
    }
    
    blacklist[0] = 1;
    
    /* mode detection (stage 2) */
    
    /* calculate syndromes of code words 0 to 5 and modes 1 to 3 */
    cw_offset = 0;
    
    for (cw_counter = 0; cw_counter < 6; cw_counter++)
    {
    codeword_length = get_codeword_length(n_codewords, n_symb, cw_counter);
    
    rs16_calculate_six_syndromes(syndromes + SYNDROME_IDX(1, cw_counter), iobuf + cw_offset,
                     codeword_length - 1);
    
    cw_offset += codeword_length;
    
    for (mode_counter = FEC_N_MODES - 1; mode_counter >= 1; mode_counter--)
    {
        for (i = 0; i < hamming_distance[mode_counter] - 1; i++)
        {
        syndromes[SYNDROME_IDX(mode_counter, cw_counter) + i] = GF16_ADD(
            syndromes[SYNDROME_IDX(1, cw_counter) + i], sig_poly_syndr[mode_counter][i]);
        }
    }
    }

    /* check for valid code words */
    for (mode_counter = 1; mode_counter < FEC_N_MODES; mode_counter++)
    {
    n_broken_cw = 0;
    for (cw_counter = 0; cw_counter < 6; cw_counter++)
    {
        broken_cw = 0;
        for (i = 0; i < hamming_distance[mode_counter] - 1; i++)
        {
        broken_cw |= syndromes[SYNDROME_IDX(mode_counter, cw_counter) + i];
        }
        if (broken_cw != 0)
        {
        n_broken_cw ++;
        }
    }
    
    if (n_broken_cw == 0)
    {
        mode      = mode_counter;
        cw_offset = 0;
        
        *epmr = cw0_get_epmr(iobuf, epmr_position);
        
        for (cw_counter = 0; cw_counter < 6; cw_counter++)
        {
        codeword_length = get_codeword_length(n_codewords, n_symb, cw_counter);
        for (i = 0; i <= EP_SIG_POLY_DEG; i++)
        {
            iobuf[cw_offset + i] = GF16_ADD(iobuf[cw_offset + i], sig_polys[mode][i]);
        }
        cw_offset += codeword_length;
        }
    }
    }
    
    if (mode < 0) /* mode hasn't been detected so far -> errors occurred in transmission */
    {
    /* calculate error locator polynomials for code words 0 to 5 */
    for (mode_counter = 1; mode_counter < FEC_N_MODES; mode_counter++)
    {
        for (cw_counter = 0; cw_counter < 6; cw_counter++)
        {
        deg_elp[DEG_ELP_IDX(mode_counter, cw_counter)] = rs16_calculate_elp(
            elp + ELP_IDX(mode_counter, cw_counter), syndromes + SYNDROME_IDX(mode_counter, cw_counter),
            t[mode_counter]);
        if (deg_elp[DEG_ELP_IDX(mode_counter, cw_counter)] > t[mode_counter])
        {
            blacklist[mode_counter] = 1;
            break;
        }
        }
    }
    
    /* risk analysis for mode candidate selection */
    for (mode_counter = 1; mode_counter < FEC_N_MODES; mode_counter++)
    {
        dec_risk_f[mode_counter].mantissa = SIMPLE_FLOAT_1_MANTISSA;
        dec_risk_f[mode_counter].exponent = 0;
        
        if (blacklist[mode_counter] == 0)
        {
        for (cw_counter = 0; cw_counter < 6; cw_counter++)
        {
            dec_risk_f[mode_counter] = simple_float_mul(
            dec_risk_f[mode_counter],
            risk_table_f[mode_counter][deg_elp[DEG_ELP_IDX(mode_counter, cw_counter)]]);
        }
        
        if (simple_float_cmp(dec_risk_f[mode_counter], ep_risk_thresh) <= 0)
        {
            mode_candidates[n_mode_candidates++] = mode_counter;
        }
        
        if (simple_float_cmp(dec_risk_f[mode_counter], risk_min_f) < 0)
        {
            risk_min_f = dec_risk_f[mode_counter];
        }
        }
    }
    assert(n_mode_candidates <= 4); // suppress false gcc warning when OPTIM=3
    
    /* sort mode candidates by risk */
    for (i = 0; i < n_mode_candidates; i++)
    {
        idx_min   = i;
        val_min_f = dec_risk_f[mode_candidates[i]];
        
        for (j = i + 1; j < n_mode_candidates; j++)
        {
        if (simple_float_cmp(dec_risk_f[mode_candidates[j]], val_min_f) < 0)
        {
            val_min_f = dec_risk_f[mode_candidates[j]];
            idx_min   = j;
        }
        }
        
        if (idx_min > i)
        {    
        tmp                      = mode_candidates[i];
        mode_candidates[i]       = mode_candidates[idx_min];
        mode_candidates[idx_min] = tmp;
        }
    }
    
    /* try out candidate modes */
    for (i = 0; i < n_mode_candidates; i++)
    {
        mode = mode_candidates[i];
        
        for (cw_counter = 0; cw_counter < 6; cw_counter++)
        {
        codeword_length = get_codeword_length(n_codewords, n_symb, cw_counter);
        
        if (deg_elp[DEG_ELP_IDX(mode, cw_counter)])
        {
            if (rs16_factorize_elp(err_pos + ERR_POS_IDX(mode, cw_counter), elp + ELP_IDX(mode, cw_counter),
                       deg_elp[DEG_ELP_IDX(mode, cw_counter)], codeword_length - 1))
            {
            /* elp did not split into distinct linear factors or error position was out of range */
            mode = -1;
            break;
            }
        }
        }
        if (mode > 0)
        {
        /* decodable mode with lowest risk has been found */
        break;
        }
    }
    
    if (mode < 0)
    {
        /* no decodable mode has been found */
        *error_report = ERROR_REPORT_BEC_MASK;
        *bfi          = 1;
        mode          = -1;
        
        *epmr = fec_estimate_epmr_from_cw0(iobuf, t, syndromes, elp, deg_elp, err_pos, err_symb, n_codewords,
                           n_symb);
        
         
        rop = mode;
        goto CLEANUP;
    }
    
    /* perform error correction */
    cw_offset     = 0;
    *error_report = 0;
    for (cw_counter = 0; cw_counter < 6; cw_counter++)
    {
        codeword_length = get_codeword_length(n_codewords, n_symb, cw_counter);
        
        if (deg_elp[DEG_ELP_IDX(mode, cw_counter)])
        {
        rs16_calculate_errors(
            err_symb + ERR_SYMB_IDX(mode, cw_counter), err_pos + ERR_POS_IDX(mode, cw_counter),
            syndromes + SYNDROME_IDX(mode, cw_counter), deg_elp[DEG_ELP_IDX(mode, cw_counter)], t[mode]);
        
        /* correct errors and sum up number of corrected bits */
        for (i = 0; i < deg_elp[DEG_ELP_IDX(mode, cw_counter)]; i++)
        {
            iobuf[err_pos[ERR_POS_IDX(mode, cw_counter) + i] + cw_offset] =
            GF16_ADD(iobuf[err_pos[ERR_POS_IDX(mode, cw_counter) + i] + cw_offset],
                 err_symb[ERR_SYMB_IDX(mode, cw_counter) + i]);
            *error_report += rs16_bit_count_table[err_symb[ERR_SYMB_IDX(mode, cw_counter) + i]];
        }

                for (i = 0; i < mode; i ++)
                {
                    if(deg_elp[DEG_ELP_IDX(mode, cw_counter)] > i)
                    {
                        mode_broken[i] = 1;
                    }
                }

        }
        
        for (i = 0; i <= EP_SIG_POLY_DEG; i++)
        {
        iobuf[cw_offset + i] = GF16_ADD(iobuf[cw_offset + i], sig_polys[mode][i]);
        }
        cw_offset += codeword_length;
    }
    
    /* set epmr according to risk value of cw0 */
    epmr_dec_fail_increment = 8;
    
    if (risk_table_f[mode][deg_elp[DEG_ELP_IDX(mode, 0)]].exponent <= -8)
    {
        epmr_dec_fail_increment -= 4;
    }
    if (risk_table_f[mode][deg_elp[DEG_ELP_IDX(mode, 0)]].exponent <= -16)
    {
        epmr_dec_fail_increment -= 4;
    }
    
    *epmr = (LC3PLUS_EpModeRequest)(cw0_get_epmr(iobuf, epmr_position) + epmr_dec_fail_increment);
    }
    
    /* mode has been successfully detected -> now check and try to correct remaining code words*/
    *n_pccw = fec_get_n_pccw(n_symb / 2, mode + 1, ccc_flag);
    if (ccc_flag == 0)
    {
    n_pccw0 = fec_get_n_pccw(n_symb / 2, mode + 1, ccc_flag);
    *n_pccw = n_pccw0;
    }
    else
    {
    n_pccw0 = 0;
    }
    
    for (cw_counter = 6; cw_counter < n_codewords; cw_counter++)
    {
    /* usual error correction scheme: syndromes -> elp's, errors, etc. */
    codeword_length                              = get_codeword_length(n_codewords, n_symb, cw_counter);
    array_of_trust[n_codewords - 1 - cw_counter] = 1;
    
    syndr_calc[t[mode] - 1](syndromes, iobuf + cw_offset, codeword_length -1);
    
    deg_elp[0] = rs16_calculate_elp(elp, syndromes, t[mode]);
    
        for (i = 0; i < mode; i ++)
        {
            if (deg_elp[0] > i)
            {
                mode_broken[i] = 1;
            }
        }

    if (deg_elp[0] > t[mode])
    {
            for (i = 0; i < 4; i ++)
            {
                mode_broken[i] = 1;
            }
        cw_offset += codeword_length;
        if (cw_counter < n_codewords - n_pccw0)
        {
            *error_report = ERROR_REPORT_BEC_MASK;
            mode          = -1;
            *bfi          = 1;
          break;
        }
        else
        {
        *bfi                                         = 2;
        array_of_trust[n_codewords - 1 - cw_counter] = 0;
        continue;
        }
    }
    
    if (deg_elp[0])
    {
        if (rs16_factorize_elp(err_pos, elp, deg_elp[0], codeword_length - 1))
        {
            cw_offset += codeword_length;
            for (i = 0; i < 4; i ++)
            {
                mode_broken[i] = 1;
            }
            if (cw_counter < n_codewords - n_pccw0)
            {
                *error_report = ERROR_REPORT_BEC_MASK;
                mode = -1;
                *bfi = 1;
                
                break;
            }
            else
            {
                *bfi                                         = 2;
                array_of_trust[n_codewords - 1 - cw_counter] = 0;
                continue;
            }
        }
        
        rs16_calculate_errors(err_symb, err_pos, syndromes, deg_elp[0], t[mode]);
        
        /* correct errors and sum up number of corrected bits */
        for (i = 0; i < deg_elp[0]; i++)
        {
        iobuf[err_pos[i] + cw_offset] = GF16_ADD(iobuf[err_pos[i] + cw_offset], err_symb[i]);
        *error_report                 += rs16_bit_count_table[err_symb[i]];
        }
    }
    cw_offset += codeword_length;
    if (risk_table_f[mode][deg_elp[0]].exponent > -16)
    {
        array_of_trust[n_codewords - 1 - cw_counter] = 0;
    }
    }
    
    *error_report &= ERROR_REPORT_BEC_MASK;
    for (i = 0; i < 4; i ++)
    {
        if (!mode_broken[i])
        {
            *error_report |= error_report_ep_ok[i];
        }
    }

    if (mode >= 0)
    {
        rop = mode + 1;
    } else {
        rop = -1;
    }
    
     

CLEANUP:
    return rop;
}

FEC_STATIC void rs16_calculate_six_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg)
{
        LC3_INT32    i;
        LC3_UINT8 buffer[15];

    assert(cw_poly_deg >= 12);

    for (i = 0; i <= cw_poly_deg; i++)
    {
        buffer[i] = cw[i];
    }

    syndromes[0] = buffer[0];
    syndromes[1] = buffer[0];
    syndromes[2] = buffer[0];
    syndromes[3] = buffer[0];
    syndromes[4] = buffer[0];
    syndromes[5] = buffer[0];

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[1], 32));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[1], 64));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[1], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[1], 48));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[1], 96));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[1], 192));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[2], 64));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[2], 48));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[2], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[2], 80));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[2], 112));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[2], 240));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[3], 128));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[3], 192));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[3], 160));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[3], 240));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[3], 16));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[3], 128));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[4], 48));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[4], 80));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[4], 240));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[4], 32));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[4], 96));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[4], 160));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[5], 96));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[5], 112));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[5], 16));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[5], 96));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[5], 112));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[5], 16));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[6], 192));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[6], 240));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[6], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[6], 160));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[6], 16));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[6], 192));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[7], 176));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[7], 144));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[7], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[7], 208));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[7], 96));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[7], 240));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[8], 80));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[8], 32));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[8], 160));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[8], 64));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[8], 112));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[8], 128));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[9], 160));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[9], 128));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[9], 240));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[9], 192));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[9], 16));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[9], 160));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[10], 112));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[10], 96));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[10], 16));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[10], 112));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[10], 96));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[10], 16));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[11], 224));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[11], 176));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[11], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[11], 144));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[11], 112));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[11], 192));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[12], 240));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[12], 160));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[12], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[12], 128));
    syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[12], 16));
    syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[12], 240));

    if (cw_poly_deg >= 13)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[13], 208));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[13], 224));
        syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[13], 160));
        syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[13], 176));
        syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[13], 96));
        syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[13], 128));
    }

    if (cw_poly_deg >= 14)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[14], 144));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[14], 208));
        syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[14], 240));
        syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[14], 224));
        syndromes[4] = GF16_ADD(syndromes[4], GF16_MUL0(buffer[14], 112));
        syndromes[5] = GF16_ADD(syndromes[5], GF16_MUL0(buffer[14], 160));
    }

     
}

FEC_STATIC void rs16_calculate_four_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg)
{
        LC3_INT32    i;
        LC3_UINT8 buffer[15];

    assert(cw_poly_deg >= 12);

    for (i = 0; i <= cw_poly_deg; i++)
    {
        buffer[i] = cw[i];
    }

    syndromes[0] = buffer[0];
    syndromes[1] = buffer[0];
    syndromes[2] = buffer[0];
    syndromes[3] = buffer[0];

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[1], 32));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[1], 64));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[1], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[1], 48));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[2], 64));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[2], 48));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[2], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[2], 80));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[3], 128));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[3], 192));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[3], 160));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[3], 240));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[4], 48));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[4], 80));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[4], 240));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[4], 32));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[5], 96));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[5], 112));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[5], 16));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[5], 96));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[6], 192));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[6], 240));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[6], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[6], 160));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[7], 176));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[7], 144));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[7], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[7], 208));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[8], 80));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[8], 32));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[8], 160));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[8], 64));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[9], 160));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[9], 128));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[9], 240));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[9], 192));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[10], 112));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[10], 96));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[10], 16));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[10], 112));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[11], 224));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[11], 176));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[11], 128));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[11], 144));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[12], 240));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[12], 160));
    syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[12], 192));
    syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[12], 128));

    if (cw_poly_deg >= 13)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[13], 208));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[13], 224));
        syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[13], 160));
        syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[13], 176));
    }

    if (cw_poly_deg >= 14)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[14], 144));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[14], 208));
        syndromes[2] = GF16_ADD(syndromes[2], GF16_MUL0(buffer[14], 240));
        syndromes[3] = GF16_ADD(syndromes[3], GF16_MUL0(buffer[14], 224));
    }

     
}

FEC_STATIC void rs16_calculate_two_syndromes(LC3_UINT8 *syndromes, LC3_UINT8 *cw, LC3_INT32 cw_poly_deg)
{
        LC3_INT32    i;
        LC3_UINT8 buffer[15];

    assert(cw_poly_deg >= 12);

    for (i = 0; i <= cw_poly_deg; i++)
    {
        buffer[i] = cw[i];
    }

    syndromes[0] = buffer[0];
    syndromes[1] = buffer[0];

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[1], 32));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[1], 64));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[2], 64));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[2], 48));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[3], 128));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[3], 192));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[4], 48));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[4], 80));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[5], 96));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[5], 112));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[6], 192));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[6], 240));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[7], 176));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[7], 144));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[8], 80));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[8], 32));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[9], 160));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[9], 128));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[10], 112));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[10], 96));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[11], 224));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[11], 176));

    syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[12], 240));
    syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[12], 160));

    if (cw_poly_deg >= 13)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[13], 208));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[13], 224));
    }

    if (cw_poly_deg >= 14)
    {
        syndromes[0] = GF16_ADD(syndromes[0], GF16_MUL0(buffer[14], 144));
        syndromes[1] = GF16_ADD(syndromes[1], GF16_MUL0(buffer[14], 208));
    }

     
}

FEC_STATIC LC3_INT8 rs16_calculate_elp(LC3_UINT8 *elp, LC3_UINT8 *syndromes, LC3_INT16 t)
/* calculates error locator polynomial vie Petterson's algorithm */
{
        LC3_INT8  ret;
        LC3_UINT8 det, det_inv, aux, all_s, *s;
        LC3_UINT8 s22, s33, s44, s13, s14, s15;
        LC3_UINT8 s23, s24, s25, s34, s35;
        LC3_UINT8 a, b, c, d, e, f;

    ret    = 0;
    all_s  = 0;
    s      = syndromes;
    elp[0] = 1;
    memset(elp + 1, 0, 3);

    switch (t)
    {
    case 3:
    {
        /* check for errors */
        all_s = s[0] | s[1] | s[2] | s[3] | s[4] | s[5];

        if (all_s == 0)
        {
            break;
        }

        /* assume 3 errors */
        s22 = GF16_MUL(s[1], s[1]);
        s33 = GF16_MUL(s[2], s[2]);
        s44 = GF16_MUL(s[3], s[3]);
        s13 = GF16_MUL(s[0], s[2]);

        det = GF16_ADD(GF16_ADD(GF16_MUL(s13, s[4]), GF16_MUL(s44, s[0])),
                       GF16_ADD(GF16_MUL(s22, s[4]), GF16_MUL(s33, s[2])));

        if (det)
        {
            det_inv = gf16_inv_table[det] << 4;

            s14 = GF16_MUL(s[0], s[3]);
            s15 = GF16_MUL(s[0], s[4]);

            s23 = GF16_MUL(s[1], s[2]);
            s24 = GF16_MUL(s[1], s[3]);
            s25 = GF16_MUL(s[1], s[4]);

            s34 = GF16_MUL(s[2], s[3]);
            s35 = GF16_MUL(s[2], s[4]);

            a = GF16_ADD(s35, s44) << 4;
            b = GF16_ADD(s15, s33) << 4;
            c = GF16_ADD(s13, s22) << 4;
            d = GF16_ADD(s34, s25) << 4;
            e = GF16_ADD(s23, s14) << 4;
            f = GF16_ADD(s24, s33) << 4;

            aux    = GF16_ADD(GF16_ADD(GF16_MUL0(a, s[3]), GF16_MUL0(d, s[4])), GF16_MUL0(f, s[5]));
            elp[3] = GF16_MUL0(aux, det_inv);

            aux    = GF16_ADD(GF16_ADD(GF16_MUL0(d, s[3]), GF16_MUL0(b, s[4])), GF16_MUL0(e, s[5]));
            elp[2] = GF16_MUL0(aux, det_inv);

            aux    = GF16_ADD(GF16_ADD(GF16_MUL0(f, s[3]), GF16_MUL0(e, s[4])), GF16_MUL0(c, s[5]));
            elp[1] = GF16_MUL0(aux, det_inv);

            if (elp[3] == 0)
            {
                ret = t+1;
            }
            else
            {
                ret = 3;
            }
            break;
        }

        /* assume two errors */
        det = GF16_ADD(GF16_MUL(syndromes[0], syndromes[2]), GF16_MUL(syndromes[1], syndromes[1]));

        if (det)
        {
            det_inv = gf16_inv_table[det] << 4;

            aux    = GF16_ADD(GF16_MUL(syndromes[1], syndromes[2]), GF16_MUL(syndromes[0], syndromes[3]));
            elp[1] = GF16_MUL0(aux, det_inv);

            aux    = GF16_ADD(GF16_MUL(syndromes[2], syndromes[2]), GF16_MUL(syndromes[1], syndromes[3]));
            elp[2] = GF16_MUL0(aux, det_inv);

            /* check remaining LSF relations */
            aux = GF16_ADD(GF16_ADD(GF16_MUL(elp[2], s[2]), GF16_MUL(elp[1], s[3])), s[4])
                    | GF16_ADD(GF16_ADD(GF16_MUL(elp[2], s[3]), GF16_MUL(elp[1], s[4])), s[5]);

            aux |= elp[2] == 0;

            if (aux != 0)
            {
                ret = t + 1;
            }
            else
            {
                ret = 2;
            }
            break;
        }

        /* assume one error */
        if (syndromes[0] != 0)
        {
            elp[1] = GF16_MUL(syndromes[1], gf16_inv_table[syndromes[0]]);

            /* check remaining LSF relations */
            aux = GF16_ADD(GF16_MUL(elp[1], s[1]), s[2])
                    | GF16_ADD(GF16_MUL(elp[1], s[2]), s[3])
                    | GF16_ADD(GF16_MUL(elp[1], s[3]), s[4])
                    | GF16_ADD(GF16_MUL(elp[1], s[4]), s[5]);

            aux |= elp[1] == 0;

            if (aux != 0)
            {
                ret = t + 1;
            }
            else
            {
                ret = 1;
            }
            break;
        }

        ret = t + 1;
        break;
    }
    case 2:
    {
        all_s = s[0] | s[1] | s[2] | s[3];

        if (all_s == 0)
        {
            break;
        }

        /* assume two errors */
        det = GF16_ADD(GF16_MUL(syndromes[0], syndromes[2]), GF16_MUL(syndromes[1], syndromes[1]));

        if (det)
        {
            det_inv = gf16_inv_table[det] << 4;

            aux    = GF16_ADD(GF16_MUL(syndromes[1], syndromes[2]), GF16_MUL(syndromes[0], syndromes[3]));
            elp[1] = GF16_MUL0(aux, det_inv);

            aux    = GF16_ADD(GF16_MUL(syndromes[2], syndromes[2]), GF16_MUL(syndromes[1], syndromes[3]));
            elp[2] = GF16_MUL0(aux, det_inv);

            if (elp[2] == 0)
            {
                ret = t + 1;
            }
            else
            {
                ret = 2;
            }
            break;
        }

        /* assume one error */
        if (syndromes[0] != 0)
        {
            elp[1] = GF16_MUL(syndromes[1], gf16_inv_table[syndromes[0]]);

            /* check remaining LSF relation */
            aux = GF16_ADD(GF16_MUL(elp[1], s[1]), s[2]) | GF16_ADD(GF16_MUL(elp[1], s[2]), s[3]);
            aux |= elp[1] == 0;

            if (aux != 0)
            {
                ret = t + 1;
            }
            else
            {
                ret = 1;
            }
            break;
        }

        ret = t + 1;
        break;
    }
    case 1:
    {
        all_s = s[0] | s[1];

        if (all_s == 0)
        {
            break;
        }

        if (syndromes[0] != 0)
        {
            elp[1] = GF16_MUL(syndromes[1], gf16_inv_table[syndromes[0]]);
            if (elp[1] == 0)
            {
                ret = t + 1;
            }
            else
            {
                ret = 1;
            }
            break;
        }

        ret = t + 1;
        break;
    }
    default: assert(0 && "calculating elp of this degree not implemented");
    }

     
    return ret;
}

FEC_STATIC LC3_INT16 rs16_factorize_elp(LC3_UINT8 *err_pos, LC3_UINT8 *elp, LC3_INT16 deg_elp, LC3_INT16 max_pos)
{
        LC3_UINT8 beta, gamma;
        LC3_INT16 zeros, err_pos0, err_pos1, err_pos2, ret;

    beta  = 0;
    gamma = 0;
    zeros = 0;
    ret   = 0;

    switch (deg_elp)
    {
    case 0: break;

    case 1:
        err_pos0 = gf16_log_g[elp[1]];
        if (err_pos0 > max_pos)
        {
            ret = 1;
            break;
        }

        err_pos[0] = (LC3_UINT8)err_pos0;
        break;

    case 2:
        zeros = rs16_elp_deg2_table[elp[1] | (elp[2] << 4)];
        if (zeros == 0)
        {  
             
            return 1;
        }

        err_pos0 = zeros & 15;
        err_pos1 = (zeros >> 4) & 15;

        if (err_pos0 > max_pos || err_pos1 > max_pos)
        {
            ret = 1;
            break;
        }

        err_pos[0] = (LC3_UINT8)err_pos0;
        err_pos[1] = (LC3_UINT8)err_pos1;
        break;

    case 3:
        /* beta = a*a + b, gamma = a*b + c */
        beta  = GF16_ADD(GF16_MUL(elp[1], elp[1]), elp[2]);
        gamma = GF16_ADD(GF16_MUL(elp[1], elp[2]), elp[3]);
        zeros = rs16_elp_deg3_table[beta | gamma << 4];

        if (zeros == 0)
        /* elp does not split over GF(16) or has multiple zeros */
        {
            ret = 1;
            break;
        }

        /* remove shift from zeros */
        err_pos0 = GF16_ADD(zeros & 15, elp[1]);
        err_pos1 = GF16_ADD((zeros >> 4) & 15, elp[1]);
        err_pos2 = GF16_ADD((zeros >> 8) & 15, elp[1]);

        if (err_pos0 == 0 || err_pos1 == 0 || err_pos2 == 0)
        {
             
            return 1;
        }

        err_pos0 = gf16_log_g[err_pos0];
        err_pos1 = gf16_log_g[err_pos1];
        err_pos2 = gf16_log_g[err_pos2];

        if (err_pos0 > max_pos || err_pos1 > max_pos || err_pos2 > max_pos)
        {
            ret = 1;
            break;
        }

        err_pos[0] = (LC3_UINT8)err_pos0;
        err_pos[1] = (LC3_UINT8)err_pos1;
        err_pos[2] = (LC3_UINT8)err_pos2;

        break;

    default: assert(0 && "invalid degree in rs16_error_locator");
    }

     
    return ret;
}

FEC_STATIC void rs16_calculate_errors(LC3_UINT8 *err_symb, LC3_UINT8 *err_pos, LC3_UINT8 *syndromes, LC3_INT8 deg_elp, LC3_INT8 t)
{
        LC3_UINT8 det_inv;
        LC3_UINT8 x0, x1, x2;
        LC3_UINT8 x0sq, x1sq, x2sq;
        LC3_UINT8 c0, c1, c2;
        LC3_UINT8 s0, s1, s2;
        LC3_UINT8 tmp;
    
    UNUSED(t);

    assert(deg_elp <= t);

    switch (deg_elp)
    {
    case 0: break;

    case 1:
        err_symb[0] = GF16_MUL(gf16_g_pow[15 - err_pos[0]], syndromes[0]);

        break;

    case 2:
        s0 = (LC3_UINT8) (syndromes[0] << 4);
        s1 = (LC3_UINT8) (syndromes[1] << 4);

        x0 = gf16_g_pow[err_pos[0]];
        x1 = gf16_g_pow[err_pos[1]];

        x0sq = GF16_MUL(x0, x0);
        x1sq = GF16_MUL(x1, x1);

        tmp     = GF16_ADD(GF16_MUL(x0sq, x1), GF16_MUL(x1sq, x0));
        det_inv = gf16_inv_table[tmp] << 4;

        tmp         = GF16_ADD(GF16_MUL0(x1sq, s0), GF16_MUL0(x1, s1));
        err_symb[0] = GF16_MUL0(tmp, det_inv);

        tmp         = GF16_ADD(GF16_MUL0(x0sq, s0), GF16_MUL0(x0, s1));
        err_symb[1] = GF16_MUL0(tmp, det_inv);

        break;

    case 3:
        s0 = syndromes[0] << 4;
        s1 = syndromes[1] << 4;
        s2 = syndromes[2] << 4;

        x0 = gf16_g_pow[err_pos[0]];
        x1 = gf16_g_pow[err_pos[1]];
        x2 = gf16_g_pow[err_pos[2]];

        x0sq = GF16_MUL(x0, x0);
        x1sq = GF16_MUL(x1, x1);
        x2sq = GF16_MUL(x2, x2);

        tmp     = GF16_MUL(GF16_ADD(x1, x0), GF16_ADD(x2, x0));
        tmp     = GF16_MUL(GF16_ADD(x2, x1), tmp);
        det_inv = gf16_inv_table[tmp] << 4;

        c0 = GF16_ADD(GF16_MUL(x1, x2sq), GF16_MUL(x2, x1sq));
        c1 = GF16_ADD(x2sq, x1sq);
        c2 = GF16_ADD(x2, x1);

        err_symb[0] = GF16_ADD(GF16_ADD(GF16_MUL0(c0, s0), GF16_MUL0(c1, s1)), GF16_MUL0(c2, s2));

        c0 = GF16_ADD(GF16_MUL(x0, x2sq), GF16_MUL(x2, x0sq));
        c1 = GF16_ADD(x2sq, x0sq);
        c2 = GF16_ADD(x2, x0);

        err_symb[1] = GF16_ADD(GF16_ADD(GF16_MUL0(c0, s0), GF16_MUL0(c1, s1)), GF16_MUL0(c2, s2));

        c0 = GF16_ADD(GF16_MUL(x0, x1sq), GF16_MUL(x1, x0sq));
        c1 = GF16_ADD(x1sq, x0sq);
        c2 = GF16_ADD(x1, x0);

        err_symb[2] = GF16_ADD(GF16_ADD(GF16_MUL0(c0, s0), GF16_MUL0(c1, s1)), GF16_MUL0(c2, s2));

        tmp         = GF16_MUL0(err_symb[0], det_inv);
        err_symb[0] = GF16_MUL(tmp, gf16_inv_table[x0]);

        tmp         = GF16_MUL0(err_symb[1], det_inv);
        err_symb[1] = GF16_MUL(tmp, gf16_inv_table[x1]);

        tmp         = GF16_MUL0(err_symb[2], det_inv);
        err_symb[2] = GF16_MUL(tmp, gf16_inv_table[x2]);

        break;

    default: assert(0 && "method not implemented\n"); break;
    }

     
}

/* hash functions for data validation */

/* hamming distance 4 */
static const LC3_UINT32 crc14_mask[16] = {0,      17989,  35978,  51919,  71956,  89937,  103838, 119771,
                                       143912, 160877, 179874, 194791, 207676, 224633, 239542, 254451};

/* hamming distance 4 */
static const LC3_UINT32 crc22_mask[16] = {0,        4788009,  9576018,  14356859, 19152036, 23933837, 28713718, 33500639,
                                       33650273, 38304072, 43214899, 47867674, 52775621, 57427436, 62346391, 67001278};

FEC_STATIC LC3_INT16 crc1(LC3_UINT8 *data, LC3_INT16 data_size, LC3_INT16 epmr, LC3_UINT8 *hash_val, LC3_INT16 hash_size, LC3_INT16 check)
{
        LC3_UINT32 const *mask;
        LC3_INT32           shift, i, fail;
        LC3_UINT32        rem;

    fail = 0;
    rem  = 0;

    assert(hash_size > 0);

    switch (hash_size)
    {
    case 2:
        shift = 14;
        mask  = crc14_mask;
        break;
    case 3:
        shift = 22;
        mask  = crc22_mask;
        break;
    default:
        shift = 0;
        mask = 0;
        assert(0 && "crc hash size not implemented");
    }

    /* data array contains 4-bit words */
    for (i = data_size - 1; i >= 0; i--)
    {
        rem = (rem << 4) ^ data[i];
        rem ^= mask[(rem >> shift) & 15];
    }

    rem = (rem << 4) ^ (epmr << 2);
    rem ^= mask[(rem >> shift) & 15];

    for (i = 0; i < 2 * hash_size - 1; i++)
    {
        rem <<= 4;
        rem ^= mask[(rem >> shift) & 15];
    }

    rem ^= ((LC3_UINT32) epmr) << shift;

    if (check)
    {
        /* test hash value */
        for (i = 0; i < 2 * hash_size; i++)
        {
            fail |= hash_val[i] ^ ((rem >> (4*i)) & 15);
        }
    }
    else
    {
        /* write hash value */
        for (i = 0; i < 2 * hash_size; i++)
        {
            hash_val[i] = (LC3_UINT8) ((rem >> (4*i)) & 15);
        }
    }

     
    return fail;
}

/* hamming distance = 4 */
static const LC3_UINT32 crc16_mask[16] = {0,      107243, 190269, 214486, 289937, 380538, 428972, 469319,
                                       579874, 621513, 671263, 761076, 832947, 857944, 938638, 1044581};

FEC_STATIC LC3_INT16 crc2(LC3_UINT8 *data, LC3_INT16 data_size, LC3_UINT8 *hash_val, LC3_INT16 hash_size, LC3_INT16 check)
{
        LC3_UINT32 const *mask;
        LC3_INT32           shift, i, fail;
        LC3_UINT32        rem;

    fail = 0;
    rem  = 0;

    assert(hash_size > 0);

    switch (hash_size)
    {
    case 2:
        shift = 16;
        mask  = crc16_mask;
        break;
    default:
        shift = 0;
        mask = 0;
        assert(0 && "crc hash size not implemented");
    }

    /* data array contains 4-bit words */
    for (i = data_size - 1; i >= 0; i--)
    {
        rem = (rem << 4) ^ data[i];
        rem ^= mask[(rem >> shift) & 15];
    }

    for (i = 0; i < 2 * hash_size; i++)
    {
        rem <<= 4;
        rem ^= mask[(rem >> shift) & 15];
    }

    if (check)
    {
        /* test hash value */
        for (i = 0; i < 2 * hash_size; i++)
        {
            fail |= hash_val[i] ^ ((rem >> (4*i)) & 15);
        }
    }
    else
    {
        /* write hash value */
        for (i = 0; i < 2 * hash_size; i++)
        {
            hash_val[i] = (LC3_UINT8) ((rem >> (4*i)) & 15);
        }
    }

     
    return fail;
}

/* simple float implementation */
FEC_STATIC simple_float simple_float_mul(simple_float op1, simple_float op2)
{
        simple_float rop;
        LC3_INT32       aux;
    
    aux          = (op1.mantissa * op2.mantissa) >> 14;
    rop.exponent = op1.exponent + op2.exponent;
    if (aux & 32768L)
    {
        aux          >>= 1;
        rop.exponent ++;
    }
    rop.mantissa = (LC3_INT16) aux;

     
    return rop;
}

/* Auxiliary */

FEC_STATIC LC3_INT16 simple_float_cmp(simple_float op1, simple_float op2)
/* returns 1 if op1 > op2, 0 if op1 = op2, and -1 if op1 < op2 */
{
        LC3_INT16 rval;
        LC3_INT16 mdiff;
        LC3_INT16 ediff;

    rval = 0;

    ediff = op1.exponent - op2.exponent;
    mdiff = (LC3_INT16) op1.mantissa - (LC3_INT16) op2.mantissa;

    if (ediff == 0)
    {
        if (mdiff > 0)
        {
            rval = 1;
        }
        if (mdiff < 0)
        {
            rval = -1;
        }
    }
    else
    {
        if (ediff > 0)
        {
            rval = 1;
        }
        if (ediff < 0)
        {
            rval = -1;
        }
    }

     
    return rval;
}

