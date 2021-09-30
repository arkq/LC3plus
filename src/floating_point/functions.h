/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/
                                                                               

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "clib.h"
#include "defines.h"
#include "float.h"
#include "lc3.h"
#include "setup_dec_lc3.h"
#include "setup_enc_lc3.h"
#include "structs.h"
#include "util.h"

/* FFT */
#include "fft/iisfft.h"

/* fft.c */
void real_fft_init(Fft* fft, LC3_INT32 length, HANDLE_IIS_FFT *handle);
void real_ifft_init(Fft* fft, LC3_INT32 length, HANDLE_IIS_FFT *handle);
void real_fft_apply(Fft* fft, const LC3_FLOAT* in,  LC3_FLOAT* out);

void fft_init(Fft* fft, LC3_INT length);
void fft_free(Fft* fft);
void real_fft_free(Fft* fft);
void fft_apply(Fft* fft, const Complex* input, Complex* output);

/* dct.c */
void dct2_init(Dct2* dct, LC3_INT length);
void dct2_free(Dct2* dct);
void dct2_apply(Dct2* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

void dct3_init(Dct3* dct, LC3_INT length);
void dct3_free(Dct3* dct);
void dct3_apply(Dct3* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

void dct4_init(Dct4* dct, LC3_INT length);
void dct4_free(Dct4* dct);
void dct4_apply(Dct4* dct, const LC3_FLOAT* input, LC3_FLOAT* output);

/* mdct.c */
void mdct_init(Mdct* mdct, LC3_INT length, LC3_INT frame_dms, LC3_INT fs_idx, LC3_INT hrmode);
void mdct_free(Mdct* mdct);
void mdct_apply(const LC3_FLOAT* input, LC3_FLOAT* output, Mdct* mdct);

#ifdef ENABLE_PADDING
LC3_INT paddingDec_fl(LC3_UINT8* bytes, LC3_INT nbbits, LC3_INT L_spec, LC3_INT bw_cutoff_bits, LC3_INT ep_enabled, LC3_INT* total_padding, LC3_INT *np_zero);
#endif

void processEncoderEntropy_fl(LC3_UINT8* bytes, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT numbytes, LC3_INT bw_cutoff_bits,
                              LC3_INT bw_cutoff_idx, LC3_INT lastnz, LC3_INT N, LC3_INT lsbMode, LC3_INT gg_idx, LC3_INT num_tns_filters,
                              LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* scf_idx, LC3_INT fac_ns_idx
                              , LC3_INT bfi_ext, LC3_INT fs_idx
                              );
void processDecoderEntropy_fl(LC3_UINT8* bytes, LC3_INT numbytes, LC3_INT* mask_side, LC3_INT* bp_side, LC3_INT N, LC3_INT fs_idx,
                              LC3_INT bw_cutoff_bits, LC3_INT* bfi, LC3_INT* gg_idx, LC3_INT* scf_idx, LC3_INT* fac_ns_idx,
                              LC3_INT* tns_numfilters, LC3_INT* tns_order, LC3_INT* ltpf_idx, LC3_INT* bw_cutoff_idx, LC3_INT* lastnz,
                              LC3_INT* lsbMode, LC3_INT frame_dms
                              );
void processQuantizeSpec_fl(LC3_FLOAT x[], LC3_FLOAT gain, LC3_INT xq[], LC3_INT nt, LC3_INT totalBits, LC3_INT* nbits, LC3_INT* nbits2, LC3_INT fs,
                            LC3_INT* lastnzout, LC3_INT* codingdata, LC3_INT* lsbMode, LC3_INT mode, LC3_INT target, LC3_INT hrmode);

void processEstimateGlobalGain_fl(LC3_FLOAT x[], LC3_INT lg, LC3_INT nbitsSQ, LC3_FLOAT* gain, LC3_INT* quantizedGain,
                                  LC3_INT* quantizedGainMin, LC3_INT quantizedGainOff, LC3_FLOAT* targetBitsOff,
                                  LC3_INT* old_targetBits, LC3_INT old_specBits, LC3_INT bq_mode
                                  , LC3_INT regBits, LC3_FLOAT frame_ms
);

void processAriDecoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT L_spec, LC3_INT fs_idx, LC3_INT enable_lpc_weighting,
                          LC3_INT tns_numfilters, LC3_INT lsbMode, LC3_INT lastnz, LC3_INT* bfi, LC3_INT* tns_order, LC3_INT fac_ns_idx,
                          LC3_INT gg_idx, uint8_t* resBits, LC3_INT* x, LC3_INT* nf_seed, LC3_INT* tns_idx, LC3_INT* zero_frame, LC3_INT numbytes,
                          LC3_INT* nbits_residual, LC3_INT* residualPresent, LC3_INT frame_dms,
              LC3_INT32 n_pc, LC3_INT32 be_bp_left, LC3_INT32 be_bp_right, LC3_INT32 enc, LC3_INT32 *b_left, LC3_INT32 *spec_inv_idx,
                          LC3_INT hrmode
                          );

void processMdctShaping_fl(LC3_FLOAT x[], LC3_FLOAT gains[], const LC3_INT bands_offset[], LC3_INT fdns_npts);

void processResidualCoding_fl(LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gain, LC3_INT L_spec, LC3_INT targetBits, LC3_INT nBits, uint8_t * resBits,
                              LC3_INT* numResBits
                              , LC3_INT hrmode
);

void processResidualDecoding_fl(LC3_INT* bitsRead, LC3_FLOAT x[], LC3_INT L_spec, uint8_t prm[], LC3_INT resQBits
                                , LC3_INT hrmode
);

void processAdjustGlobalGain_fl(LC3_INT* gg_idx, LC3_INT gg_idx_min, LC3_INT gg_idx_off, LC3_FLOAT* gain, LC3_INT target, LC3_INT nBits,
                                LC3_INT* gainChange, LC3_INT fs_idx
                                , LC3_INT16 hrmode, LC3_INT16 frame_dms
                                );

void processApplyGlobalGain_fl(LC3_FLOAT x[], LC3_INT xLen, LC3_INT global_gain_idx, LC3_INT global_gain_off);

void processNoiseFactor_fl(LC3_INT* fac_ns_idx, LC3_FLOAT x[], LC3_INT xq[], LC3_FLOAT gg, LC3_INT BW_cutoff_idx, LC3_INT frame_dms,
                           LC3_INT target_bytes
                            );

void processNoiseFilling_fl(LC3_FLOAT xq[], LC3_INT nfseed, LC3_INT fac_ns_idx, LC3_INT bw_stopband, LC3_INT frame_dms, LC3_FLOAT fac_ns_pc, LC3_INT spec_inv_idx);

void processOlpa_fl(LC3_FLOAT* s_12k8, LC3_FLOAT* mem_s12k8, LC3_FLOAT* mem_s6k4, LC3_INT* mem_old_T0, LC3_INT* T0_out,
                    LC3_FLOAT* normcorr_out, LC3_INT len, LC3_INT frame_dms);

void processTnsCoder_fl(LC3_FLOAT* x, LC3_INT bw_cutoff_idx, LC3_INT bw_fcbin, LC3_INT fs, LC3_INT N, LC3_INT frame_dms, LC3_INT nBits,
                        LC3_INT* order_out, LC3_INT* rc_idx, LC3_INT* tns_numfilters, LC3_INT* bits_out
                        , LC3_INT16 near_nyquist_flag
        );                        
void levinsonDurbin(LC3_FLOAT* r, LC3_FLOAT* out_lev, LC3_FLOAT* rc_unq, LC3_FLOAT* error, LC3_INT len);

void processTnsDecoder_fl(LC3_FLOAT* x, LC3_INT* rc_idx, LC3_INT* order, LC3_INT numfilters, LC3_INT bw_fcbin, LC3_INT N, LC3_INT fs);

void processSnsComputeScf_fl(LC3_FLOAT* x, LC3_INT tilt, LC3_INT xLen, LC3_FLOAT* gains, LC3_INT smooth, LC3_FLOAT sns_damping, LC3_FLOAT attdec_damping_factor);

void processSnsInterpolateScf_fl(LC3_FLOAT* gains, LC3_INT encoder_side, LC3_INT bands_number, LC3_FLOAT* gains_LC3_INT);

void processDetectCutoffWarped_fl(LC3_FLOAT* d2, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT* bw_idx);
void processNearNyquistdetector_fl(LC3_INT16* near_nyquist_flag, const LC3_INT fs_idx, const LC3_INT near_nyquist_index,
                                   const LC3_INT bands_number, const LC3_FLOAT* ener);

void processPerBandEnergy_fl(LC3_INT bands_number, const LC3_INT* acc_coeff_per_band, LC3_INT16 hrmode, LC3_INT16 frame_dms, LC3_FLOAT* d2, LC3_FLOAT* d);

void ProcessingIMDCT_fl(LC3_FLOAT* y, LC3_INT yLen, const LC3_FLOAT* win, LC3_INT winLen, LC3_INT last_zeros, LC3_FLOAT* mem, LC3_FLOAT* x,
                        Dct4* dct);

void ProcessingITDA_WIN_OLA_fl(LC3_FLOAT* x_tda, LC3_INT32 yLen, const LC3_FLOAT* win, LC3_INT32 winLen, LC3_INT32 last_zeros, LC3_FLOAT* mem, LC3_FLOAT* x);

void process_ltpf_coder_fl(LC3_FLOAT* xin, LC3_INT xLen, LC3_INT ltpf_enable, LC3_INT pitch_ol, LC3_FLOAT pitch_ol_norm_corr, LC3_INT frame_dms,
                           LC3_FLOAT* mem_old_x, LC3_INT memLen, LC3_FLOAT* mem_norm_corr_past, LC3_INT* mem_on, LC3_FLOAT* mem_pitch,
                           LC3_INT* param, LC3_FLOAT* mem_norm_corr_past_past, LC3_INT* bits);

void process_ltpf_decoder_fl(LC3_FLOAT* x, LC3_INT xLen, LC3_FLOAT* y, LC3_INT fs, LC3_FLOAT* mem_old_x, LC3_FLOAT* mem_old_y,
                             LC3_INT* mem_pitch_LC3_INT, LC3_INT* mem_pitch_fr, LC3_FLOAT* mem_gain, LC3_INT* mem_beta_idx, LC3_INT bfi,
                             LC3_INT* param, LC3_INT* mem_param, LC3_INT conf_beta_idx, LC3_FLOAT conf_beta, LC3_INT concealMethod, LC3_FLOAT damping
                             , LC3_INT *mem_ltpf_active
);

void process_resamp12k8_fl(LC3_FLOAT x[], LC3_INT x_len, LC3_FLOAT mem_in[], LC3_INT mem_in_len, LC3_FLOAT mem_50[], LC3_FLOAT mem_out[],
                           LC3_INT mem_out_len, LC3_FLOAT y[], LC3_INT* y_len, LC3_INT fs_idx, LC3_INT frame_dms, LC3_INT fs);

void write_bit_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT bit);
void write_uint_backward_fl(LC3_UINT8* ptr, LC3_INT* bp_side, LC3_INT* mask_side, LC3_INT val, LC3_INT numbits);

void processAriEncoder_fl(LC3_UINT8* bytes, LC3_INT bp_side, LC3_INT mask_side, LC3_INT* x, LC3_INT* tns_order, LC3_INT tns_numfilters,
                          LC3_INT* tns_idx, LC3_INT lastnz, LC3_INT* codingdata, uint8_t* res_bits, LC3_INT resBitsLen, LC3_INT lsbMode,
                          LC3_INT nbbits, LC3_INT enable_lpc_weighting);

void attack_detector_fl(LC3_FLOAT* in, LC3_INT frame_size, LC3_INT fs, LC3_INT* lastAttackPosition, LC3_FLOAT* accNrg, LC3_INT* attackFlag,
                        LC3_FLOAT* attdec_filter_mem, LC3_INT attackHandlingOn, LC3_INT attdec_nblocks, LC3_INT attdec_hangover_threshold);

void process_snsQuantizesScf_Enc(LC3_FLOAT* env, LC3_INT* index, LC3_FLOAT* envq, Dct2 dct2structSNS);

void process_snsQuantizesScf_Dec(LC3_INT* scf_idx, LC3_FLOAT* scf_q);

void processMdct_fl(LC3_FLOAT* in, LC3_FLOAT* out, Mdct* mdctStruct);

int       alloc_encoder(LC3PLUS_Enc* encoder, int channels);
void      set_enc_frame_params(LC3PLUS_Enc* encoder);
LC3PLUS_Error update_enc_bitrate(LC3PLUS_Enc* encoder, int bitrate);

LC3PLUS_Error FillEncSetup(LC3PLUS_Enc* encoder, int samplerate, int channels);

/* Setup Functions */
int       alloc_decoder(LC3PLUS_Dec* decoder, int samplerate, int channels);
void      set_dec_frame_params(LC3PLUS_Dec* decoder);
LC3PLUS_Error update_dec_bitrate(LC3PLUS_Dec* decoder, int ch, int nBytes);
LC3PLUS_Error FillDecSetup(LC3PLUS_Dec* decoder, int samplerate, int channels, LC3PLUS_PlcMode plc_mode);

int       Enc_LC3PLUS_fl(LC3PLUS_Enc* encoder, void** input, LC3_UINT8* output, int bps
, LC3_INT32 bfi_ext
);
LC3PLUS_Error Dec_LC3PLUS_fl(LC3PLUS_Dec* decoder, LC3_UINT8* input, int input_bytes, void** output, int bps, int bfi_ext);

void* balloc(void* base, size_t* base_size, size_t size);

void processPlcMain_fl(LC3_FLOAT *q_d_fl_c, LC3_FLOAT *syntM_fl_c, LC3PLUS_Dec* decoder, DecSetup* h_DecSetup, LC3_INT bfi,
               PlcAdvSetup *PlcAdvSetup, PlcSetup *PlcSetup, LC3_INT plcMeth, LC3_INT ltpf_pitch_int, LC3_INT ltpf_pitch_fr,
               LC3_INT tilt, const LC3_INT *bands_offset, LC3_INT bands_number, const LC3_INT *bands_offsetPLC,
               LC3_INT n_bandsPLC, LC3_INT16 hrmode, pcState *statePC);

void processPlcUpdate_fl(PlcAdvSetup *PlcAdvSetup, LC3_INT32 frame_length, LC3_FLOAT *syntM, LC3_FLOAT *scf_q,
             LC3_INT32 *nbLostCmpt, LC3_FLOAT *cum_alpha, LC3_INT32 bfi, LC3_INT32 *prevBfi, LC3_INT32 *prevprevBfi);

void processPlcUpdateSpec_fl(LC3_FLOAT *q_d_prev, LC3_FLOAT *q_d_fl_c, LC3_INT yLen);

void processNoiseSubstitution_fl(LC3_FLOAT* spec, LC3_FLOAT* spec_prev, LC3_INT32 yLen);

void process_cutoff_bandwidth(LC3_FLOAT* d_fl, LC3_INT len, LC3_INT bw_bin);
void update_enc_bandwidth(LC3PLUS_Enc* encoder, LC3_INT bandwidth);

/* al_fec.c */
LC3_INT16 fec_get_n_pccw(LC3_INT16 slot_bytes, LC3_INT16 fec_mode, LC3_INT16 ccc_flag);
LC3_INT16 fec_get_data_size(LC3_INT16 fec_mode, LC3_INT16 ccc_flag, LC3_INT16 slot_bytes);
LC3_INT16 fec_get_n_pc(LC3_INT16 fec_mode, LC3_INT16 n_pccw, LC3_INT16 slot_bytes);
void processReorderBitstream_fl(LC3_UINT8* bytes, LC3_INT32 n_pccw, LC3_INT32 n_pc, LC3_INT32 b_left, LC3_INT32 len);
void fec_encoder(LC3_INT16 mode, LC3_INT16 epmr, LC3_UINT8 *iobuf, LC3_INT16 data_bytes, LC3_INT16 slot_bytes, LC3_INT16 n_pccw);
LC3_INT32 fec_decoder(LC3_UINT8 *iobuf, LC3_INT16 slot_bytes, LC3_INT32 *data_bytes, LC3PLUS_EpModeRequest *epmr, LC3_INT16 ccc_flag, LC3_INT16 *n_pccw, LC3_INT32 *bfi,
                LC3_INT16 *be_bp_left, LC3_INT16 *be_bp_right, LC3_INT16 *n_pc, LC3_INT16 *m_fec);

LC3_FLOAT array_max_abs(LC3_FLOAT *in, LC3_INT32 len);

void processPcClassify_fl(LC3_INT32 pitch_present, LC3_INT32 frame_dms, LC3_FLOAT *q_d_prev, LC3_FLOAT *q_old_res, LC3_INT32 yLen, LC3_INT32 spec_inv_idx, LC3_FLOAT stab_fac, LC3_INT32 *bfi);
void processPcMain_fl(LC3_INT32 *bfi, LC3PLUS_Dec* decoder, LC3_FLOAT *sqQdec, DecSetup* h_DecSetup, LC3_INT32 pitch_present, LC3_FLOAT stab_fac, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_INT32 fac_ns_idx, pcState *statePC, LC3_INT32 spec_inv_idx, LC3_INT32 yLen);
void processPcUpdate_fl(LC3_INT32 bfi, LC3_FLOAT *q_res, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_INT32 rframe, LC3_INT32 *BW_cutoff_idx_nf, LC3_INT32 *prev_BW_cutoff_idx_nf, LC3_INT32 fac_ns_idx, LC3_FLOAT *prev_fac_ns, LC3_FLOAT *fac, LC3_FLOAT *q_old_res, LC3_FLOAT *prev_gg, LC3_INT32 spec_inv_idx, LC3_INT32 yLen);
void processPcApply_fl(LC3_FLOAT *q_res, LC3_FLOAT *q_old_res, LC3_FLOAT *q_d_prev, LC3_INT32 spec_inv_idx, LC3_INT32 yLen, LC3_INT32 gg_idx, LC3_INT32 gg_idx_off, LC3_FLOAT *prev_gg, LC3_FLOAT *fac, LC3_INT32 *pc_nbLostCmpt);

void processPlcClassify_fl(LC3_INT plcMeth, LC3_INT *concealMethod, LC3_INT32 *nbLostCmpt, LC3_INT32 bfi,
                          LC3_FLOAT *xcorr, LC3_INT32 framelength, LC3_INT32 frame_dms, LC3_INT32 pitch_int,
                          LC3_INT32 fs, const LC3_INT *band_offsets, LC3_INT32 bands_number, LC3_INT32 tilt, PlcAdvSetup *plcAd
                          , LC3_INT32 hrmode
);
void processPlcComputeStabFacMain_fl(LC3_FLOAT *scf_q, LC3_FLOAT *old_scf_q, LC3_FLOAT *old_old_scf_q, LC3_INT32 bfi, LC3_INT32 prev_bfi, LC3_INT32 prev_prev_bfi, LC3_FLOAT *stab_fac);

void processPlcDampingScramblingMain_fl(LC3_INT32 *ns_seed,
                                        LC3_INT32 *pc_seed, LC3_INT32 ns_nbLostCmpt_pc,
                                        LC3_INT32 ns_nbLostCmpt, LC3_FLOAT *stabFac, LC3_FLOAT *cum_fading_slow, LC3_FLOAT *cum_fading_fast,
                                        LC3_FLOAT *spec_prev, LC3_FLOAT *spec, LC3_INT32 spec_inv_idx, LC3_INT32 yLen, LC3_INT32 bfi,
                                        LC3_INT32 frame_dms, LC3_INT32 concealMethod, LC3_INT32 pitch_present_bfi1, LC3_INT32 pitch_present_bfi2,
                                        LC3_FLOAT *cum_fflcAtten);
void processPlcDampingScrambling_fl(LC3_FLOAT *spec, LC3_INT32 yLen, LC3_INT32 nbLostCmpt, LC3_FLOAT *stabFac, LC3_INT32 processDampScramb,
                                    LC3_FLOAT *cum_fflcAtten, LC3_INT32 pitch_present, LC3_INT32 frame_dms, LC3_FLOAT *cum_fading_slow,
                                    LC3_FLOAT *cum_fading_fast, LC3_INT32 *seed, LC3_INT32 spec_inv_idx);

void plc_phEcu_F0_refine_first(LC3_INT32 *plocs, LC3_INT32 n_plocs, LC3_FLOAT *f0est, const LC3_INT32 Xabs_len,
                                LC3_FLOAT *f0binPtr, LC3_FLOAT *f0gainPtr, const LC3_INT32 nSubm);
void processTdcLpcEstimation_fl(LC3_FLOAT *r, LC3_INT32 fs_idx, LC3_INT32 len, LC3_FLOAT *A, LC3_INT32 frame_dms);

LC3_FLOAT plc_phEcuSetF0Hz(LC3_INT32 fs, LC3_FLOAT *old_pitchPtr);

void plc_phEcu_processPLCspec2shape(LC3_INT16 prev_bfi, LC3_INT16 bfi, LC3_FLOAT q_d[], LC3_INT32 yLen, LC3_FLOAT *stPhECU_oold_grp_shape, LC3_FLOAT *stPhECU_old_grp_shape);
void plc_phEcu_LF_peak_analysis(LC3_INT32 *plocs, LC3_INT32 *n_plocs, LC3_FLOAT *f0est, const LC3_FLOAT *Xabs,
                                 LC3_FLOAT *f0binPtr,  LC3_FLOAT *f0gainPtr, const LC3_INT32 nSubm);

void plc_phEcu_F0_refine_first(LC3_INT32 *plocs, LC3_INT32 n_plocs, LC3_FLOAT *f0est, const LC3_INT32 Xabs_len,
                                LC3_FLOAT *f0binPtr, LC3_FLOAT *f0gainPtr, const LC3_INT32 nSubm);

LC3_FLOAT plc_phEcu_imax2_jacobsen_mag(const Complex *y, LC3_FLOAT *c_jacobPtr);
LC3_FLOAT plc_phEcu_interp_max(const LC3_FLOAT *y, LC3_INT32 y_len);
void      plc_phEcu_fft_spec2_sqrt_approx(const Complex* x, LC3_INT32 x_len, LC3_FLOAT* x_abs);
LC3_INT32   plc_phEcu_pitch_in_plocs(LC3_INT32* plocs, LC3_INT32 n_plocs);
void      plc_phEcu_spec_ana(LC3_FLOAT* xfp, LC3_INT32 xfp_len, const LC3_FLOAT* whr, 
              LC3_FLOAT* pfind_sensPtr, LC3_INT32* plocs,
              LC3_INT32* n_plocs, LC3_FLOAT* f0est, Complex* x, LC3_INT32* x_len, 
              LC3_FLOAT* f0hzLtpBinPtr, LC3_FLOAT* f0gainLtpPtr, LC3_INT32 bw_idx, Fft* PhEcu_Fft);
void      plc_phEcu_subst_spec(LC3_INT32* plocs, LC3_INT32 n_plocs, LC3_FLOAT* f0est, LC3_INT32 time_offs, Complex* X, LC3_INT32 X_len,
                          LC3_FLOAT* mag_chg_gr, LC3_INT32 *seed, LC3_FLOAT* alpha, LC3_FLOAT* beta, LC3_FLOAT* Xavg,
                          LC3_INT32 t_adv_in, LC3_INT32 Lprot, LC3_INT32 delta_corr, LC3_FLOAT *corr_phase_dbg,
                          LC3_FLOAT *X_i_new_re_dbg, LC3_FLOAT *X_i_new_im_dbg);
void plc_phEcu_rec_frame(Complex *X_in, LC3_INT32 xfp_len, LC3_INT32 Lecu, const LC3_FLOAT *whr, const LC3_FLOAT *winMDCT, LC3_INT32 Lprot,
               LC3_FLOAT *xfp, LC3_INT32 time_offs, LC3_FLOAT *x_out,
               Complex *full_spec_dbg, LC3_FLOAT* ifft_out_dbg, LC3_FLOAT* xsubst_dbg, 
               LC3_INT32 LA_ZEROS, LC3_INT32 LA, Fft* PhEcu_Ifft

   );
void plc_phEcu_tba_spect_Xavg(LC3_INT32 fs_idx, LC3_INT32 n_grp, LC3_FLOAT *oold_spec_shape, 
                              LC3_FLOAT *oold_EwPtr,     LC3_FLOAT *old_spec_shape,   LC3_FLOAT *old_EwPtr, 
                              LC3_FLOAT *gr_pow_left, LC3_FLOAT *gr_pow_right, LC3_FLOAT *Xavg);
void plc_phEcu_tba_per_band_gain(LC3_INT32 n_grp, LC3_FLOAT *gr_pow_left, LC3_FLOAT *gr_pow_right, LC3_FLOAT *trans, LC3_FLOAT *grp_pow_change);
void plc_phEcu_tba_trans_dect_gains(LC3_INT32 burst_len, LC3_INT32 n_grp, LC3_FLOAT *grp_pow_change,                   
                          LC3_FLOAT *stPhECU_beta_mute, LC3_FLOAT *stPhECU_mag_chg_1st,
                          LC3_FLOAT *alpha, LC3_FLOAT *beta, LC3_FLOAT *mag_chg, LC3_FLOAT *ph_dith, LC3_INT32 *tr_dec,
                          LC3_FLOAT *att_val, LC3_INT32 *attDegreeFrames, LC3_FLOAT *thresh_dbg);
void plc_phEcu_trans_burst_ana_sub(LC3_INT32 fs_idx, LC3_INT32 burst_len, LC3_INT32 n_grp, LC3_FLOAT *oold_spect_shape, 
                         LC3_FLOAT *oold_EwPtr,      LC3_FLOAT *old_spect_shape, 
                         LC3_FLOAT *old_EwPtr, LC3_FLOAT *stPhECU_beta_mute, 
                         LC3_FLOAT *stPhECU_mag_chg_1st, LC3_FLOAT *stPhECU_Xavg, LC3_FLOAT *alpha, LC3_FLOAT *beta, LC3_FLOAT *mag_chg,
                         LC3_INT32 *tr_dec_dbg, LC3_FLOAT *gpc_dbg);
void plc_phEcu_hq_ecu(
    LC3_FLOAT *f0binPtr, LC3_FLOAT *f0ltpGainPtr, 
    LC3_FLOAT *xfp, LC3_INT16 prev_bfi, LC3_INT32 *short_flag_prev, 
    LC3_INT32 fs, LC3_INT32 * time_offs,
    Complex *X_sav_m, LC3_INT32 *n_plocs, LC3_INT32 *plocs, LC3_FLOAT *f0est,  const LC3_FLOAT *mdctWin, 
    LC3_FLOAT *env_stabPtr, LC3_INT32 delta_corr, 
    LC3_FLOAT *pfind_sensPtr,
    LC3_INT32 PhECU_LA, LC3_INT32 t_adv, const LC3_FLOAT *winWhr, LC3_FLOAT *oold_grp_shape, 
    LC3_FLOAT *oold_EwPtr, LC3_FLOAT *old_grp_shape, 
    LC3_FLOAT *old_EwPtr,
    LC3_FLOAT *st_beta_mute, LC3_FLOAT *st_mag_chg_1st, LC3_FLOAT *st_Xavg, LC3_INT32 LA_ZEROS, LC3_FLOAT *x_tda, LC3_FLOAT *xsubst_dbg, Complex *X_out_m_dbg,
    LC3_INT32 *seed_dbg, LC3_FLOAT *mag_chg_dbg, LC3_INT32 *tr_dec_dbg, LC3_FLOAT *gpc_dbg, LC3_FLOAT *X_i_new_re_dbg, LC3_FLOAT *X_i_new_im_dbg, LC3_FLOAT *corr_phase_dbg
     ,Fft* PhEcu_Fft,Fft* PhEcu_Ifft
);

void processTdcPreemphasis_fl(LC3_FLOAT *in, LC3_FLOAT *pre_emph_factor, LC3_INT32 n_bands);

void processTdcTdac_fl(const LC3_FLOAT *synth_inp, const LC3_FLOAT *win, LC3_INT32 frame_length, LC3_INT32 la_zeroes, LC3_FLOAT *ola_mem);
void processTdcInverseOdft_fl(LC3_FLOAT *in, LC3_INT32 n_bands, LC3_FLOAT *out, LC3_INT32 lpc_order);

void processTdcApply_fl(const LC3_INT32 pitch_LC3_INT, const LC3_FLOAT *preemphFac, const LC3_FLOAT* A, const LC3_INT32 lpc_order, const LC3_FLOAT* pcmbufHist, const LC3_INT32 max_len_pcm_plc, const LC3_INT32 N, const LC3_INT32 frame_dms,
                                     const LC3_INT32 SampRate, const LC3_INT32 nbLostCmpt, const LC3_INT32 overlap, const LC3_FLOAT *stabFac, LC3_FLOAT harmonicBuf[MAX_PITCH], LC3_FLOAT synthHist[M],
                                     LC3_INT32* fract, LC3_INT16* seed, LC3_FLOAT* gain_c, LC3_FLOAT* alpha, LC3_FLOAT* synth);
void* balloc(void* base, size_t* base_size, size_t size);



#endif
