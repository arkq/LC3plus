/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

content
*******
README.txt
src/fec_control_unit.c
src/fec_control_unit.h
src/fec_control_unit_config.c
src/fec_control_unit_config.h
src/fec_control_unit_ferEstim.c
src/fec_control_unit_ferEstim.h
src/fec_control_unit_ferReductionEstim.c
src/fec_control_unit_ferReductionEstim.h
src/fec_control_unit_meanEstim.c
src/fec_control_unit_meanEstim.h
src/CMakeLists.txt

application
***********
- The interface is provided by means of fec_control_unit.h

- The principal usage is as follows:

    HANDLE_FEC_CONTROL_UNIT h_fecControlUnit;

    h_fecControlUnit = fecControlUnit_Create(minEpMode, maxEpMode);

    fecControlUnit_Apply(h_fecControlUnit, &epokFlags, &epmrOp, &epModeOp, &epmrTp, &epModeTp);

    fecControlUnit_Destroy(&h_fecControlUnit);

- When operated in a sink, the values epokFlags, epmrOp, epModeOp can be retrieved via the LC3plus
  API:

    epokFlags     = lc3plus_dec_get_epok_flags(decoder)
    epModeOp      = lc3plus_dec_get_m_fec(decoder);
    epmrOp        = lc3plus_dec_get_ep_mode_request(decoder);

  The return values epmrTp and epModeTp can be passed to the LC3plus encoder via:

    lc3plus_enc_set_ep_mode_request(encoder, epmrTp)
    lc3plus_enc_set_ep_mode(encoder, epModeTp)

- If the part running the fec control unit is only a source, the pointers to epokFlags and epModeOp
  should be set to NULL and the epmr will be provided in a different way. In that case, no epmrTp
  will be computed and the corresponding pointer may be set to NULL.

- If the part running the fec control unit is only a sink, the return value pointer to epModeTp
  should point to NULL and no epMode for this part will be computed. Furthermore, the value epmrTp
  will need to be transmitted to the source in a custom way.

building the lib
****************
- A cmake list CMakeLists.txt is provided together with the source code. To build the
  library run

    cmake . && make

  in the src folder. Currently, the MOS-LQO tables are specified statically at compile
  time and can be chosen by adding the option -DOPTION_SETUP=NUM with NUM ranging from
  0 to 3. The setups currently supported are

    0 : NB   and 40 byte payload
    1 : WB   and 40 byte payload
    2 : SSWB and 40 byte payload
    3 : SWB  and 80 byte payload

notes
*****
- The current  implementation is  primarily for demonstration.  It still
  uses   floating   point  data   types   and   is  greedy   in   memory
  consumption.  Furthermore,  most  parameters   of  the  algorithm  are
  currently specified at  compile time. Steps that should  be taken from
  here are:

  a) removal of float data types in mean estimator
  b) packing of BER tables
  c) making setup options dynamical
