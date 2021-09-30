/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

/*! \file fec_control_unit.h
 *  This header provides the API for the fec control unit as described in Annex D.
 *
 *  The provided implementation is mainly for demonstration and can be further
 *  optimized with respect to memory requirements.
 */

#ifndef __fec_control_unit_h__
#define __fec_control_unit_h__

#include "lc3.h" /* required for LC3PLUS_EpMode and LC3PLUS_EpModeRequest*/

typedef struct t_fecControlUnit * HANDLE_FEC_CONTROL_UNIT;

/*! creates and initializes fec control unit handle
 */


HANDLE_FEC_CONTROL_UNIT fecControlUnit_Create (LC3PLUS_EpMode         const epModeMin,    // i: min ep mode for session
                                               LC3PLUS_EpMode         const epModeMax);   // i: max ep mode for session

/*! destroys fec control unit handle and releases memory
 */

void                    fecControlUnit_Destroy(HANDLE_FEC_CONTROL_UNIT * const self);      // i


/*! Executes fec control unit algorithm for given frame data.
 *
 *  The function may be called from a sink only, or source only or source and sink device.
 *  If called from a source, it computes the EP mode for the next outgoing frame and if
 *  called from a sink, it computes the EP mode request.
 *
 *
 *  \param[in]  self            fec control unit handle initialized by fecControlUnit_Create()
 *  \param[in]  epokFlags       Pointer to epokFlags. Set to NULL if this part is a source only.
 *  \param[in]  epmrOp          Pointer to EP mode request of other part. Set to NULL if this part
 *                              is a source only.
 *  \param[in]  currentEpMode   Pointer to EP mode of last decoded frame. Set to NULL if this part
 *                              is a source only. If received frame was not decodable, as indicated
 *                              by *epokFlags == 0, the value *currentEpMode is ignored and may be set
 *                              to any value.
 *  \param[out] epmrTp          Pointer to EP mode request of this part. Set to NULL if this part is
 *                              a source only.
 *  \param[out] epModeTp        Pointer to EP mode of this part. Set to NULL if this part is a sink only.
 *  \return                     Zero on success and non-zero value if an error occured.
 */
int                     fecControlUnit_Apply  (HANDLE_FEC_CONTROL_UNIT       const self,          // i
                                               const int                   * const epokFlags,     // i
                                               const LC3PLUS_EpModeRequest * const epmrOp,        // i
                                               const int                   * const EpModeOp,      // i
                                               LC3PLUS_EpModeRequest       * const epmrTp,        // o
                                               LC3PLUS_EpMode              * const epModeTp);     // o

#endif // __fec_control_unit_h__
