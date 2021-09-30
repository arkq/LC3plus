/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __fec_control_unit_ferReductionEstim_h__
#define __fec_control_unit_ferReductionEstim_h__


typedef struct t_ferReductionEstimator * HANDLE_FER_REDUCTION_ESTIMATOR;


HANDLE_FER_REDUCTION_ESTIMATOR   ferReductionEstimator_Create (unsigned int                     const epModeNum,
                                                               unsigned int                     const length,
                                                               unsigned int                     const updateRate);

void                            ferReductionEstimator_Destroy (HANDLE_FER_REDUCTION_ESTIMATOR * const self);


void                               ferReductionEstimator_Init (HANDLE_FER_REDUCTION_ESTIMATOR   const self,
                                                               float                    const * const default_reduction_factors);

void                             ferReductionEstimator_Update (HANDLE_FER_REDUCTION_ESTIMATOR   const self,
                                                               int                              const epModesTillCurrEpMode,
                                                               float                    const * const bfis,
                                                               float                            const current_fer);


void           ferReductionEstimator_GetFerReductionEstimates (HANDLE_FER_REDUCTION_ESTIMATOR   const self,
                                                               float *                          const ferReductionValues);

#endif // __fec_control_unit_ferReductionEstim_h__
