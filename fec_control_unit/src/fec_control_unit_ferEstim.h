/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __fec_control_unit_ferEstim_h__
#define __fec_control_unit_ferEstim_h__


typedef struct t_ferEstimator * HANDLE_FER_ESTIMATOR;


HANDLE_FER_ESTIMATOR ferEstimator_Create (unsigned int           const epModeNum,
                                          unsigned int           const length,
                                          unsigned int           const updateRate);

void                ferEstimator_Destroy (HANDLE_FER_ESTIMATOR * const self);


void                   ferEstimator_Init (HANDLE_FER_ESTIMATOR   const self);

void                 ferEstimator_Update (HANDLE_FER_ESTIMATOR   const self,
                                          float          const * const bfi);


void           ferEstimator_GetFerValues (HANDLE_FER_ESTIMATOR   const self,
                                          float                * const ferValues);

void       ferEstimator_GetConsistencies (HANDLE_FER_ESTIMATOR   const self,
                                          float                * const consistency);

void       ferEstimator_SetConsistencies (HANDLE_FER_ESTIMATOR   const self,
                                          float                  const consistency);
#endif // __fec_control_unit_ferEstim_h__
