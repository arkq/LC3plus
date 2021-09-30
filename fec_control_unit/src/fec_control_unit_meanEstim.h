/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#ifndef __fec_control_unit_meanEstim_h__
#define __fec_control_unit_meanEstim_h__


typedef struct t_meanEstimator * HANDLE_MEAN_ESTIMATOR;


HANDLE_MEAN_ESTIMATOR meanEstimator_Create         (unsigned int            const length,
                                                    unsigned int            const updateRate,
                                                    unsigned int            const haveConsistency);

void                  meanEstimator_Destroy        (HANDLE_MEAN_ESTIMATOR * const self);

void                  meanEstimator_Init           (HANDLE_MEAN_ESTIMATOR   const self,
                                                    float                   const value);

int                   meanEstimator_Update         (HANDLE_MEAN_ESTIMATOR   const self,
                                                    float                   const newValue);

void                  meanEstimator_GetMeanValue   (HANDLE_MEAN_ESTIMATOR   const self,
                                                    float                 * const value);

void                  meanEstimator_GetConsistency (HANDLE_MEAN_ESTIMATOR   const self,
                                                    float                 * const consistency);

void                  meanEstimator_SetConsistency (HANDLE_MEAN_ESTIMATOR   const self,
                                                    float                   const consistency);


#endif // __fec_control_unit_meanEstim_h__
