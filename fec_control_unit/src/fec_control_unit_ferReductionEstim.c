/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "fec_control_unit_meanEstim.h"
#include "fec_control_unit_config.h"
#include "fec_control_unit_ferReductionEstim.h"

struct t_ferReductionEstimator {
    float                   m_default_fer_reduction_factors[3]; //
    float                   m_reference_fer_values[3];          //
    float                   m_current_fer;                      //
    unsigned int            m_meanEstimNum;                     // number of mean estimation buffers
    HANDLE_MEAN_ESTIMATOR * m_meanEstim;                        // m_meanEstimNum mean estimation buffer structs
}; // HANDLE_FER_REDUCTION_ESTIMATOR


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
HANDLE_FER_REDUCTION_ESTIMATOR ferReductionEstimator_Create(unsigned int const epModeNum,
                                                            unsigned int const length,
                                                            unsigned int const updateRate)
{
    HANDLE_FER_REDUCTION_ESTIMATOR self;

    assert(length % updateRate == 0);

    self = (HANDLE_FER_REDUCTION_ESTIMATOR) calloc(1, sizeof(struct t_ferReductionEstimator));
    assert(self != NULL);

    // FER reduction is estimated for each pair of consecutive EP modes.
    self->m_meanEstimNum = epModeNum - 1;

    self->m_meanEstim = (HANDLE_MEAN_ESTIMATOR*) calloc(self->m_meanEstimNum, sizeof(HANDLE_MEAN_ESTIMATOR));

    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        self->m_meanEstim[i] = meanEstimator_Create(length, updateRate, 1);
    }
    return self;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferReductionEstimator_Init(HANDLE_FER_REDUCTION_ESTIMATOR const self,
                                float                  const * const default_reduction_factors)
{
    assert(self != NULL);
    for (int i = 0; i < self->m_meanEstimNum; i ++)
    {
        self->m_default_fer_reduction_factors[i] = default_reduction_factors[i];
        self->m_reference_fer_values         [i] = 0;
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferReductionEstimator_Destroy(HANDLE_FER_REDUCTION_ESTIMATOR * const self)
{
    if (*self != NULL)
    {
        if ((*self)->m_meanEstim)
        {
            for (int i = 0; i < (*self)->m_meanEstimNum; i++)
            {
                if ((*self)->m_meanEstim[i])
                {
                    meanEstimator_Destroy(&(*self)->m_meanEstim[i]);
                    (*self)->m_meanEstim[i] = NULL;
                }
            }
            free((*self)->m_meanEstim);
            (*self)->m_meanEstim = NULL;
        }
        free (*self);
        *self = NULL;
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferReductionEstimator_Update(HANDLE_FER_REDUCTION_ESTIMATOR const self,
                                  int                            const epModesTillCurrEpMode,
                                  float                  const * const bfi,
                                  float                          const current_fer)
{
    assert(self != NULL);
    // fer reduction estimators are updated for epmodes smaller than or equal to the current ep mode
    for (int i = 0; i < epModesTillCurrEpMode; i ++)
    {
        if (bfi[i-1])
        {
            // reduction is estimated by taking the mean values over bfis in mode m for all
            // frames where bfi in mode m - 1 is 1.
            if (meanEstimator_Update(self->m_meanEstim[i], bfi[i]))
            {
                // reference value is set when mean value is updated
                self->m_reference_fer_values[i] = current_fer;
            }
        }
    }
    self->m_current_fer = current_fer;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferReductionEstimator_GetFerReductionEstimates(HANDLE_FER_REDUCTION_ESTIMATOR const self,
                                                    float *                        const ferReductionValues)
{
    assert(self != NULL);

    float consistency;

    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        meanEstimator_GetConsistency(self->m_meanEstim[i], &consistency);
        if (
            (consistency > STABILITY_THRESHOLD) &&
            (self->m_reference_fer_values[i] > self->m_current_fer / FER_REDUCTION_DELTA) &&
            (self->m_reference_fer_values[i] < self->m_current_fer * FER_REDUCTION_DELTA)
            )
        {
            // take estimated value if estimator is stable and reference fer value is close to
            // current fer value.
            meanEstimator_GetMeanValue(self->m_meanEstim[i], &(ferReductionValues[i]));
        }
        else
        {
            // take default value otherwise.
            ferReductionValues[i] = self->m_default_fer_reduction_factors[i];
        }
    }
}
