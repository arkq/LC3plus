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

#include "fec_control_unit_ferEstim.h"
#include "fec_control_unit_meanEstim.h"
#include "fec_control_unit_config.h"

struct t_ferEstimator {
    unsigned int            m_meanEstimNum;  // number of mean estimation buffers
    HANDLE_MEAN_ESTIMATOR * m_meanEstim;     // m_meanEstimNum mean estimation buffer structs
}; // HANDLE_FER_ESTIMATOR


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/



HANDLE_FER_ESTIMATOR ferEstimator_Create(unsigned int const epModeNum,
                                         unsigned int const length,
                                         unsigned int const updateRate)
{
    HANDLE_FER_ESTIMATOR self;

    assert(length % updateRate == 0);

    self = (HANDLE_FER_ESTIMATOR) calloc(1, sizeof(struct t_ferEstimator));
    assert(self != NULL);

    self->m_meanEstimNum = epModeNum;

    self->m_meanEstim = (HANDLE_MEAN_ESTIMATOR*) calloc(self->m_meanEstimNum, sizeof(HANDLE_MEAN_ESTIMATOR));
    
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        self->m_meanEstim[i] = meanEstimator_Create(length, updateRate, 1);
    }
    return self;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferEstimator_Init(HANDLE_FER_ESTIMATOR const self)
{
    assert(self != NULL);

    float fer;

    // initialization with static values
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        switch (i)
        {
        case 3: fer = FER_INIT_EPM4 ; break;
        case 2: fer = FER_INIT_EPM3 ; break;
        case 1: fer = FER_INIT_EPM2 ; break;
        case 0: fer = FER_INIT_EPM1 ; break;
        }
        meanEstimator_Init(self->m_meanEstim[i], fer);
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferEstimator_Destroy(HANDLE_FER_ESTIMATOR * const self)
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
void ferEstimator_Update(HANDLE_FER_ESTIMATOR const self,
                         float        const * const bfi)
{
    assert(self != NULL);
    
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        meanEstimator_Update(self->m_meanEstim[i],
                             bfi[i]);
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferEstimator_GetFerValues(HANDLE_FER_ESTIMATOR const self,
                               float              * const ferValues)
{
    assert(self != NULL);
    assert(ferValues != NULL);
    
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        meanEstimator_GetMeanValue(self->m_meanEstim[i], &(ferValues[i]));
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferEstimator_GetConsistencies(HANDLE_FER_ESTIMATOR const self,
                                   float              * const consistency)
{
    assert(self != NULL);
    assert(consistency != NULL);
    
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        meanEstimator_GetConsistency(self->m_meanEstim[i], &(consistency[i]));
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void ferEstimator_SetConsistencies (HANDLE_FER_ESTIMATOR   const self,
                                    float                  const consistency)
{
    assert(self != NULL);
    
    for (int i = 0; i < self->m_meanEstimNum; i++)
    {
        meanEstimator_SetConsistency(self->m_meanEstim[i], consistency);
    }
}


