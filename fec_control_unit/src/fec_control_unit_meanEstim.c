/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <stdlib.h>
#include <assert.h>

#include "fec_control_unit_meanEstim.h"

struct t_meanEstimator {
    unsigned int           m_length;                 // history length
    unsigned int           m_updateRate;             // update rate
    unsigned int           m_historyBufLen;          // length of history buffer
    unsigned int           m_frameCounter;           // frame counter (modulo update rate)
    float                  m_meanValue;              // traced value
    float                  m_consistency;            // consistency measure
    float *                m_historyBuffer;          // history buffer
}; // HANDLE_MEAN_ESTIMATOR


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
HANDLE_MEAN_ESTIMATOR meanEstimator_Create(unsigned int const length,
                                           unsigned int const updateRate,
                                           unsigned int const haveConsistency)
{
    HANDLE_MEAN_ESTIMATOR self;

    self = (HANDLE_MEAN_ESTIMATOR) calloc(1, sizeof(struct t_meanEstimator));
    assert(self != NULL);

    self->m_meanValue     = 0;
    self->m_consistency   = haveConsistency ? 0 : -1;
    self->m_length        = length;
    self->m_updateRate    = updateRate;
    self->m_historyBufLen = length / updateRate + 1;
    self->m_frameCounter  = 0;
    self->m_historyBuffer = NULL;
    if (length > 0)
    {
        self->m_historyBuffer = (float*) calloc(length, sizeof(self->m_historyBuffer));
        assert(self->m_historyBuffer != 0);
    }
    return self;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void meanEstimator_Destroy(HANDLE_MEAN_ESTIMATOR * const self)
{
    if (*self != NULL)
    {
        if ( (*self)->m_historyBuffer != NULL )
        {
            free((*self)->m_historyBuffer);
            (*self)->m_historyBuffer = NULL;
        }
        free(*self);
        *self = NULL;
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void meanEstimator_Init(HANDLE_MEAN_ESTIMATOR const self,
                        float                 const value)
{
    assert(self != NULL);
    
    for(int j = 1; j < self->m_historyBufLen; j++)
    {
        self->m_historyBuffer[j] = value * self->m_updateRate;
    }
    self->m_meanValue = value;

}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
int meanEstimator_Update(HANDLE_MEAN_ESTIMATOR const self,
                         float                 const newValue)
{
    assert(self != NULL);
    assert(self->m_frameCounter < self->m_updateRate);
    
    self->m_historyBuffer[0] += newValue;
        
    if (self->m_frameCounter == self->m_updateRate - 1)
    {
        self->m_meanValue = 0;
        for (int j = self->m_historyBufLen - 2; j >= 0; j --)
        {
            self->m_meanValue          += self->m_historyBuffer[j];
            self->m_historyBuffer[j+1]  = self->m_historyBuffer[j];
        }
        self->m_meanValue        /= self->m_length;
        self->m_historyBuffer[0]  = 0;

        if (self->m_consistency >= 0)
        {
            /* determine consistency by comparing estimated variance to expected variance */
            float var_estim     = 0;
            float p             = self->m_meanValue;
            float var_bernoulli = self->m_length * p * (1 - p);
            
            for (int j = 1; j < self->m_historyBufLen; j ++)
            {
                float tmp  = p * self->m_updateRate - self->m_historyBuffer[j];
                var_estim += tmp * tmp;
            }
            /* var_estim = 0                       , if historyBuffer contains BufLen times the same value */
            /* var_estim = var_bernoulli*updateRate, if historyBuffer contains only 0 and 1 in any distribution */
            self->m_consistency = 1;
            if (p > 0 && p < 1)
            {
                self->m_consistency = var_bernoulli / var_estim;
            }
        }
    }

    if (self->m_frameCounter++ == self->m_updateRate - 1)
    {
        self->m_frameCounter = 0;
    }

    return self->m_frameCounter == 0;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void meanEstimator_GetMeanValue (HANDLE_MEAN_ESTIMATOR   const self,
                                 float                 * const value)
{
    assert(self != NULL);
    assert(value != NULL);

    *value = self->m_meanValue;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void meanEstimator_GetConsistency (HANDLE_MEAN_ESTIMATOR const self,
                                   float               * const consistency)
{
    assert(self != NULL);
    assert(consistency != NULL);

    *consistency = self->m_consistency;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void meanEstimator_SetConsistency (HANDLE_MEAN_ESTIMATOR   const self,
                                   float                   const consistency)
{
    assert (self != NULL);

    self->m_consistency = consistency;
}


