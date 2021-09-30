/******************************************************************************
*                        ETSI TS 103 634 V1.3.1                               *
*              Low Complexity Communication Codec Plus (LC3plus)              *
*                                                                             *
* Copyright licence is solely granted through ETSI Intellectual Property      *
* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
* estoppel or otherwise.                                                      *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "fec_control_unit.h"
#include "fec_control_unit_config.h"
#include "fec_control_unit_ferEstim.h"
#include "fec_control_unit_ferReductionEstim.h"


struct t_fecControlUnit {
    LC3PLUS_EpMode                   m_epModeMin;              // min ep mode in session
    LC3PLUS_EpMode                   m_epModeMax;              // max ep mode in session
    unsigned int                     m_epModeNum;              // number of ep modes in session
    LC3PLUS_EpMode                   m_lastValidEpMode;        // ep mode of last decodable frame
    int                              m_lowConfidenceRunLength; // counter for low confidence epmrs
    HANDLE_FER_ESTIMATOR             m_ferEstim[FER_BUF_NUM];  // ptr to FER_BUF_NUM structs for fer estimation
    HANDLE_FER_REDUCTION_ESTIMATOR   m_ferReductionEstim;      // struct for fer reduction estimation
}; // HANDLE_FEC_CONTROL_UNIT



/******************************************************************************************************
 *  private functions
 *****************************************************************************************************/


/******************************************************************************/
#define SIMPLE_RAND_MAX 65535
static uint16_t seed = 12345;
static uint16_t simple_rand(void)
{
    return seed = (13849 + (seed * 31821)) & 65535;
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static void getOptEpmr(HANDLE_FEC_CONTROL_UNIT     const self,            // i
                       float               const * const FER,             // i
                       LC3PLUS_EpModeRequest     * const idxMaxMosLqo)    // o
{
    assert(self != NULL);

    float fer[4];   // [4] epmode, idx 0 -> self->m_epModeMin
    float mos[4];   // [4] epmode, idx 0 -> self->m_epModeMin

    // linear interpolation
    for (int i = 0, m = EP2EPMR(self->m_epModeMin); i < self->m_epModeNum; i ++, m++)
    {
        fer[i] = FER[i] * 100; // fer in percent
        if (fer[i] >= FER_MAX)
        {
            fer[i] = FER_MAX;
            mos[i] = mosLqo[i][FER_MAX];
        }
        else
        {
            int   fer1 = (int)fer[i];
            int   fer2 = fer1+1;
            float mos1 = mosLqo[m][fer1];
            float mos2 = mosLqo[m][fer2];
            float m  = (mos1-mos2)/(fer1-fer2);
            float n  = (fer1*mos2-fer2*mos1)/(fer1-fer2);
            mos[i] = m*fer[i]+n;
        }
    }

    
    // determine best index
    *idxMaxMosLqo = EP2EPMR(self->m_epModeMin);
    for (int m = EP2EPMR(self->m_epModeMin) + 1, i = 1; m <= EP2EPMR(self->m_epModeMax); m ++, i++)
    {
        if (mos[*idxMaxMosLqo - EP2EPMR(self->m_epModeMin)] < mos[i])
        {
            *idxMaxMosLqo = m;
        }
    }
}

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static void bfiExtrapolation(HANDLE_FEC_CONTROL_UNIT const self,              // i
                             int                     const epokFlags,         // i
                             int                     const currentEpMode,     // i
                             float                 * const bfiPerMode)        // o
{
    assert(self != NULL);

    int          epokMask;
    float        ferReductionValues[4];           // [4] epmode, idx 0 -> LC3PLUS_EP_ZERO (idx 0 not used)
    float        ferReductionCum       = 1;
    unsigned int randomNumber          = simple_rand();

    
    /***************************************************
     * set bfiPerMode for each epMode based on epokMaks
     ***************************************************/
    
    for (int m = EP2EPMR(self->m_epModeMin); m <= EP2EPMR(self->m_epModeMax); m ++)
    {
        bfiPerMode[m] = 1;
        if (epokFlags >= 0)
        {
            epokMask      = 0x0001 << m;
            bfiPerMode[m] = (epokFlags & epokMask) == 0;
        }
    }

    /*************************************************
     * get fer reduction estimates for higher modes
     *************************************************/

    ferReductionEstimator_GetFerReductionEstimates(self->m_ferReductionEstim,
                                                   &ferReductionValues[self->m_epModeMin]);
    
    /***********************************************************************
     * extrapolate bfiPerMode for higher modes, if bfi is 1 in current mode
     ***********************************************************************/
    
    for (int m = EP2EPMR(currentEpMode) + 1; m <= EP2EPMR(self->m_epModeMax); m ++)
    {
        ferReductionCum *= ferReductionValues[m];
        bfiPerMode[m]   *= ferReductionCum;

        bfiPerMode[m] = randomNumber < bfiPerMode[m] * SIMPLE_RAND_MAX;
    }
}


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static void bufferUpdate(HANDLE_FEC_CONTROL_UNIT const self,
                         int                     const currentEpMode,
                         float                   const bfiPerMode[4])
{
    assert(self != NULL);

    /***********************
     * update fer estimator
     ***********************/
    for (int i = 0; i < FER_BUF_NUM; i++)
    {
        ferEstimator_Update(self->m_ferEstim[i],
                            &bfiPerMode[EP2EPMR(self->m_epModeMin)]);
    }

    /*********************************
     * update fer reduction estimator
     *********************************/
    float ferValues[4]; // [4] epmode, idx 0 -> self->m_epModeMin

    ferEstimator_GetFerValues(self->m_ferEstim[FER_BUF_NUM - 1], ferValues);
    
    ferReductionEstimator_Update(self->m_ferReductionEstim,
                                 currentEpMode - self->m_epModeMin,
                                 &bfiPerMode[self->m_epModeMin],
                                 ferValues[0]);
}


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static LC3PLUS_EpModeRequest epmrSelection_bernoulli(HANDLE_FEC_CONTROL_UNIT const self,
                                                     float                     ferValues[3][4],
                                                     float                     consistency[3][4])
{
    assert(self != NULL);

    /***************************************************
     * determine longest stable buffer for each ep mode
     ***************************************************/
    int                   stableBuf[4];        // [4] epmode, idx 0 -> LC3PLUS_EP_ZERO
    float                 ferValueSelected[4]; // [4] epmode, idx 0 -> self->m_epModeMin
    LC3PLUS_EpModeRequest optEpmr[3];          // [stable, medium, long]
    LC3PLUS_EpModeRequest epmrTp = LC3PLUS_EPMR_HIGH;

    for (int j = 0, m = EP2EPMR(self->m_epModeMin); j < self->m_epModeNum; j ++, m++)// for each ep mode
    {
        stableBuf[m] = 0;                                               // preselect the short estimator as default
        
        if      (consistency[2][j] > STABILITY_THRESHOLD)               // Is the long estimator stable?
        {
            stableBuf[m] = 2;                                           // preselect the long estimator
            
            if      (ferValues[0][j] > FER_SM * ferValues[1][j] &&      // strong increase in short estimator?
                     ferValues[0][j] > FER_SL * ferValues[2][j])
            {
                stableBuf[m] = 0;                                       // select the short estimator
                                                                        
                ferEstimator_SetConsistencies(self->m_ferEstim[1], 0);  // mark intermediate estimator as unstable until next update
                ferEstimator_SetConsistencies(self->m_ferEstim[2], 0);  // mark long estimator as unstable until next update
            }
            else if (ferValues[1][j] > FER_ML * ferValues[2][j])        // strong increase in middle estimator?
            {
                stableBuf[m] = 1;                                       // select the middle estimator
                ferEstimator_SetConsistencies(self->m_ferEstim[2], 0);  // mark long estimator as unstable until next update
            }
        }
        else if (consistency[1][j] > STABILITY_THRESHOLD)               // Is the middle estimator stable?
        {
            stableBuf[m] = 1;                                           // preselect the middle estimator
            
            if      (ferValues[0][j] > FER_SM * ferValues[1][j])        // strong increase in short estimator?
            {
                stableBuf[m] = 0;                                       // select the short estimator
                ferEstimator_SetConsistencies(self->m_ferEstim[1], 0);  // mark intermediate estimator as unstable until next update
                ferEstimator_SetConsistencies(self->m_ferEstim[2], 0);  // mark long estimator as unstable until next update
            }
        }
        
    }

    /************************************
     * opt for lowest common denominator
     ************************************/
    int minLen = FER_BUF_NUM-1;  // corresponds to long range estimator
        
    for (int m = EP2EPMR(self->m_epModeMin); m <= EP2EPMR(self->m_epModeMax); m++)
    {
        if (stableBuf[m] < minLen)
        {
            minLen = stableBuf[m];
        }
    }

    for (int j = 0, m = EP2EPMR(self->m_epModeMin); j < self->m_epModeNum; j ++, m++) // for each ep mode
    {
        stableBuf[m]        = minLen;
        ferValueSelected[j] = ferValues[minLen][j];
    }

    
    /**********************
     * retrieve fer values
     **********************/
    getOptEpmr(self, ferValueSelected        , &(optEpmr[0])); // stable
    getOptEpmr(self, ferValues[FER_BUF_NUM-2], &(optEpmr[1])); // medium
    getOptEpmr(self, ferValues[FER_BUF_NUM-1], &(optEpmr[2])); // long

    /******************
     * posterior rules
     ******************/
    epmrTp = optEpmr[0];                                      // take stable

    if (stableBuf[optEpmr[0]] < 2 && epmrTp + 1 < optEpmr[2]) // if stable is smaller than long
    {
        epmrTp = optEpmr[2] - 1;                              // take long-1, if stable is < long-1
    }
    if (stableBuf[optEpmr[0]] < 1 && epmrTp     < optEpmr[1]) // if stable is smaller than mid
    {
        epmrTp = optEpmr[1];                                  // take mid, if stable < mid
    }
    
    
    return epmrTp;
}  


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static LC3PLUS_EpModeRequest epmrSelection(HANDLE_FEC_CONTROL_UNIT const self)
{
    assert(self != NULL);

    LC3PLUS_EpModeRequest epmrTp;
    float             ferValues  [3][4]; // [3] - buffer, [4] epmode, idx 0 -> self->m_epModeMin
    float             consistency[3][4]; // [3] - buffer, [4] epmode, idx 0 -> self->m_epModeMin

    for (int i = 0; i < FER_BUF_NUM; i++)
    {
        ferEstimator_GetFerValues    (self->m_ferEstim[i], ferValues  [i]);
        ferEstimator_GetConsistencies(self->m_ferEstim[i], consistency[i]);
    }

    epmrTp = epmrSelection_bernoulli   (self, ferValues, consistency);
    
    return epmrTp;
}


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
static void epmodeSelection(HANDLE_FEC_CONTROL_UNIT const self,        // i/o
                            LC3PLUS_EpModeRequest   const epmrOp,      // i
                            LC3PLUS_EpMode        * const epModeTp,    // o
                            LC3PLUS_EpModeRequest * const epmrTp)      // i/o
{
    assert(self != NULL);
    assert(epModeTp != NULL);

    // Set requested error protection mode in our encoder:
    // take over immediately if confidence is high or medium
    if( epmrOp < LC3PLUS_EPMR_ZERO_NC ) // high or medium confidence
    {
        *epModeTp                    = EPMR2EP(epmrOp);
        self->m_lowConfidenceRunLength = 0;
    }
    else // low confidence
    {
        // increase local encoder's protection mode if no valid (high/medium confidence)
        // EPMR has been received in the last EPMR_CONFIDENCE_LEN frames
        // and request ep_mode increase from remote encoder
        if( ++(self->m_lowConfidenceRunLength) >= EPMR_CONFIDENCE_LEN )
        {
            int increment = self->m_lowConfidenceRunLength / EPMR_CONFIDENCE_LEN;
            *epModeTp = MIN(*epModeTp + increment, self->m_epModeMax);

            if (epmrTp != NULL)
            {
                *epmrTp = MIN((*epmrTp) + increment, EP2EPMR(self->m_epModeMax));
            }
        }
    }
}


/******************************************************************************************************
 *  public functions
 *****************************************************************************************************/

/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
HANDLE_FEC_CONTROL_UNIT fecControlUnit_Create(LC3PLUS_EpMode         const epModeMin,
                                              LC3PLUS_EpMode         const epModeMax)
{
    HANDLE_FEC_CONTROL_UNIT self = NULL;

    assert (FER_BUF_NUM         == 3);
    assert (epModeMin           >= LC3PLUS_EP_ZERO);
    assert (epModeMax           <= LC3PLUS_EP_HIGH);
    assert (epModeMax-epModeMin <  4);
   
    self = (HANDLE_FEC_CONTROL_UNIT) calloc(1, sizeof(struct t_fecControlUnit));
    assert (self != NULL);
    
    self->m_epModeMin = epModeMin;
    self->m_epModeMax = epModeMax;
    self->m_epModeNum = epModeMax - epModeMin + 1;

    self->m_lastValidEpMode        = LC3PLUS_EP_ZERO;
    self->m_lowConfidenceRunLength = 0;


    for(int i = 0; i < FER_BUF_NUM; i++)
    {
        self->m_ferEstim[i] = ferEstimator_Create(self->m_epModeNum,
                                                  bufConfig[i][0],
                                                  bufConfig[i][1]);

        ferEstimator_Init(self->m_ferEstim[i]);
    }

    if ( self->m_epModeNum > 1)
    {
        self->m_ferReductionEstim = ferReductionEstimator_Create(self->m_epModeNum,
                                                                 FER_REDUTCION_LENGTH,
                                                                 FER_REDUCTION_RATE);
        
        ferReductionEstimator_Init(self->m_ferReductionEstim,
                                   &default_reduction_factors[EP2EPMR(self->m_epModeMin)]);
    }

    return self;
}


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
void fecControlUnit_Destroy(HANDLE_FEC_CONTROL_UNIT * const self)
{
    if (*self != NULL)
    {
        for (int i = 0; i < FER_BUF_NUM; i++)
        {
            if (&(*self)->m_ferEstim[i] != NULL)
            {
                ferEstimator_Destroy(&(*self)->m_ferEstim[i]);
            } 
        }

        if ((*self)->m_ferReductionEstim != NULL)
        {
            ferReductionEstimator_Destroy(&(*self)->m_ferReductionEstim);
        } 

        free (*self);
        *self = NULL;
    }
}


/*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+**+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*/
int fecControlUnit_Apply(HANDLE_FEC_CONTROL_UNIT       const self,              // i
                         const int                   * const epokFlags,         // i
                         const LC3PLUS_EpModeRequest     * const epmrOp,        // i
                         const int                   * const epModeOp,          // i
                         LC3PLUS_EpModeRequest           * const epmrTp,        // o
                         LC3PLUS_EpMode                  * const epModeTp)      // o
{
    float bfiPerMode[4]; // [4] epmode, idx 0 -> LC3PLUS_EP_ZERO
    LC3PLUS_EpMode currentEpMode;

    int sink = 0, source = 0;

    // determine whether this part is source or sink or both
    if (epokFlags != NULL && epModeOp != NULL && epmrTp != NULL)
    {
        sink = 1;
    }

    if (epmrOp != NULL && epModeTp != NULL)
    {
        source = 1;
    }

    if (!(source | sink) || self == NULL)
    {
        goto error;
    }

    if (sink)
    {
        if (!(*epokFlags == 0 || *epokFlags == 8 || *epokFlags == 12 || *epokFlags == 14 || *epokFlags == 15))
            goto error;
        if (!(*epokFlags == 0 || (*epModeOp >= LC3PLUS_EP_ZERO && *epModeOp <= LC3PLUS_EP_HIGH)))
            goto error;
    }

    if (source)
    {
        if (!(*epmrOp >= LC3PLUS_EPMR_ZERO && *epmrOp <= LC3PLUS_EPMR_HIGH_NC))
            goto error;
    }


    if (sink)
    {
        if (epokFlags == 0)
        {
            /* frame is lost: estimate currentEpMode by last valid EP mode */
            currentEpMode = self->m_lastValidEpMode;
        }
        else
        {
            currentEpMode = *epModeOp;
            self->m_lastValidEpMode = currentEpMode;
        }

        /********************************
         * determine bfi for each epMode
         ********************************/
        bfiExtrapolation(self,
                         *epokFlags,
                         currentEpMode,
                         bfiPerMode); 

        /********************
         * update estimators
         ********************/
        bufferUpdate(self,
                     currentEpMode,
                     bfiPerMode);

        /*****************
         * determine epmr
         *****************/
        *epmrTp = epmrSelection(self);
    }

    if (source)
    {
        /****************************************
         * decide on epmode for backward channel
         ****************************************/
        epmodeSelection(self,
                        *epmrOp,
                        epModeTp,
                        epmrTp);
    }
    return 0;
error:
    return -1;
}
