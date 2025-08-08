/*  *****************************************************************
 *  *****************************************************************
 *  File:    LowerBound.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    24/02/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_LOWERBOUND_H
#define DW_LOWERBOUND_H

#include "LabelingAlgorithm.h"
#include "Label.h"

namespace LowerBoundNS
{
    inline __attribute__((always_inline))
    void extendLabel(const LabelingAlgorithmNS::Label*         labelPtr,
                     LabelingAlgorithmNS::Label*               newLabel,
                     const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                     LabelingAlgorithmNS::LabelingData*        lDataPtr,
                     int                                       custJ)
    {

        int custI = labelPtr->cust;

        newLabel->vetResources[0] = labelPtr->vetResources[0] + vetMatResCost(custI, custJ, 0);
        newLabel->vetResources[1] = 0.0;
        newLabel->i = 0; //lDataPtr->getIndex(0, newLabel->vetResources[0]);
        newLabel->j = 0; //labelPtr->j;
        newLabel->cust = custJ;
        newLabel->bitSetNg = labelPtr->bitSetNg;
        newLabel->bitSetNg[custJ] = true;
        newLabel->typeLabel = LabelingAlgorithmNS::Forward;

        for(int i=0; i < labelPtr->tamRoute; ++i)
            newLabel->vetRoute[i] = labelPtr->vetRoute[i];

        newLabel->vetRoute[labelPtr->tamRoute] = custJ;
        newLabel->tamRoute = labelPtr->tamRoute + 1;

    }


    bool getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                           Eigen::VectorX<FloatType>&                vetDist,
                           int                                       dest,
                           LabelingAlgorithmNS::LabelingData*        lDataPtr);


}

#endif //DW_LOWERBOUND_H
