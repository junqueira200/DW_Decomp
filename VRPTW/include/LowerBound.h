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

namespace LowerBoundNS
{



    bool getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                           Eigen::VectorXd&                          vetDist,
                           int                                       dest,
                           LabelingAlgorithmNS::LabelingData*        lDataPtr);

    void copyDistMat(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost, Eigen::MatrixXd& distMat);

}

#endif //DW_LOWERBOUND_H
