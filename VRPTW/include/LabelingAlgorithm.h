/*  *****************************************************************
 *  *****************************************************************
 *  File:    LabelingAlgorithm.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    20/11/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_LABELINGALGORITHM_H
#define DW_LABELINGALGORITHM_H

#include "LabelingConstants.h"
#include <Eigen/Eigen>
#include <bitset>
#include "Aux.h"
#include "DW_Decomp.h"
#include <boost/array.hpp>
#include <boost/container/set.hpp>
#include <bits/stdc++.h>
#include "Label.h"
#include "NgSet.h"

//typedef double FloatType;

namespace LabelingAlgorithmNS
{
    enum LabelingTypeAlg
    {
        AlgForward,
        AlgBackward,
        AlgBidirectional
    };

    bool forwardLabelingAlgorithm(int numRes, int numCust, const Vet3D_ResCost& vetMatResCost,
                                  const MatBoundRes& vetVetBound, int dest, const NgSet& ngSet, LabelingData& lData,
                                  Eigen::MatrixXd& matColX, int& numSol, FloatType labelStart, int NumMaxLabePerBucket,
                                  bool dominaceCheck, FloatType& maxDist, Eigen::VectorX<FloatType>& vetRedCost,
                                  bool exact);

    bool bidirectionalAlgorithm(int numRes, int numCust, const Vet3D_ResCost& vetMatResCostForward,
                                const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& vetVetBound, int dest,
                                const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                                FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                                Eigen::VectorX<FloatType>& vetRedCost, bool exact, LabelingTypeAlg labelingTypeAlg);

    bool labelingAlgorithmm(int numRes, int numCust, const Vet3D_ResCost& vetMatResCostForward,
                            const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& vetVetBound, int dest,
                            const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                            FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                            Eigen::VectorX<FloatType>& vetRedCost, bool exact, LabelingTypeAlg labelingTypeAlg);


    //inline __attribute__((always_inline))
    bool checkCompleteDominance(const Label& l0, const Label& l1, int numResources);

    inline __attribute__((always_inline))
    bool checkDominanceSubSet(const Label &l0, const Label& l1)
    {
        return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
    }

    typedef Eigen::Array<FloatType, 1, NumMaxResources> VetBackwardMask;

    bool extendLabelForward(const Label& label, Label& newLabel, const Vet3D_ResCost& vetMatResCost,
                            const MatBoundRes& vetVetBound, int custI, int t, const NgSet& ngSet, int numResources);

    bool extendLabelBackward(const Label& label, Label& newLabel, const Vet3D_ResCost& vetMatResCost,
                             const MatBoundRes& vetVetBound, int custI, int t, const NgSet& ngSet,
                             int numResources, const VetBackwardMask& vetBackwardMask);

    inline __attribute__((always_inline))
    bool extendLabel(const Label& label, Label& newLabel, const Vet3D_ResCost& vetMatResCostForward,
                     const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& vetVetBound, int custI, int t,
                     const NgSet& ngSet, int numResources, const VetBackwardMask& vetBackwardMask)
    {
        if(label.typeLabel == Forward)
            return extendLabelForward(label, newLabel, vetMatResCostForward, vetVetBound, custI, t, ngSet, numResources);
        else
            return extendLabelBackward(label, newLabel, vetMatResCostBackward, vetVetBound, custI, t, ngSet,
                                       numResources, vetBackwardMask);
    }

    inline __attribute__((always_inline))
    void applyLabelMask(Label* label, int numResorces, const VetBackwardMask& vetBackwardMask)
    {
        for(int i=0; i < numResorces; ++i)
            label->vetResources[i] *= vetBackwardMask[i];
    }

    void startBackwardLabel(Label* labelPtr, VetBackwardMask& vetBackwardMask, const MatBoundRes& vetVetBound,
                            int numResources, int dest, double labelStart, int i, LabelingData& lData);

    Bucket* dominanceIntraBucketForward(int cust, Label* label, LabelingData& lData, LabelHeap& labelHeap, int numRes,
                                        int dest, int& correctPos);

    Bucket* dominanceIntraBucketBackward(int cust, Label* label, LabelingData& lData, LabelHeap& labelHeap, int numRes,
                                        int dest, int& correctPos);

    //inline __attribute__((always_inline))
    Bucket* dominanceIntraBucket(int cust, Label* label, LabelingData& lData, LabelHeap& labelHeap, int numRes,
                                 int dest, int& correctPos);

    void changeTypeAlg(LabelingTypeAlg& labelingTypeAlg);

    /*
    {
        if(label->typeLabel == Forward)
            return dominanceIntraBucketForward(cust, label, lData, labelHeap, numRes, dest, correctPos);
        else
            return dominanceIntraBucketBackward(cust, label, lData, labelHeap, numRes, 0, correctPos);
    }
    */

    Label* getLabel();
    void rmLabel(Label* label);
    void writeNgSet(Label* label, const NgSet& ngset);


}
#endif //DW_LABELINGALGORITHM_H
