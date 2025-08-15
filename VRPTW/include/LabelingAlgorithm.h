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

    inline Vector<FloatType> vetValueOfReducedCostsG;

    enum LabelingTypeAlg
    {
        AlgForward,
        AlgBackward,
        AlgBidirectional
    };

    class SortRoute
    {
    public:

        Eigen::VectorXi vetRoute;
        Eigen::VectorXi vetSortRoute;
        FloatType dist;

        void computeDistance(const EigenMatrixRowD& matDist);

    };

    class HashSortRoute
    {
    public:

        SortRoute& sortRoute;
        uint64_t hash;

        HashSortRoute(SortRoute& sortRoute_);
    };


    inline
    bool operator == (const HashSortRoute &sol0, const HashSortRoute &sol1)
    {
        if(sol0.hash != sol1.hash)
            return false;

        if(sol0.sortRoute.vetSortRoute.size() != sol1.sortRoute.vetSortRoute.size())
            return false;

        for(int i=0; i < sol0.sortRoute.vetSortRoute.size(); ++i)
        {
            if(sol0.sortRoute.vetSortRoute[i] != sol1.sortRoute.vetSortRoute[i])
                return false;
        }

        return true;
    }

    inline
    std::size_t  hash_value(const HashSortRoute& solXHash){return solXHash.hash;}


    typedef Vector<std::unique_ptr<LabelingAlgorithmNS::SortRoute>>        VetSortRoute;
    typedef boost::unordered_multiset<LabelingAlgorithmNS::HashSortRoute>  MultsetSortRoute;

    bool checkIfRouteIsInSet(Label* label, VetSortRoute& vetSortRoute, MultsetSortRoute& multsetSourRoute,
                             const EigenMatrixRowD& matDist);


    bool bidirectionalAlgorithm(int numRes, int numCust, const Vet3D_ResCost& vetMatResCostForward,
                                const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& vetVetBound, int dest,
                                const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                                FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                                Eigen::VectorX<FloatType>& vetRedCost, bool exact, LabelingTypeAlg labelingTypeAlg,
                                bool arc, Eigen::VectorX<FloatType>* vetLowerBoundRedCost);

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
    Bucket* dominanceIntraBucketSlow(int cust, Label* label, LabelingData& lData, LabelHeap* labelHeap, int numRes,
                                     int dest, int& correctPos);

    Bucket* dominanceIntraBucketFast0(int cust, Label* label, LabelingData& lData, LabelHeap* labelHeap, int numRes,
                                     int dest, int& correctPos);

    Bucket* dominanceIntraBucketFast1(int cust, Label* label, LabelingData& lData, LabelHeap* labelHeap, int numRes,
                                     int dest, int& correctPos);


    void changeTypeAlg(LabelingTypeAlg& labelingTypeAlg);

    Label* mergeForwardAndBackward(Label* forwardPtr, Label* backwardPtr, const ArrayResources& vetMaxResources,
                                   const MatBoundRes& vetVetBound, int numResorces);


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
    void resetLabelPool();
    void startPool();

    void writeNgSet(Label* label, const NgSet& ngset);
    void startGlobalMemory(const Vector<VectorI>& vetRoutes);
    void addToVetRoutesG(const VectorI& route);
    void cleanVetRouteG();
    std::string printRoute(Label* label);
    void invertRoutes(Vector<VectorI>& vetRoutes);

    bool checkIfAllLabelsInHeapHaveA_Bucket(LabelHeap& labelHeap, LabelingData& lData);
}
#endif //DW_LABELINGALGORITHM_H
