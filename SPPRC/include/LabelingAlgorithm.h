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

#include "Aux.h"
//#include "DW_Decomp.h"
#include "LabelingConstants.h"
#include <Eigen/Eigen>
#include <bitset>
#include <boost/array.hpp>
#include <boost/container/set.hpp>
#include <boost/unordered_set.hpp>
// #include <bits/stdc++.h>
#include "Label.h"
#include "NgSet.h"

// typedef double FloatType;

namespace LabelingAlgorithmNS
{

inline Vector<FloatType> vetValueOfReducedCostsG;
inline int               G_numRecDominance = 2;

enum LabelingTypeAlg
{
    AlgForward,
    AlgBackward,
    AlgBidirectional
};

class Arc
{
  public:
    int       i, j;
    FloatType redCost;

    void set(int i_, int j_, FloatType redCost_)
    {
        i = i_;
        j = j_;
        redCost = redCost_;
    }
};

INLINE
bool operator<(const Arc &arc0, const Arc &arc1) { return arc0.redCost < arc1.redCost; }

class SortRoute
{
  public:
    Eigen::VectorXi vetRoute;
    Eigen::VectorXi vetSortRoute;
    FloatType       dist;

    void computeDistance(const EigenMatrixRowD &matDist);
};

class HashSortRoute
{
  public:
    SortRoute &sortRoute;
    uint64_t   hash;

    HashSortRoute(SortRoute &sortRoute_);
};

inline bool operator==(const HashSortRoute &sol0, const HashSortRoute &sol1)
{
    if(sol0.hash != sol1.hash)
        return false;

    if(sol0.sortRoute.vetSortRoute.size() != sol1.sortRoute.vetSortRoute.size())
        return false;

    for(int i = 0; i < sol0.sortRoute.vetSortRoute.size(); ++i)
    {
        if(sol0.sortRoute.vetSortRoute[i] != sol1.sortRoute.vetSortRoute[i])
            return false;
    }

    return true;
}

inline std::size_t hash_value(const HashSortRoute &solXHash) { return solXHash.hash; }

typedef Vector<std::unique_ptr<LabelingAlgorithmNS::SortRoute>>       VetSortRoute;
typedef boost::unordered_multiset<LabelingAlgorithmNS::HashSortRoute> MultsetSortRoute;

bool checkIfRouteIsInSet(Label                 *label,
                         VetSortRoute          &vetSortRoute,
                         MultsetSortRoute      &multsetSourRoute,
                         const EigenMatrixRowD &matDist);

bool bidirectionalAlgorithm(int                        	numRes,
                            int                        	numCust,
                            const Vet3D_ResCost       	&vetMatResCostForwardConst,
                            const Vet3D_ResCost       	&vetMatResCostBackward,
                            const MatBoundRes         	&vetVetBound,
                            int                        	dest,
                            const NgSet               	&ngSet,
                            LabelingData              	&lData,
                            //Eigen::MatrixXd           &matColX,
                            int                       	&numSol,
                            //FloatType                  labelStart,
                            int                        	NumMaxLabePerBucket,
                            bool                       	dominaceCheck,
                            FloatType                 	&maxDist,
                            Eigen::VectorX<FloatType> 	&vetRedCost,
                            bool                       	exact,
                            LabelingTypeAlg            	labelingTypeAlg,
                            bool                       	arc,
                            //Eigen::VectorXd           *vetLowerBoundRedCost,
                            EigenMatrixRowD       		*ptrMatDp,
                            double                     	percentage,
                            int*					   	vetSol,
                            int					   		&sizeVetSol,
                            int							capacity_Q);

inline __attribute__((always_inline)) bool checkDominanceSubSet(const Label &l0, const Label &l1)
{
    return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
}

// inline __attribute__((always_inline))
std::pair<bool, bool> checkCompleteDominancePair(const Label &l0, const Label &l1, int numResources);


// first index i with vet[i]->vetResources[0] >= value
INLINE
int lowerBoundIdx(FloatType value, const Bucket* bucket)
{
    int lo = 0, hi = bucket->sizeVetPtrLabel;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if(bucket->vetPtrLabel[mid]->vetResources[0] < value)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

// first index i with vet[i]->vetResources[0] > value
INLINE
int upperBoundIdx(FloatType value, const Bucket* bucket)
{
    int lo = 0, hi = bucket->sizeVetPtrLabel;
    while (lo < hi)
    {
        int mid = lo + (hi - lo) / 2;
        if(!(value < bucket->vetPtrLabel[mid]->vetResources[0]))
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo;
}

INLINE
void removeBatchFromBucket(const VectorI& idx, int cnt, Bucket* bucket)
{
    if (!bucket || cnt <= 0) return;
    const int n = bucket->sizeVetPtrLabel;
    if (n <= 0 || idx[0] >= n) return;

    assert(idx[0] >= 0 && "Index cannot be negative");
    int w = idx[0];
    int k = 0;
    for (int r = idx[0]; r < n; ++r)
    {
        if (r == idx[k])
        {
            assert((k == 0 || idx[k] > idx[k - 1]) && "idx must be strictly ascending and unique");
            if (bucket->vetPtrLabel[r])
                bucket->vetPtrLabel[r]->active = false;
            ++k;
            if (k == cnt)
            {
                // All target elements removed! Bulk copy remaining pointers and exit.
                ++r;
                const int remaining = n - r;
                if (remaining > 0)
                {
                    std::memmove(&bucket->vetPtrLabel[w],
                                 &bucket->vetPtrLabel[r],
                                 remaining * sizeof(bucket->vetPtrLabel[0]));
                    w += remaining;
                }
                break;
            }
            continue;
        }
        bucket->vetPtrLabel[w++] = bucket->vetPtrLabel[r];
    }

    // Clear trailing pointers
    for (int t = w; t < n; ++t)
        bucket->vetPtrLabel[t] = nullptr;

    bucket->sizeVetPtrLabel = w;
}

INLINE
void removeFromBucket(int i, Bucket* bucket)
{
    Label* labelPtr = bucket->vetPtrLabel[i];
    labelPtr->active = false;

    for(int t=i; t < (bucket->sizeVetPtrLabel-1); ++t)
        bucket->vetPtrLabel[t] = bucket->vetPtrLabel[t+1];

    bucket->vetPtrLabel[bucket->sizeVetPtrLabel-1] = nullptr;
    bucket->sizeVetPtrLabel -= 1;
}

INLINE
bool checkCompleteDominance(const Label &l0, const Label &l1, int numResources)
{
    if(l0.cust != l1.cust)
    {
        std::cout << "ERROR, lo.cust(" << l0.cust << ") != l1.cust(" << l1.cust << ")\n\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    if(l0.typeLabel == Forward) [[likely]]
    {
        if(doubleGreater(l0.vetResources[0], l1.vetResources[0], FloatEp) ||
           doubleGreater(l0.vetResources[1], l1.vetResources[1], FloatEp))
            return false;

        for(int r=2; r < numResources; ++r)
        {
            if(doubleGreater(l0.vetResources[r], l1.vetResources[r], FloatEp))
                return false;
        }

    }
    else
    {
        if(doubleGreater(l0.vetResources[0], l1.vetResources[0], FloatEp) ||
           doubleLess(l0.vetResources[1], l1.vetResources[1], FloatEp))
            return false;

        for(int r=2; r < numResources; ++r)
        {
            if(doubleLess(l0.vetResources[r], l1.vetResources[r], FloatEp))
                return false;
        }
    }

    return checkDominanceSubSet(l0, l1);

}


typedef Eigen::Array<FloatType, 1, NumMaxResources> VetBackwardMask;

bool extendLabelForward(const Label         &label,
                        Label               &newLabel,
                        const Vet3D_ResCost &vetMatResCost,
                        const MatBoundRes   &vetVetBound,
                        int                  custI,
                        int                  t,
                        const NgSet         &ngSet,
                        int                  numResources);

bool extendLabelBackward(const Label           &label,
                         Label                 &newLabel,
                         const Vet3D_ResCost   &vetMatResCost,
                         const MatBoundRes     &vetVetBound,
                         int                    custI,
                         int                    t,
                         const NgSet           &ngSet,
                         int                    numResources,
                         const VetBackwardMask &vetBackwardMask);

inline __attribute__((always_inline)) bool extendLabel(const Label           &label,
                                                       Label                 &newLabel,
                                                       const Vet3D_ResCost   &vetMatResCostForward,
                                                       const Vet3D_ResCost   &vetMatResCostBackward,
                                                       const MatBoundRes     &vetVetBound,
                                                       int                    custI,
                                                       int                    t,
                                                       const NgSet           &ngSet,
                                                       int                    numResources,
                                                       const VetBackwardMask &vetBackwardMask)
{
    if(label.typeLabel == Forward)
        return extendLabelForward(
            label, newLabel, vetMatResCostForward, vetVetBound, custI, t, ngSet, numResources);
    else
        return extendLabelBackward(label,
                                   newLabel,
                                   vetMatResCostBackward,
                                   vetVetBound,
                                   custI,
                                   t,
                                   ngSet,
                                   numResources,
                                   vetBackwardMask);
}

inline __attribute__((always_inline)) void
applyLabelMask(Label *label, int numResorces, const VetBackwardMask &vetBackwardMask)
{
    for(int i = 0; i < numResorces; ++i)
        label->vetResources[i] *= vetBackwardMask[i];
}

void startBackwardLabel(Label             *labelPtr,
                        VetBackwardMask   &vetBackwardMask,
                        const MatBoundRes &vetVetBound,
                        int                numResources,
                        int                dest,
                        double             labelStart,
                        int                i,
                        LabelingData      &lData);

Bucket *dominanceIntraBucketForward(int           cust,
                                    Label        *label,
                                    LabelingData &lData,
                                    LabelHeap    &labelHeap,
                                    int           numRes,
                                    int           dest,
                                    int          &correctPos);

Bucket *dominanceIntraBucketBackward(int           cust,
                                     Label        *label,
                                     LabelingData &lData,
                                     LabelHeap    &labelHeap,
                                     int           numRes,
                                     int           dest,
                                     int          &correctPos);

// inline __attribute__((always_inline))
Bucket *dominanceIntraBucketSlow(int           cust,
                                 Label        *label,
                                 LabelingData &lData,
                                 LabelHeap    *labelHeap,
                                 int           numRes,
                                 int           dest,
                                 int          &correctPos);

Bucket *dominanceIntraBucketFast0(int           cust,
                                  Label        *label,
                                  LabelingData &lData,
                                  LabelHeap    *labelHeap,
                                  int           numRes,
                                  int           dest,
                                  int          &correctPos);

Bucket *dominanceIntraBucketFast1(int           cust,
                                  Label        *label,
                                  LabelingData &lData,
                                  LabelHeap    *labelHeap,
                                  int           numRes,
                                  int           dest,
                                  int          &correctPos);

bool dominanceInterBuckets(int           cust,
                           Label        *label,
                           LabelingData &lData,
                           LabelHeap    *labelHeap,
                           int           numRes,
                           int           dest,
                           int          &correctPos);

void changeTypeAlg(LabelingTypeAlg &labelingTypeAlg);

Label *mergeForwardAndBackward(Label                *forwardPtr,
                               Label                *backwardPtr,
                               const ArrayResources &vetMaxResources,
                               const MatBoundRes    &vetVetBound,
                               int                   numResorces,
                               const NgSet          &ngSet);

/*
{
    if(label->typeLabel == Forward)
        return dominanceIntraBucketForward(cust, label, lData, labelHeap, numRes, dest,
correctPos); else return dominanceIntraBucketBackward(cust, label, lData, labelHeap,
numRes, 0, correctPos);
}
*/

Label *getLabel();
void   rmLabel(Label *label);
void   resetLabelPool();
void   startPool();

void        writeNgSet(Label *label, const NgSet &ngset);
void        startGlobalMemory(const Vector<VectorI> &vetRoutes);
void        addToVetRoutesG(const VectorI &route);
void        cleanVetRouteG();
std::string printRoute(Label *label);
void        invertRoutes(Vector<VectorI> &vetRoutes);

bool checkIfAllLabelsInHeapHaveA_Bucket(LabelHeap &labelHeap, LabelingData &lData);

void selectXArcs(const Vet3D_ResCost &vetMatResCostForwardConst,
                 Vet3D_ResCost       &vetMatResCostForward,
                 double               percentage,
                 int                  numCust);
} // namespace LabelingAlgorithmNS
#endif // DW_LABELINGALGORITHM_H
