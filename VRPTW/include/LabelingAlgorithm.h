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

#include <Eigen/Eigen>
#include "Grafo.h"
#include <bitset>
#include "Aux.h"
#include "DW_Decomp.h"
#include "safe_vector.h"
#include "safe_matrix.h"
#include <boost/array.hpp>
#include <set>
#include "Grafo.h"
#include "Aux.h"
#include <boost/container/set.hpp>
#include <bits/stdc++.h>
#include "safe_3D_matrix.h"


typedef double FloatType;

namespace LabelingAlgorithmNS
{

    inline std::vector<int> vetRoteG;
    inline FloatType minDistG = std::numeric_limits<FloatType>::max();
    inline FloatType maxDistG = std::numeric_limits<FloatType>::min();
    inline bool exactLabelingG = false;

    constexpr int   NumMaxResources   = 2;
    constexpr int   NumMaxRoute       = 256;
    constexpr int   NumMaxCust        = 100;
    constexpr int   NgSetSize         = 5;
    constexpr int   NumBuckets        = 10;
    constexpr int   vetPtrLabelSize   = 10;
    constexpr bool  NullFlush         = false;
    constexpr bool  Print             = false;
    inline    int   numMaxLabelG      = 2000; // 500

    constexpr bool  DominaIterBuckets = true;

    class LabelSet;

    //typedef std::set<LabelSet>::iterator LabelIt;

    struct Bound
    {
        FloatType lowerBound = -std::numeric_limits<FloatType>::infinity();
        FloatType upperBound = std::numeric_limits<FloatType>::infinity();
    };

    std::ostream& operator<< (std::ostream& out, const Bound &bound);

    /// Access (i,j,r) where (i,j) is an arc and r a resource
    typedef Vector3D<FloatType, false> Vet3D_ResCost;

    /// Access (i,r) where i is a customer and r a resource
    typedef Eigen::Matrix<Bound, -1, -1, Eigen::RowMajor> MatBoundRes;

    class NgSet
    {
    public:

        EigenMatrixRowI matNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);
        [[nodiscard]] bool contain(int i, int j) const;
        void setNgSets(const EigenMatrixRowD &matDist);
    };

    struct Node
    {
        int       i;
        FloatType dist;

        bool operator<(const Node &node) const
        {
            return dist < node.dist;
        }
    };


    class Label
    {
    public:

        bool    active   = false;
        int     i        = -1;
        int     j        = -1;
        int     cust     = -1;
        int     pos      = -1;
        int     tamRoute = 0;
        void*   it;
        //int numResources = 1;
        Eigen::Array<FloatType, 1, NumMaxResources> vetResources;
        std::bitset<NumMaxCust>                     bitSetNg;
        boost::array<int, NumMaxRoute>              vetRoute;


        Label() = default;

    };

    class LabelCmp
    {
    public:

        bool operator()(Label *l0, Label *l1) const
        {
            //return doubleLess(l0->vetResources[0], l1->vetResources[0], std::numeric_limits<FloatType>::epsilon());
            return l0->vetResources[0] < l1->vetResources[0];// && l0->vetResources[1] < l1->vetResources[1];
        }

        bool isGreater(Label *l0, Label *l1)const
        {
            return l0->vetResources[0] > l1->vetResources[0];
        }

    };

    class LabelHeap
    {
    public:

        Vector<Label*> vet;
        int heapSize = 0;

        explicit LabelHeap(int tam){vet = Vector<Label*>(tam, nullptr);}
        int parent(int i) { return (i-1)/2;}
        int left(int i) { return (2*i + 1);}
        int right(int i) { return (2*i + 2);}
        bool empty(){return heapSize==0;}

        void heapify(int i);
        [[nodiscard]]Label* extractMin();
        void decreaseKey(int i, FloatType val);
        [[nodiscard]]Label* getMin(){return vet[0];}
        void deleteKey(int i);
        void insertKey(Label* label);
    };


    std::ostream& operator<< (std::ostream& out, const Label &label);


    class Bucket
    {
    public:
        Eigen::VectorX<Label*> vetPtrLabel;
        int                    sizeVetPtrLabel;

        Bucket()
        {
            vetPtrLabel.resize(vetPtrLabelSize);
        }
        void flush()
        {
            if(NullFlush)
            {
                for(int i = 0; i < vetPtrLabel.size(); ++i)
                    vetPtrLabel[i] = nullptr;
            }

            sizeVetPtrLabel = 0;
        }

        void addLabel(Label *labelPtr);
//        bool delLabel(Label *labelPtr);

    };

    class Step
    {
    public:

        FloatType start;
        FloatType end;
        FloatType stepSize;
    };

    class MatBucket
    {
    public:

        Eigen::Matrix<Bucket, -1, -1, Eigen::RowMajor> mat;
        MatBucket()=default;

    };

    class LabelingData
    {
    public:

        /// Access first the customer and after (i,j), where i is the component of the first resource and j
        ///     the component of the second one.
        Eigen::VectorX<MatBucket>       vetMatBucket;

        /// Access first the resorce
        Vector<Matrix<Bound, false>>    vetMatBound;

        Eigen::Vector<Step, 2>          vetStepSize;

        Eigen::Vector<int, 2>           vetNumSteps;

        /// Given a bucket index (i, j)
        EigenMatrixRowI                 matBucketIndex;

        int numMainResources;
        int numCust;
        int numMaxSteps;

        GraphNS::Graph<int> graphBucket;

        LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_, int numMainResources_, int numCust_);
        LabelingData()=default;

        void flushLabel();
        int getIndex(int resource, FloatType val);
        void removeLabel(Label *label);
        Label* getBestLabel(int cust);

        inline __attribute__((always_inline))
        int getIndexGraphBucket(int i, int j)
        {
            return i*vetNumSteps[1] + j;
        }


        // std::set<Label*, LabelCmp> &setLabel
        void dominanceInterBuckets(LabelHeap& labelHeap, int numRes, const int localNumMaxLabel);

        void setupGraphBucket();

    };

    typedef boost::container::multiset<Label*, LabelCmp>::iterator LabelSetIt;


    bool forwardLabelingAlgorithm(const int                     numRes,
                                  const int                     numCust,
                                  const Vet3D_ResCost&          vetMatResCost,
                                  const MatBoundRes&            vetVetBound,
                                  const int                     dest,
                                  const NgSet&                  ngSet,
                                  LabelingData&                 lData,
                                  Eigen::MatrixXd&              matColX,
                                  int&                          numSol,
                                  const FloatType               labelStart,
                                  int                           NumMaxLabePerBucket,
                                  bool                          dominaceCheck,
                                  FloatType&                    maxDist,
                                  Eigen::VectorX<FloatType>&    vetRedCost);

    //inline __attribute__((always_inline))
    bool checkCompleteDominance(const Label& l0, const Label& l1, int numResources);

    inline __attribute__((always_inline))
    bool checkDominanceSubSet(const Label &l0, const Label& l1)
    {
        return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
    }

    bool extendLabel(const Label&          label,
                     Label&                newLabel,
                     const Vet3D_ResCost&  vetMatResCost,
                     const MatBoundRes&    vetVetBound,
                     int                   custI,
                     int                   custJ,
                     const NgSet&          ngSet,
                     int                   numResources);

    void removeCycles(Label &label, const int numCust);
    void removeCycles2(Label &label, const int numCust);
    void updateLabelCost(Label &label, const Vet3D_ResCost &vetMatResCost, FloatType labelStart);

    // Linear index for a nxn matrix
    inline __attribute__((always_inline))
    int getIndex(int i, int j, int numClie){return i*numClie+j;}

    bool checkDistance(const Eigen::Matrix<FloatType, -1, -1, Eigen::RowMajor> &matDist);
    bool containRoute(const Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> &vetLabel, int numSol, Label* label);

    bool labelHaveRoute(std::vector<int> &vetRoute, Label *label);


    Bucket* dominanceIntraBucket(int           cust,
                                 Label*        label,
                                 LabelingData& lData,
                                 LabelHeap&    labelHeap,
                                 int           numRes,
                                 int           dest,
                                 int&          correctPos);

}
#endif //DW_LABELINGALGORITHM_H
