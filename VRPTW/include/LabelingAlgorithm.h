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
#include <boost/array.hpp>
#include <set>
#include "Grafo.h"
//#include "MemoryPool.h"
#include "Aux.h"
#include <boost/container/set.hpp>
#include <bits/stdc++.h>


typedef double FloatType;

namespace LabelingAlgorithmNS
{
    inline std::vector<int> vetRoteG;

    constexpr int   NumMaxResources   = 2;
    constexpr int   NumMaxRoute       = 300;
    constexpr int   NumMaxCust        = 100;
    constexpr int   NgSetSize         = 5;
    constexpr int   NumBuckets        = 10;
    constexpr int   vetPtrLabelSize   = 5;
    constexpr bool  NullFlush         = false;
    constexpr bool  Print             = false;
    constexpr int   NumMaxLabel       = 10000;
    constexpr bool  DominaIterBuckets = true;

    class LabelSet;

    //typedef std::set<LabelSet>::iterator LabelIt;

    struct Bound
    {
        FloatType lowerBound = -std::numeric_limits<FloatType>::infinity();
        FloatType upperBound = std::numeric_limits<FloatType>::infinity();
    };

    std::ostream& operator<< (std::ostream& out, const Bound &bound);

    // TODO Fix!
    // Fist access the resource, and then the cost of an arc (i,j)
    typedef Eigen::Array<Eigen::Matrix<FloatType, -1, -1, Eigen::RowMajor>, 1, -1> VetMatResCost;

    // Fist access the resource, and then the bound of a customer
    typedef Eigen::Vector<Eigen::Vector<Bound, -1>, -1> VetVetResBound;

    class NgSet
    {
    public:

        EigenMatrixRowI matNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);
        bool contain(int i, int j) const;
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
            return l0->vetResources[0] < l1->vetResources[0];
        }

    };


    std::ostream& operator<< (std::ostream& out, const Label &label);


    class Bucket
    {
    public:
        // Bound: [lower;upper)
        // Eigen::Vector<Bound, 2> vetBound;
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

        /** Access first the customer and after (i,j), where i is the component of the first resource
          *  and j the component of the second one.
          */
        Eigen::VectorX<MatBucket>                                       vetMatBucket;
        Eigen::Matrix<Eigen::Vector<Bound, 2>, -1, -1, Eigen::RowMajor> matBound;
        Eigen::Vector<Step, 2>                                          vetStepSize;
        Eigen::Vector<int, 2>                                           vetNumSteps;
        EigenMatrixRowI                                                 matBucketIndex; // Given a bucket index (i, j)

        int numMainResources;
        int numCust;
        int numMaxSteps;

        GraphNS::Graph<int> graphBucket;

        LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_,
                     int numMainResources_,
                     int numCust_);
        LabelingData()=default;

        void flushLabel();
        int getIndex(int resource, FloatType val);
        void removeLabel(Label *label);
        Label* getBestLabel(int cust);

        void checkMat();

        inline __attribute__((always_inline))
        int getIndexGraphBucket(int i, int j)
        {
            return i*vetNumSteps[1] + j;
        }


        // std::set<Label*, LabelCmp> &setLabel
        void dominanceInterBuckets(std::multiset<Label*, LabelCmp> &setLabel, int numRes, int localNumMaxLabel);//, MemoryPool_NS::Pool<Label> &poolTemp);
        void setupGraphBucket();

    };

    typedef std::multiset<Label*, LabelCmp>::iterator LabelSetIt;


    bool forwardLabelingAlgorithm(const int                     numRes,
                                  const int                     numCust,
                                  const VetMatResCost&          vetMatResCost,
                                  const VetVetResBound&         vetVetBound,
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

    inline __attribute__((always_inline))
    bool checkDominance(const Label& l0, const Label& l1, int numResources)
    {
        if(l0.cust != l1.cust)
        {
            std::cout<<"ERROR, lo.cust("<<l0.cust<<") != l1.cust("<<l1.cust<<")\n\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

        // Check the resources
        //#pragma GCC unroll NumMaxResources
        for(int i=0; i < NumMaxResources; ++i)
        {
            // l0.vetResources[i] > l1.vetResources[i]
            if(doubleLess(l1.vetResources[i], l0.vetResources[i], std::numeric_limits<FloatType>::epsilon()))
                return false;

            if((i+1) == numResources)
                break;
        }

        // TODO errado?
        //if(l0.bitSetNg == 0)
        //    return false;

        // Check if l0 is a subset of l1
        return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
    }

    bool extendLabel(const Label&          label,
                     Label&                newLabel,
                     const VetMatResCost&  vetMatResCost,
                     const VetVetResBound& vetVetBound,
                     int                   custI,
                     int                   custJ,
                     const NgSet&          ngSet,
                     int                   numResources);

    void removeCycles(Label &label, const int numCust);
    void removeCycles2(Label &label, const int numCust);
    void updateLabelCost(Label &label, const VetMatResCost &vetMatResCost, FloatType labelStart);

    // Linear index for a nxn matrix
    inline __attribute__((always_inline))
    int getIndex(int i, int j, int numClie){return i*numClie+j;}

    bool checkDistance(const Eigen::Matrix<FloatType, -1, -1, Eigen::RowMajor> &matDist);
    bool containRoute(const Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> &vetLabel, int numSol, Label* label);

    bool labelHaveRoute(std::vector<int> &vetRoute, Label *label);
    void checkDataStructs(Label* ptrLabel, LabelingData& lData, std::multiset<Label*, LabelCmp>& set);
    void eraseLabelFromSet(Label* ptrLabel, std::multiset<Label*, LabelCmp>& set);
    Label* dominanceIntraBucket(Label* label, Bucket &bucket, LabelSetIt& set);

}
#endif //DW_LABELINGALGORITHM_H
