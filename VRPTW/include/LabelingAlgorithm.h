//
// Created by igor on 20/11/24.
//

#ifndef DW_LABELINGALGORITHM_H
#define DW_LABELINGALGORITHM_H

#include <Eigen/Eigen>
#include "Grafo.h"
#include <bitset>
#include "Aux.h"

namespace LabelingAlgorithmNS
{

    constexpr int NumMaxResources = 1;
    constexpr int NumMaxRoute     = 300;
    constexpr int NumMaxCust      = 150;
    constexpr int NgSetSize       = 5;
    constexpr int NumBuckets      = 10;
    constexpr int vetPtrLabelSize = 5;
    constexpr bool NullFlush      = true;


    struct Bound
    {
        double lowerBound = -std::numeric_limits<double>::infinity();
        double upperBound = std::numeric_limits<double>::infinity();
    };

    // TODO Fix!
    // Fist access the resource, and then the cost of an arc (i,j)
    typedef Eigen::Array<Eigen::Array<double, NumMaxCust, NumMaxCust>, 1, NumMaxResources> VetMatResCost;

    // Fist access the resource, and then the bound of a customer
    typedef Eigen::Array<Eigen::Array<Bound, 1, NumMaxCust>, 1, NumMaxResources> VetVetResBound;

    class NgSet
    {
    public:

        EigenMatrixRow matNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);
        bool contain(int i, int j) const;
        void setNgSets(const EigenMatrixRow &matDist);
    };

    struct Node
    {
        int i;
        double dist;

        bool operator<(const Node &node) const
        {
            return dist < node.dist;
        }
    };

    class Label
    {
    public:

        //int numResources = 1;
        Eigen::Array<double, 1, NumMaxResources> vetResources;
        Eigen::VectorXi vetRoute;
        int tamRoute = 0;
        std::bitset<NumMaxCust> bitSet;

        Label() = default;

    };


    class Bucket
    {
    public:
        // Bound: [lower;upper)
        Eigen::Vector<Bound, 2> vetBound;
        Eigen::VectorX<Label*> vetPtrLabel;
        int sizeVetPtrLabel;

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
    };

    class Step
    {
    public:

        double start;
        double end;
        double stepSize;
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

        Eigen::VectorX<MatBucket> vetMatBucket;
        Eigen::Vector<Step, 2> vetStepSize;
        Eigen::Vector<int, 2> vetNumSteps;
        int numMainResources;
        int numCust;
        int numMaxSteps;

        LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_,
                     int numMainResources_,
                     int numCust_);

        void flushLabel();
    };


    void forwardLabelingAlgorithm(const int numRes,
                                  const int numCust,
                                  const VetMatResCost& vetMatResCost,
                                  const VetVetResBound& vetVetBound,
                                  int dest,
                                  const NgSet &ngSet,
                                  LabelingData &lData);

    bool checkDominance(const Label& l0, const Label& l1, int numResources);

    bool extendLabel(const Label &label,
                     Label &newLabel,
                     const VetMatResCost& vetMatResCost,
                     const VetVetResBound& vetVetBound,
                     int custI,
                     int custJ,
                     const NgSet &ngSet,
                     int numResources);

}
#endif //DW_LABELINGALGORITHM_H
