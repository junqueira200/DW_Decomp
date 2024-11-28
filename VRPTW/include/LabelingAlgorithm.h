//
// Created by igor on 20/11/24.
//

#ifndef DW_LABELINGALGORITHM_H
#define DW_LABELINGALGORITHM_H

#include <Eigen/Eigen>
#include "Grafo.h"
#include <bitset>
#include "Aux.h"
#include "DW_Decomp.h"
#include "safe_vector.h"
#include <boost/array.hpp>

namespace LabelingAlgorithmNS
{

    constexpr int NumMaxResources = 2;
    constexpr int NumMaxRoute     = 300;
    constexpr int NumMaxCust      = 20;
    constexpr int NgSetSize       = 5;
    constexpr int NumBuckets      = 10;
    constexpr int vetPtrLabelSize = 5;
    constexpr bool NullFlush      = true;

    constexpr bool Print          = false;

    struct Bound
    {
        double lowerBound = -std::numeric_limits<double>::infinity();
        double upperBound = std::numeric_limits<double>::infinity();
    };

    std::ostream& operator<< (std::ostream& out, const Bound &bound);

    // TODO Fix!
    // Fist access the resource, and then the cost of an arc (i,j)
    typedef Eigen::Array<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>, 1, -1> VetMatResCost;

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
        boost::array<int, 200> vetRoute;
        int tamRoute = 0;
        std::bitset<NumMaxCust> bitSet;

        int i    = -1;
        int j    = -1;
        int cust = -1;
        int pos  = -1;

        Label() = default;

    };


    std::ostream& operator<< (std::ostream& out, const Label &label);


    class Bucket
    {
    public:
        // Bound: [lower;upper)
        // Eigen::Vector<Bound, 2> vetBound;
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

        void addLabel(Label *labelPtr);
//        bool delLabel(Label *labelPtr);

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

        /** Access first the customer and after (i,j), where i is the component of the first resource
          *  and j the component of the second one.
          */
        Eigen::VectorX<MatBucket> vetMatBucket;
        Eigen::Matrix<Eigen::Vector<Bound, 2>, -1, -1, Eigen::RowMajor> matBound;
        Eigen::Vector<Step, 2> vetStepSize;
        Eigen::Vector<int, 2> vetNumSteps;
        int numMainResources;
        int numCust;
        int numMaxSteps;

        LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_,
                     int numMainResources_,
                     int numCust_);
        LabelingData()=default;

        void flushLabel();
        int getIndex(int resource, double val);
        void removeLabel(Label *label);

        bool isCustEmpty(int cust)
        {
//std::cout<<"isCustEmpty\n";
//std::cout<<"vetNumSteps[1]("<<vetNumSteps[1]<<")\n";

            MatBucket &matBucket = vetMatBucket[cust];

            for(int i=0; i < vetNumSteps[0]; ++i)
            {
                for(int j=0; j < vetNumSteps[1]; ++j)
                {
//std::cout<<"before\n";
                    Bucket &bucket = matBucket.mat(i, j);
//std::cout<<"got bucket\n";
//std::cout<<"bucket.sizeVetPtrLabel("<<bucket.sizeVetPtrLabel<<")\n";

                    for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
                    {
                        // TODO: FIX
                        Label* label = bucket.vetPtrLabel[t];
                        if(label->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                            return false;
                    }
                }
            }

//std::cout<<"END\n";
            return true;

        }
    };


    bool forwardLabelingAlgorithm(const int numRes, const int numCust, const VetMatResCost &vetMatResCost,
                                  const VetVetResBound &vetVetBound, const int dest, const NgSet &ngSet,
                                  LabelingData &lData, Eigen::VectorXd &vetX, const double labelStart);

    bool checkDominance(const Label& l0, const Label& l1, int numResources);

    bool extendLabel(const Label &label,
                     Label &newLabel,
                     const VetMatResCost& vetMatResCost,
                     const VetVetResBound& vetVetBound,
                     int custI,
                     int custJ,
                     const NgSet &ngSet,
                     int numResources);

    void removeCycles(Label &label, const int numCust);
    void updateLabelCost(Label &label, const VetMatResCost &vetMatResCost);


    // Linear index for a nxn matrix
    inline __attribute__((always_inline))
    int getIndex(int i, int j, int numClie){return i*numClie+j;}
}
#endif //DW_LABELINGALGORITHM_H
