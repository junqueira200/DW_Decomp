#ifndef LABEL_H
#define LABEL_H

#include "LabelingConstants.h"
#include "safe_vector.h"
#include <boost/array.hpp>
#include <Eigen/Eigen>
#include <bitset>
#include "Grafo.h"
#include "safe_matrix.h"
#include "Aux.h"

namespace LabelingAlgorithmNS
{
    enum TypeLabel
    {
        Forward,
        Backward
    };

    /// Label must be a FLAT data structure
    struct alignas(64) Label
    {
    public:


        bool        active    = false;
        TypeLabel   typeLabel = Forward;
        // First dimension for MatBucket
        int         i         = -1;
        // Second dimension for MatBucket
        int         j         = -1;
        int         cust      = -1;
        int         posHeap   = -1;
        int         tamRoute  = 0;
        // Position of the label from vetPtrLabel in Bucket class
        int         posBucket = -1;

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


    struct alignas(64)  Bucket
    {
    public:
        Eigen::VectorX<Label*> vetPtrLabel;
        int                    sizeVetPtrLabel = 0;

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
        std::string print(int numResorces);
        void removeElement(int i);
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
        Eigen::VectorX<MatBucket>       vetMatBucketForward;
        Eigen::VectorX<MatBucket>       vetMatBucketBackward;

    public:

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

        Bucket* getBucket(Label* label)
        {

            if(label->typeLabel == Forward)
                return &vetMatBucketForward[label->cust].mat(label->i, label->j);
            else
                return &vetMatBucketBackward[label->cust].mat(label->i, label->j);
        }

        // std::set<Label*, LabelCmp> &setLabel
        void dominanceInterBuckets(LabelHeap& labelHeap, int numRes, int localNumMaxLabel,
                                   Eigen::VectorX<MatBucket>& vetMatBucket, TypeLabel typeLabel);

        void setupGraphBucket();
        bool compareVetMatBucket(const ArrayResources& vetMaxResouces);

        inline __attribute__((always_inline))
        void whiteLabelIndex(Label* label)
        {
            label->i = getIndex(0, label->vetResources[0]);
            label->j = getIndex(1, label->vetResources[1]);
        }

        void checkVetMatBucketBackward();
        void checkVetMatBucketForward();
        Vector<std::pair<int,int>> getListOfIndexForMerge(const Label& label, const ArrayResources& vetMaxResouces,
                                                          int numResources);


    };

    bool searchLabel(Label* label, Bucket& bucket);

}

#endif // LABEL_H
