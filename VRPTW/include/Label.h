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
#include "NgSet.h"

namespace LabelingAlgorithmNS
{
    /// Stores the start indexes for the reduced cost and demand
    typedef Eigen::Array<int, 1, 2> IndexStart;
    /// Stores the end indexes for the reduced cost and demand
    typedef Eigen::Array<int, 1, 2> IndexEnd;

    struct Index
    {
        IndexStart start;
        IndexEnd   end;
    };

    enum TypeLabel
    {
        Forward,
        Backward
    };

    /// Label must be a FLAT data structure
    struct Label
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
        // Order of label's creation
        int64_t     index     = -1;

        //Eigen::Array<FloatType, 1, NumMaxResources> vetResources;
        std::array<FloatType, NumMaxResources>      vetResources;
        std::bitset<NumMaxCust>                     bitSetNg;
        boost::array<int, NumMaxRoute>              vetRoute;


        Label() = default;

    };

    INLINE
    void setVetResources0(Label& label)
    {
        memset(&label.vetResources[0], 0, sizeof(FloatType)*NumMaxResources);
    }

    class LabelCmp
    {
    public:

        INLINE
        bool operator()(Label *l0, Label *l1) const
        {
            //return doubleLess(l0->vetResources[0], l1->vetResources[0], std::numeric_limits<FloatType>::epsilon());
            if(MinHeap)
                return l0->HEAP_KEY < l1->HEAP_KEY;// && l0->vetResources[1] < l1->vetResources[1];
            else
                return l0->HEAP_KEY > l1->HEAP_KEY;
        }

        INLINE
        bool isGreater(Label *l0, Label *l1)const
        {
            if(MinHeap)
                return l0->HEAP_KEY > l1->HEAP_KEY;
            else
                return l0->HEAP_KEY < l1->HEAP_KEY;
        }

    };

    class LabelHeap
    {
    public:

        Vector<Label*> vet;
        int heapSize = 0;

        explicit LabelHeap(int tam){vet = Vector<Label*>(tam, nullptr);}
        INLINE int parent(int i) { return (i-1)/2;}
        INLINE int left(int i) { return (2*i + 1);}
        INLINE int right(int i) { return (2*i + 2);}
        INLINE bool empty(){return heapSize==0;}

        void heapify(int i);
        [[nodiscard]]Label* extractTop();
        void decreaseKey(int i, KEY_TYPE val);
        [[nodiscard]]Label* getTop(){return vet[0];}
        void deleteKey(int i);
        void insertKey(Label* label);
    };


    std::ostream& operator<< (std::ostream& out, const Label &label);


    struct Bucket
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
        void addElement(int pos, Label* labelPtr);
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

        /// Each cell has 2 index for start and end. The index are valid for every
        ///  label(forward) make a merge in vetMatBucketBackward
        Eigen::MatrixX<Index>           matForwardRange;

         /// Each cell has 2 index for start and end. The index are valid for every
        ///  label(backward) make a merge in vetMatBucketForward
        Eigen::MatrixX<Index>           matBackwardRange;

        Eigen::Vector<Step, 2>          vetStepSize;

        Eigen::Vector<int, 2>           vetNumSteps;

        /// Given a bucket index (i, j)
        EigenMatrixRowI                 matBucketIndex;

        int numMainResources;

        /// Have a copy of the deposit
        int numCust;
        int numMaxSteps;

        GraphNS::Graph<int> graphBucket;

        ArrayResources vetMaxResources;

        EigenMatrixRowD matDist;

        LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_, int numMainResources_, int numCust_,
                     const ArrayResources& vetMaxResources_);
        LabelingData()=default;

        void flushLabel();
        int getIndex(int resource, FloatType val);
        void removeLabel(Label *label);
        Label* getBestLabel(int cust);

        INLINE
        int getIndexGraphBucket(int i, int j)
        {
            return i*vetNumSteps[1] + j;
        }

        INLINE
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

        INLINE
        void whiteLabelIndex(Label* label)
        {
            label->i = getIndex(0, label->vetResources[0]);
            label->j = getIndex(1, label->vetResources[1]);
        }

        void checkVetMatBucketBackward();
        void checkVetMatBucketForward();
        Index getListOfIndexForMerge(const Label& label);

        int doMerge(Label* label, const ArrayResources& vetMaxResources, const MatBoundRes& vetVetBound, int numResorces,
                    const NgSet& ngSet);

        INLINE
        int getSecondDeposit(){return numCust-1;}

        INLINE
        int getNumberOfSolutions(){ return vetMatBucketForward[getSecondDeposit()].mat(0,0).sizeVetPtrLabel;}

        void checkLabels();

    };

    bool searchLabel(Label* label, Bucket& bucket);
    std::string printIndex(const Index& index);

    INLINE
    bool labelLessForward(Label* label0, Label* label1, int numResorces)
    {
        int numEqual = 0;
        for(int i=0; i < numResorces; ++i)
        {
            if(label0->vetResources[i] > label1->vetResources[i])
                return false;
            if(doubleEqual(label0->vetResources[i], label1->vetResources[i], 1E-5))
                numEqual += 1;
        }

        if(numEqual == numResorces)
            return false;

        return true;
    }

    INLINE
    bool labelEqual(Label* label0, Label* label1, int numResorces)
    {
        for(int i=0; i < numResorces; ++i)
        {
            if(!doubleEqual(label0->vetResources[i], label1->vetResources[i], 1E-5))
                return false;
        }

        return true;
    }

}

#endif // LABEL_H
