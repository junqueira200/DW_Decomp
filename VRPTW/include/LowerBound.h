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


    class FloatTypeHeap
    {
    public:

        Vector<FloatType> vet;
        int heapSize = 0;

        explicit FloatTypeHeap(int tam){vet = Vector<FloatType>(tam, 0.0);}
        int parent(int i) { return (i-1)/2;}
        int left(int i) { return (2*i + 1);}
        int right(int i) { return (2*i + 2);}
        bool empty(){return heapSize==0;}

        void heapify(int i);
        [[nodiscard]]FloatType extractMin();
        void decreaseKey(int i, FloatType val);
        [[nodiscard]]FloatType getMin(){return vet[0];}
        void deleteKey(int i);
        void insertKey(FloatType val);
    };

    [[nodiscard]]
    bool bellmanFord(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                     Eigen::VectorX<FloatType>&                vetDist,
                     int                                       src);

    [[nodiscard]]
    bool getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                           Eigen::VectorX<FloatType>&                vetDist);

}

#endif //DW_LOWERBOUND_H
