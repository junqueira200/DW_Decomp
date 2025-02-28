/*  *****************************************************************
 *  *****************************************************************
 *  File:    LowerBound.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    24/02/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "LowerBound.h"

void LowerBoundNS::FloatTypeHeap::heapify(int i)
{

}

FloatType LowerBoundNS::FloatTypeHeap::extractMin()
{

}

void LowerBoundNS::FloatTypeHeap::decreaseKey(int i, FloatType val)
{

}

void LowerBoundNS::FloatTypeHeap::deleteKey(int i)
{

}

void LowerBoundNS::FloatTypeHeap::insertKey(FloatType val)
{

    if(heapSize == (int)vet.size())
    {
        vet.resize(2*heapSize);
        for(int i=heapSize; i < (int)vet.size(); ++i)
            vet[i] = 0;
    }

    heapSize += 1;
    int i = heapSize - 1;
    int iParent = parent(i);
    vet[i] = val;

    while(i!=0 && vet[iParent] > vet[i])//vet[iParent]->vetResources[0] > vet[i]->vetResources[0])
    {
        std::swap(vet[i], vet[iParent]);
        i = iParent;
        iParent = parent(i);
    }

}
