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

bool LowerBoundNS::bellmanFord(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                               Eigen::VectorX<FloatType>&                vetDist,
                               int                                       src)
{

    constexpr FloatType InfFloatType = std::numeric_limits<FloatType>::infinity();
    vetDist.setConstant(InfFloatType);
    vetDist[src] = 0.0;
    const int size = vetDist.size();

    for(int t=0; t < size; ++t)
    {
        // Relax the edges
        for(int i=0; i < size; ++i)
        {
            for(int j=0; j < size; ++j)
            {
                if(i == j || vetMatResCost(i, j, 0) == InfFloatType || vetDist[i] == InfFloatType)
                    continue;

                FloatType newDist = vetDist[i] + vetMatResCost(i, j, 0);
                if(newDist < vetDist[j])
                    vetDist[j] = newDist;
            }
        }
    }

    // Check if there is a negative cycle by relaxing the edges for the last time
    bool negCycle = false;

    for(int i=0; i < size; ++i)
    {
        for(int j=0; j < size; ++j)
        {
            if(i == j || vetMatResCost(i, j, 0) == InfFloatType || vetDist[i] == InfFloatType)
                continue;

            FloatType newDist = vetDist[i] + vetMatResCost(i, j, 0);
            if(newDist < vetDist[j])
            {
                negCycle = true;
                break;
            }
        }

        if(negCycle)
            break;
    }

    if(negCycle)
    {
        //std::cout<<"negCycle\n";
        return false;
    }

    //std::cout<<"getLB\n";

    return true;

}

/** ******************************************************************************************
 *  ******************************************************************************************
 *
 *  @param vetMatResCost
 *  @param vetDist          Returns the LB distance from the i-th customer to the last one
 *
 *  ******************************************************************************************
 *  ******************************************************************************************
 */
bool LowerBoundNS::getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                                     Eigen::VectorX<FloatType>&                vetDist)
{

    static Eigen::VectorX<FloatType> vetDistAux(vetDist);
    vetDist.setConstant(std::numeric_limits<FloatType>::infinity());
    vetDist[vetDist.size()-1] = 0.0;

    for(int i=0; i < (vetDist.size()-1); ++i)
    {
        if(!bellmanFord(vetMatResCost, vetDistAux, i))
            return false;

        vetDist[i] = vetDistAux[vetDist.size()-1];
    }


    return true;
}