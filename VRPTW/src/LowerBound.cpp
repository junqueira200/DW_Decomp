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


using namespace LabelingAlgorithmNS;

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
                                     Eigen::VectorX<FloatType>&                vetDist,
                                     int                                       dest,
                                     LabelingAlgorithmNS::LabelingData*        lDataPtr)

{
    static LabelHeap labelHeap(1000000);
    startPool();

    int indexI = lDataPtr->getIndex(0, 0.0);
    vetDist.setConstant(MaxFloatType);

    //lDataPtr->vetMatBucketForward[0].mat(0, 1).flush();

    for(int cust=0; cust < dest; ++cust)
    {
        labelHeap.heapSize = 0;
        lDataPtr->flushLabel();
        resetLabelPool();

        Label* labelPtr = getLabel();

        labelPtr->typeLabel      = Forward;
        labelPtr->cust           = cust;
        labelPtr->tamRoute       = 1;
        labelPtr->vetRoute[0]    = cust;
        labelPtr->bitSetNg       = 0;
        labelPtr->bitSetNg[cust] = true;

        labelPtr->vetResources.setZero();
        labelPtr->i               = 0;
        labelPtr->j               = 0;
        labelPtr->posHeap         = 0;
        labelPtr->posBucket       = 0;

        Bucket* bucket = lDataPtr->getBucket(labelPtr);

        bucket->vetPtrLabel.resize(1);
        bucket->sizeVetPtrLabel = 1;
        bucket->vetPtrLabel[0]  = labelPtr;

        labelHeap.insertKey(labelPtr);

        vetDist[cust] = MaxFloatType;

        std::bitset<NumMaxCust> bitSet = 0;
        bitSet[cust] = true;

        while(!labelHeap.empty())
        {
            labelPtr = labelHeap.extractMin();
            //std::cout<<"cust: "<<labelPtr->cust<<"\n\n";
            lDataPtr->removeLabel(labelPtr);

            bitSet[labelPtr->cust] = true;

            for(int next=1; next <= dest; ++next)
            {
                if(labelPtr->bitSetNg[next] || next == labelPtr->cust || (labelPtr->cust == 0 && next == dest) ||
                        bitSet[next])
                    continue;

                Label* labelAux = getLabel();
                extendLabel(labelPtr, labelAux, vetMatResCost, lDataPtr, next);
                int correctPos = -1;
                Bucket* bucket = dominanceIntraBucket(next, labelAux, *lDataPtr, &labelHeap, 1, dest, correctPos);

                if(!bucket)
                {
                    rmLabel(labelAux);
                    continue;
                }

                if((bucket->sizeVetPtrLabel+1) > bucket->vetPtrLabel.size())
                    bucket->vetPtrLabel.conservativeResize(2*bucket->sizeVetPtrLabel);


                if(bucket->sizeVetPtrLabel != 0)
                {
                    for(int i = bucket->sizeVetPtrLabel-1; i >= correctPos; --i)
                    {
                        int index = i+1;
                        bucket->vetPtrLabel[index]            = bucket->vetPtrLabel[i];
                        bucket->vetPtrLabel[index]->posBucket = index;
                    }
                }

                bucket->vetPtrLabel[correctPos] = labelAux;
                bucket->vetPtrLabel[correctPos]->posBucket = correctPos;
                bucket->sizeVetPtrLabel += 1;

                if(next != dest)
                    labelHeap.insertKey(labelAux);            }
        }

        labelPtr = getLabel();
        labelPtr->cust = MaxInt;
        labelPtr->vetResources[0] = MaxFloatType;

        //for(int i=0; i <= lDataPtr->vetNumSteps(0); ++i)
        {
            Label* aux = getLabel();
            aux->cust = dest;
            aux->i = 0;
            aux->j = 0;
            aux->typeLabel = Forward;

            Bucket* bucket = lDataPtr->getBucket(aux);
            for(int j=0; j < bucket->sizeVetPtrLabel; ++j)
            {
                if(bucket->vetPtrLabel[j]->vetResources[0] < labelPtr->vetResources[0])
                    labelPtr = bucket->vetPtrLabel[j];
            }
        }

        vetDist[cust] = labelPtr->vetResources[0];
        //std::cout<<vetDist[cust]<<"\n";


    }

    std::cout<<"vetDist: "<<vetDist.transpose()<<"\n\n";
    return true;
}
