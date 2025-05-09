/*  *****************************************************************
 *  *****************************************************************
 *  File:    LabelingAlgorithm.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    20/11/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "LabelingAlgorithm.h"
#include "MemoryPool.h"
#include <list>
#include <algorithm>
#include "LowerBound.h"

using namespace LabelingAlgorithmNS;
using namespace LowerBoundNS;

LabelingAlgorithmNS::NgSet::NgSet()
{
    matNgSet = EigenMatrixRowI(NumMaxCust, NgSetSize);//Eigen::MatrixXi(NumMaxCust, NgSetSize);
    matNgSet.setConstant(-1);
}

// Checks if j is in the ng-set of i; TODO: busca binaria
bool LabelingAlgorithmNS::NgSet::contain(int i, int j) const
{
    if(!active)
        return true;

    for(int jj=0; jj < NgSetSize; ++jj)
    {
        if(matNgSet(i, jj) == j)
            return true;

        if((jj+1) == ngSetSize)
            break;
    }

    return false;
}

LabelingAlgorithmNS::NgSet::NgSet(int numCust_, int ngSetSize_)
{
    numCust   = numCust_;
    ngSetSize = ngSetSize_;

    matNgSet = Eigen::MatrixXi(numCust, ngSetSize);
}

void LabelingAlgorithmNS::NgSet::setNgSets(const EigenMatrixRowD &matDist)
{

    Eigen::VectorX<Node> vet(numCust-1);

    for(int i=0; i < numCust; ++i)
    {
        //std::cout<<"i("<<i<<")\n";
        int next = 0;
        for(int j=0; j < numCust; ++j)
        {
            if(i == j)
                continue;

            vet[next].dist = matDist(i, j);
            vet[next].i    = j;
            next += 1;

            //std::cout<<"("<<j<<","<<matDist(i,j)<<"); ";
        }


        std::sort(vet.begin(), vet.end());
        next = 0;

        //std::cout<<"\tsort:\n";

        for(int j=0; j < ngSetSize; ++j)
        {
            matNgSet(i, j) = vet[j].i;
            //std::cout<<"("<<vet[j].i<<","<<vet[j].dist<<"); ";
        }

        if((i+1) == numCust)
            break;
    }

}


inline MemoryPool_NS::Pool<LabelingAlgorithmNS::Label> labelPoolG;
inline bool startLabelPoolG = false;


/** **********************************************************************************************************
 *  **********************************************************************************************************
 *
 * @brief Executes the forward labeling algorithm for resource-constrained path enumeration.
 *
 * @details This function performs label generation and extension from the origin/destination to
 * the destination/origin node, applying resource constraints, dominance checks. It is
 * used in column generation methods for routing problems.
 *
 * @param numRes                Number of resources (e.g., time, capacity).
 * @param numCust               Number of customers or nodes in the problem.
 * @param vetMatResCost         3D matrix of resource costs between nodes.
 * @param vetVetBound           Matrix of upper/lower bounds for resource usage.
 * @param dest                  Index of the destination node.
 * @param ngSet                 Ng set to enforce route feasibility constraints.
 * @param lData                 Reference to data structures used in labeling (e.g., buckets, label pools).
 * @param matColX               Output matrix to store generated columns (paths).
 * @param numSol                Reference to an integer tracking the number of generated solutions.
 * @param labelStart            Initial cost or reduced cost value to begin labeling with.
 * @param NumMaxLabePerBucket   Maximum number of labels allowed per bucket.
 * @param dominaceCheck         Whether to apply dominance checks during labeling.
 * @param maxDist               Reference to a float storing the maximum distance (or cost) observed.
 * @param vetRedCost            Vector of reduced costs per customer, used for pruning.
 * @param exact                 If true, applies exact pricing (e.g., for exact solution of the subproblem).
 *
 * @return                      Returns \c true if at least one feasible solution (column) was found;
 *                                  otherwise, returns \c false.
 *
 * **********************************************************************************************************
 * **********************************************************************************************************
 */
bool
LabelingAlgorithmNS::forwardLabelingAlgorithm(const int                     numRes,
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
                                              Eigen::VectorX<FloatType>&    vetRedCost,
                                              bool                          exact)
{
    //static Eigen::VectorX<FloatType> vetLowerBoundDist(numCust);
    //const bool vetLowerBoundDistValid =  getDistLowerBound(vetMatResCost, vetLowerBoundDist);

    matColX.setZero();
    if(vetRoteG.size() > 0)
    {
        std::cout<<"vetRoteG: ";
        for(int i:vetRoteG)
            std::cout<<i<<" ";
        std::cout<<"\n";

        std::cout<<"*********************LABELING*********************\n\n";
    }

    //checkDistance(vetMatResCost[0]);

    //dominaceCheck = false;

    if(NumMaxLabePerBucket == -1)
        NumMaxLabePerBucket = std::numeric_limits<int>::max();

    //std::cout<<"NumMaxLabePerBucket: "<<NumMaxLabePerBucket<<"\n\n";

    numSol = 0;
    static Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> vetLabel;
    vetLabel.setZero();

    maxDist = -std::numeric_limits<FloatType>::max();
    //redCost = std::numeric_limits<FloatType>::max();

    if(Print)
    {
        std::cout << "************************************************************\n";
        std::cout << "*****************FORWARD LABELING ALGORITHM*****************\n\n";
        std::cout << "numCust: " << numCust << "\n";
    }
    //static bool labelPoolStart = false;//MemoryPool_NS::Pool<Label> labelPool(44, 400);
    if(!startLabelPoolG)
    {
        startLabelPoolG = true;
        labelPoolG.startPool(5, 100);
    }

    labelPoolG.resetPool(false);
    lData.flushLabel();

    static LabelHeap labelHeap(10000);
    labelHeap.vet.setAll(nullptr);
    labelHeap.heapSize = 0;



    Label *labelPtr = labelPoolG.getT();
    if(labelPtr == nullptr)
    {
        std::cout << "Pool.getT return nullptr!";
        PRINT_DEBUG("", "");
        throw "ERRO";
    }

    //labelPtr->vetRoute.resize(10);
    labelPtr->tamRoute = 1;

    labelPtr->vetRoute[0] = 0;
    labelPtr->vetResources.setZero();
    labelPtr->vetResources[0] = labelStart;
    labelPtr->bitSetNg = 0;

    if(Print)
    {   std::cout << "\nsetup label0\n";
        std::cout << *labelPtr << "\n";
    }

    int i = lData.getIndex(0, labelStart);
    int j = lData.getIndex(1, 0.0);

    labelPtr->i = i;
    labelPtr->j = j;
    labelPtr->cust = 0;
    labelPtr->active = true;
    labelPtr->posBucket = 0;

    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel[0]  = labelPtr;
    lData.vetMatBucketForward[0].mat(i, j).sizeVetPtrLabel = 1;

    //listLabel.push_back(labelPtr);

    //(*(LabelSetIt*)labelPtr->it) = setLabel.insert(labelPtr);
    labelHeap.insertKey(labelPtr);
    labelPtr = nullptr;

    int numIt = 0;

    Label *labelPtrBest = nullptr;
    int maxSize = 0;
    int maxSizeVetPtrLabel = 0;
    int localNumMaxLabel = numMaxLabelG;

    Label* ptrLabelTarget = nullptr;
    bool printSize = false;

    while(!labelHeap.empty() &&
         ((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel < DW_DecompNS::NumMaxSolSubProb) || exact))
    {
        if((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel >= 10 && labelHeap.heapSize >= 25 * numCust &&
            !exact))
            break;

        // TODO remover
        //lData.checkMat();

        if(ptrLabelTarget)
        {
            std::cout<<"ptrLabelTarget("<<ptrLabelTarget->active<<"): "<<*ptrLabelTarget<<"\n\n";
        }

        if(labelHeap.heapSize > localNumMaxLabel && DominaIterBuckets)
        {
            lData.dominanceInterBuckets(labelHeap, numRes, localNumMaxLabel, lData.vetMatBucketForward);


            if(labelHeap.heapSize > localNumMaxLabel)
            {
                int max = std::max(labelHeap.heapSize, localNumMaxLabel);
                localNumMaxLabel =  (int)(max * 1.05);
            }


            if(exactLabelingG)
            {
                std::cout << "heapSize: " << labelHeap.heapSize << "; maxHeapSize:" << localNumMaxLabel << "; maxDist:"
                          << maxDistG << "; minDist: " << minDistG << "\n";
            }

        }


        maxSize = std::max(maxSize, labelHeap.heapSize);


        if(Print)
            std::cout << "numIt: " << numIt << "\n";

        labelPtr = labelHeap.extractMin();
        if(labelPtr == nullptr)
        {

            std::cout<<"listLabel have a nullptr\n";
            PRINT_DEBUG("", "");
            throw "ERRO";
        }

        if(Print || labelHaveRoute(vetRoteG, labelPtr) || exactLabelingG)
        {
            std::cout << "labelPtr: " << labelPtr << "\n";
            std::cout << *labelPtr << "\n\n";
        }


        //std::cout<<"Extract label("<<labelPtr->vetRoute<<")\n";
        int lastCust = labelPtr->cust;

        if(!labelPtr->active)
        {
            PRINT_DEBUG("", "");
            throw "ERROR";

            labelPoolG.delT(labelPtr);
            continue;
        }

        if(lastCust == dest)
        {
            if(Print)
                std::cout<<"lastCust==dest("<<dest<<")\n";
            continue;
        }

        // Remove labelPtr from vetMatBucket
        lData.removeLabel(labelPtr);

        // Extend label
        for(int t=1; t < numCust; ++t)
        {
            //std::cout<<"t: "<<t<<"\n";
            if(t == lastCust || vetMatResCost(labelPtr->cust, t, 0) == std::numeric_limits<FloatType>::infinity())
                continue;

            if(labelPtr->bitSetNg[t] == 1)
                continue;


            if(Print)
                std::cout<<"\tt("<<t<<")\n";

            Label *labelPtrAux = labelPoolG.getT();
            if(labelPtrAux == nullptr)
            {
                std::cout<<"Pool.getT return nullptr!";
                PRINT_DEBUG("", "");
                throw "ERRO";
            }


            /*
            if(labelPtrAux->it == nullptr)
            {
                LabelSetIt* labelSetIt = new LabelSetIt;
                labelPtrAux->it = (void*)labelSetIt;
            }
            */


            if(extendLabel(*labelPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelPtr->cust, t, ngSet, numRes))
            {
                int correctPos = 0;

                if(t == dest)
                {
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    removeCycles2(*labelPtrAux, numCust);
                    updateLabelCost(*labelPtrAux, vetMatResCost, labelStart);


                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }
                }

                Bucket* bucket = dominanceIntraBucket(t, labelPtrAux, lData, labelHeap, numRes, dest, correctPos);

                if(!bucket)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if((bucket->sizeVetPtrLabel+1) > NumMaxLabePerBucket && t != dest)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if(t == dest && labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                {

                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                int size = bucket->sizeVetPtrLabel+1;
                if(size > bucket->vetPtrLabel.size())
                {
                    bucket->vetPtrLabel.conservativeResize(2*size);
                }

                if(bucket->sizeVetPtrLabel != 0)
                {
                    for(int i = bucket->sizeVetPtrLabel-1; i >= correctPos; --i)
                    {
                        int index = i+1;
                        bucket->vetPtrLabel[index]            = bucket->vetPtrLabel[i];
                        bucket->vetPtrLabel[index]->posBucket = index;
                    }

                }

                bucket->vetPtrLabel[correctPos]            = labelPtrAux;
                bucket->vetPtrLabel[correctPos]->posBucket = correctPos;
                bucket->sizeVetPtrLabel                   += 1;

                if(t != dest)
                    labelHeap.insertKey(labelPtrAux);


            }
            else
                labelPoolG.delT(labelPtrAux);

        }

//std::cout<<"END Extend label\n";

        if(Print)
            std::cout<<"\n########################################################################################\n\n";

        labelPoolG.delT(labelPtr);
        numIt += 1;
    }


    if(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel > 0)
    {
        numSol = std::min(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel, DW_DecompNS::NumMaxSolSubProb);
        for(int l=0; l < numSol; ++l)
        {
            Label* label = lData.vetMatBucketForward[dest].mat(0, 0).vetPtrLabel[l];
            label->vetRoute[label->tamRoute-1] = label->vetRoute[0];
            auto &vetRoute = label->vetRoute;


            for(int i = 0; i < (label->tamRoute - 1); ++i)
            {
                int index = getIndex(vetRoute[i], vetRoute[i+1], numCust-1);
                matColX(index, l) = 1.0;
            }

            vetRedCost[l] = label->vetResources[0];
        }

        return true;
    }
    else
    {
        return false;
    }
}


/** **********************************************************************************************************
 *  **********************************************************************************************************
 *
 * @brief Executes the forward labeling algorithm for resource-constrained path enumeration.
 *
 * @details This function performs label generation and extension from the origin to the destination
 * node, applying resource constraints, dominance checks. It is
 * used in column generation methods for routing problems.
 *
 * @param numRes                Number of resources (e.g., time, capacity).
 * @param numCust               Number of customers or nodes in the problem.
 * @param vetMatResCost         3D matrix of resource costs between nodes.
 * @param vetVetBound           Matrix of upper/lower bounds for resource usage.
 * @param dest                  Index of the destination node.
 * @param ngSet                 Ng set to enforce route feasibility constraints.
 * @param lData                 Reference to data structures used in labeling (e.g., buckets, label pools).
 * @param matColX               Output matrix to store generated columns (paths).
 * @param numSol                Reference to an integer tracking the number of generated solutions.
 * @param labelStart            Initial cost or reduced cost value to begin labeling with.
 * @param NumMaxLabePerBucket   Maximum number of labels allowed per bucket.
 * @param dominaceCheck         Whether to apply dominance checks during labeling.
 * @param maxDist               Reference to a float storing the maximum distance (or cost) observed.
 * @param vetRedCost            Vector of reduced costs per customer, used for pruning.
 * @param exact                 If true, applies exact pricing (e.g., for exact solution of the subproblem).
 *
 * @return                      Returns \c true if at least one feasible solution (column) was found;
 *                                  otherwise, returns \c false.
 *
 * **********************************************************************************************************
 * **********************************************************************************************************
 */
bool
LabelingAlgorithmNS::bidirectionalAlgorithm(const int                     numRes,
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
                                            Eigen::VectorX<FloatType>&    vetRedCost,
                                            bool                          exact)
{
    //static Eigen::VectorX<FloatType> vetLowerBoundDist(numCust);
    //const bool vetLowerBoundDistValid =  getDistLowerBound(vetMatResCost, vetLowerBoundDist);

    matColX.setZero();
    if(vetRoteG.size() > 0)
    {
        std::cout<<"vetRoteG: ";
        for(int i:vetRoteG)
            std::cout<<i<<" ";
        std::cout<<"\n";

        std::cout<<"*********************LABELING*********************\n\n";
    }

    //checkDistance(vetMatResCost[0]);

    //dominaceCheck = false;

    if(NumMaxLabePerBucket == -1)
        NumMaxLabePerBucket = std::numeric_limits<int>::max();

    //std::cout<<"NumMaxLabePerBucket: "<<NumMaxLabePerBucket<<"\n\n";

    numSol = 0;
    static Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> vetLabel;
    vetLabel.setZero();

    maxDist = -std::numeric_limits<FloatType>::max();
    //redCost = std::numeric_limits<FloatType>::max();

    if(Print)
    {
        std::cout << "************************************************************\n";
        std::cout << "*****************FORWARD LABELING ALGORITHM*****************\n\n";
        std::cout << "numCust: " << numCust << "\n";
    }

    if(!startLabelPoolG)
    {
        startLabelPoolG = true;
        labelPoolG.startPool(5, 100);
    }

    labelPoolG.resetPool(false);
    lData.flushLabel();

    static LabelHeap labelHeap(10000);
    labelHeap.vet.setAll(nullptr);
    labelHeap.heapSize = 0;

    Label* labelForwardPtr  = labelPoolG.getT();
    Label* labelBackwardPtr = labelPoolG.getT();

    //labelPtr->vetRoute.resize(10);
    labelForwardPtr->tamRoute = 1;
    labelForwardPtr->vetRoute[0] = 0;
    labelForwardPtr->vetResources.setZero();
    labelForwardPtr->vetResources[0] = labelStart;
    labelForwardPtr->bitSetNg = 0;
    labelForwardPtr->bitSetNg[0] = true;

    int i = lData.getIndex(0, labelStart);
    int j = lData.getIndex(1, 0.0);

    labelForwardPtr->i         = i;
    labelForwardPtr->j         = j;
    labelForwardPtr->cust      = 0;
    labelForwardPtr->active    = true;
    labelForwardPtr->posBucket = 0;

    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel[0]  = labelForwardPtr;
    lData.vetMatBucketForward[0].mat(i, j).sizeVetPtrLabel = 1;

    //labelBackwardPtr

    labelBackwardPtr->tamRoute        = 1;
    labelBackwardPtr->vetRoute[0]     = dest;
    labelBackwardPtr->typeLabel       = backward;
    labelBackwardPtr->vetResources.setZero();
    labelBackwardPtr->vetResources[0] = labelStart;
    labelBackwardPtr->bitSetNg        = 0;
    labelBackwardPtr->bitSetNg[dest]  = true;
    labelBackwardPtr->i               = i;
    labelBackwardPtr->j               = j;
    labelBackwardPtr->cust            = dest;
    labelBackwardPtr->active          = true;
    labelBackwardPtr->posBucket       = 0;


    /*
    lData.vetMatBucketBackward[dest].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucketBackward[dest].mat(i, j).vetPtrLabel[0]  = labelBackwardPtr;
    lData.vetMatBucketBackward[dest].mat(i, j).sizeVetPtrLabel = 1;
    labelHeap.insertKey(labelBackwardPtr);
    */


    labelHeap.insertKey(labelForwardPtr);
    labelForwardPtr = nullptr;
    labelBackwardPtr = nullptr;

    int numIt = 0;

    Label *labelPtrBest = nullptr;
    int maxSize = 0;
    int maxSizeVetPtrLabel = 0;
    int localNumMaxLabel = numMaxLabelG;

    Label* ptrLabelTarget = nullptr;
    bool printSize = false;

    while(!labelHeap.empty() &&
          ((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel < DW_DecompNS::NumMaxSolSubProb) || exact))
    {
        if((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel >= 10 && labelHeap.heapSize >= 25 * numCust &&
            !exact))
            break;

        if(ptrLabelTarget)
        {
            std::cout<<"ptrLabelTarget("<<ptrLabelTarget->active<<"): "<<*ptrLabelTarget<<"\n\n";
        }

        if(labelHeap.heapSize > localNumMaxLabel && DominaIterBuckets)
        {
            lData.dominanceInterBuckets(labelHeap, numRes, localNumMaxLabel, lData.vetMatBucketForward);



            if(labelHeap.heapSize > localNumMaxLabel)
            {
                int max = std::max(labelHeap.heapSize, localNumMaxLabel);
                localNumMaxLabel =  (int)(max * 1.05);
            }


            if(exactLabelingG)
            {
                std::cout << "heapSize: " << labelHeap.heapSize << "; maxHeapSize:" << localNumMaxLabel << "; maxDist:"
                          << maxDistG << "; minDist: " << minDistG << "\n";
            }

        }


        maxSize = std::max(maxSize, labelHeap.heapSize);


        if(Print)
            std::cout << "numIt: " << numIt << "\n";

        labelForwardPtr = labelHeap.extractMin();
        if(labelForwardPtr == nullptr)
        {

            std::cout<<"listLabel have a nullptr\n";
            PRINT_DEBUG("", "");
            throw "ERRO";
        }

        if(Print || labelHaveRoute(vetRoteG, labelForwardPtr) || exactLabelingG)
        {
            std::cout << "labelPtr: " << labelForwardPtr << "\n";
            std::cout << *labelForwardPtr << "\n\n";
        }


        //std::cout<<"Extract label("<<labelPtr->vetRoute<<")\n";
        int lastCust = labelForwardPtr->cust;

        if(!labelForwardPtr->active)
        {
            PRINT_DEBUG("", "");
            throw "ERROR";

            labelPoolG.delT(labelForwardPtr);
            continue;
        }

        if(lastCust == dest && labelForwardPtr->typeLabel == forward)
        {
            if(Print)
                std::cout<<"lastCust==dest("<<dest<<")\n";
            continue;
        }

        if(lastCust == 0 && labelForwardPtr->typeLabel == backward)
            continue;

        // Remove labelPtr from vetMatBucket
        lData.removeLabel(labelForwardPtr);

        // Extend label
        for(int t=0; t < numCust; ++t)
        {


            int tAux = t;
            //std::cout<<"t: "<<t<<"\n";
            if(labelForwardPtr->typeLabel  == forward &&
               vetMatResCost(labelForwardPtr->cust, t, 0) == std::numeric_limits<FloatType>::infinity())
                continue;

            else if(labelForwardPtr->typeLabel  == backward &&
               vetMatResCost(t, labelForwardPtr->cust, 0) == std::numeric_limits<FloatType>::infinity())
                continue;


            if(labelForwardPtr->bitSetNg[t] == 1)
                continue;

            if(t == labelForwardPtr->cust)
                continue;


            if(Print)
                std::cout<<"\tt("<<t<<")\n";

            Label *labelPtrAux = labelPoolG.getT();
            if(labelPtrAux == nullptr)
            {
                std::cout<<"Pool.getT return nullptr!";
                PRINT_DEBUG("", "");
                throw "ERRO";
            }


            /*
            if(labelPtrAux->it == nullptr)
            {
                LabelSetIt* labelSetIt = new LabelSetIt;
                labelPtrAux->it = (void*)labelSetIt;
            }
            */

            if(extendLabel(*labelForwardPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelForwardPtr->cust, t, ngSet,
                           numRes))
            {
                int correctPos = 0;

                if(t == dest && labelPtrAux->typeLabel == forward)
                {
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    removeCycles2(*labelPtrAux, numCust);
                    updateLabelCost(*labelPtrAux, vetMatResCost, labelStart);


                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }
                }


                if(t == 0 && labelPtrAux->typeLabel == backward)
                {
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {
                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    removeCycles2(*labelPtrAux, numCust);

                    // Reverse the array
                    std::reverse(labelPtrAux->vetRoute.begin(),
                                 labelPtrAux->vetRoute.begin()+labelPtrAux->tamRoute);


                    updateLabelCost(*labelPtrAux, vetMatResCost, labelStart);


                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    //change the backward type to forward
                    labelPtrAux->typeLabel = forward;
                    labelPtrAux->cust      = dest;
                    tAux = dest;

                }

                Bucket* bucket = dominanceIntraBucket(tAux, labelPtrAux, lData, labelHeap, numRes, dest, correctPos);

                if(!bucket)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if((bucket->sizeVetPtrLabel+1) > NumMaxLabePerBucket && t != dest && t != 0)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }


                int size = bucket->sizeVetPtrLabel+1;
                if(size > bucket->vetPtrLabel.size())
                {
                    bucket->vetPtrLabel.conservativeResize(2*size);
                }

                if(bucket->sizeVetPtrLabel != 0)
                {
                    for(int i = bucket->sizeVetPtrLabel-1; i >= correctPos; --i)
                    {
                        int index = i+1;
                        bucket->vetPtrLabel[index]            = bucket->vetPtrLabel[i];
                        bucket->vetPtrLabel[index]->posBucket = index;
                    }

                }

                bucket->vetPtrLabel[correctPos]            = labelPtrAux;
                bucket->vetPtrLabel[correctPos]->posBucket = correctPos;
                bucket->sizeVetPtrLabel                   += 1;

                if(t != dest)
                    labelHeap.insertKey(labelPtrAux);


            }
            else
                labelPoolG.delT(labelPtrAux);

        }

//std::cout<<"END Extend label\n";

        if(Print)
            std::cout<<"\n########################################################################################\n\n";

        labelPoolG.delT(labelForwardPtr);
        numIt += 1;
    }


    if(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel > 0)
    {
        numSol = std::min(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel, DW_DecompNS::NumMaxSolSubProb);
        for(int l=0; l < numSol; ++l)
        {
            Label* label = lData.vetMatBucketForward[dest].mat(0, 0).vetPtrLabel[l];
            label->vetRoute[label->tamRoute-1] = label->vetRoute[0];
            auto &vetRoute = label->vetRoute;


            for(int i = 0; i < (label->tamRoute - 1); ++i)
            {
                int index = getIndex(vetRoute[i], vetRoute[i+1], numCust-1);
                matColX(index, l) = 1.0;
            }

            vetRedCost[l] = label->vetResources[0];
        }

        return true;
    }
    else
    {
        return false;
    }
}

/** *********************************************************************************************
 *  *********************************************************************************************
 *
 * @brief Extends an existing label to create a new one.
 *
 *
 * @details This function attempts to extend an existing label by computing the new label's
 * resource consumption based on given resource cost matrices, bounds, and
 * customer indices. The ng-set ...
 *
 * @param label             The current label to be extended.
 * @param newLabel          Reference to the label that will store the result of the extension.
 * @param vetMatResCost     The 3D matrix containing resource costs.
 * @param vetVetBound       The matrix containing resource bounds.
 * @param custI             The index of the current customer.
 * @param custJ             The index of the next customer.
 * @param ngSet             The set representing non-neighbor restrictions.
 * @param numResources      The number of resources considered.
 *
 * @return Returns \c true if the label was successfully extended, otherwise returns \c false.
 *
 * *********************************************************************************************
 * *********************************************************************************************
 */
bool LabelingAlgorithmNS::extendLabel(const Label&          label,
                                      Label&                newLabel,
                                      const Vet3D_ResCost&  vetMatResCost,
                                      const MatBoundRes&    vetVetBound,
                                      int                   custI,
                                      int                   custJ,
                                      const NgSet&          ngSet,
                                      const int             numResources)
{


    // Goes through resources
    bool boundOk = true;

    if(label.typeLabel == backward)
        std::swap(custI, custJ);

    newLabel.typeLabel = label.typeLabel;

    for(int i=0; i < numResources; ++i)//NumMaxResources; ++i)
    {

        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] + vetMatResCost(custI, custJ, i);
        boundOk = boundOk && (newLabel.vetResources[i] <= vetVetBound(custJ, i).upperBound);
    }

    if(!boundOk)
        return false;

    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ its in custI ngSet
    if(ngSet.contain(custJ, custI))
    {
        newLabel.bitSetNg[custJ] = true;
    }


    for(int i=0; i < label.tamRoute; ++i)
        newLabel.vetRoute[i] = label.vetRoute[i];


    if(Print)
        std::cout<<"APOS COPY ROUTE\n";

    newLabel.vetRoute[label.tamRoute] = custJ;
    newLabel.tamRoute = label.tamRoute+1;

    return true;

}

/** *********************************************************************************************
 *  *********************************************************************************************
 *
 * @param vetStepSize_
 * @param numMainResources_
 * @param numCust_              Number of customers of the instance plus one
 *
 *  *********************************************************************************************
 *  *********************************************************************************************
 */
LabelingAlgorithmNS::LabelingData::LabelingData(const Eigen::Vector<Step, 2>& vetStepSize_,
                                                int                           numMainResources_,
                                                int                           numCust_)
{
    vetStepSize      = vetStepSize_;
    numMainResources = std::min(numMainResources_, 2);
    numCust          = numCust_;
    vetNumSteps      = {1, 1};

    // Calculates the number of steps of each resource
    for(int r=0; r < numMainResources; ++r)
    {
        FloatType resource = vetStepSize[r].start;
        int step = 1;

        while(resource <= vetStepSize[r].end)
        {
            resource += vetStepSize[r].stepSize;
            step += 1;
        }

        if(step > 1)
            step -= 1;

        vetNumSteps[r] = step;

//std::cout<<"Resource("<<r<<") have ("<<step<<") steps\n";

    }

    numMaxSteps = std::max(vetNumSteps[0], vetNumSteps[1]);

    vetMatBucketForward = Eigen::VectorX<MatBucket>(numCust);//(numMaxSteps, numMaxSteps);
    vetMatBucketBackward = Eigen::VectorX<MatBucket>(numCust);
    for(int i=0; i < numCust; ++i)
    {
        vetMatBucketForward[i].mat = Eigen::Matrix<Bucket, -1, -1, Eigen::RowMajor>(vetNumSteps[0], vetNumSteps[1]);
        vetMatBucketBackward[i].mat = Eigen::Matrix<Bucket, -1, -1, Eigen::RowMajor>(vetNumSteps[0], vetNumSteps[1]);
    }


    vetMatBound = Vector<Matrix<Bound, false>>(2);

    // Generates matBound
    //matBound = Eigen::Matrix<Eigen::Vector<Bound, 2>, -1, -1, Eigen::RowMajor>(vetNumSteps[0], vetNumSteps[1]);
    vetMatBound[0] = Matrix<Bound, false>(vetNumSteps[0], vetNumSteps[1]);
    vetMatBound[1] = Matrix<Bound, false>(vetNumSteps[0], vetNumSteps[1]);

    for(int i=0; i < vetNumSteps[0]; ++i)
    {
        FloatType r00 = vetStepSize[0].start + i*vetStepSize[0].stepSize;

        FloatType r01 = r00 + vetStepSize[0].stepSize;
        for(int j=0; j < vetNumSteps[1]; ++j)
        {
            FloatType r10 = vetStepSize[1].start + j*vetStepSize[1].stepSize;
            FloatType r11 = r10 + vetStepSize[1].stepSize;

            if(i != 0)
                //matBound(i, j)[0].lowerBound = r00;
                vetMatBound[0].get(i, j).lowerBound = r00;
            else
                //matBound(i, j)[0].lowerBound = -std::numeric_limits<FloatType>::infinity();
                vetMatBound[0].get(i, j).lowerBound = -std::numeric_limits<FloatType>::infinity();

            if(i != (vetNumSteps[0]-1))
                //matBound(i, j)[0].upperBound = r01;
                vetMatBound[0].get(i, j).upperBound = r01;
            else
                //matBound(i, j)[0].upperBound = std::numeric_limits<FloatType>::infinity();
                vetMatBound[0].get(i, j).upperBound = std::numeric_limits<FloatType>::infinity();

            //matBound(i, j)[1].lowerBound = r10;
            vetMatBound[1].get(i, j).lowerBound = r10;

            //matBound(i, j)[1].upperBound = r11;
            vetMatBound[1].get(i, j).upperBound = r11;
        }
    }

    //setupGraphBucket();

//    matBound(vetNumSteps[0]-1, 0)[0].upperBound = std::numeric_limits<double>::infinity();
//    matBound(0, 0)[0].lowerBound = -std::numeric_limits<double>::infinity();

    //for(int r=0; r < 2; ++r)
/*    {
        //std::cout<<"r("<<r<<"\n";

        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            double r00 = vetStepSize[0].start + i * vetStepSize[0].stepSize;
            double r01 = r00 + vetStepSize[0].stepSize;
            for(int j = 0; j < vetNumSteps[1]; ++j)
                std::cout<<"<"<<matBound(i, j)[0]<<" | "<<matBound(i, j)[1]<<"> ";

            std::cout<<"\n";
        }

        std::cout<<"\n\n";
    }*/

}

int LabelingAlgorithmNS::LabelingData::getIndex(int resource, FloatType val)
{

    //if(numMainResources == 1 && resource == 1)
    //    return 0;


    const FloatType start    = vetStepSize[resource].start;
    const FloatType end      = vetStepSize[resource].end;
    const FloatType stepSize = vetStepSize[resource].stepSize;

    const bool valGreterStart = val > start;
    const bool valLessEnd     = val < end;
    int index                 = 0;

    if(valGreterStart && valLessEnd)
        index = (int)((val-start)/stepSize);

    else if(valGreterStart)
        index = vetNumSteps[resource]-1;


    return index;

}


void LabelingAlgorithmNS::LabelingData::flushLabel()
{

    for(MatBucket& mat:vetMatBucketForward)
    {

        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[1]; ++j)
            {
//                mat.mat(i, j).vetPtrLabel.conservativeResize()
                mat.mat(i, j).flush();
            }

            if(numMainResources == 1)
                break;
        }
    }
}

std::ostream& LabelingAlgorithmNS::operator<< (std::ostream& out, const Bound &bound)
{
    out<<"["<<bound.lowerBound<<"; "<<bound.upperBound<<")";
    return out;
}

void LabelingAlgorithmNS::LabelingData::removeLabel(Label *label)
{
    if(Print)
        std::cout<<"Begin removeLabel\n";

    int cust = label->cust;
    int i    = label->i;
    int j    = label->j;

    int pos  = label->posBucket;
    int size = 0;

    Bucket *bucket = nullptr;

    if(label->typeLabel == forward)
        bucket = &vetMatBucketForward[cust].mat(i, j);

    else
        bucket = &vetMatBucketBackward[cust].mat(i,j);


    size = bucket->sizeVetPtrLabel;

    if(bucket->vetPtrLabel[pos] != label)
    {
        std::cout<<"ERROR in posBucket!\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    if(size > 1)
    {
        for(int t=pos; t < (bucket->sizeVetPtrLabel-1); ++t)
        {
            bucket->vetPtrLabel[t] = bucket->vetPtrLabel[t+1];
            bucket->vetPtrLabel[t]->posBucket = t;
        }
    }

    bucket->sizeVetPtrLabel -= 1;
}

void LabelingAlgorithmNS::Bucket::addLabel(LabelingAlgorithmNS::Label *labelPtr)
{

    if(sizeVetPtrLabel == vetPtrLabel.size())
    {
        vetPtrLabel.conservativeResize(2 * (sizeVetPtrLabel + 1));
        for(int i=sizeVetPtrLabel; i < 2*(sizeVetPtrLabel+1); ++i)
            vetPtrLabel[i] = nullptr;
    }

    vetPtrLabel[sizeVetPtrLabel] = labelPtr;
    labelPtr->posHeap = sizeVetPtrLabel;
    sizeVetPtrLabel += 1;
}


std::ostream&  LabelingAlgorithmNS::operator<<(std::ostream& out, const Label &label)
{

    out<<"Resouces(";
    for(int r=0; r < NumMaxResources; ++r)
        out<<label.vetResources[r]<<" ";

    out<<"); Route(";
    for(int i=0; i < label.tamRoute; ++i)
        out<<label.vetRoute[i]<<" ";

    out << ")";//bitSet(" << label.bitSetNg << ")";

    return out;
}

void LabelingAlgorithmNS::removeCycles(Label &label, const int numCust)
{
    static Eigen::VectorXi vetCust(numCust);
    vetCust.setConstant(-1);


    int i=0;

    if(Print)
        std::cout<<"Label antes: "<<label<<"\n\n";

    while(i < label.tamRoute)
    {
        if(vetCust[label.vetRoute[i]] != -1)
        {

            for(int k=i; k < (label.tamRoute-1); ++k)
            {
                label.vetRoute[k] = label.vetRoute[k+1];

            }

            label.tamRoute -= 1;

            if(Print)
                std::cout<<label<<"\n";

        }
        else
        {
            vetCust[label.vetRoute[i]] = i;
            i += 1;
        }
    }

}


void LabelingAlgorithmNS::removeCycles2(Label &label, const int numCust)
{

    static Eigen::VectorXi vetCust(numCust);
    vetCust.setConstant(-1);

    static Eigen::VectorXi vetRoute(NumMaxRoute);
    int num = 0;
    vetRoute.fill(0);

    for(int i=0; i < label.tamRoute; ++i)
    {
        if(vetCust[label.vetRoute[i]] == -1)
        {
            vetCust[label.vetRoute[i]] = i;
            vetRoute[num] = label.vetRoute[i];
            num += 1;
        }
    }

    for(int i=0; i < num; ++i)
        label.vetRoute[i] = vetRoute[i];

    label.tamRoute = num;


}

void LabelingAlgorithmNS::updateLabelCost(Label &label, const Vet3D_ResCost &vetMatResCost, FloatType labelStart)
{

    label.vetResources[0] = labelStart;
    for(int i=0; i < (label.tamRoute-1); ++i)
    {
        label.vetResources[0] += vetMatResCost(label.vetRoute[i], label.vetRoute[i+1], 0);
    }

}


LabelingAlgorithmNS::Label* LabelingAlgorithmNS::LabelingData::getBestLabel(int cust)
{
    MatBucket &matBucket = vetMatBucketForward[cust];
    Label *labelBest = nullptr;

    for(int i=0; i < vetNumSteps[0]; ++i)
    {
        for(int j=0; j < vetNumSteps[1]; ++j)
        {
            Bucket &bucket = matBucket.mat(i, j);

            for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
            {
                // TODO: FIX
                Label* label = bucket.vetPtrLabel[t];
                if(label->vetResources[0] < -DW_DecompNS::TolObjSubProb && !labelBest)
                {
                    labelBest = label;
                }
                else if(labelBest)
                {
                    if(label->vetResources[0] < labelBest->vetResources[0])
                        labelBest = label;
                }
            }
        }
    }

    return labelBest;

}

void LabelingAlgorithmNS::LabelingData::setupGraphBucket()
{
    graphBucket = GraphNS::Graph<int>(vetNumSteps[0]*vetNumSteps[1]);
    std::cout<<"Num nodes: "<<vetNumSteps[0]*vetNumSteps[1]<<"\n";
    std::cout<<"Num steps: "<<vetNumSteps[0]<<"; "<<vetNumSteps[1]<<"\n\n";

    //for(int cust=1; cust < numCust; ++cust)
    {
        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[0]; ++j)
            {
                int nodeId0 = getIndexGraphBucket(i, j);

                for(int ii=i; ii < vetNumSteps[0]; ++ii)
                {
                    for(int jj=j+1; jj < vetNumSteps[1]; ++jj)
                    {
                        int nodeId1 = getIndexGraphBucket(ii, jj);
                        graphBucket.addArc(nodeId0, nodeId1, 1.0);
                        //std::cout<<"Add arc: ("<<nodeId0<<", "<<nodeId1<<")\t"<<"i("<<i<<"), j("<<j<<"), ii("<<ii<<"), jj("<<jj<<")\n";
                    }
                }

            }
        }
    }

    // Bucked (i, j) pode dominar bucked (i', j') onde i' >= i e j' > j;
}


void LabelingAlgorithmNS::LabelingData::dominanceInterBuckets(LabelHeap&                 labelHeap,
                                                              int                        numRes,
                                                              const int                  localNumMaxLabel,
                                                              Eigen::VectorX<MatBucket>& vetMatBucket)
{
    if(Print)
        std::cout<<"dominanceInterBuckets\n\n";
    // Bucked (i, j) pode dominar bucked (i', j') onde i' > i e j' > j;

    int numDel = 0;

    for(int cust=1; cust < numCust; ++cust)
    {
        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[1]; ++j)
            {
                int nodeId0 = getIndexGraphBucket(i, j);
                Bucket &b0 = vetMatBucket[cust].mat(i, j);
                if(b0.sizeVetPtrLabel == 0)
                    continue;

                for(int ii=i; ii < vetNumSteps[0]; ++ii)
                {
                    for(int jj=j; jj < vetNumSteps[1]; ++jj)
                    {
                        if(ii == i && jj == j)
                            continue;

                        Bucket &b1 = vetMatBucket[cust].mat(ii, jj);
                        if(b1.sizeVetPtrLabel == 0)
                            continue;

                        int t1;

                        bool completeDominance = true;
                        if(ii > i && jj > j)
                            completeDominance = false;

                        for(int t0=0; t0 < b0.sizeVetPtrLabel; ++t0)
                        {
                            t1 = 0;
                            Label *label0 = b0.vetPtrLabel[t0];
                            while(t1 < b1.sizeVetPtrLabel)
                            {

                                Label *label1 = b1.vetPtrLabel[t1];
                                bool dominance = false;

                                if(completeDominance)
                                    dominance = checkCompleteDominance(*label0, *label1, numRes);
                                else
                                    dominance = checkDominanceSubSet(*label0, *label1);

                                //if(checkCompleteDominance(*label0, *label1, numRes))// && !checkDominance(*label1, *label0, numRes))
                                if(dominance)
                                {
                                    //std::cout<<"Domina\n";
                                    // Rm label1
                                    //eraseLabelFromSet(b1.vetPtrLabel[t1], setLabel);

                                    Label* label = b1.vetPtrLabel[t1];
                                    if(label != labelHeap.vet[label->posHeap])
                                    {
                                        std::cout << "ERROR, label(" << label << ") != labelHeap.vet[label->pos](" << labelHeap.vet[label->posHeap] << "\n\n";
                                        PRINT_DEBUG("", "");
                                        throw "ERROR";
                                    }

                                    labelHeap.deleteKey(b1.vetPtrLabel[t1]->posHeap);
                                    //bucket.vetPtrLabel[k]->active = false;
                                    labelPoolG.delT(b1.vetPtrLabel[t1]);

                                    if(t1 == (b1.sizeVetPtrLabel-1))
                                    {
                                        //labelPool.delT(bucket.vetPtrLabel[k]);
                                        b1.vetPtrLabel[t1] = nullptr;
                                    }
                                    else
                                    {
                                        //std::swap(b1.vetPtrLabel[t1], b1.vetPtrLabel[b1.sizeVetPtrLabel-1]);
                                        for(int i=t1; i < (b1.sizeVetPtrLabel-1); ++i)
                                        {
                                            b1.vetPtrLabel[i] = b1.vetPtrLabel[i+1];
                                            b1.vetPtrLabel[i]->posBucket = i;
                                        }

                                        //labelPool.delT(bucket.vetPtrLabel[bucket.sizeVetPtrLabel-1]);
                                        b1.vetPtrLabel[b1.sizeVetPtrLabel-1] = nullptr;
                                    }

                                    b1.sizeVetPtrLabel -= 1;

                                    numDel += 1;
                                    if(labelHeap.heapSize < int(0.7*localNumMaxLabel))
                                    {
                                        if(Print)
                                            std::cout<<"\t"<<"FORAM DELETADOS: "<<numDel<<" LABELS\n\n";
                                        return;
                                    }

                                    continue;
                                }
                                else
                                    t1 += 1;
                            }

                        }
                    }
                }
            }
        }
    }

    if(Print)
        std::cout<<"\t"<<"FORAM DELETADOS: "<<numDel<<" LABELS\n\n";

}

bool LabelingAlgorithmNS::checkDistance(const Eigen::Matrix<FloatType, -1, -1, Eigen::RowMajor> &matDist)
{


    for(int i=0; i < matDist.rows(); ++i)
    {
        for(int j=0; j < matDist.cols(); ++j)
        {
            if(matDist(i, j) < 0.0)
                return false;
        }
    }

    return true;

}

bool LabelingAlgorithmNS::containRoute(const Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> &vetLabel,
                                       int                                                           numSol,
                                       Label*                                                        label)
{
    if(numSol == 0)
        return false;

    for(int i=0; i < numSol; ++i)
    {
        if(vetLabel[i]->tamRoute != label->tamRoute)
            continue;

        bool equal = true;
        for(int ii=0; ii < label->tamRoute; ++ii)
        {
            if(vetLabel[i]->vetRoute[ii] != label->vetRoute[ii])
            {
                equal = false;
                break;
            }
        }

        if(equal)
            return true;
    }

    return false;
}

bool LabelingAlgorithmNS::labelHaveRoute(std::vector<int> &vetRoute, Label *label)
{
    if(vetRoute.empty())
        return false;
    int min = std::min((int)vetRoute.size(), label->tamRoute);

    for(int i=0; i < min; ++i)
    {
        if(vetRoute[i] != label->vetRoute[i] && vetRoute[i] != 0)
            return false;
    }

    return true;
}


/** **********************************************************************************************
 *  **********************************************************************************************
 *
 * @brief Performs intra-bucket dominance check for a given label
 * @details
 * Performs an intra-bucket dominance check for a given label. Since the labels are
 *      sorted by reduced cost, this function is responsible for finding the correct
 *      position of (labelPtrAux) in the bucket. It first locates the correct
 *      position while testing the dominance of (<label>, labelPtrAux). After
 *      finding the correct position, it then tests the dominance of the others labels
 *      ( labelPtrAux, <label>), where <label> are labels which are greater then
 *      labelPtrAux in the reduced cost.
 *
 * @param cust          Customer node
 * @param labelPtrAux   Pointer to the label being checked
 * @param lData         Reference to the labeling data structure
 * @param labelHeap     Reference to the heap of labels
 * @param numRes        Number of resources being considered
 * @param dest          Destination node
 * @param correctPos    Reference to an integer storing the correct position
 *                          for insertion or modification
 * @return              Pointer to the modified or resulting bucket after dominance check.
 *
 *  **********************************************************************************************
 *  **********************************************************************************************
 */
Bucket* LabelingAlgorithmNS::dominanceIntraBucket(int           cust,
                                                  Label*        labelPtrAux,
                                                  LabelingData& lData,
                                                  LabelHeap&    labelHeap,
                                                  int           numRes,
                                                  int           dest,
                                                  int&          correctPos)
{
    const bool rmFromHeap = (cust != dest);
    Bucket* bucketPtr = nullptr;

    if(cust != dest)
    {
        // Find index from resources
        int i = lData.getIndex(0, labelPtrAux->vetResources[0]);
        int j = 0;
        if(lData.numMainResources > 1)
            j = lData.getIndex(1, labelPtrAux->vetResources[1]);

        labelPtrAux->i = i;
        labelPtrAux->j = j;
        labelPtrAux->cust = cust;

        if(Print)
            std::cout << "\t\ti(" << i << "); j(" << j << ")\n";

        if(labelPtrAux->typeLabel == forward)
            bucketPtr = &lData.vetMatBucketForward[cust].mat(i, j);
        else
            bucketPtr = &lData.vetMatBucketBackward[cust].mat(i, j);
    }
    else
    {

        labelPtrAux->i = 0;
        labelPtrAux->j = 0;
        labelPtrAux->cust = cust;

        bucketPtr = &lData.vetMatBucketForward[cust].mat(0, 0);
    }

    int k=0;

    int dom = 0;

    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(bucketPtr->vetPtrLabel[correctPos]->vetResources[0] <= labelPtrAux->vetResources[0])
            correctPos += 1;
        else
            break;
    }

    // Labels before correctPos can dominate labelPtrAux;
    for(int ii=0; ii < correctPos; ++ii)
    {
        bool canDomindate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(bucketPtr->vetPtrLabel[ii]->vetResources[r] > labelPtrAux->vetResources[r])
            {
                canDomindate = false;
                break;
            }
        }

        if(canDomindate)
        {
            if(checkDominanceSubSet(*bucketPtr->vetPtrLabel[ii], *labelPtrAux))
            {

                //labelPoolG.delT(labelPtrAux);
                labelPtrAux = nullptr;
                return nullptr;
            }
        }
    }

    // Search for labels with have equal resources
    /*
    for(int ii=correctPos-1; ii >= 0; --ii)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[ii]->vetResources[0], labelPtrAux->vetResources[0], 1E-8))
        {
            bool canDominate = true;
            for(int r=1; r < numRes; ++r)
            {
                if(!doubleEqual(bucketPtr->vetPtrLabel[ii]->vetResources[r], labelPtrAux->vetResources[r], 1E-8))
                {
                    canDominate = false;
                    break;
                }
            }

            if(canDominate)
            {
                if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[ii]))
                {

                    if(rmFromHeap)
                        labelHeap.deleteKey(bucketPtr->vetPtrLabel[ii]->posHeap);

                    labelPoolG.delT(bucketPtr->vetPtrLabel[ii]);

                    for(int t=ii; t < (bucketPtr->sizeVetPtrLabel-1); ++t)
                    {
                        bucketPtr->vetPtrLabel[t] = bucketPtr->vetPtrLabel[t+1];
                        bucketPtr->vetPtrLabel[t]->posBucket = t;
                    }

                    bucketPtr->sizeVetPtrLabel -= 1;
                    correctPos -= 1;
                }
            }
        }
        else
            break;
    }
    */

    // Labels after correctPos (inclusive) can be dominated by labelPtrAux
    int bucketSizeVetPtrLabel = bucketPtr->sizeVetPtrLabel;
    int ii = correctPos;

    while(ii < bucketSizeVetPtrLabel)
    {
        bool canDominate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(labelPtrAux->vetResources[r] > bucketPtr->vetPtrLabel[ii]->vetResources[r])
            {
                canDominate = false;
                break;
            }
        }

        if(canDominate)
        {
            if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[ii]))
            {
                // Remove bucket.vetPtrLabel[ii]
                if(rmFromHeap)
                    labelHeap.deleteKey(bucketPtr->vetPtrLabel[ii]->posHeap);

                labelPoolG.delT(bucketPtr->vetPtrLabel[ii]);

                for(int t=ii; t < (bucketPtr->sizeVetPtrLabel-1); ++t)
                {
                    bucketPtr->vetPtrLabel[t] = bucketPtr->vetPtrLabel[t+1];
                    bucketPtr->vetPtrLabel[t]->posBucket = t;
                }


                bucketPtr->sizeVetPtrLabel -= 1;
                bucketSizeVetPtrLabel  -= 1;
                continue; // while ii
            }
        }

        ii += 1;
    }

    if(labelPtrAux == nullptr)
        return nullptr;

    maxDistG = std::max(maxDistG, labelPtrAux->vetResources[0]);
    minDistG = std::min(minDistG, labelPtrAux->vetResources[0]);

    labelPtrAux->active = true;

    return bucketPtr;

}

bool LabelingAlgorithmNS::checkCompleteDominance(const Label& l0, const Label& l1, int numResources)
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
        //if(doubleLess(l1.vetResources[i], l0.vetResources[i], std::numeric_limits<FloatType>::epsilon()))
        if(l0.vetResources[i] > l1.vetResources[i])
            return false;

        if((i+1) == numResources)
            break;
    }

    // Check if l0 is a subset of l1
    //return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
    return checkDominanceSubSet(l0, l1);
}

static LabelCmp labelCmp;

Label* LabelHeap::extractMin()
{
    if(heapSize == 0)
        return nullptr;
    if(heapSize == 1)
    {
        heapSize += -1;
        Label* root = vet[0];
        vet[0] = nullptr;
        return root;
    }

    Label* root = vet[0];
    vet[0] = vet[heapSize-1];
    vet[0]->posHeap = 0;
    heapSize += -1;
    heapify(0);
    return root;

}

void LabelHeap::decreaseKey(int i, FloatType val)
{
    vet[i]->vetResources[0] = val;
    int iParent = parent(i);

    while(i!=0 && vet[iParent]->vetResources[0] > vet[i]->vetResources[0])
    {
        std::swap(vet[i], vet[iParent]);
        std::swap(vet[i]->posHeap, vet[iParent]->posHeap);
        i = iParent;
        iParent = parent(i);
    }

}

void LabelHeap::deleteKey(int i)
{
    decreaseKey(i, -std::numeric_limits<FloatType>::max());
    Label* label = extractMin();

}

void LabelHeap::heapify(int i)
{
    int l, r, smallert;
    bool cond;

    do
    {

        l = left(i);
        r = right(i);
        smallert = i;

        if(l < heapSize && labelCmp(vet[l], vet[i]))//vet[l]->vetResources[0] < vet[i]->vetResources[i])
            smallert = l;

        if(r < heapSize && labelCmp(vet[r], vet[smallert]))//vet[r]->vetResources[0] < vet[smallert]->vetResources[0])
            smallert = r;

        cond = (smallert != i);
        if(cond)
        {
            std::swap(vet[i], vet[smallert]);
            std::swap(vet[i]->posHeap, vet[smallert]->posHeap);
        }

    }
    while(cond);

}

void LabelHeap::insertKey(Label *label)
{
    if(heapSize == (int)vet.size())
    {
        vet.resize(2*heapSize);
        for(int i=heapSize; i < (int)vet.size(); ++i)
            vet[i] = nullptr;
    }

    heapSize += 1;
    int i = heapSize - 1;
    int iParent = parent(i);
    vet[i] = label;
    label->posHeap = i;

    while(i!=0 && labelCmp.isGreater(vet[iParent], vet[i]))//!labelCmp(vet[iParent], vet[i]))
    {
        std::swap(vet[i], vet[iParent]);
        std::swap(vet[i]->posHeap, vet[iParent]->posHeap);
        i = iParent;
        iParent = parent(i);
    }
}
