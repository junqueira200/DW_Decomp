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
//#include "NgSet.h"
#include "LabelingUtils.h"

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

    int maxSize = 0;
    int localNumMaxLabel = numMaxLabelG;

    Label* ptrLabelTarget = nullptr;

    while(!labelHeap.empty() &&
         ((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel < DW_DecompNS::NumMaxSolSubProb) || exact))
    {
        if((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel >= 10 && labelHeap.heapSize >= 25 * numCust &&
            !exact))
            break;

        // TODO remover
        //lData.checkMat();

        lData.checkVetMatBucketForward();

        if(ptrLabelTarget)
        {
            std::cout<<"ptrLabelTarget("<<ptrLabelTarget->active<<"): "<<*ptrLabelTarget<<"\n\n";
        }

        if(labelHeap.heapSize > localNumMaxLabel && DominaIterBuckets)
        {
            lData.dominanceInterBuckets(labelHeap, numRes, localNumMaxLabel, lData.vetMatBucketForward, Forward);


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


            if(extendLabelForward(*labelPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelPtr->cust, t, ngSet, numRes))
            {
                if(Print)
                    std::cout<<"\t"<<*labelPtrAux<<"\n\n";

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


                    //std::cout<<"\nPrimeira Solucao:\n"<<*labelPtrAux<<"\n"<<labelPtrAux->bitSetNg<<"\n";
                    //exit(-1);

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

        /*
        if(numIt == 5)
        {
            std::cout<<"Error!\n";
            PRINT_DEBUG("", "");
            exit(-1);
        }
        */
    }


    // TODO Remover
    //std::cout<<"Primeiras Solucoes: \n";

    if(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel > 0)
    {
        numSol = std::min(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel, DW_DecompNS::NumMaxSolSubProb);
        for(int l=0; l < numSol; ++l)
        {
            Label* label = lData.vetMatBucketForward[dest].mat(0, 0).vetPtrLabel[l];
            //std::cout<<*label<<"\n";
            label->vetRoute[label->tamRoute-1] = label->vetRoute[0];
            auto &vetRoute = label->vetRoute;


            for(int i = 0; i < (label->tamRoute - 1); ++i)
            {
                int index = getIndex(vetRoute[i], vetRoute[i+1], numCust-1);
                matColX(index, l) = 1.0;
            }

            vetRedCost[l] = label->vetResources[0];
        }

        //exit(-1);
        return true;
    }
    else
    {
        return false;
    }
}


void LabelingAlgorithmNS::checkHeap(LabelHeap& heap, LabelingData& lData)
{
    for(int k=0; k < heap.heapSize; ++k)
    {
        Label* label = heap.vet[k];

        int posBucket = label->posBucket;
        int posHeap   = label->posHeap;

        if(k != posHeap)
        {
            std::cout<<"posHeap esta errada\n";
            std::cout<<"posHeap("<<posHeap<<"), k("<<k<<")\n";
            std::cout<<*label<<"\n\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

        Bucket* bucket = lData.getBucket(label);
        if(bucket->vetPtrLabel[posBucket] != label)
        {

            std::cout<<"\nposBucket esta errada\n";

            std::cout<<"posBucket("<<posBucket<<")\n";
            std::cout<<*label<<"\n\n";

            PRINT_DEBUG("", "");
            throw "ERROR";

        }

    }
}


/** **********************************************************************************************************
 *  **********************************************************************************************************
 *
 * @brief Executes the Forward labeling algorithm for resource-constrained path enumeration.
 *
 * @details This function performs label generation and extension from the origin to the destination
 * node, applying resource constraints, dominance checks. It is
 * used in column generation methods for routing problems.
 *
 * @param numRes                Number of resources (e.g., time, capacity).
 * @param numCust               Number of customers or nodes in the problem.
 * @param vetMatResCost         3D matrix of resource costs between nodes.
 * @param matBoundRes           Matrix of upper/lower bounds for resource usage.
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
bool LabelingAlgorithmNS::
     bidirectionalAlgorithm(const int numRes, const int numCust, const Vet3D_ResCost& vetMatResCostForward,
                            const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& matBoundRes, const int dest,
                            const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                            const FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                            Eigen::VectorX<FloatType>& vetRedCost, bool exact)
{
    //vetRoteG = {19, 18, 5, 10, 1, 12, 9, 17, 9, 17, 9, 17};

    //static Eigen::VectorX<FloatType> vetLowerBoundDist(numCust);
    //const bool vetLowerBoundDistValid =  getDistLowerBound(vetMatResCost, vetLowerBoundDist);
    VetBackwardMask vetBackwardMask;
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
        std::cout << "******************************************************************\n";
        std::cout << "*****************BIDIRECTIONAL LABELING ALGORITHM*****************\n\n";
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

    Label* labelPtr  = labelPoolG.getT();
    Label* labelBackwardPtr = labelPoolG.getT();


    // TODO create a start label
    //labelPtr->vetRoute.resize(10);
    labelPtr->tamRoute = 1;
    labelPtr->vetRoute[0] = 0;
    labelPtr->vetResources.setZero();
    labelPtr->vetResources[0] = labelStart;
    labelPtr->bitSetNg = 0;
    labelPtr->bitSetNg[0] = true;

    int i = lData.getIndex(0, labelStart);
    int j = lData.getIndex(1, 0.0);

    labelPtr->i         = i;
    labelPtr->j         = j;
    labelPtr->cust      = 0;
    labelPtr->active    = true;
    labelPtr->posBucket = 0;

    /*
    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel[0]  = labelForwardPtr;
    lData.vetMatBucketForward[0].mat(i, j).sizeVetPtrLabel = 1;
    */

    //labelBackwardPtr

    startBackwardLabel(labelBackwardPtr, vetBackwardMask, matBoundRes, numRes, dest, labelStart, i, lData);
    labelHeap.insertKey(labelBackwardPtr);



    //labelHeap.insertKey(labelForwardPtr);
    labelPtr = nullptr;
    labelBackwardPtr = nullptr;

    int numIt = 0;

    Label *labelPtrBest = nullptr;
    int maxSize = 0;
    int maxSizeVetPtrLabel = 0;
    int localNumMaxLabel = numMaxLabelG;

    Label* ptrLabelTarget = nullptr;
    bool printSize = false;


    ArrayResources vetMaxResources;
    vetMaxResources.setConstant(-std::numeric_limits<FloatType>::max());
    vetMaxResources[0] = (FloatType)0.0;

    for(int r=1; r < numRes; ++r)
    {
        for(int i=0; i < numCust; ++i)
        {
            if(matBoundRes(i, r).upperBound > vetMaxResources[r])
                vetMaxResources[r] = matBoundRes(i, r).upperBound;
        }
    }

    //std::cout<<"Max Resources: "<<vetMaxResources<<"\n";

    while(!labelHeap.empty() &&
          ((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel < DW_DecompNS::NumMaxSolSubProb) || exact))
    {
        if((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel >= 10 && labelHeap.heapSize >= 25 * numCust &&
            !exact))
            break;

        //checkHeap(labelHeap, lData);
        lData.checkVetMatBucketBackward();

        if(ptrLabelTarget)
        {
            std::cout<<"ptrLabelTarget("<<ptrLabelTarget->active<<"): "<<*ptrLabelTarget<<"\n\n";
        }

        if(labelHeap.heapSize > localNumMaxLabel && DominaIterBuckets)
        {
            lData.dominanceInterBuckets(labelHeap, numRes, localNumMaxLabel, lData.vetMatBucketForward, Backward);

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

        if(labelHaveRoute(vetRoteG, labelPtr))
        {
            std::cout<<"labelHaveRoute\n"<<*labelPtr<<"\ni("<<labelPtr->i<<"), j("<<j<<")\n";
        }


        //std::cout<<"Extract label("<<labelPtr->vetRoute<<")\n";
        int lastCust = labelPtr->cust;

        if(!labelPtr->active)
        {
            std::cout<<"labelPtr is not active!\n"<<*labelPtr<<"\n";
            PRINT_DEBUG("", "");
            throw "ERROR";

            labelPoolG.delT(labelPtr);
            continue;
        }

        if(lastCust == dest && labelPtr->typeLabel == Forward)
        {
            if(Print)
                std::cout<<"lastCust==dest("<<dest<<")\n";
            continue;
        }

        if(lastCust == 0 && labelPtr->typeLabel == Backward)
            continue;

        // Remove labelPtr from vetMatBucket
        lData.removeLabel(labelPtr);

        // Extend label
        for(int t=0; t < numCust; ++t)
        {


            int tAux = t;
            //std::cout<<"t: "<<t<<"\n";
            if(labelPtr->typeLabel  == Forward &&
               vetMatResCostForward(labelPtr->cust, t, 0) == std::numeric_limits<FloatType>::infinity())
                continue;

            else if(labelPtr->typeLabel  == Backward &&
               vetMatResCostForward(t, labelPtr->cust, 0) == std::numeric_limits<FloatType>::infinity())
                continue;


            if(labelPtr->bitSetNg[t] == 1)
                continue;

            if(t == labelPtr->cust)
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

            if(extendLabel(*labelPtr, *labelPtrAux, vetMatResCostForward, vetMatResCostBackward, matBoundRes,
                           labelPtr->cust, t, ngSet, numRes, vetBackwardMask))
            {
                if(Print)
                    std::cout<<"\t"<<*labelPtrAux<<"\n";

                int correctPos = 0;

                if(t == dest && labelPtrAux->typeLabel == Forward)
                {
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    removeCycles2(*labelPtrAux, numCust);
                    updateLabelCost(*labelPtrAux, vetMatResCostForward, labelStart);


                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                }


                if(t == 0 && labelPtrAux->typeLabel == Backward)
                {

                    //applyLabelMask(labelPtrAux, numRes, vetBackwardMask);
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {
                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    removeCycles2(*labelPtrAux, numCust);

                    // Reverse the array
                    std::reverse(labelPtrAux->vetRoute.begin(),
                                 labelPtrAux->vetRoute.begin()+labelPtrAux->tamRoute);


                    updateLabelCost(*labelPtrAux, vetMatResCostForward, labelStart);


                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    //change the backward type to forward
                    convertLabelBackwardToForward(labelPtrAux, vetMaxResources, numRes);
                    labelPtrAux->cust      = dest;
                    labelPtrAux->i = 0;
                    labelPtrAux->j = 0;
                    tAux = dest;

                    /*
                    static int numP = 0;
                    std::cout<<*labelPtrAux<<"\n\n";
                    if(numP == 5)
                        exit(-1);

                    numP += 1;
                    */

                }

                Bucket* bucket = dominanceIntraBucket(tAux, labelPtrAux, lData, labelHeap, numRes, dest, correctPos);

                if(!bucket)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if((bucket->sizeVetPtrLabel+1) > NumMaxLabePerBucket && tAux != dest)
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

                if(tAux != dest && labelPtrAux->typeLabel == Forward)
                    labelHeap.insertKey(labelPtrAux);

                else if(labelPtrAux->typeLabel == Backward)
                    labelHeap.insertKey(labelPtrAux);

                //checkHeap(labelHeap, lData);
            }
            else
                labelPoolG.delT(labelPtrAux);

        }

//std::cout<<"END Extend label\n";

        if(Print)
            std::cout<<"\n########################################################################################\n\n";

        labelPoolG.delT(labelPtr);
        numIt += 1;

        /*
        if(numIt == 5)
        {
            std::cout<<"Error!\n";
            PRINT_DEBUG("", "");
            exit(-1);
        }
        */
    }

    //std::cout<<"Solucoes: \n";

    if(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel > 0)
    {
        numSol = std::min(lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel, DW_DecompNS::NumMaxSolSubProb);
        for(int l=0; l < numSol; ++l)
        {
            Label* label = lData.vetMatBucketForward[dest].mat(0, 0).vetPtrLabel[l];
            //std::cout<<*label<<"\n";
            label->vetRoute[label->tamRoute-1] = label->vetRoute[0];
            auto &vetRoute = label->vetRoute;


            for(int i = 0; i < (label->tamRoute - 1); ++i)
            {
                int index = getIndex(vetRoute[i], vetRoute[i+1], numCust-1);
                matColX(index, l) = 1.0;
            }

            vetRedCost[l] = label->vetResources[0];
        }

        // TODO remover
        //exit(-1);
        //std::cout<<"*************************************************************************************\n\n";
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
 * @param t                 The index of the next customer.
 * @param ngSet             The set representing non-neighbor restrictions.
 * @param numResources      The number of resources considered.
 *
 * @return Returns \c true if the label was successfully extended, otherwise returns \c false.
 *
 * *********************************************************************************************
 * *********************************************************************************************
 */
bool LabelingAlgorithmNS::
     extendLabelForward(const Label& label, Label& newLabel, const Vet3D_ResCost& vetMatResCost,
                        const MatBoundRes& vetVetBound, int custI, int t, const NgSet& ngSet, int numResources)
{


    // Goes through resources
    bool boundOk = true;

    newLabel.typeLabel = label.typeLabel;

    for(int i=0; i < numResources; ++i)//NumMaxResources; ++i)
    {

        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] + vetMatResCost(custI, t, i);
        boundOk = boundOk && (newLabel.vetResources[i] <= vetVetBound(t, i).upperBound);
    }

    if(!boundOk)
        return false;

    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ(t) its in custI ngSet
    if(ngSet.contain(t, custI))
    {
        newLabel.bitSetNg[t] = true;
    }


    for(int i=0; i < label.tamRoute; ++i)
        newLabel.vetRoute[i] = label.vetRoute[i];

    newLabel.vetRoute[label.tamRoute] = t;
    newLabel.tamRoute = label.tamRoute+1;
    newLabel.active = true;
    newLabel.typeLabel = Forward;

    return true;

}


bool LabelingAlgorithmNS::
     extendLabelBackward(const Label& label, Label& newLabel, const Vet3D_ResCost& vetMatResCostBackward,
                         const MatBoundRes& vetVetBound, int custI, int t, const NgSet& ngSet, int numResources,
                         const Eigen::Array<FloatType, 1, NumMaxResources>& vetBackwardMask)
{


    // Goes through resources
    bool boundOk = true;
    newLabel.typeLabel = label.typeLabel;
    newLabel.vetResources[0] = label.vetResources[0] + vetMatResCostBackward(t, custI, 0);

    if(Print)
        std::cout<<"\t\t("<<t<<","<<custI<<"): "<<vetMatResCostBackward(t, custI, 0)<<"\n";

    for(int i=1; i < numResources; ++i)//NumMaxResources; ++i)
    {
        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] - vetMatResCostBackward(t, custI, i);
        boundOk = boundOk && ((vetBackwardMask[i]*newLabel.vetResources[i]) <= vetVetBound(t, i).upperBound);
        boundOk = boundOk && ((vetBackwardMask[i]*newLabel.vetResources[i]) >= vetVetBound(t, i).lowerBound);
    }

    if(!boundOk)
        return false;

    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ(t) its in custI ngSet
    if(ngSet.contain(custI, t))
    {
        newLabel.bitSetNg[custI] = true;
    }


    for(int i=0; i < label.tamRoute; ++i)
        newLabel.vetRoute[i] = label.vetRoute[i];

    newLabel.vetRoute[label.tamRoute] = t;
    newLabel.tamRoute = label.tamRoute+1;
    newLabel.active = true;
    newLabel.typeLabel = Backward;

    return true;
}



void LabelingAlgorithmNS::
     startBackwardLabel(Label* labelPtr, VetBackwardMask& vetBackwardMask, const MatBoundRes& vetVetBound,
                        int numResources, int dest, double labelStart, int i, LabelingData& lData)
{
    labelPtr->tamRoute        = 1;
    labelPtr->vetRoute[0]     = dest;
    labelPtr->typeLabel       = Backward;
    labelPtr->vetResources.setZero();
    labelPtr->vetResources[0] = labelStart;
    labelPtr->bitSetNg        = 0;
    labelPtr->bitSetNg[dest]  = true;

    labelPtr->i               = i;
    labelPtr->cust            = dest;
    labelPtr->active          = true;
    labelPtr->posBucket       = 0;

    vetBackwardMask.setOnes();

    for(int r=1; r < numResources; ++r)
    {
        FloatType max = -std::numeric_limits<FloatType>::infinity();
        for(int i=0; i <= dest; ++i)
            max = std::max(vetVetBound(i, r).upperBound, max);

        if(max == std::numeric_limits<FloatType>::infinity())
            vetBackwardMask[r] = -1;
        else
            labelPtr->vetResources[r] = max;

    }

    labelPtr->j = lData.getIndex(1, labelPtr->vetResources[1]);
    int j = labelPtr->j;

    lData.vetMatBucketBackward[dest].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucketBackward[dest].mat(i, j).vetPtrLabel[0]  = labelPtr;
    lData.vetMatBucketBackward[dest].mat(i, j).sizeVetPtrLabel = 1;

}


/*
 *
 * 50
 *
 * custI    = 50
 * custJ(t) = 5
 *
 * (5 50)
 *
 */

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
LabelingAlgorithmNS::LabelingData::LabelingData(const Eigen::Vector<Step, 2>& vetStepSize_, int numMainResources_,
                                                int numCust_)
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


    for(MatBucket& mat:vetMatBucketBackward)
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
    {
        std::cout<<"Begin removeLabel\n";
        std::cout<<*label<<"\n";
    }

    int cust = label->cust;
    int i    = label->i;
    int j    = label->j;

    int pos  = label->posBucket;
    int size = 0;

    Bucket *bucket = nullptr;

    if(label->typeLabel == Forward)
        bucket = &vetMatBucketForward[cust].mat(i, j);

    else
        bucket = &vetMatBucketBackward[cust].mat(i, j);


    size = bucket->sizeVetPtrLabel;

    // TODO ERRO OCORRE AQUI
    if(bucket->vetPtrLabel[pos] != label)
    {
        std::cout<<"ERROR in posBucket!\n";
        std::cout<<"cust("<<cust<<"), pos("<<pos<<")\n";

        std::cout<<"Backward: "<<(label->typeLabel == Backward)<<"\n";

        for(int i=0; i < bucket->sizeVetPtrLabel; ++i)
        {
            if(label == bucket->vetPtrLabel[i])
                std::cout<<"Achou na posicao errada!\n";
        }

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


void LabelingAlgorithmNS::LabelingData::
     dominanceInterBuckets(LabelHeap& labelHeap, int numRes, const int localNumMaxLabel,
                           Eigen::VectorX<MatBucket>& vetMatBucket, TypeLabel typeLabel)
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
                        // TODO: FIX
                        bool completeDominance = true;//true;

                        if(typeLabel == Forward && ii > i && jj > j)
                            completeDominance = false;
                        else if(typeLabel == Forward && ii > i && jj < j)
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
                                        std::cout << "ERROR, label(" << label << ") != labelHeap.vet[label->pos](" <<
                                                     labelHeap.vet[label->posHeap] << "\n\n";
                                        PRINT_DEBUG("", "");
                                        throw "ERROR";
                                    }
                                    \
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
                                       int numSol, Label* label)
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
    if(label->tamRoute != (int)vetRoute.size())
        return false;

    int min = std::min((int)vetRoute.size(), label->tamRoute);

    for(int i=0; i < min; ++i)
    {
        if(vetRoute[i] != label->vetRoute[i] && vetRoute[i] != 0)
            return false;
    }

    return true;
}


/** ********************************************************************************************************************
 *  ********************************************************************************************************************
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
 *  ********************************************************************************************************************
 *  ********************************************************************************************************************
 */
Bucket* LabelingAlgorithmNS::dominanceIntraBucketForward(int cust, Label* labelPtrAux, LabelingData& lData,
                                                         LabelHeap& labelHeap, int numRes, int dest, int& correctPos)
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

        bucketPtr = &lData.vetMatBucketForward[cust].mat(i, j);
    }
    else
    {

        labelPtrAux->i = 0;
        labelPtrAux->j = 0;
        labelPtrAux->cust = cust;
        bucketPtr = &lData.vetMatBucketForward[cust].mat(0, 0);
    }

    // Get the correct possition
    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleGreaterEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0],
                              FloatEp))
            break;
        correctPos += 1;
    }

    int before = correctPos;

    // Search for labels with have equal resources
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], FloatEp))
        {
            bool equal = true;
            for(int r=1; r < numRes; ++r)
            {
                if(!doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[r], labelPtrAux->vetResources[r],
                               FloatEp))
                {
                    equal = false;
                    break;
                }
            }

            if(equal)
            {
                if(checkDominanceSubSet(*bucketPtr->vetPtrLabel[correctPos], *labelPtrAux))
                {
                    labelPtrAux = nullptr;
                    return nullptr;
                }
                else if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[correctPos]))
                {
                    if(rmFromHeap)
                        labelHeap.deleteKey(bucketPtr->vetPtrLabel[correctPos]->posHeap);

                    labelPoolG.delT(bucketPtr->vetPtrLabel[correctPos]);

                    for(int ii=correctPos; ii < (bucketPtr->sizeVetPtrLabel-1); ++ii)
                    {
                        bucketPtr->vetPtrLabel[ii] = bucketPtr->vetPtrLabel[ii+1];
                        bucketPtr->vetPtrLabel[ii]->posBucket = ii;
                    }

                    bucketPtr->sizeVetPtrLabel -= 1;

                }
                else
                    correctPos += 1;
            }
            else
                correctPos += 1;

        }
        else
            break;
    }

    //std::cout<<"end while\n\n";

    // Labels before correctPos can dominate labelPtrAux;
    for(int ii=0; ii < before; ++ii)
    {
        //if(posEqual > 0 && ii == posEqual)
        //    break;

        bool canDomindate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(doubleGreater(bucketPtr->vetPtrLabel[ii]->vetResources[r], labelPtrAux->vetResources[r], FloatEp))
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

    // Labels after correctPos (inclusive) can be dominated by labelPtrAux
    int bucketSizeVetPtrLabel = bucketPtr->sizeVetPtrLabel;
    int ii = correctPos;

    while(ii < bucketSizeVetPtrLabel)
    {
        bool canDominate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(doubleGreater(labelPtrAux->vetResources[r], bucketPtr->vetPtrLabel[ii]->vetResources[r], FloatEp))
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

/*

Bucket* LabelingAlgorithmNS::dominanceIntraBucketForward(int           cust,
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

        if(labelPtrAux->typeLabel == Forward)
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

*/


/*
 *
 *
Bucket* LabelingAlgorithmNS::dominanceIntraBucketBackward(int           			cust,
                                                          Label*        			labelPtrAux,
                                                          LabelingData& 			lData,
                                                          LabelHeap&    			labelHeap,
                                                          int           			numRes,
                                                          int           			dest,
                                                          int&          			correctPos)
{
    // TODO Altera vetMatBucketBackward

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

        bucketPtr = &lData.vetMatBucketBackward[cust].mat(i, j);
    }
    else
    {

        PRINT_DEBUG("", "");
        std::cout<<"cust("<<cust<<"), dest("<<dest<<")\n";
        throw "ERROR";

        labelPtrAux->i = 0;
        labelPtrAux->j = 0;
        labelPtrAux->cust = cust;
        bucketPtr = &lData.vetMatBucketForward[cust].mat(0, 0);
    }

    if(bucketPtr->sizeVetPtrLabel == 0)
    {
        correctPos = 0;
        return bucketPtr;
    }

    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleLess(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], FloatEp))
            correctPos += 1;
        else
            break;
    }

    // TODO check
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], FloatEp))
        {
            if(doubleLessEqual(labelPtrAux->vetResources[1], bucketPtr->vetPtrLabel[correctPos]->vetResources[1],
                               FloatEp))
            correctPos += 1;
        else
            break;
        }
        else
            break;
    }


    // Labels before correctPos can dominate labelPtrAux;
    for(int ii=0; ii < correctPos; ++ii)
    {
        bool canDomindate = true;
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

    // Labels after correctPos (inclusive) can be dominated by labelPtrAux
    int bucketSizeVetPtrLabel = bucketPtr->sizeVetPtrLabel;
    int ii = correctPos;

    while(ii < bucketSizeVetPtrLabel)
    {
        bool canDominate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(labelPtrAux->vetResources[r] < bucketPtr->vetPtrLabel[ii]->vetResources[r])
            {
                canDominate = false;
                break;
            }
        }

        if(canDominate)
        {
            if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[ii]))
            {

                if(Print && labelHaveRoute(vetRoteG, bucketPtr->vetPtrLabel[ii]))
                {
                    std::cout<<"AQUI!";
                    PRINT_DEBUG("", "");
                    throw "ERROR";
                }


                // Remove bucket.vetPtrLabel[ii]
                if(rmFromHeap)
                {
                    labelHeap.deleteKey(bucketPtr->vetPtrLabel[ii]->posHeap);
                }

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

 *
 * */

Bucket* LabelingAlgorithmNS::
        dominanceIntraBucketBackward(int cust, Label* labelPtrAux, LabelingData& lData, LabelHeap& labelHeap, int numRes,
                                     int dest, int& correctPos)
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

        bucketPtr = &lData.vetMatBucketBackward[cust].mat(i, j);
    }
    else
    {

        labelPtrAux->i = 0;
        labelPtrAux->j = 0;
        labelPtrAux->cust = cust;
        bucketPtr = &lData.vetMatBucketForward[cust].mat(0, 0);
    }

    // Get the correct possition
    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleGreaterEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0],
                              FloatEp))
            break;
        correctPos += 1;
    }

    int before = correctPos;

    // Search for labels with have equal resources
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], FloatEp))
        {
            bool equal = true;
            for(int r=1; r < numRes; ++r)
            {
                if(!doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[r], labelPtrAux->vetResources[r],
                               FloatEp))
                {
                    equal = false;
                    break;
                }
            }

            if(equal)
            {
                if(checkDominanceSubSet(*bucketPtr->vetPtrLabel[correctPos], *labelPtrAux))
                {
                    labelPtrAux = nullptr;
                    return nullptr;
                }
                else if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[correctPos]))
                {
                    if(rmFromHeap)
                        labelHeap.deleteKey(bucketPtr->vetPtrLabel[correctPos]->posHeap);

                    labelPoolG.delT(bucketPtr->vetPtrLabel[correctPos]);

                    for(int ii=correctPos; ii < (bucketPtr->sizeVetPtrLabel-1); ++ii)
                    {
                        bucketPtr->vetPtrLabel[ii] = bucketPtr->vetPtrLabel[ii+1];
                        bucketPtr->vetPtrLabel[ii]->posBucket = ii;
                    }

                    bucketPtr->sizeVetPtrLabel -= 1;

                }
                else
                    correctPos += 1;
            }
            else
                correctPos += 1;

        }
        else
            break;
    }

    //std::cout<<"end while\n\n";

    // Labels before correctPos can dominate labelPtrAux;
    for(int ii=0; ii < before; ++ii)
    {
        //if(posEqual > 0 && ii == posEqual)
        //    break;

        bool canDomindate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(doubleLess(bucketPtr->vetPtrLabel[ii]->vetResources[r], labelPtrAux->vetResources[r], FloatEp))
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

    // Labels after correctPos (inclusive) can be dominated by labelPtrAux
    int bucketSizeVetPtrLabel = bucketPtr->sizeVetPtrLabel;
    int ii = correctPos;

    while(ii < bucketSizeVetPtrLabel)
    {
        bool canDominate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(doubleLess(labelPtrAux->vetResources[r], bucketPtr->vetPtrLabel[ii]->vetResources[r], FloatEp))
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

    if(l0.typeLabel == Forward)
    {
        // Check the resources
        //#pragma GCC unroll NumMaxResources
        for(int i=0; i < numResources; ++i)
        {
            //if(l0.vetResources[i] > l1.vetResources[i])
            if(doubleGreater(l0.vetResources[i], l1.vetResources[i], FloatEp) &&
               !doubleEqual(l0.vetResources[i], l1.vetResources[i], FloatEp))
                return false;
        }
    }
    else
    {
        //if(l0.vetResources[0] > l1.vetResources[0])
        if(doubleGreater(l0.vetResources[0], l1.vetResources[0], FloatEp) &&
           !doubleEqual(l0.vetResources[0], l1.vetResources[0], FloatEp))
            return false;

        // Check the resources
        //#pragma GCC unroll NumMaxResources
        for(int i=1; i < numResources; ++i)
        {
            //if(l0.vetResources[i] < l1.vetResources[i])
            if(doubleLess(l0.vetResources[i], l1.vetResources[i], FloatEp) &&
               !doubleEqual(l0.vetResources[i], l1.vetResources[i], FloatEp))
                return false;
        }
    }

    // Check if l0 is a subset of l1
    //return (l0.bitSetNg & l1.bitSetNg) == l0.bitSetNg;
    return checkDominanceSubSet(l0, l1);
}

void LabelingAlgorithmNS::convertLabelBackwardToForward(Label* label, const ArrayResources &vetMaxResources,
                                                        int numResources)
{
    if(label->typeLabel == Forward)
    {
        std::cout<<"Label("<<label<<") is alrey of the Forward type!\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    for(int r=1; r < numResources; ++r)
        label->vetResources[r] = vetMaxResources[r] - label->vetResources[r];

    label->typeLabel = Forward;
}

bool LabelingAlgorithmNS::LabelingData::compareVetMatBucket(const ArrayResources& vetMaxResouces)
{

    Label* labelPtr = labelPoolG.getT();

    for(int cust=1; cust < (numCust-1); ++cust)
    {
        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[1]; ++j)
            {
                Bucket& bucket = vetMatBucketBackward[cust].mat(i, j);
                for(int k=0; k < bucket.sizeVetPtrLabel; ++k)
                {
                    Label* labelBackward = bucket.vetPtrLabel[k];
                    copyLabel(*labelBackward, *labelPtr, 2);
                    convertLabelBackwardToForward(labelPtr, vetMaxResouces, 2);

                    labelPtr->i = getIndex(0, labelPtr->vetResources[0]);
                    labelPtr->j = getIndex(1, labelPtr->vetResources[1]);

                    Bucket& bucketForward = vetMatBucketForward[cust].mat(labelPtr->i, labelPtr->j);

                    if(!searchLabel(labelPtr, bucketForward))
                    {
                        std::cout<<"ERROR\nLabel: "<<*labelPtr<<", nao foi achado em vetMatBucketForwar\n";
                        PRINT_DEBUG("", "");
                        exit(-1);
                    }
                }
            }
        }
    }

    return true;
}

/// Search the label in the bucket. The variable label(pointer) is NOT in the bucket, it is search using the
/// 	resorces and the route
bool LabelingAlgorithmNS::searchLabel(Label* label, Bucket& bucket)
{
    for(int i=0; i < bucket.sizeVetPtrLabel; ++i)
    {
        Label* labelAux = bucket.vetPtrLabel[i];
        if(doubleEqual(label->vetResources[0], labelAux->vetResources[0]) &&
           doubleEqual(label->vetResources[1], labelAux->vetResources[1]))
        {
            if(labelAux->tamRoute != label->tamRoute)
                return false;

            for(int j=0; j < label->tamRoute; ++j)
            {
                if(label->vetRoute[i] != labelAux->vetRoute[i])
                    return false;
            }

            return true;
        }
    }

    return false;
}

bool LabelingAlgorithmNS::
     labelingAlgorithmm(int numRes, int numCust, const Vet3D_ResCost& vetMatResCostForward,
                        const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& vetVetBound, int dest,
                        const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                        FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                        Eigen::VectorX<FloatType>& vetRedCost, bool exact, LabelingTypeAlg labelingTypeAlg)
{

    if(labelingTypeAlg == AlgForward)
    {
        static bool print = false;
        if(!print)
        {
            std::cout<<"Forward!\n\n";
            print = true;
        }
        return forwardLabelingAlgorithm(numRes, numCust, vetMatResCostForward, vetVetBound, dest, ngSet, lData, matColX,
                                        numSol, labelStart, NumMaxLabePerBucket, dominaceCheck, maxDist, vetRedCost,
                                        exact);
    }
    else
    {
        static bool print = false;
        if(!print)
        {
            std::cout<<"BIDIRECTIONAL!\n\n";
            print = true;
        }

        return bidirectionalAlgorithm(numRes, numCust, vetMatResCostForward, vetMatResCostBackward, vetVetBound, dest,
                                      ngSet, lData, matColX, numSol, labelStart, NumMaxLabePerBucket, dominaceCheck,
                                      maxDist, vetRedCost, exact);

    }
}
