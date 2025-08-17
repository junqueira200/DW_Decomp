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

// These routes will be tracked during the labeling algorithm
inline Vector<VectorI> vetRoutesG;
inline MemoryPool_NS::Pool<LabelingAlgorithmNS::Label> labelPoolG;
inline bool startLabelPoolG = false;

void LabelingAlgorithmNS::startGlobalMemory(const Vector<VectorI>& vetRoutes)
{
    vetRoutesG = Vector(vetRoutes);
    std::println("startGlobalMemory");

    for(const VectorI& vet:vetRoutesG)
    {
        std::cout<<vet<<"\n";
    }

    std::println();

}
void LabelingAlgorithmNS::invertRoutes(Vector<VectorI>& vetRoutes)
{

    for(VectorI& route:vetRoutes)
    {
        std::reverse(route.begin(), route.end());
    }

}

void LabelingAlgorithmNS::addToVetRoutesG(const VectorI& route)
{
    vetRoutesG.push_back(route);
}

void LabelingAlgorithmNS::cleanVetRouteG()
{
    vetRoutesG = Vector<VectorI>();
}

LabelingAlgorithmNS::NgSet::NgSet()
{
    //vetNgSet  = Eigen::VectorX<std::bitset<NumMaxCust>>(NumMaxCust);
    vetNgSet = VetBitSet(NumMaxCust);
}


LabelingAlgorithmNS::NgSet::NgSet(int numCust_, int ngSetSize_)
{
    numCust   = numCust_;
    ngSetSize = ngSetSize_;
    vetNgSet  = VetBitSet(numCust);
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
        vetNgSet[i] = 0;
        for(int j=0; j < ngSetSize; ++j)
        {
            vetNgSet[i][vet[j].i] = true;
        }

        if((i+1) == numCust)
            break;
    }

}

Label* LabelingAlgorithmNS::getLabel()
{
    return labelPoolG.getT();
}

void LabelingAlgorithmNS::rmLabel(Label* label)
{
    labelPoolG.delT(label);
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
bool LabelingAlgorithmNS::bidirectionalAlgorithm(const int numRes, const int numCust, const Vet3D_ResCost& vetMatResCostForward,
                          const Vet3D_ResCost& vetMatResCostBackward, const MatBoundRes& matBoundRes, const int dest,
                          const NgSet& ngSet, LabelingData& lData, Eigen::MatrixXd& matColX, int& numSol,
                          const FloatType labelStart, int NumMaxLabePerBucket, bool dominaceCheck, FloatType& maxDist,
                          Eigen::VectorX<FloatType>& vetRedCost, bool exact, LabelingTypeAlg labelingTypeAlg, bool arc,
                          Eigen::VectorX<FloatType>* vetLowerBoundRedCost)
{
    //vetRoteG = {19, 18, 5, 10, 1, 12, 9, 17, 9, 17, 9, 17};


    //static Eigen::VectorX<FloatType> vetLowerBoundDist(numCust);
    //const bool vetLowerBoundDistValid =  getDistLowerBound(vetMatResCost, vetLowerBoundDist);
    VetBackwardMask vetBackwardMask;
    //matColX.setZero();
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

    if(Print && PrintG)
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

    static LabelHeap labelHeap(100000);
    //labelHeap.vet.setAll(nullptr);
    labelHeap.heapSize = 0;

    Label* labelPtr  = labelPoolG.getT();
    Label* labelBackwardPtr = labelPoolG.getT();

    labelPtr->vetResources.setZero();
    labelBackwardPtr->vetResources.setZero();

    int i = lData.getIndex(0, labelStart);
    int j = lData.getIndex(1, 0.0);

    int nextIndex = 0;


    if(labelingTypeAlg == AlgForward || labelingTypeAlg == AlgBidirectional)
    {

        // TODO create a start label
        labelPtr->tamRoute = 1;
        labelPtr->vetRoute[0] = 0;
        labelPtr->vetResources.setZero();
        labelPtr->vetResources[0] = labelStart;
        labelPtr->bitSetNg = 0;
        labelPtr->bitSetNg[0] = true;

        labelPtr->i         = i;
        labelPtr->j         = j;
        labelPtr->cust      = 0;
        labelPtr->active    = true;
        labelPtr->posBucket = 0;
        labelPtr->typeLabel = Forward;
        labelPtr->index = nextIndex++;

        lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel.resize(1);
        lData.vetMatBucketForward[0].mat(i, j).vetPtrLabel[0]  = labelPtr;
        lData.vetMatBucketForward[0].mat(i, j).sizeVetPtrLabel = 1;

        labelHeap.insertKey(labelPtr);

    }

    if(labelingTypeAlg == AlgBackward || labelingTypeAlg == AlgBidirectional)
    {
        startBackwardLabel(labelBackwardPtr, vetBackwardMask, matBoundRes, numRes, dest, labelStart, i, lData);
        labelHeap.insertKey(labelBackwardPtr);
        invertRoutes(vetRoutesG);
        labelBackwardPtr->index = nextIndex++;

        for(VectorI& route:vetRoutesG)
        {
            route[0] = dest;
        }
    }



    //labelHeap.insertKey(labelForwardPtr);
    labelPtr = nullptr;
    labelBackwardPtr = nullptr;
    int numIt = 0;
    int maxSize = 0;
    int localNumMaxLabel = numMaxLabelG;
    Label* ptrLabelTarget = nullptr;


    static ArrayResources vetMaxResources;
    static bool vetMaxResourcesSet = false;

    if(!vetMaxResourcesSet)
    {
        vetMaxResources.setConstant(-MaxFloatType);
        vetMaxResources[0] = (FloatType)0.0;

        for(int r=1; r < numRes; ++r)
        {
            for(int i=0; i < numCust; ++i)
            {
                if(matBoundRes(i, r).upperBound > vetMaxResources[r])
                    vetMaxResources[r] = matBoundRes(i, r).upperBound;
            }
        }

        vetMaxResourcesSet = true;
        lData.vetMatBucketForward[dest].mat(0, 0).vetPtrLabel.resize(DW_DecompNS::NumMaxSolSubProb+5);
    }

    //std::cout<<"Max Resources: "<<vetMaxResources<<"\n";

    while(!labelHeap.empty() &&
          ((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel < DW_DecompNS::NumMaxSolSubProb) || exact))
    {
//std::cout<<"While\n";
        // 																 10                         25


        if((lData.vetMatBucketForward[dest].mat(0, 0).sizeVetPtrLabel >= 5 && labelHeap.heapSize >= 50 * numCust &&
            !exact))
        {
            break;
        }


        //checkHeap(labelHeap, lData);
        //lData.checkVetMatBucketBackward();

        if(ptrLabelTarget)
        {
            std::cout<<"ptrLabelTarget("<<ptrLabelTarget->active<<"): "<<*ptrLabelTarget<<"\n\n";
        }


        /*
        if(labelHeap.heapSize > localNumMaxLabel && DominaIterBuckets)
        {
            lData.dominanceInterBuckets(labelHeap, numRes, localNumMaxLabel, lData.vetMatBucketForward, Backward);

            if(labelHeap.heapSize > localNumMaxLabel)
            {
                int max = std::max(labelHeap.heapSize, localNumMaxLabel);
                localNumMaxLabel =  (int)(max * 1.1);
            }


            if(exactLabelingG)
            {
                std::cout << "heapSize: " << labelHeap.heapSize << "; maxHeapSize:" << localNumMaxLabel << "; maxDist:"
                          << maxDistG << "; minDist: " << minDistG << "\n";
            }

        }
        */

        //lData.checkLabels();

        maxSize = std::max(maxSize, labelHeap.heapSize);


        if(Print && PrintG)
            std::cout << "numIt: " << numIt << "\n";

        labelPtr = labelHeap.extractTop();
        if(labelPtr == nullptr)
        {

            std::cout<<"listLabel have a nullptr\n";
            PRINT_DEBUG("", "");
            throw "ERRO";
        }

        if(!labelPtr->active)
        {
            std::cout<<"All labels shuld be active!\n";
            PRINT_EXIT();
        }

        //bool haveRoute = labelHaveRoute(vetRoutesG, labelPtr);

        if((PrintG && Print) || exactLabelingG)
        {
            std::cout << "labelPtr: " << labelPtr << "\n";
            std::cout << *labelPtr << "\n\n";
        }

        //std::cout<<"Extract label("<<labelPtr->vetRoute<<")\n";
        int lastCust = labelPtr->cust;

        /*
        if(!labelPtr->active)
        {
            std::cout<<"labelPtr is not active!\n"<<*labelPtr<<"\n";
            PRINT_DEBUG("", "");
            throw "ERROR";

            labelPoolG.delT(labelPtr);
            continue;
        }
        */

        if(lastCust == dest && labelPtr->typeLabel == Forward)
        {
            if(Print && PrintG)
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
               vetMatResCostForward(labelPtr->cust, t, 0) == InfFloatType)
                continue;

            else if(labelPtr->typeLabel  == Backward &&
               vetMatResCostBackward(t, labelPtr->cust, 0) == InfFloatType)
                continue;


            if(labelPtr->bitSetNg[t] == 1 || t == labelPtr->cust)
                continue;

            if(Print && PrintG)
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

                labelPtrAux->index = nextIndex++;
                bool haveRoute = false;
                if(TrackingRoutes)
                {
                    haveRoute = labelHaveRoute(vetRoutesG, labelPtrAux);
                    if(haveRoute)
                        std::cout<<"\tFind route("<<printRoute(labelPtrAux)<<") in vector of routes; Parcial route created!\n";

                }

                if(Print && PrintG)
                    std::cout<<"\t"<<*labelPtrAux<<"\n";

                int correctPos = 0;

                if(vetLowerBoundRedCost)
                {
                    if((labelPtrAux->vetResources[0] + (*vetLowerBoundRedCost)[t]) >= -DW_DecompNS::TolObjSubProb)
                    {
                        static bool print = false;
                        if(!print)
                        {
                            std::cout<<"Funcionou\n";
                            print = true;
                        }

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }
                }

                if(t == dest && labelPtrAux->typeLabel == Forward)
                {
                    if(labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                    {

                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }

                    //removeCycles2(*labelPtrAux, numCust);
                    //updateLabelCost(*labelPtrAux, vetMatResCostForward, labelStart);


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

                    //removeCycles2(*labelPtrAux, numCust);

                    // Reverse the array
                    std::reverse(labelPtrAux->vetRoute.begin(),
                                 labelPtrAux->vetRoute.begin()+labelPtrAux->tamRoute);


                    //updateLabelCost(*labelPtrAux, vetMatResCostForward, labelStart);


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

                /*
                if(tAux == dest)
                {
                    if(checkIfRouteIsInSet(labelPtrAux, vetSortRoute, multsetSourRoute, lData.matDist))
                    {
                        labelPoolG.delT(labelPtrAux);
                        continue;
                    }
                }
                */



                Bucket* bucket = nullptr;
                if(SortLabels)
                    bucket = dominanceIntraBucketFast1(tAux, labelPtrAux, lData, &labelHeap, numRes, dest, correctPos);
                else
                    bucket = dominanceIntraBucketSlow(tAux, labelPtrAux, lData, &labelHeap, numRes, dest, correctPos);

                if(!bucket)
                {
//std::cout<<"bucket is null\n";
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if((vetValueOfReducedCostsG.size() + 1) < vetValueOfReducedCostsG.capacity())
                    vetValueOfReducedCostsG.push_back(labelPtrAux->vetResources[0]);

                if((bucket->sizeVetPtrLabel+1) > NumMaxLabePerBucket && tAux != dest)
                {
                    if(SortLabels)
                    {
                        if(correctPos == (bucket->sizeVetPtrLabel))
                        {

                            labelPoolG.delT(labelPtrAux);
                            continue;
                        }

                        Label* del = bucket->vetPtrLabel[bucket->sizeVetPtrLabel-1];
                        bucket->removeElement(bucket->sizeVetPtrLabel-1);
                        labelHeap.deleteKey(del->posHeap);
                        labelPoolG.delT(del);

                    }
                    else
                    {
                        //labelPoolG.delT(labelPtrAux);
                        //continue;

                        int indexGreater = 0;
                        FloatType val    = bucket->vetPtrLabel[0]->vetResources[0];

                        for(int i=1; i < bucket->sizeVetPtrLabel; ++i)
                        {
                            if(bucket->vetPtrLabel[i]->vetResources[0] > val)
                            {
                                val = bucket->vetPtrLabel[i]->vetResources[0];
                                indexGreater = i;
                            }
                        }

                        if(labelPtrAux->vetResources[0] > val)
                        {
                            labelPoolG.delT(labelPtrAux);
                            continue;
                        }
                        else
                        {

                            Label* del = bucket->vetPtrLabel[indexGreater];
                            bucket->removeElement(indexGreater);

                            if(del != labelHeap.vet[del->posHeap])
                            {
                                std::cout<<"ERROR. labelHeap is Wrong";
                                PRINT_EXIT();
                            }

                            labelHeap.deleteKey(del->posHeap);
                            memset((void*)del, 0, sizeof(Label));
                            labelPoolG.delT(del);
                            //std::cout<<"Delete label greater the labelAuxPtr\n";
                        }

                        if(correctPos > bucket->sizeVetPtrLabel)
                            correctPos = bucket->sizeVetPtrLabel;

                    }
                }

                bucket->addElement(correctPos, labelPtrAux);

                //std::cout<<"Insert in position: "<<correctPos<<"\n";

//std::cout<<"\tInsert in bucket\n";

                if(tAux != dest && labelPtrAux->typeLabel == Forward)
                {
                    labelHeap.insertKey(labelPtrAux);
//std::cout<<"Insert\n";

                }
                else if(labelPtrAux->typeLabel == Backward)
                    labelHeap.insertKey(labelPtrAux);

                if((labelPtrAux->cust != 0 || labelPtrAux->cust != dest) && labelingTypeAlg == AlgBidirectional)
                {
                    int numLabels = lData.doMerge(labelPtrAux, vetMaxResources, matBoundRes, numRes);
                    if(numLabels >= DW_DecompNS::NumMaxSolSubProb && !exact)
                        break;
                }

                /*
                if(!checkIfAllLabelsInHeapHaveA_Bucket(labelHeap, lData))
                {
                    std::cout<<"ERROR in processing label: \n"<<*labelPtrAux<<"\n";
                    PRINT_EXIT();
                }
                */

                //checkHeap(labelHeap, lData);
            }
            else
            {
                labelPoolG.delT(labelPtrAux);
                if(Print && PrintG)
                    std::printf("\tFailed\n");
            }

        }

//std::cout<<"END Extend label\n";

        if(Print && PrintG)
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

//std::cout<<"Is heap empty: "<<labelHeap.empty()<<"\n";
    }

//std::cout<<"End while\n";

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
        std::cout<<numSol<<"\n";

        // TODO remover
        //exit(-1);
        //std::cout<<"*************************************************************************************\n\n";
        return true;
    }
    else
    {
        std::cout<<"0\n";
        return false;
    }
}

void LabelingAlgorithmNS::startPool()
{

    if(!startLabelPoolG)
    {
        startLabelPoolG = true;
        labelPoolG.startPool(5, 100);
    }
}

std::string LabelingAlgorithmNS::printRoute(Label* label)
{
    std::string str;

    for(int i=0; i < label->tamRoute; ++i)
        str += std::to_string(label->vetRoute[i]) + " ";

    return str;
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

    //newLabel.typeLabel = label.typeLabel;

    for(int i=0; i < numResources; ++i)//NumMaxResources; ++i)
    {

        // Extend the iº resource
        //std::cout<<"\t\tresource("<<label.vetResources[i]<<") + "<<vetMatResCost(custI, t, i)<<"\n";
        newLabel.vetResources[i] = label.vetResources[i] + vetMatResCost(custI, t, i);
        boundOk = boundOk && (newLabel.vetResources[i] <= vetVetBound(t, i).upperBound);
    }

    //std::cout<<"\n";

    if(!boundOk)
    {
        return false;
    }

    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ(t) its in custI ngSet
    newLabel.bitSetNg[t] = ngSet.contain(t, custI);



//    for(int i=0; i < label.tamRoute; ++i)
//        newLabel.vetRoute[i] = label.vetRoute[i];

    memcpy(&newLabel.vetRoute[0], &label.vetRoute[0], label.tamRoute*sizeof(int));

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
    //newLabel.typeLabel = label.typeLabel;
    newLabel.vetResources[0] = label.vetResources[0] + vetMatResCostBackward(t, custI, 0);

    if(Print && PrintG)
        std::cout<<"\t\t("<<t<<","<<custI<<"): "<<vetMatResCostBackward(t, custI, 0)<<"\n";


    for(int i=1; i < numResources; ++i)//NumMaxResources; ++i)
    {
        if(Print && PrintG)
            std::printf("\t\tcost: %.1f\n\t\tUP: %.1f\n\t\tLB: %.1f\n\n",
                        vetMatResCostBackward(t, custI, i), vetVetBound(t, i).upperBound,
                        vetVetBound(t, i).lowerBound);

        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] - vetMatResCostBackward(t, custI, i);
        boundOk = boundOk && (newLabel.vetResources[i] <= vetVetBound(t, i).upperBound);
        boundOk = boundOk && (newLabel.vetResources[i] >= vetVetBound(t, i).lowerBound);
    }

    if(!boundOk)
        return false;

    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ(t) its in custI ngSet
    newLabel.bitSetNg[custI] = ngSet.contain(custI, t);




    //for(int i=0; i < label.tamRoute; ++i)
    //    newLabel.vetRoute[i] = label.vetRoute[i];

    memcpy(&newLabel.vetRoute[0], &label.vetRoute[0], label.tamRoute*sizeof(int));

    newLabel.vetRoute[label.tamRoute] = t;
    newLabel.tamRoute  = label.tamRoute+1;
    newLabel.active    = true;
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

    lData.vetMatBucketBackward[dest].mat(i, j).vetPtrLabel.resize(1);
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
LabelingAlgorithmNS::LabelingData::LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_, int numMainResources_, int numCust_,
                     const ArrayResources& vetMaxResources_)
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

    vetMaxResources = vetMaxResources_;

    // Create the index for matForwardRange and matBackwardRang

    Label label;
    label.cust = 1;
    label.active = true;

    matBackwardRange.resize(vetNumSteps[0], vetNumSteps[1]);
    matBackwardRange.setConstant({{-1,-1}, {-1, -1}});

    matForwardRange.resize(vetNumSteps[0], vetNumSteps[1]);
    matForwardRange.setConstant({{-1,-1}, {-1, -1}});

    for(int i=0; i < vetNumSteps(0); ++i)
    {
        for(int j=0; j < vetNumSteps(1); ++j)
        {
            label.vetResources[0] = vetMatBound[0](i, j).upperBound;
            label.vetResources[1] = vetMatBound[1](i, j).upperBound;

            label.typeLabel = Forward;
            matForwardRange(i, j) = getListOfIndexForMerge(label);
            //std::cout<<printIndex(matForwardRange(i, j))<<"\n\n";

            label.typeLabel = Backward;
            matBackwardRange(i, j) = getListOfIndexForMerge(label);


        }
    }

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
    if(Print && PrintG)
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
        bool find = false;

        for(int i=0; i < bucket->sizeVetPtrLabel; ++i)
        {
            if(label == bucket->vetPtrLabel[i])
            {
                std::cout<<"Achou na posicao errada!\n";
                find = true;
                break;
            }
        }

        if(!find)
        {
            std::cout<<"Label IS NOT in this Bucket!\n\n";
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
    labelPtr->posBucket = sizeVetPtrLabel;
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

    if(Print && PrintG)
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

            if(Print && PrintG)
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
    return;
    if(Print && PrintG)
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

                        //if(typeLabel == Forward && ii > i && jj > j)
                        //    completeDominance = false;
                        //else if(typeLabel == Forward && ii > i && jj < j)
                        //    completeDominance = false;

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
                                        if(Print && PrintG)
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

    if(Print && PrintG)
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

bool LabelingAlgorithmNS::labelHaveRoute(Vector<VectorI>& vetRoute, Label *label)
{
    if(vetRoute.empty() || !TrackingRoutes)
        return false;

    for(const VectorI& route:vetRoute)
    {
        int min = std::min((int)route.size(), label->tamRoute);
        for(int i=0; i < min; ++i)
        {
            if(route[i] != label->vetRoute[i] && route[i] != 0)
                return false;
        }
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

        if(Print && PrintG)
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

/**
 * @brief LabelingAlgorithmNS::dominanceIntraBucketBackward
 * @param cust
 * @param labelPtrAux
 * @param lData
 * @param labelHeap
 * @param numRes
 * @param dest
 * @param correctPos
 * @return
 */
Bucket* LabelingAlgorithmNS::
        dominanceIntraBucketBackward(int cust, Label* labelPtrAux, LabelingData& lData, LabelHeap& labelHeap, int numRes,
                                     int dest, int& correctPos)
{

    //std::cout<<"dominanceIntraBucketBackward\n\n";

    const bool rmFromHeap = (cust != dest);
    Bucket* bucketPtr = nullptr;

    // RM ini

    /*
    Label l0;
    Label l1;

    l0.active          = true;
    l0.vetResources[0] = -14007.0;
    l0.vetResources[1] = 600.0;

    l1.active          = true;
    l1.vetResources[0] = -14008.0;
    l1.vetResources[1] = 600.0;

    std::cout<<"Aqui\n";

    bucketPtr = new Bucket();
    bucketPtr->addLabel(&l0);
    bucketPtr->sizeVetPtrLabel = 1;
    labelPtrAux = &l1;

    std::cout<<bucketPtr->print(2)<<"\n";
    */
    // RM end

    if(cust != dest || true)
    {
        // Find index from resources
        int i = lData.getIndex(0, labelPtrAux->vetResources[0]);
        int j = 0;
        if(lData.numMainResources > 1)
            j = lData.getIndex(1, labelPtrAux->vetResources[1]);

        labelPtrAux->i = i;
        labelPtrAux->j = j;
        labelPtrAux->cust = cust;

        if(Print && PrintG)
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


    //std::cout<<"\t"<<bucketPtr->print(2)<<"\n";

    // Get the correct possition
    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleGreaterEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0],
                              FloatEp))
            break;
        correctPos += 1;
    }

    //std::cout<<"correctPos: "<<correctPos<<"\n";

    int before = correctPos;

    // Search for labels with have equal resources
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], FloatEp))
        {
            //FloatType a = bucketPtr->vetPtrLabel[correctPos]->vetResources[0];
            //FloatType b = labelPtrAux->vetResources[0];

            //std::cout<<"a == b("<<(a == b)<<"float: "<<doubleEqual(a, b, FloatEp)<<")\n";
            //std::cout<<"a = "<<a<<"; b = "<<b<<"\n\n";

            //std::cout<<"equal\n";
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
                    //std::cout<<"\nRet null\n";
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

    //std::cout<<"correctPos: "<<correctPos<<"\n";

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
                //std::cout<<"\nRet null\n";
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

    //std::cout<<"Correct Position is: "<<correctPos<<"\n\n";
    //PRINT_EXIT();

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
            if(l0.vetResources[i] > l1.vetResources[i])
            //if(doubleGreater(l0.vetResources[i], l1.vetResources[i], FloatEp))
                return false;
        }
    }
    else
    {
        if(l0.vetResources[0] > l1.vetResources[0])
        //if(doubleGreater(l0.vetResources[0], l1.vetResources[0], FloatEp))
            return false;

        // Check the resources
        //#pragma GCC unroll NumMaxResources
        for(int i=1; i < numResources; ++i)
        {
            if(l0.vetResources[i] < l1.vetResources[i])
            //if(doubleLess(l0.vetResources[i], l1.vetResources[i], FloatEp))
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

Bucket* LabelingAlgorithmNS::dominanceIntraBucketSlow(int cust, Label* label, LabelingData& lData, LabelHeap* labelHeap,
                                                      int numRes, int dest, int& correctPos)
{
    //std::cout<<"\t"<<lData.vetMatBucketBackward[15].mat(0,0).print(2)<<"\n\n";

    bool rmFromHeap = true;
    Bucket* bucketPtr = nullptr;

    if((label->typeLabel == Forward && cust == dest) || !labelHeap)
        rmFromHeap = false;
    else if(label->typeLabel == Backward && cust == 0)
        rmFromHeap = false;

    if(cust != dest)
    {
        // Find index from resources
        int i = lData.getIndex(0, label->vetResources[0]);
        int j = 0;
        if(lData.numMainResources > 1)
            j = lData.getIndex(1, label->vetResources[1]);

        label->i = i;
        label->j = j;
        label->cust = cust;

        if(Print && PrintG)
            std::cout << "\t\ti(" << i << "); j(" << j << ")\n";

        if(label->typeLabel == Backward)
            bucketPtr = &lData.vetMatBucketBackward[cust].mat(i, j);
        else if(label->typeLabel == Forward)
            bucketPtr = &lData.vetMatBucketForward[cust].mat(i, j);
        else
        {
            std::printf("Error, typeLabel is not Backward or Forward!\n");
            PRINT_EXIT();
        }
    }
    else
    {

        label->i = 0;
        label->j = 0;
        label->cust = cust;
        bucketPtr = &lData.vetMatBucketForward[cust].mat(0, 0);

        if(label->typeLabel == Backward)
        {
            std::cout<<"ERROR?\nlabel type equal to Backward, and cust("<<cust<<") equal to dest("<<dest<<")\n\n";
            PRINT_EXIT();
        }

    }

    //std::cout<<"\t\t* "<<bucketPtr->print(2)<<"\n\n";

    for(int i=0; i < bucketPtr->sizeVetPtrLabel; ++i)
    {
        /*
        if(bucketPtr->vetPtrLabel[i]->cust != cust)
        {
            std::cout<<"CUST DIFFERENTS!";
            PRINT_EXIT();
        }
        */

        if(checkCompleteDominance(*bucketPtr->vetPtrLabel[i], *label, numRes))
        {
            correctPos = -1;
            return nullptr;
        }
    }

    int i=0;
    while(i < bucketPtr->sizeVetPtrLabel)
    {
        if(checkCompleteDominance(*label, *bucketPtr->vetPtrLabel[i], numRes))
        {
            if(TrackingRoutes)
            {
                if(labelHaveRoute(vetRoutesG, bucketPtr->vetPtrLabel[i]))
                {
                    std::println("The route(%s) dominates the route(%s)(in vetRoutesG)", printRoute(label),
                                 printRoute(bucketPtr->vetPtrLabel[i]));
                    PRINT_EXIT();
                }
            }

            if(rmFromHeap)
                labelHeap->deleteKey(bucketPtr->vetPtrLabel[i]->posHeap);

            labelPoolG.delT(bucketPtr->vetPtrLabel[i]);

            std::swap(bucketPtr->vetPtrLabel[i], bucketPtr->vetPtrLabel[bucketPtr->sizeVetPtrLabel-1]);

            /*
            for(int t=i; t < (bucketPtr->sizeVetPtrLabel-1); ++t)
            {
                bucketPtr->vetPtrLabel[t] = bucketPtr->vetPtrLabel[t+1];
                bucketPtr->vetPtrLabel[t]->posBucket = t;
            }
            */

            bucketPtr->sizeVetPtrLabel -= 1;

            if(i <= (bucketPtr->sizeVetPtrLabel-1) && i >= 0 && bucketPtr->sizeVetPtrLabel >= 0)
                bucketPtr->vetPtrLabel[i]->posBucket = i;


            //std::cout<<"rm; size("<<bucketPtr->sizeVetPtrLabel<<"); i("<<i<<")\n";
        }
        else
            i += 1;
    }


    correctPos = bucketPtr->sizeVetPtrLabel;
    return bucketPtr;

}


Bucket* LabelingAlgorithmNS::dominanceIntraBucketFast1(int cust, Label* labelPtrAux, LabelingData& lData,
                                                       LabelHeap* labelHeap, int numRes, int dest, int& correctPos)
{


    bool rmFromHeap = true;
    Bucket* bucketPtr = nullptr;

    if((labelPtrAux->typeLabel == Forward && cust == dest) || !labelHeap)
        rmFromHeap = false;
    else if(labelPtrAux->typeLabel == Backward && cust == 0)
        rmFromHeap = false;


    if(cust != dest)
    {
        // Find index from resources
        int i = lData.getIndex(0, labelPtrAux->vetResources[0]);
        int j = lData.getIndex(1, labelPtrAux->vetResources[1]);

        labelPtrAux->i = i;
        labelPtrAux->j = j;
        labelPtrAux->cust = cust;

        if(Print)
            std::cout << "\t\ti(" << i << "); j(" << j << ")\n";

        bucketPtr = lData.getBucket(labelPtrAux);
    }
    else
    {

        labelPtrAux->i = 0;
        labelPtrAux->j = 0;
        labelPtrAux->cust = cust;
        bucketPtr = lData.getBucket(labelPtrAux);
    }

    if(bucketPtr->sizeVetPtrLabel == 0)
    {
        correctPos = 0;
        return bucketPtr;
    }

    correctPos = 0;
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleLess(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], 1E-5))
            correctPos += 1;
        else
            break;
    }


    /*
    while(correctPos < bucketPtr->sizeVetPtrLabel)
    {
        if(doubleEqual(bucketPtr->vetPtrLabel[correctPos]->vetResources[0], labelPtrAux->vetResources[0], 1E-5))
        {
            if(doubleLessEqual(labelPtrAux->vetResources[1], bucketPtr->vetPtrLabel[correctPos]->vetResources[1],
                               1E-5))
            correctPos += 1;
        else
            break;
        }
        else
            break;
    }
    */

    // Labels before correctPos can dominate labelPtrAux;
    for(int ii=0; ii < correctPos; ++ii)
    {
        bool canDomindate = true;

        for(int r=1; r < numRes; ++r)
        {
            if(!doubleLess(bucketPtr->vetPtrLabel[ii]->vetResources[r], labelPtrAux->vetResources[r], 1E-5))
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
                correctPos = -1;
                return nullptr;
            }
        }
    }

    // Search for labels with have equal resources
    int pos = correctPos;

    if(pos >= bucketPtr->sizeVetPtrLabel)
    {
        return bucketPtr;
    }

    while(labelEqual(labelPtrAux, bucketPtr->vetPtrLabel[pos], numRes))
    {
        Label* bucketLabel = bucketPtr->vetPtrLabel[pos];

        //auto end = labelPtrAux->bitSetNg & bucketLabel->bitSetNg;
        // Check if labelPtrAux dominate the label in bucketPtr
        //if(end == labelPtrAux->bitSetNg)
        if(checkDominanceSubSet(*labelPtrAux, *bucketLabel))
        {
            bucketLabel->i = -1;
            bucketLabel->j = -1;

            if(rmFromHeap)
                labelHeap->deleteKey(bucketLabel->posHeap);

            labelPoolG.delT(bucketLabel);
            bucketPtr->removeElement(pos);


        }
        // Chck if the label in bucketPtr dominate labelPtrAux
        //else if(end == bucketLabel->bitSetNg)
        else if(checkDominanceSubSet(*bucketLabel, *labelPtrAux))
        {
            correctPos = -1;
            return nullptr;
        }
        else
            pos += 1;

        if(pos >= bucketPtr->sizeVetPtrLabel)
            break;

    }


    // Labels after correctPos (inclusive) can be dominated by labelPtrAux
    int ii = pos;

    while(ii < bucketPtr->sizeVetPtrLabel)
    {
        bool canDominate = true;
        for(int r=1; r < numRes; ++r)
        {
            if(doubleLess(labelPtrAux->vetResources[r], bucketPtr->vetPtrLabel[ii]->vetResources[r], 1E-5))
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

                bucketPtr->vetPtrLabel[ii]->i = -1;
                bucketPtr->vetPtrLabel[ii]->j = -1;

                if(rmFromHeap)
                    labelHeap->deleteKey(bucketPtr->vetPtrLabel[ii]->posHeap);

                labelPoolG.delT(bucketPtr->vetPtrLabel[ii]);
                bucketPtr->removeElement(ii);

                continue; // while ii
            }
        }

        ii += 1;
    }

    if(labelPtrAux == nullptr)
        return nullptr;

    labelPtrAux->active = true;

    return bucketPtr;

}

Bucket* LabelingAlgorithmNS::dominanceIntraBucketFast0(int cust, Label* labelPtrAux, LabelingData& lData,
                                                       LabelHeap* labelHeap, int numRes, int dest, int& correctPos)
{

    const bool rmFromHeap = (cust != dest) && labelHeap;
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

        if(Print && PrintG)
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

    int pos = 0;

    if(bucketPtr->sizeVetPtrLabel == 0)
    {
        correctPos = 0;
        return bucketPtr;
    }

    while(labelLessForward(bucketPtr->vetPtrLabel[pos], labelPtrAux, numRes))
    {
        if(checkDominanceSubSet(*bucketPtr->vetPtrLabel[pos], *labelPtrAux))
            return nullptr;

        if(pos == (bucketPtr->sizeVetPtrLabel-1))
        {
            pos += 1;
            break;
        }

        pos += 1;
    }

    correctPos = pos;

    if(pos >= (bucketPtr->sizeVetPtrLabel-1))
        return bucketPtr;


    while(labelEqual(bucketPtr->vetPtrLabel[pos], labelPtrAux, numRes))
    {
        auto resultAnd = (bucketPtr->vetPtrLabel[pos]->bitSetNg & labelPtrAux->bitSetNg);

        // Checks if label in bucketPtr dominate labelPtrAux
        if(resultAnd == bucketPtr->vetPtrLabel[pos]->bitSetNg)
        {
            correctPos = -1;
            return nullptr;
        }
        // Checks if labelPtrAux dominate label in bucketPtr
        else if(resultAnd == labelPtrAux->bitSetNg)
        {
            Label* labelI = bucketPtr->vetPtrLabel[pos];
            bucketPtr->removeElement(pos);

            if(rmFromHeap)
                labelHeap->deleteKey(labelI->posHeap);

            labelPoolG.delT(labelI);
        }
        else
            pos += 1;

        if(pos == bucketPtr->sizeVetPtrLabel)
            break;

    }

    if(pos >= (bucketPtr->sizeVetPtrLabel-1))
        return bucketPtr;

    while(pos < bucketPtr->sizeVetPtrLabel)
    {
        if(labelLessForward(labelPtrAux, bucketPtr->vetPtrLabel[pos], numRes))
        {
            if(checkDominanceSubSet(*labelPtrAux, *bucketPtr->vetPtrLabel[pos]))
            {
                if(rmFromHeap)
                    labelHeap->deleteKey(bucketPtr->vetPtrLabel[pos]->posHeap);

                labelPoolG.delT(bucketPtr->vetPtrLabel[pos]);
                bucketPtr->removeElement(pos);
            }
            else
                pos += 1;
        }
        else
            pos += 1;
    }

    return bucketPtr;
}


void LabelingAlgorithmNS::changeTypeAlg(LabelingTypeAlg& labelingTypeAlg)
{

    if(labelingTypeAlg == AlgForward)
        labelingTypeAlg = AlgBackward;
    else if(labelingTypeAlg == AlgBackward)
        labelingTypeAlg = AlgForward;

}

void LabelingAlgorithmNS::writeNgSet(Label* label, const NgSet& ngSet)
{
    label->bitSetNg = 0;
    for(int i=0; i < (label->tamRoute-1); ++i)
    {
        int cliI = label->vetRoute[i];
        int cliJ = label->vetRoute[i+1];

        if(ngSet.contain(cliJ, cliI))
            label->bitSetNg[cliJ] = true;
    }
}

Label* LabelingAlgorithmNS::
       mergeForwardAndBackward(Label* forwardPtr, Label* backwardPtr, const ArrayResources& vetMaxResources,
                               const MatBoundRes& vetVetBound, int numResorces)
{
    if(forwardPtr->typeLabel == Backward)
        std::swap(forwardPtr, backwardPtr);

    Label& forward  = *forwardPtr;
    Label& backward = *backwardPtr;

    if(forward.typeLabel != Forward)
    {
        std::printf("Error\n foreard label is not of the correct type!\n");
        PRINT_THROW();
    }

    if(backward.typeLabel != Backward)
    {
        std::printf("Error\n backward label is not of the correct type!\n");
        PRINT_THROW();
    }

    if(forward.cust != backward.cust)
    {
        std::printf("");
        PRINT_THROW();
    }

    Label* result = labelPoolG.getT();
    result->typeLabel = Forward;
    result->vetResources[0] = forward.vetResources[0] + backward.vetResources[0];

    if(result->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
    {
        labelPoolG.delT(result);
        return nullptr;
    }

    bool resourcesOk = true;
    for(int i=1; i < numResorces; ++i)
    {
        result->vetResources[i] = forward.vetResources[i] + (vetMaxResources[i]-backward.vetResources[i]);
        resourcesOk = resourcesOk && result->vetResources[i] <= vetVetBound(backward.cust, i).upperBound;
        resourcesOk = resourcesOk && result->vetResources[i] >= vetVetBound(backward.cust, i).lowerBound;
    }

    if(!resourcesOk)
    {
        labelPoolG.delT(result);
        return nullptr;
    }

    // Copy route
    for(int i=0; i < forward.tamRoute; ++i)
        result->vetRoute[i] = forward.vetRoute[i];

    int tam = forward.tamRoute;
    for(int i=(backward.tamRoute-2); i >= 0; --i)
    {
        result->vetRoute[tam] = backward.vetRoute[i];
        tam += 1;
    }

    result->tamRoute = tam;
    result->active   = true;
    result->bitSetNg = forward.bitSetNg | backward.bitSetNg;
    result->i = 0;
    result->j = 0;
    result->cust = result->vetRoute[tam-1];
    result->posHeap = -1;

    //removeCycles2(*result, numCust);
    //updateLabelCost(*result, vetMatResCostForward, startLabel);

    if(result->vetResources[0] >= DW_DecompNS::TolObjSubProb)
    {
        labelPoolG.delT(result);
        result = nullptr;
    }

    return result;
}



void LabelingAlgorithmNS::resetLabelPool()
{
    labelPoolG.resetPool(false);
}

HashSortRoute::HashSortRoute(SortRoute& sortRoute_):sortRoute(sortRoute_)
{
    // https://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-16.html
    // ELF Hash algorithm
    hash = 0;
    for(int i=0; i < sortRoute.vetSortRoute.size(); ++i)
    {   //         valHash * 16
        hash = (hash<<4) + sortRoute.vetRoute[i];
        uint64_t g = hash & 0xF0000000L;

        if(g != 0)
            hash ^= g >> 24;
        hash &= ~g;
    }


}

bool LabelingAlgorithmNS::checkIfRouteIsInSet(Label* label, VetSortRoute& vetSortRoute,
                                              MultsetSortRoute& multsetSourRoute, const EigenMatrixRowD& matDist)
{
    vetSortRoute.push_back(std::make_unique<SortRoute>());

    LabelingAlgorithmNS::SortRoute& sortRoute = *vetSortRoute[vetSortRoute.size()-1];
    sortRoute.vetSortRoute.resize(label->tamRoute);
    sortRoute.vetRoute.resize(label->tamRoute);

    for(int i=0; i < label->tamRoute; ++i)
    {
        int cust = label->vetRoute[i];
        if(i == (label->tamRoute) || i == 0)
            cust = 0;

        sortRoute.vetRoute[i]     = cust;
        sortRoute.vetSortRoute[i] = cust;

    }

    sortRoute.computeDistance(matDist);

    std::sort(sortRoute.vetSortRoute.begin(), sortRoute.vetSortRoute.end());
    HashSortRoute hashSortRoute(sortRoute);

    if(multsetSourRoute.count(hashSortRoute) >= 1)
    {
        auto it = multsetSourRoute.find(hashSortRoute);

        bool equalDist = false;

        while(it != multsetSourRoute.end())
        {

            if(doubleLessEqual((*it).sortRoute.dist, sortRoute.dist), 1E-5)
            {
                vetSortRoute.pop_back();
                return true;
            }

            ++it;
        }

        std::cout<<"Route Repeted\n";
    }

    multsetSourRoute.insert(hashSortRoute);
    return false;
}

void  LabelingAlgorithmNS::SortRoute::computeDistance(const EigenMatrixRowD& matDist)
{
    dist = 0.0;
    for(int i=0; i < (vetRoute.size()-1); ++i)
    {
        dist += matDist(vetRoute[i], vetRoute[i+1]);
    }
}


bool  LabelingAlgorithmNS::checkIfAllLabelsInHeapHaveA_Bucket(LabelHeap& labelHeap, LabelingData& lData)
{

    for(int i=0; i < labelHeap.heapSize; ++i)
    {
        Label* label = labelHeap.vet[i];
        Bucket* bucket = lData.getBucket(label);

        bool error = false;

        if(label->posBucket >= bucket->sizeVetPtrLabel)
        {
            error = true;
            std::cout<<"pos >= size\n";
            std::cout<<"pos: "<<label->posBucket<<"\n";
            std::cout<<"size: "<<bucket->sizeVetPtrLabel<<"\n\n";
        }
        else
        {
            if(bucket->vetPtrLabel[label->posBucket] != label)
                error = true;
        }

        if(error)
        {
            std::cout<<"Label was not find in bucket!\n";
            std::cout<<*label<<"\n\n";
            //PRINT_EXIT();
            return false;
        }


    }

    return true;
}





