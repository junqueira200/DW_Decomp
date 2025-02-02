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
#include <boost/container/set.hpp>

LabelingAlgorithmNS::NgSet::NgSet()
{
    matNgSet = Eigen::MatrixXi(NumMaxCust, NgSetSize);
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

bool
LabelingAlgorithmNS::forwardLabelingAlgorithm(const int                     numRes,
                                              const int                     numCust,
                                              const VetMatResCost&          vetMatResCost,
                                              const VetVetResBound&         vetVetBound,
                                              const int                     dest,
                                              const NgSet&                  ngSet,
                                              LabelingData&                 lData,
                                              Eigen::MatrixXd&              matColX,
                                              int&                          numSol,
                                              const FloatType               labelStart,
                                              int                           NumMaxLabePerBucket,
                                              bool                          dominaceCheck,
                                              FloatType&                    maxDist,
                                              Eigen::VectorX<FloatType>&    vetRedCost)
{
    //dominaceCheck = false;

    if(NumMaxLabePerBucket == -1)
        NumMaxLabePerBucket = std::numeric_limits<int>::max();

    //std::cout<<"NumMaxLabePerBucket: "<<NumMaxLabePerBucket<<"\n\n";

    if(checkDistance(vetMatResCost[0]))
    {
        std::cout<<"dist > 0\n";
        numSol = 0;
        return false;

    }

    numSol = 0;
    Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> vetLabel;
    vetLabel.setZero();

    maxDist = -std::numeric_limits<FloatType>::max();
    //redCost = std::numeric_limits<FloatType>::max();

    if(Print)
    {
        std::cout << "************************************************************\n";
        std::cout << "*****************FORWARD LABELING ALGORITHM*****************\n\n";
        std::cout << "numCust: " << numCust << "\n";
    }
    static bool labelPoolStart = false;//MemoryPool_NS::Pool<Label> labelPool(44, 400);
    if(!labelPoolStart)
    {
        labelPoolStart = true;
        labelPoolG.startPool(5, 100);
    }

    labelPoolG.resetPool(false);
    lData.flushLabel();

    //std::list<Label*> listLabel;
    // std::set
    std::multiset<Label*, LabelCmp> setLabel;

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

    lData.vetMatBucket[0].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucket[0].mat(i, j).vetPtrLabel[0] = labelPtr;
    lData.vetMatBucket[0].mat(i, j).sizeVetPtrLabel = 1;

    //listLabel.push_back(labelPtr);
    setLabel.insert(labelPtr);
    labelPtr = nullptr;

    int numIt = 0;

    Label *labelPtrBest = nullptr;
    int maxSize = 0;
    int maxSizeVetPtrLabel = 0;
    int localNumMaxLabel = NumMaxLabel;

    //while(!listLabel.empty() && !labelPtrBest)
    while(!setLabel.empty() && numSol < DW_DecompNS::NumMaxSolSubProb)
    {
        if(int(setLabel.size()) > localNumMaxLabel && DominaIterBuckets && dominaceCheck)
        {
            lData.dominanceInterBuckets(setLabel, numRes, localNumMaxLabel);

            while(int(setLabel.size()) > localNumMaxLabel)
            {
                //localNumMaxLabel += NumMaxLabel;
                localNumMaxLabel *= 2;
            }
        }

        if(numSol >= 1 && setLabel.size() > NumMaxLabel)
            break;

        maxSize = std::max(maxSize, int(setLabel.size()));

        if(Print)
            std::cout << "numIt: " << numIt << "\n";
        //labelPtr = listLabel.back();
        //listLabel.pop_back();
        labelPtr = (*setLabel.begin());
        setLabel.erase(setLabel.begin());

        if(labelPtr == nullptr)
        {

            std::cout<<"listLabel have a nullptr\n";
            PRINT_DEBUG("", "");
            throw "ERRO";
        }

        if(Print)
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
            if(t == lastCust || vetMatResCost[0](labelPtr->cust, t) == std::numeric_limits<FloatType>::infinity())
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

            if(extendLabel(*labelPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelPtr->cust, t, ngSet, numRes))
            {

                maxDist = std::max(maxDist, labelPtrAux->vetResources[0]);

                if(Print)
                    std::cout<<"\t\textendLabel "<<labelPtrAux<<": "<<*labelPtrAux<<"\n";

                // Find index from resources
                i = lData.getIndex(0, labelPtrAux->vetResources[0]);
                j = 0;
                if(lData.numMainResources > 1)
                    j = lData.getIndex(1, labelPtrAux->vetResources[1]);

                labelPtrAux->i = i;
                labelPtrAux->j = j;
                labelPtrAux->cust = t;

                if(Print)
                    std::cout<<"\t\ti("<<i<<"); j("<<j<<")\n";

                Bucket &bucket = lData.vetMatBucket[t].mat(i, j);
                int k=0;


                if(dominaceCheck)
                {

                    while(k < bucket.sizeVetPtrLabel)
                    {
                        if(Print)
                            std::cout << "\t\t\tcheckDominance " << bucket.vetPtrLabel[k] << ": "
                                      << *bucket.vetPtrLabel[k] << "\n";

                        if(checkDominance(*labelPtrAux, *bucket.vetPtrLabel[k], numRes) && !checkDominance(*bucket.vetPtrLabel[k], *labelPtrAux, numRes))
                        {
                            if(bucket.vetPtrLabel[k] == labelPtrBest)
                                labelPtrBest = nullptr;

                            setLabel.erase(bucket.vetPtrLabel[k]);
                            bucket.vetPtrLabel[k]->active = false;
                            labelPoolG.delT(bucket.vetPtrLabel[k]);

                            if(k == (bucket.sizeVetPtrLabel - 1))
                            {
                                //labelPool.delT(bucket.vetPtrLabel[k]);
                                bucket.vetPtrLabel[k] = nullptr;
                            } else
                            {
                                std::swap(bucket.vetPtrLabel[k], bucket.vetPtrLabel[bucket.sizeVetPtrLabel-1]);
                                //labelPool.delT(bucket.vetPtrLabel[bucket.sizeVetPtrLabel-1]);
                                bucket.vetPtrLabel[bucket.sizeVetPtrLabel-1] = nullptr;
                            }

                            bucket.sizeVetPtrLabel -= 1;
                            continue;
                        }
                        else if(checkDominance(*bucket.vetPtrLabel[k], *labelPtrAux, numRes))
                        {
                            if(Print)
                                std::cout << "\t\t\t\t<" << bucket.vetPtrLabel[k] << ">> domina <<" << labelPtrAux
                                          << ">>\n";

                            labelPoolG.delT(labelPtrAux);
                            labelPtrAux = nullptr;
                            break;
                        }

                        k += 1;
                    }
                }

                if(!labelPtrAux)
                    continue;

                if((bucket.sizeVetPtrLabel+1) > NumMaxLabePerBucket && t != dest)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if(t == dest && labelPtrAux->vetResources[0] >= -DW_DecompNS::TolObjSubProb)
                {

                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                if(t != dest)
                    maxSizeVetPtrLabel = std::max(maxSizeVetPtrLabel, bucket.sizeVetPtrLabel);


                labelPtrAux->active = true;
                //listLabel.push_back(labelPtrAux);

                if(labelPtrAux->cust == dest && labelPtrAux->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                {

                    // TODO print
                    if(!dominaceCheck)
                        std::cout<<*labelPtrAux<<"\n";

                    removeCycles(*labelPtrAux, numCust);
                    updateLabelCost(*labelPtrAux, vetMatResCost, labelStart);


                    // TODO print
                    //if(Print)
                    if(!dominaceCheck)
                        std::cout<<"*"<<*labelPtrAux<<"\n\n";

                    if(labelPtrAux->vetResources[0] < -DW_DecompNS::TolObjSubProb && !containRoute(vetLabel, numSol, labelPtrAux))
                    {
                        //std::cout<<"*"<<*labelPtrAux<<"\n\n";
                        vetLabel[numSol] = labelPtrAux;
                        vetRedCost[numSol] = labelPtrAux->vetResources[0];
                        numSol += 1;
                    }
                    else
                        labelPoolG.delT(labelPtrAux);
                }
                else
                {
                    bucket.addLabel(labelPtrAux);
                    setLabel.insert(labelPtrAux);
                }



/*                else if(labelPtrAux->cust == dest)
                {
                    std::cout<<"*"<<*labelPtrAux<<"\n\n";
                }*/

//std::cout<<"extendLabel to t("<<t<<")\n";

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

    // TODO print
    //std::cout<<"\n########################################################################################\n\n";

    //std::cout<<"Max dist: "<<maxDist<<"\n";
    //std::cout<<"maxSizeVetPtrLabel: "<<maxSizeVetPtrLabel<<"\n";

    if(numSol > 0)
    {
        for(int l=0; l < numSol; ++l)
        {
            vetLabel[l]->vetRoute[vetLabel[l]->tamRoute-1] = vetLabel[l]->vetRoute[0];

            //std::cout << "BEST LABEL: "<<vetLabel[l]<<" "<< *vetLabel[l] << "\n";

            auto &vetRoute = vetLabel[l]->vetRoute;

            for(int i = 0; i < (vetLabel[l]->tamRoute - 1); ++i)
            {
                matColX((getIndex(vetRoute[i], vetRoute[i+1], numCust-1)), l) = 1.0;
            }
        }

        //std::cout<<"Dest: "<<dest<<"\n\n";
        //std::cout<<"\n\n";

        return true;
    }
    else
    {
        return false;
    }

}


bool LabelingAlgorithmNS::extendLabel(const Label &label,
                                      Label &newLabel,
                                      const VetMatResCost& vetMatResCost,
                                      const VetVetResBound& vetVetBound,
                                      const int custI,
                                      const int custJ,
                                      const NgSet &ngSet,
                                      const int numResources)
{


    // Goes through resources
    for(int i=0; i < NumMaxResources; ++i)
    {

        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] + vetMatResCost[i](custI, custJ);


        // Check the bound
        if(newLabel.vetResources[i] > vetVetBound[i][custJ].upperBound)
        {
            if(Print)
                std::cout<<"\tUpperBound\n";
            return false;
        }

        else if(newLabel.vetResources[i] < vetVetBound[i][custJ].lowerBound)
        {
            if(Print)
                std::cout<<"\tLowerBound\n";
            return false;
            newLabel.vetResources[i] = vetVetBound[i][custJ].lowerBound;
        }

        if((i+1) == numResources)
            break;
    }

    newLabel.bitSetNg   = 0;
    newLabel.bitSetNg   = label.bitSetNg;

    // Checks if custJ its in custI ngSet
    if(ngSet.contain(custJ, custI))
    {
        newLabel.bitSetNg[custJ] = true;
    }



    if(((int)newLabel.vetRoute.size()) < label.tamRoute+1)
    {   std::cout<<"ini resize\n";
        //newLabel.vetRoute.resize(label.vetRoute.size() + 1);
        std::cout<<"end resize\n";
        throw "ERROR, OUT OF MEMORY";
    }

    // Copy route
    for(int i=0; i < label.tamRoute; ++i)
    {
        newLabel.vetRoute[i] = label.vetRoute[i];
    }


    if(Print)
        std::cout<<"APOS COPY ROUTE\n";

    newLabel.vetRoute[label.tamRoute] = custJ;
    newLabel.tamRoute = label.tamRoute+1;

    return true;

}

/**
 * @param vetStepSize_
 * @param numMainResources_
 * @param numCust_ Number of customers of the instance plus one
 */
LabelingAlgorithmNS::LabelingData::LabelingData(const Eigen::Vector<Step, 2> &vetStepSize_,
                                                int numMainResources_,
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

    vetMatBucket = Eigen::VectorX<MatBucket>(numCust);//(numMaxSteps, numMaxSteps);
    for(int i=0; i < numCust; ++i)
        vetMatBucket[i].mat = Eigen::Matrix<Bucket, -1, -1, Eigen::RowMajor>(vetNumSteps[0], vetNumSteps[1]);

    // Generates matBound
    matBound = Eigen::Matrix<Eigen::Vector<Bound, 2>, -1, -1, Eigen::RowMajor>(vetNumSteps[0], vetNumSteps[1]);

    for(int i=0; i < vetNumSteps[0]; ++i)
    {
        FloatType r00 = vetStepSize[0].start + i*vetStepSize[0].stepSize;

        FloatType r01 = r00 + vetStepSize[0].stepSize;
        for(int j=0; j < vetNumSteps[1]; ++j)
        {
            FloatType r10 = vetStepSize[1].start + j*vetStepSize[1].stepSize;
            FloatType r11 = r10 + vetStepSize[1].stepSize;

            if(i != 0)
                matBound(i, j)[0].lowerBound = r00;
            else
                matBound(i, j)[0].lowerBound = -std::numeric_limits<FloatType>::infinity();

            if(i != (vetNumSteps[0]-1))
                matBound(i, j)[0].upperBound = r01;
            else
                matBound(i, j)[0].upperBound = std::numeric_limits<FloatType>::infinity();

            matBound(i, j)[1].lowerBound = r10;
            matBound(i, j)[1].upperBound = r11;

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
    if(resource != 0 && resource != 1)
    {
        PRINT_DEBUG("", "Wrong resource(" << resource << ")");
        throw "ERROR";
    }

    if(numMainResources == 1 && resource == 1)
        return 0;

    int start = 0;
    int index = vetNumSteps[resource]/2;
    int end   = vetNumSteps[resource]-1;

    //int indexI, indexJ;
    int vetIndex[2];

    if(resource == 0)
    {
        vetIndex[0] = index;
        vetIndex[1] = 0;
    }
    else
    {
        vetIndex[0] = 0;
        vetIndex[1] = index;
    }

    while(!(val >= matBound(vetIndex[0], vetIndex[1])[resource].lowerBound &&
            val < matBound(vetIndex[0], vetIndex[1])[resource].upperBound))
    {
        //std::cout<<"index: "<<index<<"\n";

        if(val < matBound(vetIndex[0], vetIndex[1])[resource].lowerBound)
        {
            end = vetIndex[resource]-1;
            vetIndex[resource] = start + (end-start)/2;
        }
        else
        {
            start = vetIndex[resource]+1;
            vetIndex[resource] = start + (end-start)/2;
        }
    }

    //if(!(val >= matBound(vetIndex[0], vetIndex[1])[resource].lowerBound && val < matBound(vetIndex[0], vetIndex[1])[resource].upperBound))
    //    return -1;

    return vetIndex[resource];

}

void LabelingAlgorithmNS::LabelingData::flushLabel()
{

    for(MatBucket& mat:vetMatBucket)
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

    int pos = -1;

    for(int t=0; t < vetMatBucket[cust].mat(i, j).sizeVetPtrLabel; ++t)
    {
        if(vetMatBucket[cust].mat(i, j).vetPtrLabel[t] == label)
        {
            pos = t;
            break;
        }
    }

    if(pos == -1)
    {
        PRINT_DEBUG("", "Label not found");
        throw "ERRO";
    }

    int size = vetMatBucket[cust].mat(i, j).sizeVetPtrLabel;

    if(size > 1)
        std::swap(vetMatBucket[cust].mat(i, j).vetPtrLabel[pos],
                  vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1]);

    vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1] = nullptr;
    vetMatBucket[cust].mat(i, j).sizeVetPtrLabel -= 1;

    if(Print)
        std::cout<<"End removeLabel\n";
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
    labelPtr->pos = sizeVetPtrLabel;
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

    out << ")bitSet(" << label.bitSetNg << ")";

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

}

void LabelingAlgorithmNS::updateLabelCost(Label &label, const VetMatResCost &vetMatResCost, FloatType labelStart)
{

    label.vetResources[0] = labelStart;
    for(int i=0; i < (label.tamRoute-1); ++i)
    {
        label.vetResources[0] += vetMatResCost[0](label.vetRoute[i], label.vetRoute[i+1]);
    }

}


LabelingAlgorithmNS::Label* LabelingAlgorithmNS::LabelingData::getBestLabel(int cust)
{
    MatBucket &matBucket = vetMatBucket[cust];
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


void LabelingAlgorithmNS::LabelingData::dominanceInterBuckets(std::multiset<Label*, LabelCmp> &setLabel, int numRes, const int localNumMaxLabel)
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
                    for(int jj=j+1; jj < vetNumSteps[1]; ++jj)
                    {
                        Bucket &b1 = vetMatBucket[cust].mat(ii, jj);
                        if(b1.sizeVetPtrLabel == 0)
                            continue;

                        int t1;

                        for(int t0=0; t0 < b0.sizeVetPtrLabel; ++t0)
                        {
                            t1 = 0;
                            Label *label0 = b0.vetPtrLabel[t0];
                            while(t1 < b1.sizeVetPtrLabel)
                            {

                                Label *label1 = b1.vetPtrLabel[t1];

                                if(checkDominance(*label0, *label1, numRes) && !checkDominance(&label1, &label0, numRes))
                                {
                                    //std::cout<<"Domina\n";
                                    // Rm label1
                                    setLabel.erase(b1.vetPtrLabel[t1]);
                                    //bucket.vetPtrLabel[k]->active = false;
                                    labelPoolG.delT(b1.vetPtrLabel[t1]);

                                    if(t1 == (b1.sizeVetPtrLabel-1))
                                    {
                                        //labelPool.delT(bucket.vetPtrLabel[k]);
                                        b1.vetPtrLabel[t1] = nullptr;
                                    }
                                    else
                                    {
                                        std::swap(b1.vetPtrLabel[t1], b1.vetPtrLabel[b1.sizeVetPtrLabel-1]);
                                        //labelPool.delT(bucket.vetPtrLabel[bucket.sizeVetPtrLabel-1]);
                                        b1.vetPtrLabel[b1.sizeVetPtrLabel-1] = nullptr;
                                    }

                                    b1.sizeVetPtrLabel -= 1;

                                    numDel += 1;
                                    if((int)setLabel.size() < localNumMaxLabel/2)
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
            if(matDist(i, j) < -DW_DecompNS::TolObjSubProb)
                return false;
        }
    }

    return true;

}

bool LabelingAlgorithmNS::containRoute(const Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> &vetLabel, int numSol, Label* label)
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