//
// Created by igor on 20/11/24.

// Alteracoes: checkDominance was temporary removed
//
#include "LabelingAlgorithm.h"
#include "MemoryPool.h"
#include <list>

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
LabelingAlgorithmNS::forwardLabelingAlgorithm(const int numRes,
                                              const int numCust,
                                              const VetMatResCost &vetMatResCost,
                                              const VetVetResBound &vetVetBound,
                                              const int dest,
                                              const NgSet &ngSet,
                                              LabelingData &lData,
                                              Eigen::VectorXd &vetX,
                                              const double labelStart,
                                              int NumMaxLabePerBucket)
{

    if(NumMaxLabePerBucket == -1)
        NumMaxLabePerBucket = std::numeric_limits<int>::max();

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
        labelPoolG.startPool(44, 400);
    }

    labelPoolG.resetPool(false);
    lData.flushLabel();

    //std::list<Label*> listLabel;
    std::set<Label*, LabelCmp> setLabel;

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
    labelPtr->bitSet = 0;


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

    //while(!listLabel.empty() && !labelPtrBest)
    while(!setLabel.empty() && !labelPtrBest)
    {
        if(setLabel.size() > NumMaxLabel && DominaIterBuckets)
            lData.dominanceInterBuckets(setLabel, numRes);

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

//std::cout<<"Extend label\n";
        // Extend label
        for(int t=1; t < numCust; ++t)
        {
            if(t == lastCust || labelPtr->bitSet[t] == 1)
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

                bool dominate = false;

                int k=0;



                while(k < bucket.sizeVetPtrLabel)
                {
                    if(Print)
                        std::cout<<"\t\t\tcheckDominance "<<bucket.vetPtrLabel[k]<<": "<<*bucket.vetPtrLabel[k]<<"\n";

                    if(checkDominance(*labelPtrAux, *bucket.vetPtrLabel[k], numRes))
                    {
                        if(bucket.vetPtrLabel[k] == labelPtrBest)
                            labelPtrBest = nullptr;

                        setLabel.erase(bucket.vetPtrLabel[k]);
                        //bucket.vetPtrLabel[k]->active = false;
                        labelPoolG.delT(bucket.vetPtrLabel[k]);

                        if(Print)
                            std::cout<<"\t\t\t\t("<<labelPtrAux<<") domina ("<<bucket.vetPtrLabel[k]<<")\n";

                        if(k == (bucket.sizeVetPtrLabel-1))
                        {
                            //labelPool.delT(bucket.vetPtrLabel[k]);
                            bucket.vetPtrLabel[k] = nullptr;
                        }
                        else
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
                            std::cout<<"\t\t\t\t<"<<bucket.vetPtrLabel[k]<<">> domina <<"<<labelPtrAux<<">>\n";

                        labelPoolG.delT(labelPtrAux);
                        labelPtrAux = nullptr;
                        break;
                    }

                    k += 1;
                }




                if(!labelPtrAux)
                    continue;

                if((bucket.sizeVetPtrLabel+1) > NumMaxLabePerBucket && t != dest)
                {
                    labelPoolG.delT(labelPtrAux);
                    continue;
                }

                maxSizeVetPtrLabel = std::max(maxSizeVetPtrLabel, bucket.sizeVetPtrLabel);

                bucket.addLabel(labelPtrAux);
                labelPtrAux->active = true;
                //listLabel.push_back(labelPtrAux);

                if(labelPtrAux->cust == dest && labelPtrAux->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                {

//std::cout<<*labelPtrAux<<"\n";

                    removeCycles(*labelPtrAux, numCust);
                    updateLabelCost(*labelPtrAux, vetMatResCost);

                    //if(Print)
                        //std::cout<<"*"<<*labelPtrAux<<"\n\n";

                    if(labelPtrAux->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                    {
                        std::cout<<"*"<<*labelPtrAux<<"\n\n";
                        if(labelPtrBest)
                        {
                            if(labelPtrAux->vetResources[0] < labelPtrBest->vetResources[0])
                                labelPtrBest = labelPtrAux;
                        } else
                            labelPtrBest = labelPtrAux;

                        if(labelPtrAux->vetResources[0] <= RedCostCut)
                            break;   // for(int t=1)

                    }
                }

                setLabel.insert(labelPtrAux);

/*                else if(labelPtrAux->cust == dest)
                {
                    std::cout<<"*"<<*labelPtrAux<<"\n\n";
                }*/

//std::cout<<"extendLabel to t("<<t<<")\n";

            }
            else
                labelPoolG.delT(labelPtrAux);

        }

        if(labelPtrBest)
        {
            if(labelPtrBest->vetResources[0] <= RedCostCut)
                break;
        }

//std::cout<<"END Extend label\n";

        if(Print)
            std::cout<<"\n########################################################################################\n\n";

        labelPoolG.delT(labelPtr);
        numIt += 1;
    }

    std::cout<<"MAX SIZE: "<<maxSize<<"\n";
    std::cout<<"maxSizeVetPtrLabel: "<<maxSizeVetPtrLabel<<"\n\n";

    if(labelPtrBest)
    {
        //std::cout << "BEST LABEL: " << *labelPtrBest << "\n";

        //removeCycles(*labelPtrBest, numCust);
        //updateLabelCost(*labelPtrBest, vetMatResCost);
        labelPtrBest->vetRoute[labelPtrBest->tamRoute-1] = labelPtrBest->vetRoute[0];

        std::cout << "BEST LABEL: " << *labelPtrBest << "\n";


        if(labelPtrBest->vetResources[0] > -DW_DecompNS::TolObjSubProb)
            return false;


        std::cout << "listLabel.size: " << setLabel.size() << "\n";

        auto &vetRoute = labelPtrBest->vetRoute;

        for(int i = 0; i < (labelPtrBest->tamRoute - 1); ++i)
        {
            vetX(getIndex(vetRoute[i], vetRoute[i+1], numCust-1)) = 1.0;
        }


        //std::cout<<"X: "<<vetX.transpose()<<"\n";


/*
std::cout<<"*****************END FORWARD LABELING ALGORITHM*****************\n";
std::cout<<"****************************************************************\n\n";

*/
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

    newLabel.bitSet = 0;
    newLabel.bitSet = label.bitSet;

    // Checks if custJ its in custI ngSet
    if(custJ != 0 && ngSet.contain(custJ, custI))
        newLabel.bitSet[custJ] = true;

    if((newLabel.vetRoute.size()) < label.tamRoute+1)
    {   std::cout<<"ini resize\n";
        //newLabel.vetRoute.resize(label.vetRoute.size() + 1);
        std::cout<<"end resize\n";
        throw "ERRO, OUT OF MEMORY";
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
        double resource = vetStepSize[r].start;
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
        double r00 = vetStepSize[0].start + i*vetStepSize[0].stepSize;

        double r01 = r00 + vetStepSize[0].stepSize;
        for(int j=0; j < vetNumSteps[1]; ++j)
        {
            double r10 = vetStepSize[1].start + j*vetStepSize[1].stepSize;
            double r11 = r10 + vetStepSize[1].stepSize;

            if(i != 0)
                matBound(i, j)[0].lowerBound = r00;
            else
                matBound(i, j)[0].lowerBound = -std::numeric_limits<double>::infinity();

            if(i != (vetNumSteps[0]-1))
                matBound(i, j)[0].upperBound = r01;
            else
                matBound(i, j)[0].upperBound = std::numeric_limits<double>::infinity();

            matBound(i, j)[1].lowerBound = r10;
            matBound(i, j)[1].upperBound = r11;

        }
    }

    setupGraphBucket();

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

int LabelingAlgorithmNS::LabelingData::getIndex(int resource, double val)
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

    while(!(val >= matBound(vetIndex[0], vetIndex[1])[resource].lowerBound && val < matBound(vetIndex[0], vetIndex[1])[resource].upperBound))
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

        if(start == index && start == end)
        {
            if(!(val >= matBound(vetIndex[0], vetIndex[1])[resource].lowerBound && val < matBound(vetIndex[0], vetIndex[1])[resource].upperBound))
                return -1;

        }
    }

    return vetIndex[resource];

/*    if(resource == 0)
    {
        while(!(val >= matBound(index, 0)[0].lowerBound && val < matBound(index, 0)[0].upperBound))
        {
            //std::cout<<"index: "<<index<<"\n";

            if(val < matBound(index, 0)[0].lowerBound)
            {
                end = index-1;
                index = start + (end-start)/2;
            }
            else
            {
                start = index+1;
                index = start + (end-start)/2;
            }

            if(start == index && start == end)
            {
                if(!(val >= matBound(index, 0)[0].lowerBound && val < matBound(index, 0)[0].upperBound))
                    return -1;

            }
        }

        return index;
    }
    else
    {

        while(!(val >= matBound(0, index)[1].lowerBound && val < matBound(0, index)[1].upperBound))
        {
            //std::cout<<"index: "<<index<<"\n";

            if(val < matBound(0, index)[1].lowerBound)
            {
                end = index-1;
                index = start + (end-start)/2;
            }
            else
            {
                start = index+1;
                index = start + (end-start)/2;
            }

            if(start == index && start == end)
            {
                if(!(val >= matBound(0, index)[1].lowerBound && val < matBound(0, index)[1].upperBound))
                    return -1;

            }
        }

        return index;
    }*/



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
        std::swap(vetMatBucket[cust].mat(i, j).vetPtrLabel[pos], vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1]);

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
        for(int i=sizeVetPtrLabel; i < 2*sizeVetPtrLabel; ++i)
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

    out<<")bitSet("<<label.bitSet<<")";

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
//                prox += 1;
                //label.tamRoute -= 1;

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

void LabelingAlgorithmNS::updateLabelCost(Label &label, const VetMatResCost &vetMatResCost)
{

    label.vetResources[0] = 0.0;
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


void LabelingAlgorithmNS::LabelingData::dominanceInterBuckets(std::set<Label*, LabelCmp> &setLabel, int numRes)
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

                                if(checkDominance(*label0, *label1, numRes))
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
                                    if(setLabel.size() < NumMaxLabel/2)
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