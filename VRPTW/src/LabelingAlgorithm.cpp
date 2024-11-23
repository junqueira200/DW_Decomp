//
// Created by igor on 20/11/24.
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
        std::cout<<"i("<<i<<")\n";
        int next = 0;
        for(int j=0; j < numCust; ++j)
        {
            if(i == j)
                continue;

            vet[next].dist = matDist(i, j);
            vet[next].i    = j;
            next += 1;

            std::cout<<"("<<j<<","<<matDist(i,j)<<"); ";
        }


        std::sort(vet.begin(), vet.end());
        next = 0;

        std::cout<<"\tsort:\n";

        for(int j=0; j < ngSetSize; ++j)
        {
            matNgSet(i, j) = vet[j].i;
            std::cout<<"("<<vet[j].i<<","<<vet[j].dist<<"); ";
        }

        if((i+1) == numCust)
            break;
    }

}



void LabelingAlgorithmNS::forwardLabelingAlgorithm(const int numRes,
                                                   const int numCust,
                                                   const VetMatResCost& vetMatResCost,
                                                   const VetVetResBound& vetVetBound,
                                                   const int dest,
                                                   const NgSet &ngSet,
                                                   LabelingData &lData)
{

    static MemoryPool_NS::Pool<Label> labelPool(4, 8);
    labelPool.resetPool(false);
    lData.flushLabel();

    std::list<Label*> listLabel;

    Label *labelPtr = labelPool.getT();

    labelPtr->vetRoute.resize(10);
    labelPtr->tamRoute = 1;
    labelPtr->vetRoute[0] = 0;
    labelPtr->vetResources[0] = 0.0;

std::cout<<"\nsetup label0\n";

    int i = lData.getIndex(0, 0.0);
    int j = 0; // lData.getIndex(0, 0.0);

    labelPtr->i = i;
    labelPtr->j = j;
    labelPtr->cust = 0;

    lData.vetMatBucket[0].mat(i, j).vetPtrLabel.resize(10);
    lData.vetMatBucket[0].mat(i, j).vetPtrLabel[0] = labelPtr;
    lData.vetMatBucket[0].mat(i, j).sizeVetPtrLabel = 1;

    listLabel.push_back(labelPtr);
    labelPtr = nullptr;


    while(!listLabel.empty())
    {
        labelPtr = listLabel.front();
        listLabel.pop_front();

        std::cout<<labelPtr->vetRoute[0]<<"\n";

        // Remove labelPtr from vetMatBucket
        lData.removeLabel(labelPtr);

        // Extend label
        for(int t=1; i < numCust; ++t)
        {
std::cout<<"\tfor t("<<t<<")\n";

            Label *labelPtrAux = labelPool.getT();
            if(extendLabel(*labelPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelPtr->cust, t, ngSet, numRes))
            {
std::cout<<"\t\textendLabel\n";
                // Find index from resources
                i = lData.getIndex(0, labelPtrAux->vetResources[0]);
std::cout<<"\t\ti("<<i<<")\n";
                j = 0;
                if(lData.numMainResources > 1)
                    j = lData.getIndex(1, labelPtrAux->vetResources[1]);

std::cout<<"\t\ti("<<i<<"); j("<<j<<")\n";

                Bucket &bucket = lData.vetMatBucket[t].mat(i, j);
                bucket.addLabel(labelPtrAux);
                listLabel.push_back(labelPtrAux);

std::cout<<"extendLabel to t("<<t<<")\n";

            }
            else
                labelPool.delT(labelPtrAux);

        }


        throw "ERRO";
    }


}

bool LabelingAlgorithmNS::checkDominance(const Label& l0, const Label& l1, const int numResources)
{

    // Check the resources
    for(int i=0; i < NumMaxResources; ++i)
    {
        if(l0.vetResources[i] > l1.vetResources[i])
            return false;

        if((i+1) == numResources)
            break;
    }

    // Check if l0 is a subset of l1
    return ((l0.bitSet&l1.bitSet)==l0.bitSet);


}

// TODO Fix
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
            return false;

        else if(newLabel.vetResources[i] < vetVetBound[i][custJ].lowerBound)
        {
            return false;
            newLabel.vetResources[i] = vetVetBound[i][custJ].lowerBound;
        }

        if((i+1) == numResources)
            break;
    }

    // Checks if custJ its in custI ngSet
    if(custJ != 0 && ngSet.contain(custJ, custI))
        newLabel.bitSet[custJ] = true;


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
std::cout<<"Resource("<<r<<") have ("<<step<<") steps\n";

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

            matBound(i, j)[0].lowerBound = r00;
            matBound(i, j)[0].upperBound = r01;

            matBound(i, j)[1].lowerBound = r10;
            matBound(i, j)[1].upperBound = r11;

        }
    }

    //for(int r=0; r < 2; ++r)
    {
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
    }

}

int LabelingAlgorithmNS::LabelingData::getIndex(int resource, double val)
{
    if(resource != 0 && resource != 1)
        PRINT_DEBUG("", "Wrong resource("<<resource<<")");

    if(numMainResources == 1 && resource == 1)
        return 0;

    int start = 0;
    int index = vetNumSteps[resource]/2;
    int end   = vetNumSteps[resource]-1;

    if(resource == 0)
    {
        while(!(val >= matBound(index, 0)[0].lowerBound && val < matBound(index, 0)[0].upperBound))
        {
            std::cout<<"index: "<<index<<"\n";

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
            std::cout<<"index: "<<index<<"\n";

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
    }



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
        PRINT_DEBUG("", "");
    int size = vetMatBucket[cust].mat(i, j).sizeVetPtrLabel;
    std::swap(vetMatBucket[cust].mat(i, j).vetPtrLabel[pos], vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1]);
    vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1] = nullptr;

    vetMatBucket[cust].mat(i, j).sizeVetPtrLabel -= 1;
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
    sizeVetPtrLabel += 1;
}

bool LabelingAlgorithmNS::Bucket::delLabel(LabelingAlgorithmNS::Label *labelPtr)
{
    return false;
}
