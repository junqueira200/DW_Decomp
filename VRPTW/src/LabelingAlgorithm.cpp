//
// Created by igor on 20/11/24.
//
#include "LabelingAlgorithm.h"
#include "MemoryPool.h"

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

void LabelingAlgorithmNS::NgSet::setNgSets(const EigenMatrixRow &matDist)
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
    vetNumSteps      = {0, 0};

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
        vetMatBucket[i].mat = Eigen::Matrix<Bucket, -1, -1, Eigen::RowMajor>(numMaxSteps, numMaxSteps);


}

void LabelingAlgorithmNS::LabelingData::flushLabel()
{

    for(MatBucket& mat:vetMatBucket)
    {

        for(int i=0; i < numMaxSteps; ++i)
        {
            for(int j=0; j < numMaxSteps; ++j)
            {
//                mat.mat(i, j).vetPtrLabel.conservativeResize()
                mat.mat(i, j).flush();
            }

            if(numMainResources == 1)
                break;
        }

    }

}
