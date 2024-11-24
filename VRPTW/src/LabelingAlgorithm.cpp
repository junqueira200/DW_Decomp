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

std::cout<<"numCust: "<<numCust<<"\n";

    static MemoryPool_NS::Pool<Label> labelPool(4, 8);
    labelPool.resetPool(false);
    lData.flushLabel();

    std::list<Label*> listLabel;

    Label *labelPtr = labelPool.getT();

    labelPtr->vetRoute.resize(10);
    labelPtr->tamRoute = 1;

    labelPtr->vetRoute[0] = 0;
//    labelPtr->vetResources[0] = 0.0;
    labelPtr->vetResources.setZero();

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

    int numIt = 0;

    Label *labelPtrBest = nullptr;

    while(!listLabel.empty())// && !labelPtrBest)
    {
//std::cout<<"before front\n";
        labelPtr = listLabel.back();
std::cout<<"labelPtr: "<<labelPtr<<"\n";
std::cout<<*labelPtr<<"\n\n";

        listLabel.pop_back();


        //std::cout<<"Extract label("<<labelPtr->vetRoute<<")\n";
        int lastCust = labelPtr->cust;

        if(lastCust == dest)
        {
            std::cout<<"lastCust==dest("<<dest<<")\n";
            continue;
        }


        // Remove labelPtr from vetMatBucket
        lData.removeLabel(labelPtr);

        // Extend label
        for(int t=1; t < numCust; ++t)
        {
            if(t == lastCust || labelPtr->bitSet[t] == 1)
                continue;


std::cout<<"\tt("<<t<<")\n";

            Label *labelPtrAux = labelPool.getT();

            if(extendLabel(*labelPtr, *labelPtrAux, vetMatResCost, vetVetBound, labelPtr->cust, t, ngSet, numRes))
            {
std::cout<<"\t\textendLabel;"<<*labelPtrAux<<"\n";
                // Find index from resources
                i = lData.getIndex(0, labelPtrAux->vetResources[0]);
                j = 0;
                if(lData.numMainResources > 1)
                    j = lData.getIndex(1, labelPtrAux->vetResources[1]);

                labelPtrAux->i = i;
                labelPtrAux->j = j;
                labelPtrAux->cust = t;

std::cout<<"\t\ti("<<i<<"); j("<<j<<")\n";

                Bucket &bucket = lData.vetMatBucket[t].mat(i, j);

                for(int k=0; k < bucket.sizeVetPtrLabel; ++k)
                {
                    std::cout<<"\t\t\tcheckDominance: "<<*bucket.vetPtrLabel[k]<<"\n";
                    if(checkDominance(*labelPtrAux, *bucket.vetPtrLabel[k], numRes))
                    {
std::cout<<"\t\t\t\t("<<labelPtrAux<<") domina ("<<bucket.vetPtrLabel[k]<<")\n";
//std::cout<<"\t\t\t\t"<<*bucket.vetPtrLabel[k]<<"\n\n";
                    }
                }


//std::cout<<"\t\tDepois\n";

                bucket.addLabel(labelPtrAux);
                listLabel.push_back(labelPtrAux);

                if(labelPtrAux->cust == dest && labelPtrAux->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                    labelPtrBest = labelPtrAux;

//std::cout<<"extendLabel to t("<<t<<")\n";

            }
            else
                labelPool.delT(labelPtrAux);



        }

std::cout<<"\n########################################################################################\n\n";

        labelPool.delT(labelPtr);

        if(numIt == 15)
            break;
            //throw "ERRO";

        numIt += 1;
    }

    if(labelPtrBest)
    {
        std::cout << "BEST LABEL: " << *labelPtrBest << "\n";

        removeCycles(*labelPtrBest, numCust);
        updateLabelCost(*labelPtrBest, vetMatResCost);
        labelPtrBest->vetRoute[labelPtrBest->tamRoute-1] = labelPtrBest->vetRoute[0];

        std::cout << "BEST LABEL: " << *labelPtrBest << "\n";
    }

    return;

    MatBucket &matBucket = lData.vetMatBucket[dest];
    Label *labelDest = nullptr;

    for(int i=0; i < lData.vetNumSteps[0]; ++i)
    {
        for(int j=0; j < lData.vetNumSteps[1]; ++j)
        {
//std::cout<<"before\n";
            Bucket &bucket = matBucket.mat(i, j);
//std::cout<<"got bucket\n";
//std::cout<<"bucket.sizeVetPtrLabel("<<bucket.sizeVetPtrLabel<<")\n";

            for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
            {
                // TODO: FIX
                Label* label = bucket.vetPtrLabel[t];

                if(label->vetResources[0] < -DW_DecompNS::TolObjSubProb)
                {
                    if(labelDest)
                    {
                        if(label->vetResources[0] >= labelDest->vetResources[0])
                            continue;
                    }


                    labelDest = label;
                    //break;
                }
            }

            //if(labelDest)
                //break;
        }

        //if(labelDest)
          //  break;
    }

    if(labelDest != nullptr)
        std::cout<<"Dest: "<<(*labelDest)<<"\n\n";

    std::cout<<"listLabel.empty: "<<listLabel.empty()<<"\n";

    removeCycles(*labelDest, numCust);
    updateLabelCost(*labelDest, vetMatResCost);

    labelDest->vetRoute[labelDest->tamRoute-1] = labelDest->vetRoute[0];

    std::cout<<*labelDest<<"\n";

    throw "ERRO";

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
    return ((l1.bitSet&l0.bitSet)==l1.bitSet);


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
std::cout<<"extendLabel\n";

    // Goes through resources
    for(int i=0; i < NumMaxResources; ++i)
    {
/*
std::cout<<"\t\t\tResourcee("<<i<<"); "<<custI<<", "<<custJ<<"\n";
std::cout<<"\t\t\t\tvetMatResCost("<<i<<")("<<custI<<","<<custJ<<") = ";
std::cout<<vetMatResCost[i](custI, custJ)<<"\n";
*/
        // Extend the iº resource
        newLabel.vetResources[i] = label.vetResources[i] + vetMatResCost[i](custI, custJ);

//std::cout<<"\t\t\t\t"<<newLabel.vetResources[i]<<"\n";

        // Check the bound
        if(newLabel.vetResources[i] > vetVetBound[i][custJ].upperBound)
        {
std::cout<<"\tUpperBound\n";
            return false;
        }

        else if(newLabel.vetResources[i] < vetVetBound[i][custJ].lowerBound)
        {
std::cout<<"\tLowerBound\n";
            return false;
            newLabel.vetResources[i] = vetVetBound[i][custJ].lowerBound;
        }

        if((i+1) == numResources)
            break;
    }

    newLabel.bitSet = label.bitSet;

    // Checks if custJ its in custI ngSet
    if(custJ != 0 && ngSet.contain(custJ, custI))
        newLabel.bitSet[custJ] = true;

    if((newLabel.vetRoute.size()) < label.tamRoute+1)
        newLabel.vetRoute.resize(label.vetRoute.size()+1);

std::cout<<"copy route\n";
std::cout<<"label.tamRoute("<<label.tamRoute<<")\n";
std::cout<<"newLabel.vetRoute.size("<<newLabel.vetRoute.size()<<")\n";

    // Copy route
    for(int i=0; i < label.tamRoute; ++i)
    {
        newLabel.vetRoute[i] = label.vetRoute[i];
    }

    newLabel.vetRoute[label.tamRoute] = custJ;
    newLabel.tamRoute = label.tamRoute+1;

std::cout<<"END\n";

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

//    matBound(vetNumSteps[0]-1, 0)[0].upperBound = std::numeric_limits<double>::infinity();
//    matBound(0, 0)[0].lowerBound = -std::numeric_limits<double>::infinity();

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
    std::swap(vetMatBucket[cust].mat(i, j).vetPtrLabel[pos], vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1]);
    vetMatBucket[cust].mat(i, j).vetPtrLabel[size-1] = nullptr;

    vetMatBucket[cust].mat(i, j).sizeVetPtrLabel -= 1;


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

    while(i < label.tamRoute)
    {
        if(vetCust[label.vetRoute[i]] != -1)
        {
            int start = vetCust[label.vetRoute[i]]+1;
            int end   = i+1;
            int prox  = end;

std::cout<<"start("<<start<<"); end("<<end<<")\n";

            for(int k=start; k < end; ++k)
            {
                label.vetRoute[k] = label.vetRoute[prox];
                prox += 1;
                //label.tamRoute -= 1;

            }

            for(int k=end; k < label.tamRoute; ++k)
            {
                label.vetRoute[k] = label.vetRoute[prox];
                prox += 1;
            }

            label.tamRoute -= end-start;
            i = start;
            std::cout<<label<<"\n";

        }
        else
            vetCust[label.vetRoute[i]] = i;

        i += 1;
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
