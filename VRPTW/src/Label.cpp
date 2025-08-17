#include "Label.h"
#include "DW_Decomp.h"
#include "LabelingAlgorithm.h"

using namespace LabelingAlgorithmNS;

void LabelingAlgorithmNS::LabelingData::checkVetMatBucketBackward()
{
    // TODO remover
    return;
    for(int k=1; k < vetMatBucketBackward.size(); ++k)
    {
        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[1]; ++j)
            {
                Bucket& bucket = vetMatBucketBackward[k].mat(i, j);

                for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
                {
                    for(int tt=(t+1); tt < bucket.sizeVetPtrLabel; ++tt)
                    {
                        FloatType a = bucket.vetPtrLabel[t]->vetResources[0];
                        FloatType b = bucket.vetPtrLabel[tt]->vetResources[0];

                        if(doubleGreater(a, b, (FloatType)1E-3))
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<a<<"\n"<<b<<"\n";
                            std::cout<<"doubleGreater("<<doubleGreater(a, b, (FloatType)1E-3)<<"); >"<<(a>b)<<"\n";

                            std::cout<<bucket.print(2)<<"\n";

                            PRINT_EXIT();
                        }
                    }
                }
            }
        }
    }
}

void LabelingAlgorithmNS::LabelingData::checkVetMatBucketForward()
{
    // TODO remover
    return;
    for(int k=1; k < vetMatBucketBackward.size(); ++k)
    {
        for(int i=0; i < vetNumSteps[0]; ++i)
        {
            for(int j=0; j < vetNumSteps[1]; ++j)
            {
                Bucket& bucket = vetMatBucketBackward[k].mat(i, j);

                for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
                {
                    for(int tt=(t+1); tt < bucket.sizeVetPtrLabel; ++tt)
                    {
                        if(doubleGreater(bucket.vetPtrLabel[t]->vetResources[0],
                                         bucket.vetPtrLabel[tt]->vetResources[0], (FloatType)1E-3))
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<
                                       bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                            PRINT_EXIT();
                        }
                    }
                }
            }
        }
    }
}

static LabelCmp labelCmp;

Label* LabelHeap::extractTop()
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

void LabelHeap::decreaseKey(int i, KEY_TYPE val)
{
    if(vet[i] == nullptr)
    {
        std::cout<<"heapSize("<<heapSize<<")\n";
        std::cout<<"vet["<<i<<"] is equal to null\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    vet[i]->HEAP_KEY = val;
    int iParent = parent(i);

    while(i!=0 && labelCmp.isGreater(vet[iParent], vet[i]))//vet[iParent]->vetResources[0] > vet[i]->vetResources[0])
    {
        std::swap(vet[i], vet[iParent]);
        std::swap(vet[i]->posHeap, vet[iParent]->posHeap);
        i = iParent;
        iParent = parent(i);
    }

}

void LabelHeap::deleteKey(int i)
{
    if(i > heapSize)
    {
        std::cout<<"i("<<i<<"); size("<<heapSize<<")\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    Label* labelI = vet[i];
    decreaseKey(i, KEY_EXTREME);


    Label* label = extractTop();

    if(labelI != label)
    {
        std::cout<<"labelI("<<labelI<<") diffrent from top label ("<<label<<")\n\n";
        std::cout<<"Size: "<<heapSize<<"\n\n";
        PRINT_EXIT();
    }

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

std::string LabelingAlgorithmNS::Bucket::print(int numResorces)
{
    std::string str;

    for(int i=0; i < sizeVetPtrLabel; ++i)
    {
        Label* label = vetPtrLabel[i];
        if(label)
        {

            str += "("+std::to_string(label->cust)+"[";
            for(int t=0; t < numResorces; ++t)
                str += std::format("{:.5f}", label->vetResources[t]) + ", ";
            str += "]; ";
            for(int t=0; t < label->tamRoute; ++t)
                str += std::to_string(label->vetRoute[t]) + " ";
            str += ");  ";
        }
    }

    return str;
}

void LabelingAlgorithmNS::Bucket::removeElement(int i)
{
    for(int t=i; t < (sizeVetPtrLabel-1); ++t)
    {
        vetPtrLabel[t] = vetPtrLabel[t+1];
        vetPtrLabel[t]->posBucket = t;
    }

    sizeVetPtrLabel += -1;
}

void LabelingAlgorithmNS::Bucket::addElement(int pos, Label* labelPtr)
{

    if(sizeVetPtrLabel > 0)
    {

        int size = sizeVetPtrLabel+1;
        if(size > vetPtrLabel.size())
            vetPtrLabel.conservativeResize(2*size);

        for(int i=sizeVetPtrLabel-1; i >= pos; --i)
        {
            int index = i+1;
            vetPtrLabel[index] = vetPtrLabel[i];
            vetPtrLabel[index]->posBucket = index;
        }
    }

    vetPtrLabel[pos] = labelPtr;
    labelPtr->posBucket = pos;

    sizeVetPtrLabel += 1;
}

Index LabelingAlgorithmNS::LabelingData::getListOfIndexForMerge(const Label& label)
{

    Eigen::Array<double, 1, 2> arrayResorces;
    Eigen::Array<int, 1, 2> indexStart, indexEnd;
    arrayResorces.setZero();
    indexStart.setZero();

    indexEnd(0) = vetNumSteps(0) - 1;
    indexEnd(1) = vetNumSteps(1) - 1;

    Vector<std::pair<int, int>> vetPoints;
    vetPoints.reserve(10);


    for(int i=(vetNumSteps(0)-1); i >= 0; --i)
    {
        double lb = vetMatBound[0](i, 0).lowerBound;
        //std::printf("lb(%.2f), res(%.2f)\n", lb, label.vetResources[0]);
        if( lb == -std::numeric_limits<double>::infinity() || (lb+label.vetResources[0]) < -DW_DecompNS::TolObjSubProb)
        {
            indexEnd(0) = i;
            break;
        }
    }


    if(label.typeLabel == Forward)
    {
        for(int i=0; i < vetNumSteps(1); ++i)
        {
            double ub = vetMatBound[1](0, i).upperBound;

            if(label.vetResources[1] < ub)
            {
                indexStart(1) = i;
                break;
            }
        }

    }
    else
    {
        double forwardDemand = vetMaxResources[1] - label.vetResources[1];
        for(int i=(vetNumSteps(1)-1); i >= 0; --i)
        {
            double lb = vetMatBound[1](0, i).lowerBound;

            if((forwardDemand+lb) <= vetMaxResources[1])
            {
                indexEnd(1) = i;
                break;
            }
        }
    }

    return {indexStart, indexEnd};
}


std::string LabelingAlgorithmNS::printIndex(const Index& index)
{
    std::string str;

    str += "["+ std::to_string(index.start(0)) + "; " + std::to_string(index.end(0)) + "]\n";
    str += "["+ std::to_string(index.start(1)) + "; " + std::to_string(index.end(1)) + "]\n";

    return str;
}

int LabelingAlgorithmNS::LabelingData::doMerge(Label* label, const ArrayResources& vetMaxResources,
                                               const MatBoundRes& vetVetBound, int numResorces)
{

    //checkLabels();

    Index* index 		 = nullptr;
    MatBucket* matBucket = nullptr;

    if(label->typeLabel == Forward)
    {
        index = &matForwardRange(label->i, label->j);
        matBucket = &vetMatBucketBackward[label->cust];
    }
    else
    {
        index = &matBackwardRange(label->i, label->j);
        matBucket = &vetMatBucketForward[label->cust];
    }


    // Go through the first index, aka reduced cost
    for(int i=index->start(0); i <= index->end(0); ++i)
    {
        // Go through the second index, aka demand
        for(int j=index->start(1); j <= index->end(1); ++j)
        {
            Bucket& bucket = matBucket->mat(i, j);

            // Go through labels
            for(int t=0; t < bucket.sizeVetPtrLabel; ++t)
            {
                Label* labelAux = bucket.vetPtrLabel[t];

                if(label->typeLabel == Forward)
                {
                    if(labelAux->typeLabel == Forward)
                    {
                        std::printf("Error!, label have type equal to forward, and also labelAux!\n");
                        PRINT_EXIT();
                    }
                }
                else // Backward
                {
                    if(labelAux->typeLabel == Backward)
                    {
                        std::printf("Error!, label have type equal to backward, and also labelAux!\n");
                        PRINT_EXIT();
                    }
                }

                Label* result   = mergeForwardAndBackward(label, labelAux, vetMaxResources, vetVetBound, numResorces);

                if(!result)
                {
                    continue;
                }

                int correctPos = -1;
                Bucket* bucket = dominanceIntraBucketSlow(result->cust, result, *this, nullptr, numResorces, result->cust, correctPos);

                // TODO menory leek! :(
                if(!bucket)
                    continue;

                if((bucket->sizeVetPtrLabel+1) > bucket->vetPtrLabel.size())
                    bucket->vetPtrLabel.conservativeResize(2*bucket->sizeVetPtrLabel);

                for(int p=bucket->sizeVetPtrLabel-1; p >= correctPos; --p)
                {
                    int index = (p+1);
                    bucket->vetPtrLabel[index] = bucket->vetPtrLabel[p];
                    bucket->vetPtrLabel[index]->posBucket = index;
                }

                bucket->vetPtrLabel[correctPos] = result;
                result->posBucket = correctPos;
                bucket->sizeVetPtrLabel += 1;

                //std::cout<<"MERGE:\n"<<*label<<"\n"<<*labelAux<<"\n"<<*result<<"\n\n";
                //PRINT_EXIT();
            }
        }
    }

    return getNumberOfSolutions();

}

void LabelingAlgorithmNS::LabelingData::checkLabels()
{

    for(int t=0; t < 2; ++t)
    {
        Eigen::VectorX<MatBucket>* vetMatBucket = nullptr;
        TypeLabel type;

        if(t == 0)
        {
            vetMatBucket = &vetMatBucketForward;
            type = Forward;
        }
        else
        {
            vetMatBucket = &vetMatBucketBackward;
            type = Backward;
        }

        for(int k=0; k < numCust; ++k)
        {
            MatBucket& matBucket = (*vetMatBucket)(k);
            for(int i=0; i < vetNumSteps[0]; ++i)
            {
                for(int j=0; j < vetNumSteps[1]; ++j)
                {
                    Bucket& bucket = matBucket.mat(i, j);

                    for(int l=0; l < bucket.sizeVetPtrLabel; ++l)
                    {
                        Label* label = bucket.vetPtrLabel[l];
                        if(label->typeLabel != type)
                        {
                            std::cout<<"Error in label of ("<<k<<"), pos: ("<<i<<", "<<j<<")\n";
                            std::cout<<"Label Type("<<label->typeLabel<<") shuld be: "<<type<<"\n\n";
                            std::cout<<*label<<"\n\n";
                            PRINT_EXIT();

                        }

                        if(label->cust != k)
                        {
                            std::cout<<"Error in label of ("<<k<<"), pos: ("<<i<<", "<<j<<")\n";
                            std::cout<<"Cust("<<label->cust<<") shuld be "<<k<<"\n";
                            std::cout<<*label<<"\n\n";
                            PRINT_EXIT();
                        }
                    }
                }
            }
        }
    }
}







