#include "Label.h"

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
    if(vet[i] == nullptr)
    {
        std::cout<<"heapSize("<<heapSize<<")\n";
        std::cout<<"vet["<<i<<"] is equal to null\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

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
    if(i > heapSize)
    {
        std::cout<<"i("<<i<<"); size("<<heapSize<<")\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    decreaseKey(i, -std::numeric_limits<FloatType>::max());
    Label* label = extractMin();

    memset((void*)label, 0, sizeof(Label));

    label->i = -1;
    label->j = -1;
    label->posBucket = -1;
    label->posHeap   = -1;
    label->cust      = -1;

/*
 *  TODO deleteKey
    if(labelHaveRoute(vetRoteG, label))
    {
        std::cout<<"AQUI!";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }
*/

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
