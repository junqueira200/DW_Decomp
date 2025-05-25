#include "Label.h"


void LabelingAlgorithmNS::LabelingData::checkVetMatBucketBackward()
{
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
                                              bucket.vetPtrLabel[tt]->vetResources[0], FloatEp))
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
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
                                         bucket.vetPtrLabel[tt]->vetResources[0], FloatEp))
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                            PRINT_EXIT();
                        }
                    }
                }
            }
        }
    }
}


