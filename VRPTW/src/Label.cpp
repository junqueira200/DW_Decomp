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
                    bool equal = true;
                    for(int tt=(t+1); tt < bucket.sizeVetPtrLabel; ++tt)
                    {
                        if(doubleEqual(bucket.vetPtrLabel[t]->vetResources[0], bucket.vetPtrLabel[tt]->vetResources[0],
                                       FloatEp) && !equal)
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<"Reduced Costs are equal, but equal(false)\n";
                            PRINT_EXIT();
                        }

                        if(doubleEqual(bucket.vetPtrLabel[t]->vetResources[0], bucket.vetPtrLabel[tt]->vetResources[0],
                                       FloatEp))
                        {
                            if(!(doubleGreaterEqual(bucket.vetPtrLabel[t]->vetResources[1],
                                                    bucket.vetPtrLabel[tt]->vetResources[1], FloatEp)))
                            {
                                std::cout<<"ERROR:\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                                std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                                std::cout<<(bucket.vetPtrLabel[t]->vetResources[0] >= bucket.vetPtrLabel[tt]->vetResources[0])<<"\n";
                                PRINT_EXIT();
                            }
                        }
                        else if(doubleGreater(bucket.vetPtrLabel[t]->vetResources[0],
                                              bucket.vetPtrLabel[tt]->vetResources[0], FloatEp))
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                            PRINT_EXIT();
                        }
                        else
                            equal = false;
                    }
                }
            }
        }
    }
}

void LabelingAlgorithmNS::LabelingData::checkVetMatBucketForward()
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
                    bool equal = true;
                    for(int tt=(t+1); tt < bucket.sizeVetPtrLabel; ++tt)
                    {
                        if(bucket.vetPtrLabel[t]->vetResources[0] == bucket.vetPtrLabel[tt]->vetResources[0] && !equal)
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<"Reduced Costs are equal, but equal(false)\n";
                            PRINT_EXIT();

                        }

                        if(bucket.vetPtrLabel[t]->vetResources[0] == bucket.vetPtrLabel[tt]->vetResources[0])
                        {
                            if(!(bucket.vetPtrLabel[t]->vetResources[1] <= bucket.vetPtrLabel[tt]->vetResources[1]))
                            {
                                std::cout<<"ERROR:\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                                std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                                std::cout<<(bucket.vetPtrLabel[t]->vetResources[0] >= bucket.vetPtrLabel[tt]->vetResources[0])<<"\n";
                                PRINT_EXIT();
                            }
                        }
                        else if(bucket.vetPtrLabel[t]->vetResources[0] > bucket.vetPtrLabel[tt]->vetResources[0])
                        {
                            std::cout<<"ERROR\n"<<*bucket.vetPtrLabel[t]<<"\n"<<*bucket.vetPtrLabel[tt]<<"\n";
                            std::cout<<bucket.vetPtrLabel[t]->vetResources[0]<<"\n"<<bucket.vetPtrLabel[tt]->vetResources[0]<<"\n";
                            PRINT_EXIT();
                        }
                        else
                            equal = false;
                    }
                }
            }
        }
    }


}


