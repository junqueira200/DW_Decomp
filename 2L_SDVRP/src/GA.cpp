#include "GA.h"
#include "InputOutput.h"

using namespace GA_NS;
using namespace ParseInputNS;
using namespace RandomKeyNS;

bool GA_NS::ga(SolucaoNS::Bin  &bin,
               SolucaoNS::Rota &route,
               VectorI* 		ptrVetItems,
               int 				sizeVetItems)
{
    static Vector<RandomKey> vetRandKeyAux(input.gaPopulationSize);
    static VectorI vetItems(InstanceNS::instanciaG.numItens);
    static Vector<RandomKeyPair> vetRandKey(input.gaPopulationSize);
    static RandomKey randKey, bestRandKey;
    static SolucaoNS::Penalty penalty;
    double bestRandKeyPenalty = std::numeric_limits<double>::max();

    const int numberOfElite   = input.gaElitePercetence*input.gaPopulationSize;
    const int numberOfDiscart = input.gaDiscartPercetence*input.gaPopulationSize;
    const int numberOfRest	  = input.gaPopulationSize - numberOfDiscart - numberOfElite;

    const std::pair<int, int> indexElite(0, numberOfElite-1);
    const std::pair<int, int> indexRest(numberOfElite, numberOfElite+numberOfRest-1);
    const std::pair<int, int> indexDiscart(numberOfElite+numberOfRest,
                                           input.gaPopulationSize-1);

    //std::printf("Elite:\t %d:%d\n", indexElite.first, indexElite.second);
    //std::printf("Rest:\t %d:%d\n", indexRest.first, indexRest.second);
    //std::printf("Discart:\t %d:%d\n", indexDiscart.first, indexDiscart.second);

    static Vector<RandomKey> vetRandKeyAux2(numberOfRest);
    static Vector<RandomKey*> vetRandKeyAux2Ptr(numberOfRest);

    if(numberOfRest <= 0 || numberOfDiscart <= 0 || numberOfElite <= 0)
    {
        std::printf("Error, numberOfElite(%d), numberOfDiscart(%d), numberOfRest(%d)<=0\n",
                    numberOfElite, numberOfDiscart, numberOfRest);
        PRINT_THROW();
    }

    //std::printf("numberOfElite: \t %d\nnumberOfDiscart: \t %d\nnumberOfRest: \t %d\n",
    //            numberOfElite, numberOfDiscart, numberOfRest);

    if(!ptrVetItems)
    {
        sizeVetItems = InstanceNS::copiaItensClientes(route.vetRota,
                                                      route.numPos,
                                                      vetItems);
    }
    else
        copyVet(*ptrVetItems, vetItems, sizeVetItems);

    std::sort(vetItems.begin(), vetItems.begin()+sizeVetItems);

    for(int i=0; i < input.gaPopulationSize; ++i)
        vetRandKey[i].second = &vetRandKeyAux[i];

    for(int i=0; i < numberOfRest; ++i)
        vetRandKeyAux2Ptr[i] = &vetRandKeyAux2[i];


    startPopulation(vetRandKey, sizeVetItems);

    for(int i=0; i < input.gaNumberOfGenerations; ++i)
    {
        // Evaluate the population
        for(auto& it:vetRandKey)
        {
            randKey.vetRandomKey = it.second->vetRandomKey;
            if(decoderRandomKey(randKey, bin, vetItems, sizeVetItems, route, penalty))
                return true;

            it.first = penalty.getValue();
        }

        // Sort the population using the penalty value
        std::sort(vetRandKey.begin(), vetRandKey.end(),
                  [](RandomKeyPair& p0, RandomKeyPair& p1) -> bool
                  {
                        return p0.first < p1.first;
                  });

        // Copy the best
        if(doubleLess(vetRandKey[0].first, bestRandKeyPenalty))
        {
            bestRandKey.vetRandomKey = vetRandKey[0].second->vetRandomKey;
            bestRandKeyPenalty = vetRandKey[0].first;
        }


        // Remove the worse randomKeys
        int quant = 0;

        for(int i=(input.gaPopulationSize-1); quant < numberOfDiscart && i >= 0; --i)
        {
            vetRandKey[i].second->generateRandomKey(sizeVetItems);
            quant += 1;
        }

        // Crossover
        int nextSon = 0;

        for(int j=0; j < numberOfRest; ++j)
        {
            int index  = RandNs::getRandInt(indexElite.first, indexElite.second);
            int index2 = RandNs::getRandInt(indexRest.first, indexRest.second);

            RandomKey* randKeyElite = vetRandKey[index].second;
            RandomKey* randKeyRest  = vetRandKey[index2].second;
            RandomKey* randKeySon   = vetRandKeyAux2Ptr[nextSon];

            for(int t=0; t < sizeVetItems; ++t)
            {
                double randVal = RandNs::getRandDouble();
                if(randVal < input.gaProbabilityToInheritFromElite)
                    randKeySon->vetRandomKey[t] = randKeyElite->vetRandomKey[t];
                else
                    randKeySon->vetRandomKey[t] = randKeyRest->vetRandomKey[t];

            }

            nextSon += 1;
        }

        // Swap the generated sons to the population
        for(int j=indexRest.first; j <= indexRest.second; ++j)
            std::swap(vetRandKey[j].second, vetRandKeyAux2Ptr[j-indexRest.first]);


    }

    return false;

}

void GA_NS::startPopulation(Vector<RandomKeyPair> &vetRandKey, int vetItemsSize)
{
    for(auto& it:vetRandKey)
        it.second->generateRandomKey(vetItemsSize);
}


