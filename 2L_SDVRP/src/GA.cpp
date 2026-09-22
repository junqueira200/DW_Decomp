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

    if(!ptrVetItems)
    {
        sizeVetItems = InstanceNS::copiaItensClientes(route.vetRota,
                                                      route.numPos,
                                                      vetItems);
    }
    else
        copyVet(*ptrVetItems, vetItems, sizeVetItems);

    for(int i=0; i < input.gaPopulationSize; ++i)
        vetRandKey[i].second = &vetRandKeyAux[i];


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
        bestRandKey.vetRandomKey = vetRandKey[0].second->vetRandomKey;


        // Remove the worse randomKeys

        int quant = 0;
        int quantMax = input.gaDiscartPercetence*input.gaPopulationSize;

        for(int i=(input.gaPopulationSize-1); i >= 0; --i)
        {
            vetRandKey[i].second->generateRandomKey(sizeVetItems);
            quant += 1;

            if(quant >= quantMax)
                break;
        }

        // Crossover


        // Mutation
    }

    return false;

}

void GA_NS::startPopulation(Vector<RandomKeyPair> &vetRandKey, int vetItemsSize)
{
    for(auto& it:vetRandKey)
        it.second->generateRandomKey(vetItemsSize);
}


