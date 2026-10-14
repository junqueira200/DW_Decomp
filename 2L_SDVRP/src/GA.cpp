/* ****************************************
 * ****************************************
 *  Data:    22/09/26
 *  Arquivo: GA.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#include "GA.h"
#include "InputOutput.h"
#include <set>

using namespace GA_NS;
using namespace ParseInputNS;
using namespace RandomKeyNS;

constexpr bool PrintGA = false;

bool GA_NS::ga(SolucaoNS::Bin      &bin,
               SolucaoNS::Rota     &route,
               const VectorI *const ptrVetItems,
               int                  sizeVetItems)
{
    if constexpr(PrintGA)
        std::printf("GA INI\n");

    static Vector<RandomKey>     vetRandKeyAux(input.gaPopulationSize);
    static VectorI               vetItems(InstanceNS::instanciaG.numItens);
    static Vector<RandomKeyPair> vetRandKey(input.gaPopulationSize);
    static RandomKey             randKey, bestRandKey;
    static SolucaoNS::Penalty    penalty;
    double                       bestRandKeyPenalty = std::numeric_limits<double>::max();

    const int subPopulationSize = input.gaPopulationSize / input.gaNumberOfIslands;
    std::printf("PopulationSize: %d\n\n", input.gaPopulationSize);
    std::printf("subPopulationSize: %d\n\n", subPopulationSize);
    if(input.gaPopulationSize % input.gaNumberOfIslands != 0)
    {
        std::printf("Error, gaPopulationSize(%d) dont divide gaNumberOfIslands(%d)!\n",
                    input.gaPopulationSize,
                    input.gaNumberOfIslands);
        PRINT_THROW();
    }

    const int numberOfElite = subPopulationSize * input.gaElitePercetence;
    const int numberOfDiscart = subPopulationSize * input.gaDiscartPercetence;
    const int numberOfRest = subPopulationSize - numberOfDiscart - numberOfElite;

    const int numberOfMigrationIndividuals = std::ceil(subPopulationSize * input.gaMigrationRate);

    if constexpr(PrintGA)
        std::printf("numberOfMigrationIndividuals: %d\n", numberOfMigrationIndividuals);

    if(numberOfMigrationIndividuals > numberOfDiscart)
    {
        std::printf("Error: numberOfMigrationIndividuals (%d) > numberOfDiscart (%d)!\n",
                    numberOfMigrationIndividuals,
                    numberOfDiscart);
        PRINT_THROW();
    }

    if(numberOfMigrationIndividuals > numberOfElite)
    {
        std::printf("Error: numberOfMigrationIndividuals (%d) > numberOfElite (%d)!\n",
                    numberOfMigrationIndividuals,
                    numberOfElite);
        PRINT_THROW();
    }

    // const std::pair<int, int> indexElite(0, numberOfElite-1);
    // const std::pair<int, int> indexRest(numberOfElite, numberOfElite+numberOfRest-1);
    // const std::pair<int, int> indexDiscart(numberOfElite+numberOfRest,
    //                                        subPopulationSize-1);

    Vector<std::pair<int, int>> vetIndexElite;
    Vector<std::pair<int, int>> vetIndexRest;
    Vector<std::pair<int, int>> vetIndexDiscart;

    for(int i = 0; i < input.gaNumberOfIslands; ++i)
    {
        if constexpr(PrintGA)
            std::printf("Island: %d\n\n", i);

        vetIndexElite.push_back(
            std::make_pair(i * subPopulationSize + 0, i * subPopulationSize + numberOfElite - 1));

        vetIndexRest.push_back(
            std::make_pair(i * subPopulationSize + numberOfElite,
                           i * subPopulationSize + numberOfElite + numberOfRest - 1));

        vetIndexDiscart.push_back(
            std::make_pair(i * subPopulationSize + numberOfElite + numberOfRest,
                           i * subPopulationSize + subPopulationSize - 1));

        if constexpr(PrintGA)
        {
            std::printf("\tElite:\t\t %d:%d\n", vetIndexElite[i].first, vetIndexElite[i].second);

            std::printf("\tRest:\t\t %d:%d\n", vetIndexRest[i].first, vetIndexRest[i].second);

            std::printf(
                "\tDiscart:\t %d:%d\n", vetIndexDiscart[i].first, vetIndexDiscart[i].second);
        }
    }

    static Vector<RandomKey>   vetRandKeyAux2(numberOfRest);
    static Vector<RandomKey *> vetRandKeyAux2Ptr(numberOfRest);

    if(numberOfRest <= 0 || numberOfDiscart <= 0 || numberOfElite <= 0)
    {
        std::printf("Error, numberOfElite(%d), numberOfDiscart(%d), numberOfRest(%d)<=0\n",
                    numberOfElite,
                    numberOfDiscart,
                    numberOfRest);
        PRINT_THROW();
    }

    // std::printf("numberOfElite: \t %d\nnumberOfDiscart: \t %d\nnumberOfRest: \t %d\n",
    //             numberOfElite, numberOfDiscart, numberOfRest);

    if(!ptrVetItems)
    {
        sizeVetItems = InstanceNS::copiaItensClientes(route.vetRota, route.numPos, vetItems);
    }
    else
        copyVet(*ptrVetItems, vetItems, sizeVetItems);

    std::sort(vetItems.begin(), vetItems.begin() + sizeVetItems);

    for(int i = 0; i < input.gaPopulationSize; ++i)
        vetRandKey[i].second = &vetRandKeyAux[i];

    for(int i = 0; i < numberOfRest; ++i)
        vetRandKeyAux2Ptr[i] = &vetRandKeyAux2[i];

    startPopulation(vetRandKey, sizeVetItems);

    // PRINT_THROW();

    for(int i = 0; i < input.gaNumberOfGenerations; ++i)
    {
        if constexpr(PrintGA)
            std::printf("Iteration: %d\n", i);

        // Evaluate the population
        if constexpr(PrintGA)
            std::printf("\tINI decoderRandomKey\n");
        for(auto &it : vetRandKey)
        {
            randKey.vetRandomKey = it.second->vetRandomKey;
            if(decoderRandomKey(randKey, bin, vetItems, sizeVetItems, route, penalty))
            {
                if constexpr(PrintGA)
                    std::printf("GA END\n");
                return true;
            }

            it.first = penalty.getValue();
        }

        if constexpr(PrintGA)
            std::printf("\tEND decoderRandomKey\n");

        for(int island = 0; island < input.gaNumberOfIslands; ++island)
        {
            // Sort the population using the penalty value

            std::sort(vetRandKey.begin() + vetIndexElite[island].first,
                      vetRandKey.begin() + vetIndexElite[island].first + subPopulationSize,
                      [](RandomKeyPair &p0, RandomKeyPair &p1) -> bool
                      { return p0.first < p1.first; });

            int indexBest = vetIndexElite[island].first;
            // Copy the best
            if(doubleLess(vetRandKey[indexBest].first, bestRandKeyPenalty))
            {
                bestRandKey.vetRandomKey = vetRandKey[indexBest].second->vetRandomKey;
                bestRandKeyPenalty = vetRandKey[indexBest].first;
            }

            // Remove the worse randomKeys
            int quant = 0;
            if constexpr(PrintGA)
                std::printf("\tINI Remove\n");

            for(int i = vetIndexDiscart[island].second;
                quant < numberOfDiscart && i >= vetIndexDiscart[island].first;
                --i)
            {
                vetRandKey[i].second->generateRandomKey(sizeVetItems);
                quant += 1;
            }

            // Crossover
            int nextSon = 0;

            if constexpr(PrintGA)
                std::printf("\tINI Crossover\n");
            for(int j = 0; j < numberOfRest; ++j)
            {
                int index =
                    RandNs::getRandInt(vetIndexElite[island].first, vetIndexElite[island].second);

                int index2 =
                    RandNs::getRandInt(vetIndexRest[island].first, vetIndexRest[island].second);

                RandomKey *randKeyElite = vetRandKey[index].second;
                RandomKey *randKeyRest = vetRandKey[index2].second;
                RandomKey *randKeySon = vetRandKeyAux2Ptr[nextSon];

                for(int t = 0; t < sizeVetItems; ++t)
                {
                    double randVal = RandNs::getRandDouble();
                    if(randVal < input.gaProbToInheritFromElite)
                        randKeySon->vetRandomKey[t] = randKeyElite->vetRandomKey[t];
                    else
                        randKeySon->vetRandomKey[t] = randKeyRest->vetRandomKey[t];
                }

                nextSon += 1;
            }

            nextSon = 0;
            // Swap the generated sons to the population
            for(int j = vetIndexRest[island].first; j <= vetIndexRest[island].second; ++j)
            {
                std::swap(vetRandKey[j].second, vetRandKeyAux2Ptr[nextSon]);
                nextSon += 1;
            }
        }

        // Migration
        if(input.gaNumberOfIslands > 1 && i > 0 && (i % input.gaMigrationInterval) == 0)
        {
            // std::printf("Migration\n");
            for(int islandSend = 0; islandSend < input.gaNumberOfIslands; ++islandSend)
            {

                int islandRecive = (islandSend + 1) % input.gaNumberOfIslands;

                for(int j = 0; j < numberOfMigrationIndividuals; ++j)
                {
                    int selecedInd = vetIndexElite[islandSend].first + j;
                    int reciveInd = vetIndexDiscart[islandRecive].second - j;

                    for(int t = 0; t < sizeVetItems; ++t)
                    {
                        vetRandKey[reciveInd].second->vetRandomKey[t] =
                            vetRandKey[selecedInd].second->vetRandomKey[t];
                    }
                }
            }
        }
    }

    std::printf("GA END\n");
    return false;
}

void GA_NS::startPopulation(Vector<RandomKeyPair> &vetRandKey, int vetItemsSize)
{
    for(auto &it : vetRandKey)
        it.second->generateRandomKey(vetItemsSize);
}
