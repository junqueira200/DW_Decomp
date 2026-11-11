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
#include "ConstrutivoBin2.h"
#include "AxleWeights.h"
#include "StopData.h"
#include <set>
#include <print>

using namespace GA_NS;
using namespace ParseInputNS;
using namespace RandomKeyNS;
using namespace ConstrutivoBin2NS;
using namespace AxleWeightsNS;
using namespace SolucaoNS;
using namespace StopDataNS;

constexpr bool PrintGA = false;
constexpr bool PrintGAFail = false;   // print best penalty breakdown when GA fails

using namespace SolucaoNS;

bool GA_NS::ga(SolucaoNS::Bin      		&bin,
               SolucaoNS::Rota     		&route,
               const VectorI *const 	ptrVetItems,
               int                  	sizeVetItems,
               Vector<SolucaoNS::Bin>	*ptrVetBin)
{
    int size = 0;

    if(ptrVetBin)
    {
        for(Bin& bin:*ptrVetBin)
            bin.reset();

        size = ptrVetBin->size();
        std::print("Update ");
    }
    else
        bin.reset();

    int nextBin = 0;

    if constexpr(PrintGA)
        std::printf("GA INI\n");

    if(bin.vazio())
        bin.addEp(PontoZero);

    //static Vector<double*>		 vetPtrRandomKey(size);
    static Bin					 binAux;
    static Vector<RandomKey>     vetRandKeyAux(input.gaPopulationSize);
    static VectorI               vetItems(InstanceNS::instanciaG.numItens);
    static Vector<RandomKeyPair> vetRandKey(input.gaPopulationSize);
    static RandomKey             randKey, bestRandKey;
    static SolucaoNS::Penalty    penalty;
    static StopData				 stopData;
    double                       bestRandKeyPenalty = std::numeric_limits<double>::max();
    static SolucaoNS::Penalty    bestPenalty;
    bestPenalty.set0();

    const int subPopulationSize = input.gaPopulationSize / input.gaNumberOfIslands;

    if constexpr(PrintGA)
    {
        std::printf("PopulationSize: %d\n\n", input.gaPopulationSize);
        std::printf("subPopulationSize: %d\n\n", subPopulationSize);
    }
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

    int numberOfMigrationIndividuals = std::ceil(subPopulationSize * input.gaMigrationRate);

    if constexpr(PrintGA)
        std::printf("numberOfMigrationIndividuals: %d\n", numberOfMigrationIndividuals);

    if(numberOfMigrationIndividuals > numberOfDiscart && input.gaNumberOfIslands >= 2)
    {
        std::printf("Error: numberOfMigrationIndividuals (%d) > numberOfDiscart (%d)!\n",
                    numberOfMigrationIndividuals,
                    numberOfDiscart);
        numberOfMigrationIndividuals = std::max(numberOfDiscart-1, 1);
        //PRINT_THROW();
    }

    if(numberOfMigrationIndividuals > numberOfElite && input.gaNumberOfIslands >= 2)
    {
        std::printf("Error: numberOfMigrationIndividuals (%d) > numberOfElite (%d)!\n",
                    numberOfMigrationIndividuals,
                    numberOfElite);
        numberOfMigrationIndividuals = std::max(numberOfElite-1, 1);
        //PRINT_THROW();
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

    const double totalVolume = getTotalVolume(bin, vetItems, sizeVetItems);
    const double areaLeft    = getTotalLeftArea(bin, vetItems, sizeVetItems);
    static const double maxF = 2*semiTrailer.computeMaxF();
    const double areaSupp    = getTotalAreaSupp(bin, vetItems, sizeVetItems);
    const double totalWeight = getTotalWeight(bin, vetItems, sizeVetItems);

    stopData.compute(route, vetItems, sizeVetItems);

    /*
    std::printf(
        "Total Volume: %.6f\n"
//        "Area Left:    %.6f\n"
        "Max FK:       %.6f\n"
        "Area Supp:    %.6f\n"
        "Total Weight: %.6f\n\n",
        totalVolume,
//        areaLeft,
        maxF,
        areaSupp,
        totalWeight
        );
    */

    //PRINT_THROW();
    penalty.setNorm(Penalty::TypePenalty::AxleWeights, maxF);
    penalty.setNorm(Penalty::TypePenalty::Compactness, areaLeft);   // compactness penalty is an area (cm2)
    penalty.setNorm(Penalty::TypePenalty::Lifo, totalVolume);
    penalty.setNorm(Penalty::TypePenalty::LoadBalancing, totalWeight);
    penalty.setNorm(Penalty::TypePenalty::Support, totalVolume);
    penalty.setNorm(Penalty::TypePenalty::VolItemsNotPacked, totalVolume);
    //std::cout<<"vetNorm: "<<penalty.vetNorm<<"\n";

    //PRINT_THROW();

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
            if(it.second->evaluated)
            {
                it.first = it.second->fitness;
                continue;
            }

            const int chromSize = sizeVetItems + InstanceNS::instanciaG.numItens;
            for(int i = 0; i < chromSize; ++i)
                randKey.vetRandomKey[i] = it.second->vetRandomKey[i];

            copiaBin(bin, binAux);

            if(decoderRandomKey(randKey, binAux, vetItems, sizeVetItems, route, penalty, stopData))
            {
                if constexpr(PrintGA)
                    std::printf("GA END\n");

                bool doContinue = true;

                if(ptrVetBin)
                {
                    for(int i=0; i < nextBin; ++i)
                    {
                        if(binAreEqual((*ptrVetBin)[i], binAux))
                        {
                            doContinue = false;
                            break;
                        }
                    }

                    if(doContinue)
                    {
                        copiaBin(binAux, (*ptrVetBin)[nextBin]);
                        nextBin += 1;

                        if(nextBin == ptrVetBin->size())
                           return true;
                    }
                }
                else
                {
                    copiaBin(binAux, bin);
                    return true;
                }
            }

            it.first = penalty.getValue();
            if(it.first < bestRandKeyPenalty)
            {
                bestRandKeyPenalty = it.first;
                bestPenalty = penalty;          // raw values of the best chromosome
                bestRandKey.vetRandomKey = it.second->vetRandomKey;
            }
            it.second->fitness = it.first;
            it.second->evaluated = true;

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

                randKeySon->evaluated = false;

                for(int t = 0; t < sizeVetItems; ++t)
                {
                    double randVal = RandNs::getRandDouble();
                    if(randVal < input.gaProbToInheritFromElite)
                    {
                        randKeySon->vetRandomKey[t] = randKeyElite->vetRandomKey[t];
                        int index = sizeVetItems+vetItems[t];
                        randKeySon->vetRandomKey[index] = randKeyElite->vetRandomKey[index];
                    }
                    else
                    {
                        randKeySon->vetRandomKey[t] = randKeyRest->vetRandomKey[t];
                        int index = sizeVetItems+vetItems[t];
                        randKeySon->vetRandomKey[index] = randKeyRest->vetRandomKey[index];
                    }
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

                    for(int t = 0; t < sizeVetItems+2*InstanceNS::instanciaG.numItens; ++t)
                    {
                        vetRandKey[reciveInd].second->vetRandomKey[t] =
                            vetRandKey[selecedInd].second->vetRandomKey[t];
                    }

                    vetRandKey[reciveInd].second->fitness = vetRandKey[selecedInd].second->fitness;
                    vetRandKey[reciveInd].second->evaluated = vetRandKey[selecedInd].second->evaluated;
                }
            }
        }
    }

    if constexpr(PrintGAFail)
    {
        std::printf("GA failed. Best penalty %.4f:\n", bestRandKeyPenalty);
        bestPenalty.print();

        // Decode the best chromosome again (the decoder is deterministic) to list the items left out
        static const char *errNames[] = {"Geometric", "Lifo", "Support", "Axle", "Balance",
                                         "Compact", "Fragility", "Boundary"};
        randKey.vetRandomKey = bestRandKey.vetRandomKey;
        copiaBin(bin, binAux);
        decoderRandomKey(randKey, binAux, vetItems, sizeVetItems, route, penalty, stopData);

        std::printf("Not packed (%d of %d items):\n", (int)vetLastNotPacked.size(), sizeVetItems);
        std::printf("  %5s %5s %4s %18s %8s %10s  %s\n",
                    "item", "cust", "stop", "dims (x y z)", "weight", "volume", "blocked by");
        for(auto [itemId, err] : vetLastNotPacked)
        {
            InstanceNS::Item &it = InstanceNS::instanciaG.vetItens[itemId];
            std::printf("  %5d %5d %4d %6.0f %5.0f %5.0f %8.1f %10.0f  %s\n",
                        itemId, it.customer, findPos(route, itemId),
                        it.vetDim[0], it.vetDim[1], it.vetDim[2],
                        it.weight, it.volume, errNames[err]);
        }
        if(!vetLastSoftPlaced.empty())
        {
            std::printf("Placed with violation (%d items):\n", (int)vetLastSoftPlaced.size());
            for(auto [itemId, err] : vetLastSoftPlaced)
            {
                InstanceNS::Item &it = InstanceNS::instanciaG.vetItens[itemId];
                std::printf("  %5d %5d %4d %6.0f %5.0f %5.0f %8.1f %10.0f  %s\n",
                            itemId, it.customer, findPos(route, itemId),
                            it.vetDim[0], it.vetDim[1], it.vetDim[2],
                            it.weight, it.volume, errNames[err]);
            }
        }
        std::printf("Bin fill: %d items, %.1f%% of volume\n\n", binAux.numItens,
                    100.0 * (totalVolume - bestPenalty[Penalty::TypePenalty::VolItemsNotPacked]) /
                            (binAux.binDim[0] * binAux.binDim[1] * binAux.binDim[2]));
    }

    return false;
}

void GA_NS::startPopulation(Vector<RandomKeyPair> &vetRandKey, int vetItemsSize)
{
    for(auto &it : vetRandKey)
        it.second->generateRandomKey(vetItemsSize);
}
