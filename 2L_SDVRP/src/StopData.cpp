#include "Solucao.h"

using namespace InstanceNS;

namespace StopDataNS
{

StopData::StopData():	sizeJ(100*instanciaG.numClientes),
                        matStopItem(10, sizeJ),
                        sizeMatStopItem(10),
                        vetCustStop(instanciaG.numClientes),
                        vetStopsSet(10)
{}

void StopData::reset()
{
    sizeMatStopItem.setAll(0);
    vetCustStop.setAll(-1);
    for(auto& it:vetStopsSet)
        it.erase(it.begin(), it.end());

}

void StopData::compute(SolucaoNS::Rota &route, VectorI &vetItems, int sizeVetItems)
{
    reset();
    for(int i=1; i < (route.numPos-1); ++i)
        vetCustStop[route.vetRota[i]] = i;

    for(int i=0; i < sizeVetItems; ++i)
    {
        int item = vetItems[i];

        int cust = instanciaG.vetItens[item].customer;
        int stop = vetCustStop[cust];

        if(stop <= 0 || stop >= 10)
        {
            std::printf("ERROR, STOP(%d) <= 0\ncust: %d\nItem: %d\n", stop, cust, item);
            PRINT_THROW();
        }

        if(sizeMatStopItem[stop] >= sizeJ)
        {
            std::printf("Erro, j acess in matStopItem\n");
            std::printf("j: %d", sizeMatStopItem[stop]);
            std::printf("max: %d\n", sizeJ);
            PRINT_THROW();
        }

        matStopItem(stop, sizeMatStopItem[stop]) = item;
        sizeMatStopItem[stop] += 1;
    }


    for(int stop=(route.numPos-2); stop > 0; --stop)
    {
        // Copy items in stop to stop - 1
        for(int i=0; i < sizeMatStopItem[stop]; ++i)
        {
            if(sizeMatStopItem[stop-1] >= sizeJ)
            {
                std::printf("Error: Overflow during copy at stop %d (j: %d, max: %d)\n",
                            stop - 1, sizeMatStopItem[stop - 1], sizeJ);
                PRINT_THROW();
            }

            matStopItem(stop-1, sizeMatStopItem[stop-1]) = matStopItem(stop, i);
            sizeMatStopItem[stop-1] += 1;
        }

    }

    for(int stop=0; stop <= (route.numPos-2); ++stop)
    {
        for(int i=0; i < sizeMatStopItem[stop]; ++i)
            vetStopsSet[stop].insert(matStopItem(stop, i));
    }

}

}
