/* ****************************************
 * ****************************************
 *  Data:    20/09/26
 *  Arquivo: RandomKey.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/


#include "RandomKey.h"
#include "ConstrutivoBin2.h"
#include "rand.h"


#include <set>

using namespace ConstrutivoBin2NS;
using namespace InstanceNS;
using namespace RandomKeyNS;

void RandomKeyNS::RandomKey::generateRandomKey(int numItems)
{
    assertm(numItems > MaxNumItemsBin, "Error, MaxNumItemsBin is less them numItems");

    for(int i = 0; i < numItems+2*instanciaG.numItens; ++i)
        vetRandomKey[i] = RandNs::getRandDouble();

    this->evaluated = false;
}

bool RandomKeyNS::decoderRandomKey(RandomKey          	&randKey,
                                   SolucaoNS::Bin     	&bin,
                                   const VectorI      	&vetItems,
                                   int                 	vetItemsSize,
                                   SolucaoNS::Rota    	&route,
                                   SolucaoNS::Penalty 	&penalty,
                                   StopDataNS::StopData	&stopData)
{
    static VectorI vetItemsTemp(MaxNumItemsBin);
    copyVet(vetItems, vetItemsTemp, vetItemsSize);

    //bin.reset();
    sort_two_vectors(randKey.vetRandomKey.data(), vetItemsTemp, vetItemsSize);
    //grupByCustomer(vetItemsTemp, vetItemsSize, route);
    penalty.set0();

    randKey.evaluated = true;

    return packItemsIntoBin(bin, route, vetItemsTemp, vetItemsSize, false, false, &penalty, false,
                            &stopData, &randKey.vetRandomKey[0]);
}

void RandomKeyNS::grupByCustomer(VectorI &vetItems, int sizeVetItems, SolucaoNS::Rota &route)
{
    if(sizeVetItems <= 0 || route.numPos <= 2)
        return;

    // Use a copy to safely hold the rearranged items
    static VectorI temp(instanciaG.numItens);
    int idx = 0;

    // Iterate through customers in reverse delivery order
    // (from the last customer at route.numPos - 2 down to the first at index 1)
    for(int stop = route.numPos - 2; stop >= 1; --stop)
    {
        int cust = route.vetRota[stop];
        // Retain the pre-computed relative order of items for this customer
        for(int i = 0; i < sizeVetItems; ++i)
        {
            int item = vetItems[i];
            if(InstanceNS::instanciaG.vetItens[item].customer == cust)
            {
                temp[idx++] = item;
            }
        }
    }

    // Safety check: ensure all items belonged to customers in this route
    if(idx != sizeVetItems)
    {
        std::printf("ERROR in grupByCustomer: item count mismatch (assigned: %d, expected: %d)\n",
                    idx, sizeVetItems);
        PRINT_THROW();
    }

    // Copy the rearranged items back into vetItems
    for(int i = 0; i < sizeVetItems; ++i)
        vetItems[i] = temp[i];

}
