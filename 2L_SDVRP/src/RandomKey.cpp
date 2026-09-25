#include "RandomKey.h"
#include "rand.h"
#include "ConstrutivoBin2.h"
#include <set>

using namespace ConstrutivoBin2NS;

void RandomKeyNS::RandomKey::generateRandomKey(int numItems)
{
    assertm(numItems > MaxNumItemsBin, "Error, MaxNumItemsBin is less them numItems");

    for(int i=0; i < numItems; ++i)
        vetRandomKey[i] = RandNs::getRandDouble();
}



bool RandomKeyNS::decoderRandomKey(RandomKey&			randKey,
                                   SolucaoNS::Bin&		bin,
                                   const VectorI&		vetItems,
                                   int 					vetItemsSize,
                                   SolucaoNS::Rota& 	route,
                                   SolucaoNS::Penalty& 	penalty)
{
    static VectorI vetItemsTemp(MaxNumItemsBin);
    copyVet(vetItems, vetItemsTemp, vetItemsSize);

    bin.reset();
    sort_two_vectors(randKey.vetRandomKey.data(), vetItemsTemp, vetItemsSize);
    penalty.set0();

    return
        packItemsIntoBin(bin, route, vetItemsTemp, vetItemsSize, false, false, &penalty);

}
