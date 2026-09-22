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
                                   VectorI&				vetItems,
                                   int 					vetItemsSize,
                                   SolucaoNS::Rota& 	route,
                                   SolucaoNS::Penalty& 	penalty)
{
    /*
    std::set<int> set;
    for(int i=0; i < vetItemsSize; ++i)
        set.insert(vetItems[i]);
    */
    std::sort(vetItems.begin(), vetItems.begin()+vetItemsSize);

    bin.reset();
    sort_two_vectors(randKey.vetRandomKey.data(), vetItems, vetItemsSize);

    /*
    for(int i=0; i < vetItemsSize; ++i)
    {
        if(!set.contains(vetItems[i]))
        {
            std::printf("Error, item(%d) is not in the set\n", vetItems[i]);
            PRINT_THROW();
        }
    }

    set = std::set<int>();
    for(int i=0; i < vetItemsSize; ++i)
    {
        if(!set.contains(vetItems[i]))
            set.insert(vetItems[i]);
        else
        {
            std::printf("Error, item(%d) is repetead\n", vetItems[i]);
            PRINT_THROW();
        }
    }

    */
    penalty.set0();

    return packItemsIntoBin(bin, route, vetItems, vetItemsSize, false, false, &penalty);
}
