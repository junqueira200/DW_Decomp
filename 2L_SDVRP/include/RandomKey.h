#ifndef RANDOMKEY_H
#define RANDOMKEY_H

#include "Instancia.h"
#include "Solucao.h"

#include <numeric>
#include <algorithm>
#include <functional>

namespace RandomKeyNS
{

    class RandomKey
    {
      public:

        Array<double, MaxNumItemsBin> vetRandomKey;

        RandomKey(){};
        void generateRandomKey(int numItems);

    };

    bool decoderRandomKey(RandomKey& 			randKey,
                          SolucaoNS::Bin&	 	bin,
                          VectorI&	 			vetItems,
                          int	 				vetItemsSize,
                          SolucaoNS::Rota&		route,
                          SolucaoNS::Penalty& 	penalty);

    template <typename Key, typename Value, typename Compare = std::less<Key>>
    void sort_two_vectors(Key* 				keys,
                          Vector<Value>& 	values,
                          int 				vecSize,
                          Compare 			comp = Compare{})
    {
        // Validate bounds
        if (vecSize > MaxNumItemsBin || vecSize > MaxNumItemsBin)
        {
            PRINT_DEBUGG("", "vecSize exceeds vector capacity.")
            PRINT_THROW();
        }


        static Vector<size_t> p(MaxNumItemsBin);
        std::iota(p.begin(), p.begin()+vecSize, 0);

        std::sort(p.begin(), p.begin()+vecSize, [&](size_t i, size_t j) {
            return comp(keys[i], keys[j]);
        });

        static Vector<Key> temp_keys(MaxNumItemsBin);
        static std::vector<Value> temp_values(MaxNumItemsBin);

        for (size_t i = 0; i < vecSize; ++i)
        {
            temp_keys[i] 	= keys[p[i]];
            temp_values[i] 	= values[p[i]];
        }


        for (size_t i = 0; i < vecSize; ++i)
        {
            keys[i] 	= temp_keys[i];
            values[i]	= temp_values[i];
        }
    }


}

#endif // RANDOMKEY_H
