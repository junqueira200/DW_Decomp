#ifndef GA_H
#define GA_H

#include "Solucao.h"
#include "RandomKey.h"


using  RandomKeyPair = std::pair<double, RandomKeyNS::RandomKey*>;

namespace GA_NS
{



    bool ga(SolucaoNS::Bin& bin,
            SolucaoNS::Rota &route,
            VectorI* 		vetItems	=nullptr,
            int 			sizeVetItems=0);
    void startPopulation(Vector<RandomKeyPair>& vetRandKey, int vetItemsSize);
}

#endif // GA_H
