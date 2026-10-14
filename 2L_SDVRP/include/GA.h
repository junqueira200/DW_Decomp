/* ****************************************
 * ****************************************
 *  Data:    22/09/26
 *  Arquivo: GA.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef GA_H
#define GA_H

#include "RandomKey.h"
#include "Solucao.h"

using RandomKeyPair = std::pair<double, RandomKeyNS::RandomKey *>;

namespace GA_NS
{

bool ga(SolucaoNS::Bin      &bin,
        SolucaoNS::Rota     &route,
        const VectorI *const vetItems = nullptr,
        int                  sizeVetItems = 0);
void startPopulation(Vector<RandomKeyPair> &vetRandKey, int vetItemsSize);
} // namespace GA_NS

#endif // GA_H
