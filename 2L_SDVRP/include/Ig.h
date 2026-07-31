/* ****************************************
 * ****************************************
 *  Data:    11/12/24
 *  File:    IG.cpp
 *  Author:  Igor de Andrade Junqueira
 *  Project: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_IG_H
#define INC_2L_SDVRP_IG_H

#include "Solucao.h"

namespace IgNs
{
    bool metaheuristicaIg(SolucaoNS::Solucao &solucao);

    INLINE
    void rmRota(SolucaoNS::Solucao &sol, int r)
    {
        sol.distTotal -= sol.vetRota[r].distTotal;
        sol.vetRota[r].reset();
        sol.vetBin[r].reset();
    }

} // namespace IgNs

#endif // INC_2L_SDVRP_IG_H
