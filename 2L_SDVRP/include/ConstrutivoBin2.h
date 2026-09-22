/* ****************************************
 * ****************************************
 *  Data:    14/09/26
 *  Arquivo: ConstrutivoBin.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L/3L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_CONSTRUTIVOBIN2_H
#define INC_2L_SDVRP_CONSTRUTIVOBIN2_H

#include "Instancia.h"
#include "Solucao.h"
#include <vector>

namespace ConstrutivoBin2NS
{

enum class SortStrategy
{
    ReverseRoute_AreaDesc = 0,
    ReverseRoute_VolumeDesc,
    ReverseRoute_WeightDesc
};

/**
 * @brief Main heuristic entry point. Packs all items of the given route into its bin.
 *
 * @param rota The vehicle route containing the customer sequence and bin pointer.
 * @param print If true, outputs diagnostic information.
 * @return true if all items were packed feasibly.
 * @return false if the items could not be packed feasibly.
 */
bool packRoute(SolucaoNS::Rota &rota, VectorI _itemIds, bool print = false);

/**
 * @brief Attempts to pack a specific list of items into a bin using a given sorting strategy.
 */
bool packItemsIntoBin(SolucaoNS::Bin          &bin,
                      SolucaoNS::Rota         &rota,
                      const VectorI           &itemIds,
                      int					   numItems,
                      bool					   sortItemIds 	= true,
                      bool                     print 		= false,
                      SolucaoNS::Penalty*	   penalty		= nullptr);

/**
 * @brief Collects all item IDs associated with the customers along the route.
 */
void getRouteItems(const SolucaoNS::Rota &rota, VectorI &outItemIds);

} // namespace ConstrutivoBinNS

#endif // INC_2L_SDVRP_CONSTRUTIVOBIN_H
