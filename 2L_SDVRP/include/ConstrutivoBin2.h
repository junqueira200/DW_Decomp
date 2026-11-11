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
#include "StopData.h"

namespace ConstrutivoBin2NS
{
// Items not packed in the last packItemsIntoBin pass: (itemId, deepest failed check)
inline std::vector<std::pair<int, int>> vetLastNotPacked;


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
bool packItemsIntoBin(SolucaoNS::Bin     	&bin,
                      SolucaoNS::Rota    	&rota,
                      const VectorI      	&itemIds,
                      int                 	numItems,
                      bool                	sortItemIds = true,
                      bool                	print 		= false,
                      SolucaoNS::Penalty 	*penalty 	= nullptr,
                      bool				  	resetBin 	= false,
                      StopDataNS::StopData	*stopData 	= nullptr,
                      double				*randomKey  = nullptr);

/**
 * @brief Collects all item IDs associated with the customers along the route.
 */
void getRouteItems(const SolucaoNS::Rota &rota, VectorI &outItemIds);

double getTotalVolume(const SolucaoNS::Bin& bin, const VectorI& itemIds, int numItems);
double getTotalLeftArea(const SolucaoNS::Bin& bin, const VectorI& itemIds, int numItems);
double getTotalAreaSupp(const SolucaoNS::Bin& bin, const VectorI& itemIds, int numItems);
double getTotalWeight(const SolucaoNS::Bin& bin, const VectorI& itemIds, int numItems);

void penalizeSolution(SolucaoNS::Rota 		&route,
                      SolucaoNS::Bin 		&bin,
                      SolucaoNS::Penalty 	&penalty,
                      StopDataNS::StopData 	&stopData);

} // namespace ConstrutivoBin2NS

#endif // INC_2L_SDVRP_CONSTRUTIVOBIN_H
