/* ****************************************
 * ****************************************
 *  Data:    10/02/25
 *  Arquivo: BinPackingCP.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_BINPACKINGCP_H
#define INC_2L_SDVRP_BINPACKINGCP_H

#include "Solucao.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"

using namespace operations_research;
using namespace operations_research::sat;

namespace BinPackingCP_NS
{

    bool cpSatBinPacking(SolucaoNS::Bin &binResult, VectorI &vetItens, int tamVet);
    void testaCpSatBinPacking(int numItens);
    void criaEPs(SolucaoNS::Bin &bin);

}

#endif //INC_2L_SDVRP_BINPACKINGCP_H
