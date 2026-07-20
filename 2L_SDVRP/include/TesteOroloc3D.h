#ifndef TESTEOROLOC3D_H
#define TESTEOROLOC3D_H

#include "Instancia.h"
#include "Solucao.h"

#include "BCRoutingParams.h"
#include "LoadingChecker.h"
#include "ProblemParameters.h"

namespace TesteOroloc3D_NS
{
    enum StatusOroloc3D
    {
        INFEASIBLE = 0,
        TIME_LIMIT = 1,
        FEASIBLE = 2
    };

    void testeOroloc3D();
    void testeOroloc3D_2();

    void convertVectorOfItensToVectorOfCuboids(
        const VectorI                         &vetItens,
        std::vector<ContainerLoading::Cuboid> &vetCuboids,
        int                                    numItems,
        SolucaoNS::Rota                       &rota);

    void appendToFile(const std::string &fileName, const std::string &content);
    void readSolOroloc3D_2(SolucaoNS::Solucao &sol);
    void printSol(SolucaoNS::Solucao &sol);
    void writeToFile(const std::string& str, const std::string& strFile);

    void writeSoltionOroloc2(const SolucaoNS::Bin& bin, const SolucaoNS::Rota& route);


} // namespace TesteOroloc3D_NS

#endif // TESTEOROLOC3D_H
