#ifndef TESTEOROLOC3D_H
#define TESTEOROLOC3D_H

#include "Instancia.h"
#include "Solucao.h"

#include "ProblemParameters.h"
#include "BCRoutingParams.h"
#include "LoadingChecker.h"

namespace TesteOroloc3D_NS
{
    enum StatusOroloc3D
    {
        INFEASIBLE = 0,
        TIME_LIMIT = 1,
        FEASIBLE   = 2
    };

    void testeOroloc3D();
    void testeOroloc3D_2();

    void convertVectorOfItensToVectorOfCuboids(const VectorI& vetItens, std::vector<ContainerLoading::Cuboid>& vetCuboids,
                                               int numItems);
    void appendToFile(const std::string& fileName, const std::string& content);
    void readSolOroloc3D_2(SolucaoNS::Solucao& sol);
    void printSol(SolucaoNS::Solucao& sol);

}

#endif // TESTEOROLOC3D_H
