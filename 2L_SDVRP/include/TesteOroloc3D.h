#ifndef TESTEOROLOC3D_H
#define TESTEOROLOC3D_H

#include "Instancia.h"

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
    void convertVectorOfItensToVectorOfCuboids(const VectorI& vetItens, std::vector<ContainerLoading::Cuboid>& vetCuboids);
    void appendToFile(const std::string& fileName, const std::string& content);

}

#endif // TESTEOROLOC3D_H
