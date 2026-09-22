#ifndef C_API_H
#define C_API_H

#include "Instancia.h"
#include "InputOutput.h"
#include <unordered_set>

using namespace ParseInputNS;

typedef std::unordered_set<SolucaoNS::Rota, SolucaoNS::HashRoute> HASH_ROUTE;

struct RouteData
{
    HASH_ROUTE routeSetFeasible;
    HASH_ROUTE routeSetInfeasible;
};

inline RouteData routeData;
constexpr bool   useHash = true;


extern "C"
{
    void ini_3D_Packing(char* strInst_c, int oroloc3D);
    int testRoute(int* vet_c, int vetSize, int onlyHeuristic=0, int doInverseRoute=0);
    int testRouteCapVol(int* vet_c, int vetSize);
    double getDistanceRoute(int* vet_c, int vetSize);
    void roundDistances();

    int routeIsInFeasibleSet(int* vet_c, int vetSize);
    int routeIsInNotfeasibleSet(int* vet_c, int vetSize);
    int heuristicPacking(int* vet_c, int vetSize);
    int exactPacking(int* vet_c, int vetSize);

    void setClassical3DPackingProblem(int doPrint);
    void setLoadingOnlyProblem(int doPrint);
    void setOroloc3DProblem(int doPrint);
    void saveProblem(Input& inputTemp);
    void setProblem(Input& inputTemp);

    int 	getNumberOfCustoms();
    int 	getNumberOfTrucks();
    int 	getDemandFromCustomr(int i);
    double 	getVolumeFromCustomr(int i);
    double 	getDistance(int i, int j);
    int 	getVehicleCapacity();
    double 	getVehicleVolume();

    void setDualFeasibleFunction(double ep0, double ep2);
    double getDualVolume(int cust);
    double getDualVolumeTotal();

    inline bool doBreakTestRoute = false;
    constexpr bool useDualFunction = false;

    // TODO: Criar uma estrutura para as rotas;
    // TODO: Criar uma funcao para recuperar as solucoes do binpacking
}

#endif // C_API_H
