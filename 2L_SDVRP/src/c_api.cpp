#include "c_api.h"
#include "ConstrutivoBin.h"
#include "Instancia.h"
#include "TesteOroloc3D.h"
#include <unordered_set>
#include <omp.h>

using namespace ConstrutivoBinNS;
using namespace InstanceNS;
using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace TesteOroloc3D_NS;

void ini_3D_Packing(char *strInst_c)
{

    static bool doInit = true;

    if(doInit)
    {
        std::string strInst(strInst_c);
        input.strInstCompleto = strInst;
        InstanceNS::readOroloc3D2(input.strInstCompleto);

        RandNs::startEngine(output.semente, false);
        input.strInst = getNomeInstancia(input.strInstCompleto);

        startConstGlobalVaribles();

        output.setup();
        std::cout << "INST: " << input.strInst << " SEMENTE: " << RandNs::estado_ << " "
                  << output.data << "";

        doInit = false;
    }

}

bool testRoute(int *vet_c, int vetSize)
{
    static std::unordered_set<SolucaoNS::Rota, SolucaoNS::HashRoute> routeSetFeasible;
    static std::unordered_set<SolucaoNS::Rota, SolucaoNS::HashRoute> routeSetInfeasible;
    static SolucaoNS::Rota 					   route;
    static SolucaoNS::Bin  					   bin;
    static VectorI							   vetItems(instanciaG.numItens);
    static bool 							   doLink = true;

    if(doLink)
    {
        route.binPtr = &bin;
        doLink = false;
    }

    bin.reset();
    route.reset();

    for(int i=0; i < vetSize; ++i)
        route.vetRota[i] = vet_c[i];

    route.numPos = vetSize;
    route.computeDistance();

    if(routeSetFeasible.contains(route))
        return true;

    if(routeSetInfeasible.contains(route))
        return false;

    int numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);
    std::reverse(vetItems.begin(), vetItems.begin() + bin.numItens);

    bool feasible =
    ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
                                            std::numeric_limits<int>::max(), &route);

    if(feasible)
    {
        routeSetFeasible.insert(route);
        return true;
    }


    std::vector<Cuboid>   vetCuboids;
    Collections::IdVector stopIds;

    convertVectorOfItensToVectorOfCuboids(
        vetItems, vetCuboids, bin.numItens, route);
    // int lastCustomerId =
    // instanciaG.vetItens[bin.vetItens[bin.numItens-1]].customer;

    for(int i = 1; i < route.numPos - 1; ++i)
    {
        stopIds.push_back(route.vetRota[i]);
        // if(lastCustomerId == sol.vetRota[veic].vetRota[i])
        //      break;
    }

    static InputParameters inputParam;
    inputParam.ContainerLoading.LoadingProblem.Variant =
        LoadingProblemParams::VariantType::AllConstraints;
    inputParam.SetLoadingFlags();

    static LoadingChecker loadingChecker(inputParam.ContainerLoading);
    static Container      container((int)instanciaG.vetDimVeiculo[0],
                        (int)instanciaG.vetDimVeiculo[1],
                        (int)instanciaG.vetDimVeiculo[2],
                        (int)instanciaG.maxPayload);

           //
    PackingType    lastType;
    StatusOroloc3D statusOroc3D;
    double         tempoCpu;
    int totalTry = 0;

    for(int i=0; i < 5; ++i)
    {
        std::vector<Array<int, 4>> vetArray;
        // std::cout<<"n: "<<n<<"\n";
        double ompStart = omp_get_wtime();
        auto status = loadingChecker.ConstraintProgrammingSolver(
        PackingType::Complete, container, stopIds, vetCuboids, input.cpSatTime, vetArray);
        // std::cout<<"ret\n";
        double ompEnd = omp_get_wtime();

        if(status == LoadingStatus::Infeasible || status == LoadingStatus::Invalid)
        {
            routeSetInfeasible.insert(route);
            return false;
        }
        else if(status == LoadingStatus::FeasOpt)
        {
            routeSetFeasible.insert(route);
            return true;
        }
    }

    std::printf("ERROR, STATUS FROM OR TOOLS IS TIME LIMIT\n");
    PRINT_THROW();

    return false;
}
