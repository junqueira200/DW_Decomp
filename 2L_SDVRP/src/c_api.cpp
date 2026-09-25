#include "c_api.h"
#include "ConstrutivoBin.h"
#include "Instancia.h"
#include "TesteOroloc3D.h"
#include <unordered_set>
#include <omp.h>
#include "DualFeasibleFunctions.h"
#include "GA.h"

using namespace ConstrutivoBinNS;
using namespace InstanceNS;
using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace TesteOroloc3D_NS;
using namespace DualFeasibleFunctionsNS;
using namespace GA_NS;

void ini_3D_Packing(char *strInst_c, int oroloc3D)
{

    //static bool doInit = true;
    //std::printf("Inst: %s\n", strInst_c);
    //if(doInit)

    //std::printf("oroloc3D: %d\n", oroloc3D);
    //PRINT_THROW()

    if(oroloc3D == 0)
        setClassical3DPackingProblem(true);
    else
        setOroloc3DProblem(true);

    std::string strInst(strInst_c);
    input.strInstCompleto = strInst;
    if(oroloc3D)
        InstanceNS::readOroloc3D2(input.strInstCompleto);
    else
        InstanceNS::read3dInstance(input.strInstCompleto);

    output.semente =  RandNs::startEngine(output.semente, false);
    input.strInst = getNomeInstancia(input.strInstCompleto);

    startConstGlobalVaribles();

    output.setup();
    std::cout << "INST: " << input.strInst << " SEMENTE: " << RandNs::estado_ << " "
              << output.data << "";

}

int routeIsInFeasibleSet(int* vet_c, int vetSize)
{
    static SolucaoNS::Rota 	route;

    route.reset();
    for(int i=0; i < vetSize; ++i)
        route.vetRota[i] = vet_c[i];

    route.numPos = vetSize;
    route.computeDistance();

    return (int)routeData.routeSetFeasible.contains(route);

}

int routeIsInNotfeasibleSet(int* vet_c, int vetSize)
{
    static SolucaoNS::Rota 	route;

    route.reset();
    for(int i=0; i < vetSize; ++i)
        route.vetRota[i] = vet_c[i];

    route.numPos = vetSize;
    route.computeDistance();

    return (int)routeData.routeSetInfeasible.contains(route);
}

int heuristicPacking(int* vet_c, int vetSize)
{
    static SolucaoNS::Rota 	route;
    static SolucaoNS::Bin  	bin;
    static VectorI			vetItems(instanciaG.numItens);

    route.binPtr = &bin;
    route.reset();
    bin.reset();

    for(int i=0; i < vetSize; ++i)
        route.vetRota[i] = vet_c[i];

    route.numPos = vetSize;
    route.computeDistance();

    int numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);

    double totalVolume = 0.0;
    int totalDemand    = 0;
    for(int i=0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        totalVolume += item.vetDim[0]*item.vetDim[1]*item.vetDim[2];
        totalDemand += instanciaG.vetItens[i].weight;
    }

    //for(int i=0; i < vetSize; ++i)
    //    totalDemand += instanciaG.vetDemandaCliente[route.vetRota[i]];

    static double volumeOfVeicule = getVehicleVolume();
    if(totalVolume > volumeOfVeicule)
        return 0;

    if(totalDemand > (int)instanciaG.maxPayload)
        return 0;

    if(useHash)
    {
        if(routeData.routeSetFeasible.contains(route))
            return 1;

        if(routeData.routeSetInfeasible.contains(route))
            return 0;
    }

    std::reverse(vetItems.begin(), vetItems.begin() + numItems);

    bool feasible = ga(bin, route, &vetItems, numItems);
         //ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
         //                                        std::numeric_limits<int64_t>::max(), &route);
    if(feasible)
    {
        if(useHash)
            routeData.routeSetFeasible.insert(route);

        return 1;
    }
    else
        return 0;

}

int exactPacking(int* vet_c, int vetSize)
{
    //std::printf("exactPacking\n");
    SolucaoNS::Rota 	route;
    SolucaoNS::Bin  	bin;
    VectorI			vetItems(instanciaG.numItens);

    route.binPtr = &bin;
    route.reset();
    bin.reset();

    for(int i=0; i < vetSize; ++i)
    {
        //std::printf("%d ", vet_c[i]);
        route.vetRota[i] = vet_c[i];
    }

    //std::printf("\n");

    route.numPos = vetSize;
    route.computeDistance();

    int numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);

    //std::printf("vetItems: ");
    //for(int i=0; i < numItems; ++i)
    //    std::printf("%d ", vetItems[i]);

    //std::printf("\n\n");

    double totalVolume = 0.0;
    double totalDemand = 0.0;
    for(int i=0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        totalVolume += item.vetDim[0]*item.vetDim[1]*item.vetDim[2];
        //totalDemand += instanciaG.vetItens[i].weight;
    }

    for(int i=0; i < vetSize; ++i)
        totalDemand += instanciaG.vetDemandaCliente[route.vetRota[i]];

    static double volumeOfVeicule = getVehicleVolume();
    if(totalVolume > volumeOfVeicule)
    {
        std::printf("max Vol; totalVolume: %.2f; volumeOfVeicule: %.2f\n\n", totalVolume,
                     volumeOfVeicule);
        PRINT_THROW();
        return 0;
    }

    if(totalDemand > instanciaG.maxPayload)
    {
        std::printf("maxPayload\ntotalDemand: %f; maxPayload: %f", totalDemand,
                    instanciaG.maxPayload);
        PRINT_THROW();
        return 0;
    }

    if(useHash)
    {
        if(routeData.routeSetFeasible.contains(route))
            return 1;

        if(routeData.routeSetInfeasible.contains(route))
            return 0;
    }

    std::reverse(vetItems.begin(), vetItems.begin() + numItems);

    //bool feasible =
    //    ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
    //                                            50, &route);

    numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);
    std::reverse(vetItems.begin(), vetItems.begin() + numItems);

    bool dualFeasible = check(vetItems);
    if(!dualFeasible)
    {
        std::printf("dual infeasible!\n");
    }



    std::vector<Cuboid>   vetCuboids;
    Collections::IdVector stopIds;

    convertVectorOfItensToVectorOfCuboids(
        vetItems, vetCuboids, numItems, route);
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
    static Container    container((int)instanciaG.vetDimVeiculo[0],
                               (int)instanciaG.vetDimVeiculo[1],
                               (int)instanciaG.vetDimVeiculo[2],
                               (int)instanciaG.maxPayload);

    PackingType    lastType;
    StatusOroloc3D statusOroc3D;
    double         tempoCpu;
    int totalTry = 0;

    std::printf("\n\n************************************\n");
    std::printf(    "**************INI CP-SAT************\n\n");

    double int_fk, int_fFA, int_fRA, int_fTA;

    for(int i=0; i < 1; ++i)
    {
        bin.reset();

        std::vector<Array<int, 4>> vetArray;
        // std::cout<<"n: "<<n<<"\n";
        //double ompStart = omp_get_wtime();
        auto status = loadingChecker.ConstraintProgrammingSolver(
            PackingType::Complete, container, stopIds, vetCuboids, input.cpSatTime, vetArray,
            int_fk, int_fFA, int_fRA, int_fTA);
        // std::cout<<"ret\n";
        //double ompEnd = omp_get_wtime();

        std::printf("************************************\n");
        std::printf("**************END CP-SAT************\n\n");

        if(status == LoadingStatus::Infeasible)
        {

            if(!dualFeasible)
            {
                std::printf("Packing is NOT dualFeasible\nIt worked!\n");
                //PRINT_THROW();
            }

            //if(feasible)
            //{
            //    std::printf("ERROR, LoadingStatus::Infeasible, and heuristic solved\n\n");
            //    PRINT_THROW();
            //}
            if(useHash)
                routeData.routeSetInfeasible.insert(route);

            return 0;
        }
        else if(status == LoadingStatus::Invalid)
        {
            std::printf("ERROR, status: Invalid\n");
            PRINT_THROW();
        }
        else if(status == LoadingStatus::Unknown)
        {
            //std::printf("Status: Unknown, Time Limit?");
            std::printf("ERROR, status: Unknown\n");
            PRINT_THROW();
            continue;
        }
        else if(status == LoadingStatus::FeasOpt)
        {
            int item = 0;
            for(Array<int, 4> &array : vetArray)
            {
                bin.vetItemId[item] = vetItems[item];
                bin.vetPosItem[item].set(array[0], array[1], array[2]);
                bin.vetRotacao[item] = (InstanceNS::Rotation)array[3];

                item += 1;
            }

            bin.numItens = numItems;
            bin.computeLoadingBalancing();

            if(!bin.checkFeasibility(&route, true, true))
            {
                std::printf("ERROR, feasible solution from CP model is not feasible!\n");
                std::cout<<bin.printPlot()<<"\n";
                PRINT_THROW();
                continue;
            }

            if(useHash)
                routeData.routeSetFeasible.insert(route);

            if(!checkIfPackedAllTheItems(bin, vetItems, numItems))
            {
                PRINT_THROW();
            }

            if(!dualFeasible)
            {
                std::printf("Error, packing is dual infeasible! CP generated a feasible Solution\n");
                PRINT_THROW();
            }

            return 1;
        }
        else
        {
            std::printf("Status Unknown: %d\n", status);
            PRINT_THROW();
        }

    }

    PRINT_THROW();
    return 0;
}

int testRouteCapVol(int* vet_c, int vetSize)
{
    double totalVolume = 0.0;
    double totalDemand = 0.0;

    SolucaoNS::Rota 	route;
    for(int i=0; i < vetSize; ++i)
        route.vetRota[i] = vet_c[i];

    route.numPos = vetSize;
    double totalVolItems = 0;
    VectorI vetItems;

    int numItems = copiaItensClientes(route.vetRota, vetSize, vetItems, true);

    for(int i=0; i < numItems; ++i)
        totalVolItems += instanciaG.vetItens[vetItems[i]].volume;

    for(int i=0; i < vetSize; ++i)
    {
        //Item& item = instanciaG.vetItens[vetItems[i]];
        totalVolume += instanciaG.vetVolumeCliente[vet_c[i]];
        totalDemand += instanciaG.vetDemandaCliente[vet_c[i]];
        //totalDemand += instanciaG.vetItens[i].weight;
    }

    if(!doubleEqual(totalVolItems, totalVolume, 1E-5))
    {
        std::printf("Error, totalVolItems(%f) != totalVolume(%f)", totalVolItems,
                                                                   totalVolume);
        PRINT_THROW();
    }

    static double volumeOfVeicule = getVehicleVolume();
    if(totalVolume > volumeOfVeicule)
        return 0;

    if(totalDemand > instanciaG.maxPayload)
    {
        std::printf("Demand\n");
        return 0;
    }

    return 1;
}

double getDistanceRoute(int* vet_c, int vetSize)
{
    double dist = 0.0;
    for(int i=0; i < (vetSize-1); ++i)
        dist += getDistance(vet_c[i], vet_c[i+1]);

    return dist;
}

void roundDistances()
{
    for(int i=0; i < instanciaG.numClientes; ++i)
    {
        for(int j=0; j < instanciaG.numClientes; ++j)
            instanciaG.matDist.get(i, j) = std::round(getDistance(i, j));
    }
}

int testRoute(int *vet_c, int vetSize, int onlyHeuristic, int doInverseRoute)
{

    //std::printf("testRoute\n");
    //onlyHeuristic = false;
    //std::printf("onlyHeuristic: %d\n", onlyHeuristic);

    HASH_ROUTE& 	routeSetFeasible   = routeData.routeSetFeasible;
    HASH_ROUTE& 	routeSetInfeasible = routeData.routeSetInfeasible;
    SolucaoNS::Rota route;
    SolucaoNS::Rota	routeInverse;
    SolucaoNS::Bin  bin;
    VectorI			vetItems(instanciaG.numItens);

    route.binPtr = &bin;

    bin.reset();
    route.reset();
    routeInverse.reset();

    //std::printf("vetSize: %d\n", vetSize);

    for(int i=0; i < vetSize; ++i)
    {
        route.vetRota[i] = vet_c[i];
        //std::printf("%d ", vet_c[i]);
    }
    //std::printf("\n");
    route.numPos = vetSize;
    route.computeDistance();

    if(doInverseRoute)
    {
        SolucaoNS::copiaRota(route, routeInverse);
        std::reverse(routeInverse.vetRota.begin(),
                     routeInverse.vetRota.begin()+routeInverse.numPos);
        routeInverse.computeDistance();

        //std::cout<<"Inverse: "<<routeInverse.printRota(false)<<"\n";
    }

    //std::cout<<"Testing route: "<<route.printRota()<<"\n";

    // 0 all infesible
    // 1 only forwerd route is fesible
    // 2 Both routes are feasible
    // 3 only backward route is feasible

    if(useHash)
    {
        if(routeSetFeasible.contains(route) && routeSetFeasible.contains(routeInverse) &&
           doInverseRoute)
            return 2;
        else if(routeSetFeasible.contains(route))
            return 1;
        else if(routeSetFeasible.contains(routeInverse) && doInverseRoute)
            return 3;

        if(routeSetInfeasible.contains(route) && routeSetInfeasible.contains(routeInverse) &&
           doInverseRoute)
            return 0;

        else if(routeSetInfeasible.contains(route))
            return 0;
    }


    int numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);

    //std::printf("vetItems: ");
    //for(int i=0; i < numItems; ++i)
    //    std::printf("%d ", vetItems[i]);

    //std::printf("\n\n");

    double totalVolume = 0.0;
    double totalDemand = 0.0;
    for(int i=0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        totalVolume += item.vetDim[0]*item.vetDim[1]*item.vetDim[2];
        //totalDemand += instanciaG.vetItens[i].weight;
    }

    for(int i=0; i < vetSize; ++i)
        totalDemand += instanciaG.vetDemandaCliente[route.vetRota[i]];

    static double volumeOfVeicule = getVehicleVolume();
    if(totalVolume > volumeOfVeicule)
    {
        //std::printf("max Vol; totalVolume: %.2f; volumeOfVeicule: %.2f\n\n", totalVolume,
        //             volumeOfVeicule);
        //PRINT_THROW();
        return 0;
    }

    if(totalDemand > instanciaG.maxPayload)
    {
        //std::printf("maxPayload\ntotalDemand: %d; maxPayload: %f", totalDemand,
        //            instanciaG.maxPayload);
        //PRINT_THROW();
        return 0;
    }

    bool dualFeasible = check(vetItems, numItems);

    //std::printf("numItems: %d\n", numItems);

    //std::printf("vetItems: ");
    //for(int i=0; i < numItems; ++i)
    //    std::printf("%d ", vetItems[i]);

    //std::printf("\n");

    std::reverse(vetItems.begin(), vetItems.begin() + numItems);

    //std::printf("vetItems: ");
    //for(int i=0; i < numItems; ++i)
    //    std::printf("%d ", vetItems[i]);


    /*
    if(onlyHeuristic == 0)
    {
        std::printf("\n\n*****************************************\n");
        std::printf("**************INI CONSTRUTIVO************\n\n");
    }
    */

    bool feasible = false;

    if(useDualFunction)
        check(vetItems);


    if(onlyHeuristic >= 1)
    {
        feasible = ga(bin, route, &vetItems, numItems);
        //ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
        //                                        25, &route);
        if(feasible)
        {
            //std::printf("GA work!\n");
            if(!dualFeasible)
            {
                std::printf("Error, packing is dual infeasible!\n");
                PRINT_THROW();
            }

            if(useHash)
                routeSetFeasible.insert(route);


            /*
            if(useDualFunction && doBreakTestRoute)
            {
                std::printf("Heuristic is feasible ??\n");
                goto jmpCp;
                PRINT_THROW();
            }
            */

        }
        //else
        //    std::printf("GA dident work!\n");
        //if(feasible)
        //    return true;
        return feasible;

        //input.cpSatTime = 0.000001;
        //onlyHeuristic = 0;
    }
    else
    {
        feasible = ga(bin, route, &vetItems, numItems);
            //ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
            //                                        50, &route);

        if(!dualFeasible && feasible)
        {
            std::printf("Error, packing is dual infeasible!\n");
            PRINT_THROW();
        }
    }

    //std::printf("Construtivo: %d\n", feasible);

    //std::printf("**************END CONSTRUTIVO************\n");
    //std::printf("*****************************************\n\n");

    numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);

    if(feasible)
    {

        if(!checkIfPackedAllTheItems(bin, vetItems, numItems))
        {
            PRINT_THROW();
        }

        //std::printf("Plot: \n%s\n", bin.printPlot().c_str());

        if(useDualFunction && doBreakTestRoute)
        {
            std::printf("Heuristic is feasible ??\n");
            goto jmpCp;
            PRINT_THROW();
        }

        if(useHash)
            routeSetFeasible.insert(route);

        if(!doInverseRoute)
            return 1;
        else
        //(doInverseRoute)
        {
            bool inverse = testRoute(&routeInverse.vetRota[0], routeInverse.numPos, 0, 0);
            if(inverse)
                return 2;
            else
                return 1;
        }
    }

    //if(onlyHeuristic >= 1)
    //    return 0;

    jmpCp:

    std::vector<Cuboid>   vetCuboids;
    Collections::IdVector stopIds;

    convertVectorOfItensToVectorOfCuboids(
        vetItems, vetCuboids, numItems, route);
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
    static Container    container((int)instanciaG.vetDimVeiculo[0],
                        (int)instanciaG.vetDimVeiculo[1],
                        (int)instanciaG.vetDimVeiculo[2],
                        (int)instanciaG.maxPayload);

           //
    PackingType    lastType;
    StatusOroloc3D statusOroc3D;
    double         tempoCpu;
    int totalTry = 0;

    std::printf("\n\n************************************\n");
    std::printf(    "**************INI CP-SAT************\n\n");

    double int_fk, int_fFA, int_fRA, int_fTA;

    for(int i=0; i < 1; ++i)
    {
        bin.reset();

        std::vector<Array<int, 4>> vetArray;
        // std::cout<<"n: "<<n<<"\n";
        //double ompStart = omp_get_wtime();
        auto status = loadingChecker.ConstraintProgrammingSolver(
        PackingType::Complete, container, stopIds, vetCuboids, input.cpSatTime, vetArray,
            int_fk, int_fFA, int_fRA, int_fTA);
        // std::cout<<"ret\n";
        //double ompEnd = omp_get_wtime();

        std::printf("************************************\n");
        std::printf("**************END CP-SAT************\n\n");

        if(status == LoadingStatus::Infeasible)
        {
            if(feasible)
            {
                std::printf("Error, heuristic is feasible\n");
                PRINT_THROW();
            }

            if(useHash)
                routeSetInfeasible.insert(route);

            /*if(feasibleRev)
            {
                std::printf("A rota inversa eh viavel!\n");
                PRINT_THROW();
            }
            */
            if(doInverseRoute)
            {
                bool inverse = testRoute(&routeInverse.vetRota[0], routeInverse.numPos,
                                         0, 0);
                if(inverse)
                    return 3;
                else
                    return 0;

            }

            if(useDualFunction && doBreakTestRoute)
            {
                std::printf("\nInfeasible\n");
                PRINT_THROW();
            }


            return 0;
        }
        else if(status == LoadingStatus::Invalid)
        {
            std::printf("ERROR, status: Invalid\n");
            PRINT_THROW();
        }
        else if(status == LoadingStatus::Unknown)
        {
            //std::printf("Status: Unknown, Time Limit?");
            return 0;
            PRINT_THROW();
            continue;
        }
        else if(status == LoadingStatus::FeasOpt)
        {
            int item = 0;
            for(Array<int, 4> &array : vetArray)
            {
                bin.vetItemId[item] = vetItems[item];
                bin.vetPosItem[item].set(array[0], array[1], array[2]);
                bin.vetRotacao[item] = (InstanceNS::Rotation)array[3];

                item += 1;
            }

            bin.numItens = numItems;
            bin.computeLoadingBalancing();

            if(!bin.checkFeasibility(&route, true, true))
            {
                std::printf("ERROR, feasible solution from CP model is not feasible!\n");
                std::cout<<bin.printPlot()<<"\n";
                PRINT_THROW();
                continue;
            }

            if(useHash)
                routeSetFeasible.insert(route);

            if(!checkIfPackedAllTheItems(bin, vetItems, numItems))
            {
                PRINT_THROW();
            }

            if(doInverseRoute)
            {
                bool inverse = testRoute(&routeInverse.vetRota[0], routeInverse.numPos,
                                         0, 0);
                if(inverse)
                    return 2;
                else
                    return 1;

            }


            if(useDualFunction && doBreakTestRoute)
            {
                std::printf("\nFeasOpt\n");
                PRINT_THROW();
            }


            return 1;
        }

    }

    //std::printf("ERROR, STATUS FROM OR TOOLS IS TIME LIMIT\n");
    //PRINT_THROW();

    return false;
}

int getNumberOfCustoms()
{

    return instanciaG.numClientes;
}

int getNumberOfTrucks()
{
    return instanciaG.numVeiculos;
}

int getDemandFromCustomr(int i)
{
    return instanciaG.vetDemandaCliente[i];
}

double getDistance(int i, int j)
{
    if(input.instOroloc3D_2)
        return instanciaG.matDist(i, j) * 0.000001;
    else
        return instanciaG.matDist(i, j);
}

int getVehicleCapacity()
{
    return instanciaG.maxPayload;
}

double getVolumeFromCustomr(int cust)
{

    if(cust == 0)
        return 0.0;

    int ini = instanciaG.matCliItensIniFim(cust, 0);
    int end = instanciaG.matCliItensIniFim(cust, 1);

    double vol = 0.0;

    if(input.instOroloc3D_2)
    {
        for(int i=ini; i <= end; ++i)
        {
            const Item& item = instanciaG.vetItens[i];
            vol += (item.vetDim[0])*(item.vetDim[1])*(item.vetDim[2]);
            //std::printf("Item: %d; volume: %f\n", i, instanciaG.vetItens[i].volume);
        }
    }
    else
    {
        for(int i=ini; i <= end; ++i)
        {
            const Item& item = instanciaG.vetItens[i];
            vol += item.vetDim[0]*item.vetDim[1]*item.vetDim[2];
            //std::printf("Item: %d; volume: %f\n", i, instanciaG.vetItens[i].volume);
        }
    }

    return vol;
}

double getVehicleVolume()
{

    if(input.instOroloc3D_2)
    {    return instanciaG.vetDimVeiculo[0]*0.001*instanciaG.vetDimVeiculo[1]*0.001*
                instanciaG.vetDimVeiculo[2]*0.001;
    }
    else
    {
        return instanciaG.vetDimVeiculo[0]*instanciaG.vetDimVeiculo[1]*
               instanciaG.vetDimVeiculo[2];
    }

}

void setClassical3DPackingProblem(int doPrint)
{
    if(doPrint)
        std::printf("Seting parameters for classical 3D loading; NO LIFO\n");

    input.instOroloc3D_2  		= false;
    input.axleWights      		= false;
    input.balancedLoading 		= false;
    input.compactness     		= false;
    input.lifo			  		= false;
    input.mlifo			  		= false;
    input.removeFromShortSide	= true;
    input.fragility				= false;
    input.support				= false;

    if(input.mlifo && !input.lifo)
    {
        std::printf("Error, lifo needs to be active when mlifo is!\n");
        PRINT_THROW();
    }

    std::printf("Solving Problem:\n\taxleWights: \t\t %d \n", input.axleWights);
    std::printf("\tbalancedLoading: \t %d\n", input.balancedLoading);
    std::printf("\tcompactness: \t\t %d\n", input.compactness);
    std::printf("\tlifo: \t\t\t %d\n", input.lifo);
    std::printf("\tmlifo: \t\t\t %d\n", input.mlifo);
}

void setLoadingOnlyProblem(int doPrint)
{
    if(doPrint)
        std::printf("Seting parameters for classical loading only\n");

    input.instOroloc3D_2  		= false;
    input.axleWights      		= false;
    input.balancedLoading 		= false;
    input.compactness     		= false;
    input.lifo			  		= false;
    input.mlifo			  		= false;
    input.removeFromShortSide	= true;
    input.fragility				= false;
    input.support				= false;

    if(input.mlifo && !input.lifo)
    {
        std::printf("Error, lifo needs to be active when mlifo is!\n");
        PRINT_THROW();
    }

    std::printf("Solving Problem:\n\taxleWights: \t\t %d \n", input.axleWights);
    std::printf("\tbalancedLoading: \t %d\n", input.balancedLoading);
    std::printf("\tcompactness: \t\t %d\n", input.compactness);
    std::printf("\tlifo: \t\t\t %d\n", input.lifo);
    std::printf("\tmlifo: \t\t\t %d\n", input.mlifo);
}

void setOroloc3DProblem(int doPrint)
{
    if(doPrint)
        std::printf("Seting parameters for Oroloc3D loading\n");



    input.instOroloc3D_2  		= true;
    input.axleWights      		= true;
    input.balancedLoading 		= true;
    input.compactness     		= true;
    input.lifo			  		= true;
    input.mlifo			  		= true;
    input.removeFromShortSide	= false;

    std::printf("Solving Problem:\n\taxleWights: \t\t %d \n", input.axleWights);
    std::printf("\tbalancedLoading: \t %d\n", input.balancedLoading);
    std::printf("\tcompactness: \t\t %d\n", input.compactness);
    std::printf("\tlifo: \t\t\t %d\n", input.lifo);
    std::printf("\tmlifo: \t\t\t %d\n", input.mlifo);

    if(input.mlifo && !input.lifo)
    {
        std::printf("Error, lifo needs to be active when mlifo is!\n");
        PRINT_THROW();
    }
}

void saveProblem(Input& inputTemp)
{
    inputTemp.instOroloc3D_2 		= input.instOroloc3D_2;
    inputTemp.axleWights     		= input.axleWights;
    inputTemp.balancedLoading		= input.balancedLoading;
    inputTemp.compactness			= input.compactness;
    inputTemp.lifo					= input.lifo;
    inputTemp.mlifo					= input.mlifo;
    inputTemp.removeFromShortSide	= input.removeFromShortSide;
}

void setProblem(Input& inputTemp)
{

    input.instOroloc3D_2 		= inputTemp.instOroloc3D_2;
    input.axleWights     		= inputTemp.axleWights;
    input.balancedLoading		= inputTemp.balancedLoading;
    input.compactness			= inputTemp.compactness;
    input.lifo					= inputTemp.lifo;
    input.mlifo					= inputTemp.mlifo;
    input.removeFromShortSide	= inputTemp.removeFromShortSide;
}

void setDualFeasibleFunction(double ep0, double ep2)
{
    instanciaG.vetVolNormDualCust.setAll(0.0);

    for(int i=0; i < instanciaG.numItens; ++i)
    {
        Item& item = instanciaG.vetItens[i];
        double d0 = funcionU(ep0, item.vetDimNor[0]);
        double d1 = funcionU(ep0, item.vetDimNor[1]);
        double d2 = funcionU(ep2, item.vetDimNor[2]);
        item.volNormDual = d0*d1*d2;

        instanciaG.vetVolNormDualCust[item.customer] += item.volNormDual;
    }

    /*
    for(int i=0; i < instanciaG.numClientes; ++i)
    {
        std::printf("%d: %.2f\n", i, instanciaG.vetVolNormDualCust[i]);
    }

    std::printf("\n\n");
    */
}

double getDualVolume(int cust)
{
    return instanciaG.vetVolNormDualCust[cust];
}

double getDualVolumeTotal()
{

    return instanciaG.volNormal;
}
