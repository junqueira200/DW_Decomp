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

void ini_3D_Packing(char *strInst_c, int oroloc3D)
{

    //static bool doInit = true;
    //std::printf("Inst: %s\n", strInst_c);
    //if(doInit)

    //std::printf("oroloc3D: %d\n", oroloc3D);
    //PRINT_THROW()

    if(oroloc3D == 0)
        setClassical3DPackingProblem();
    else
        setOroloc3DProblem();

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

int testRoute(int *vet_c, int vetSize, int onlyHeuristic, int doInverseRoute)
{
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
        route.vetRota[i] = vet_c[i];

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

    int numItems = copiaItensClientes(route.vetRota, route.numPos, vetItems);

    double totalVolume = 0.0;
    int totalDemand    = 0;
    for(int i=0; i < numItems; ++i)
    {
        totalVolume += instanciaG.vetItens[i].volume;
        //totalDemand += instanciaG.vetItens[i].weight;
    }

    for(int i=0; i < vetSize; ++i)
        totalDemand += instanciaG.vetDemandaCliente[route.vetRota[i]];

    static double volumeOfVeicule = getVehicleVolume();
    if(totalVolume > volumeOfVeicule)
    {
        //std::printf("maxPayload\n");
        return 0;
    }

    if(totalDemand > (int)instanciaG.maxPayload)
    {
        //std::printf("maxPayload\ntotalDemand: %d; maxPayload: %f", totalDemand,
        //            instanciaG.maxPayload);
        return 0;
    }

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

    if(onlyHeuristic >= 1)
    {
        feasible =
        ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
                                                25, &route);
        if(feasible)
            routeSetFeasible.insert(route);
        //if(feasible)
        //    return true;
        return feasible;

        //input.cpSatTime = 0.000001;
        //onlyHeuristic = 0;
    }
    else
    {
        feasible =
            ConstrutivoBinNS::construtivoBinPacking(bin, vetItems, numItems, input.aphaBin,
                                                std::numeric_limits<int>::max()-1, &route);
    }

    //std::printf("Construtivo: %d\n", feasible);

    //std::printf("**************END CONSTRUTIVO************\n");
    //std::printf("*****************************************\n\n");

    if(feasible)
    {
        routeSetFeasible.insert(route);
        if(!doInverseRoute)
            return 1;
        else
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

    //std::printf("\n\n************************************\n");
    //std::printf(    "**************INI CP-Sat************\n\n");

    int int_fk, int_fFA, int_fRA, int_fTA;

    for(int i=0; i < 1; ++i)
    {
        bin.reset();

        std::vector<Array<int, 4>> vetArray;
        // std::cout<<"n: "<<n<<"\n";
        double ompStart = omp_get_wtime();
        auto status = loadingChecker.ConstraintProgrammingSolver(
        PackingType::Complete, container, stopIds, vetCuboids, input.cpSatTime, vetArray,
            int_fk, int_fFA, int_fRA, int_fTA);
        // std::cout<<"ret\n";
        double ompEnd = omp_get_wtime();

        //std::printf("************************************\n");
        //std::printf("**************END CP-Sat************\n\n");

        if(status == LoadingStatus::Infeasible)
        {
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
                continue;
            }

            routeSetFeasible.insert(route);

            if(doInverseRoute)
            {
                bool inverse = testRoute(&routeInverse.vetRota[0], routeInverse.numPos,
                                         0, 0);
                if(inverse)
                    return 2;
                else
                    return 1;

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
            vol += (item.vetDim[0]*0.001)*(item.vetDim[1]*0.001)*(item.vetDim[2]*0.001);
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

void setClassical3DPackingProblem()
{
    std::printf("Seting parameters for classical 3D loading; NO LIFO\n");
    input.instOroloc3D_2  		= false;
    input.axleWights      		= false;
    input.balancedLoading 		= false;
    input.compactness     		= false;
    input.lifo			  		= false;
    input.mlifo			  		= false;
    input.removeFromShortSide	= true;
    input.fragility				= true;
}

void setOroloc3DProblem()
{
    std::printf("Seting parameters for Oroloc3D loading\n");
    input.instOroloc3D_2  		= true;
    input.axleWights      		= true;
    input.balancedLoading 		= true;
    input.compactness     		= true;
    input.lifo			  		= false;
    input.mlifo			  		= true;
    input.removeFromShortSide	= false;
}
