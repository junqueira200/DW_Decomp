#include "TesteOroloc3D.h"
#include "InputOutput.h"
#include <fstream>
#include <omp.h>
#include "MILP.h"
#include "Solucao.h"
#include "AxleWeights.h"

using namespace TesteOroloc3D_NS;
using namespace InstanceNS;
using namespace ParseInputNS;

using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace MILP_NS;
using namespace SolucaoNS;
using namespace AxleWeightsNS;

void TesteOroloc3D_NS::testeOroloc3D()
{
    //std::printf("testeOroloc3D\n");
    //std::printf("%s\n", ParseInputNS::input.strInstCompleto.c_str());
     //= ParseInputNS::input.strInstCompleto.replace(".oroloc3D", ".csv")
    std::string strFileSol = ParseInputNS::input.strInstCompleto;
    std::string strFileSolItems = ParseInputNS::input.strInstCompleto;

    strFileSol.replace(strFileSol.find_last_of("."), std::string::npos, "_sol_orders.csv");
    strFileSolItems.replace(strFileSolItems.find_last_of("."), std::string::npos, "_sol_items.csv");

    std::printf("Reading %s\n", strFileSol.c_str());
    std::printf("Reading %s\n", strFileSolItems.c_str());


    std::ifstream file(strFileSol), fileItems(strFileSolItems);
    assertm(!file.is_open(), "Cant open the file: : "<<strFileSol);
    assertm(!fileItems.is_open(), "Cant open the file: : "<<strFileSolItems);

    std::string trash;
    std::getline(file, trash);
    std::getline(fileItems, trash);

    Vector<VectorI> vetOrders;
    Vector<InstanceNS::Rotation> vetRot;
    Vector<VectorI> vetVetRoutes;
    Vector<Bin> vetBin(instanciaG.numVeiculos);


    VectorI vetAux;
    int truckIdOld  = 1;
    int itemId      = -1;
    int truckId     = -1;
    int orderId     = -1;
    int unloadOrder = -1;
    int rot         = -1;
    int unload      = -1;
    Array<int, 7> vetTrash;

    std::map<int, InstanceNS::Rotation> mapRotation_axisToRot;
    mapRotation_axisToRot[-1] = InstanceNS::Rot0;
    mapRotation_axisToRot[2]  = InstanceNS::Rot1;
    mapRotation_axisToRot[0]  = InstanceNS::Rot2;

    int px, py, pz;
    int nextItemId = 0;


    /*
    for(Bin& bin:vetBin)
    {
        bin.vetItemId = VectorI();
        bin.vetPosItem = Vector<Ponto>();
        bin.numItens = 0;
    }
    */

    while(true)
    {
        fileItems>>itemId>>orderId>>truckId>>unload>>px>>py>>pz>>rot;
        if(truckId == -1)
            break;

        //std::printf("%d %d %d\n", px, py, pz);
        Bin& bin = vetBin[truckId-1];
        int nextItem = bin.numItens;
        bin.vetPosItem[nextItem] = Ponto(px, py, pz);
        bin.vetRotacao[nextItem] = mapRotation_axisToRot[rot];
        bin.numItens += 1;
    }



    //std::cout<<"map size: "<<instanciaG.mapOrderIdItem.size()<<"\n";
    for(int i=0; i < instanciaG.mapOrderIdItem.size(); ++i)
    {
        file>>orderId>>truckId>>unloadOrder;

        /*
        for(int j=0; j < 4; ++j)
        {
            fileItems>>vetTrash[j];
            //std::printf("%d, ", vetTrash[j]);
        }
        */

        //fileItems>>px>>py>>pz>>rot;
        //std::printf("%d %d %d\n", px, py, pz);


        //vetRot.push_back(mapRotation_axisToRot[rot]);
        //std::printf("\n%d\n", (int)mapRotation_axisToRot[rot]);

        if(truckId == truckIdOld)
            vetAux.push_back(orderId);
        else
        {


            //std::cout<<vetAux<<"\n";
            vetOrders.push_back(vetAux);

            if(truckId == -1)
                break;

            vetAux = VectorI();
            vetAux.push_back(orderId);
            truckIdOld = truckId;
        }
    }


    truckId = 0;
    for(auto& vet:vetOrders)
    {
        //std::cout<<"Truck "<<truckId<<": "<<vet<<"\n";
        VectorI vetRoute;
        for(int order:vet)
        {
            int cust = instanciaG.mapOrderIdCust[order];
            if(!vetRoute.empty())
            {
                if(vetRoute[vetRoute.size()-1] != cust)
                    vetRoute.push_back(cust);
            }
            else
                vetRoute.push_back(cust);
        }

        vetVetRoutes.push_back(vetRoute);

        truckId += 1;
    }

    //std::cout<<"*************\n";


    truckId = 0;

    auto vetPackingType = {PackingType::Complete};//, PackingType::NoSupport, PackingType::LoadingOnly};
    std::map<PackingType, std::string> mapPackingTypeToString;

    mapPackingTypeToString[PackingType::Complete] = "Complete";
    mapPackingTypeToString[PackingType::NoSupport] = "NoSupport";
    mapPackingTypeToString[PackingType::LoadingOnly] = "LoadingOnly";
//    mapPackingTypeToString[PackingType::NoSupportNoFragility] = "NoSupport";

    std::map<StatusOroloc3D, std::string> mapStatusOroloc3D_ToString;
    mapStatusOroloc3D_ToString[INFEASIBLE] = "INFEASIBLE";
    mapStatusOroloc3D_ToString[TIME_LIMIT] = "TIME_LIMIT";
    mapStatusOroloc3D_ToString[FEASIBLE]   = "FEASIBLE";

    GRBEnv env;
    //exit(-1);
    for(auto& vet:vetVetRoutes)
    {

        //if(truckId == 2)
        {
        std::string output = std::format("{}; {}; ",  input.strInst , truckId);
        //std::cout<<"Truck "<<truckId<<": "<<vet<<"; \n";
        VectorI vetItems;
        copiaItensClientes(vet, vet.size(), vetItems, true);

        Bin& bin = vetBin[truckId];
        bin.vetItemId = vetItems;
        bin.numItens  = (int)vetItems.size();

        std::printf("Truck %d: ", truckId);
        //semiTrailer.checkAxleWeights(bin);

        //std::cout<<vetItems<<"\n";


        std::vector<Cuboid> vetCuboids;
        Collections::IdVector stopIds;

        convertVectorOfItensToVectorOfCuboids(vetItems, vetCuboids);
        for(int cust:vet)
            stopIds.push_back(cust);

        InputParameters inputParam;
        inputParam.ContainerLoading.LoadingProblem.Variant = LoadingProblemParams::VariantType::AllConstraints;
        inputParam.SetLoadingFlags();

        LoadingChecker loadingChecker(inputParam.ContainerLoading);
        Container container((int)instanciaG.vetDimVeiculo[0], (int)instanciaG.vetDimVeiculo[1],
                            (int)instanciaG.vetDimVeiculo[2], (int)instanciaG.maxPayload);

        //
        PackingType lastType;
        StatusOroloc3D statusOroc3D;
        double tempoCpu;


        /*
        GRBModel model(env);
        model.set(GRB_IntParam_Threads, 8);
        model.set(GRB_IntParam_SolutionLimit, 1);
        model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_DoubleParam_TimeLimit, 60*10);
        model.set(GRB_DoubleParam_Heuristics, 1);
        model.set(GRB_IntParam_PumpPasses, 10000);
        model.set(GRB_IntParam_ZeroObjNodes, 100000);

        Variables variables(model, vetItems, vetItems.size());
        Bin bin;

        bin.numItens = vetItems.size();
        bin.vetItemId = vetItems;

        addBasicConstraints(model, variables, bin);
        std::printf("Truck %d: \n\n", truckId);
        model.optimize();

        */


        for(PackingType type:vetPackingType)
        {

            double ompStart = omp_get_wtime();
            //PackingType::LoadingOnly

            auto status = loadingChecker.ConstraintProgrammingSolver(type, container, stopIds, vetCuboids, 10*60);

            double ompEnd = omp_get_wtime();

            tempoCpu = ompEnd-ompStart;
            std::printf("\t%2.f S; ", tempoCpu);

            lastType = type;

            bool doBreak = false;

            if(status != LoadingStatus::FeasOpt)
            {
                //std::printf("Error in ConstraintProgrammingSolver: ");
                if(status == LoadingStatus::Infeasible)
                {
                    std::printf("INFEASIBLE\n");
                    statusOroc3D = INFEASIBLE;
                    output += "INFEASIBLE; ";
                    //break;
                }
                else
                {
                    std::printf("Time Limit\n");
                    statusOroc3D = TIME_LIMIT;
                    output += "TIME_LIMIT; ";
                }
            }
            else
            {
                std::printf("Loading was successful!\n");
                statusOroc3D = FEASIBLE;
                output += "FEASIBLE; ";
                doBreak = true;
            }

            output += std::format("{:.1f}; ", tempoCpu);

            if(doBreak)
                break;
        }

        switch (lastType)
        {
        case PackingType::Complete:
            output += "NO_RUN; Inf; NO_RUN; Inf";
            break;

        case PackingType::NoSupport:
            output += "NO_RUN; Inf";
        default:
            break;
        }

        //output += std::format("{:.1f}; {}; {}", tempoCpu, mapPackingTypeToString[lastType], mapStatusOroloc3D_ToString[statusOroc3D]);
        std::cout<<output<<"\n";
        appendToFile("../oroloc3D.csv", output);


        }


        truckId += 1;
//        break;



    }




    file.close();
    fileItems.close();

}


void TesteOroloc3D_NS::convertVectorOfItensToVectorOfCuboids(const VectorI& vetItens, std::vector<Cuboid>& vetCuboids)
{
    vetCuboids = std::vector<Cuboid>(vetItens.size());

    int pos = 0;

    for(int i=0; i < (int)vetItens.size(); ++i)
    {
        Cuboid& cuboid = vetCuboids[i];
        Item& item = instanciaG.vetItens[vetItens[i]];

        if(i > 0)
        {
            Item& item_1 = instanciaG.vetItens[vetItens[i]];
            if(item.customer != item_1.customer)
                pos += 1;


        }

        cuboid = Cuboid((size_t)i, (size_t)vetItens[i], item.vetDim[0], item.vetDim[1], item.vetDim[2], true,
                        Fragility::None, pos, item.weight);

    }
}


void TesteOroloc3D_NS::appendToFile(const std::string& fileName, const std::string& content)
{
    std::ofstream file(fileName, std::ios::app);  // Append

    if (!file.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << fileName << std::endl;
        PRINT_DEBUGG("", "");
        exit(-1);
    }

    file<<content<<"\n";

    file.close();
}

