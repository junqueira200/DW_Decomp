#include "TesteOroloc3D.h"
#include "InputOutput.h"
#include <fstream>
#include <omp.h>
#include "MILP.h"
#include "Solucao.h"
#include "AxleWeights.h"
#include "ConstrutivoBin.h"

using namespace TesteOroloc3D_NS;
using namespace InstanceNS;
using namespace ParseInputNS;

using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace MILP_NS;
using namespace SolucaoNS;
using namespace AxleWeightsNS;
using namespace ConstrutivoBinNS;


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

    std::multimap<int, PontoRot> mapItem_idToPontRot;

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

    Vector<PontoRot> vetPontoRot(instanciaG.numItens);
    //truckId = 0;

    while(true)
    {
        fileItems>>itemId>>orderId>>truckId>>unload>>px>>py>>pz>>rot;
        if(truckId == -1)
        {
            break;
        }

        vetPontoRot[instanciaG.mapItem_IdItem[itemId]].set(px, py, pz, mapRotation_axisToRot[rot]);

        mapItem_idToPontRot.insert({itemId, PontoRot(px, py, pz, mapRotation_axisToRot[rot])});

        //std::printf("%d %d %d %d\n",truckId, px, py, pz);
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
        //std::printf("%d %d\n", orderId, truckId);

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

    Vector<VectorI> vetVetItems(instanciaG.numVeiculos);
    for(int i=0; i < instanciaG.numVeiculos; ++i)
    {
        vetBin[i].numItens = 0;
    }

    truckId = 0;
    for(auto& vet:vetOrders)
    {
        //std::cout<<"Truck "<<truckId<<": "<<vet<<"\n";
        VectorI vetRoute;
        for(int order:vet)
        {
            const VectorI& vetItens = instanciaG.mapOrderIdItem[order];
            for(int itemId:vetItens)
                vetVetItems[truckId].push_back(itemId);

            int cust = instanciaG.mapOrderIdCust[order];
            if(!vetRoute.empty())
            {
                if(vetRoute[vetRoute.size()-1] != cust)
                    vetRoute.push_back(cust);
            }
            else
                vetRoute.push_back(cust);
        }

        //std::cout<<"Truck "<<truckId<<": "<<vetVetItems[truckId]<<"\n";
        for(int i=0; i < vetVetItems[truckId].size(); ++i)
        {
            int item = vetVetItems[truckId][i];
            PontoRot& pontoRot = vetPontoRot[item];

            int& next = vetBin[truckId].numItens;

            vetBin[truckId].vetItemId[next] = item;
            vetBin[truckId].vetPosItem[next].vetDim = pontoRot.ponto.vetDim;
            vetBin[truckId].vetRotacao[next] = pontoRot.rot;

            next += 1;

        }

        vetVetRoutes.push_back(vetRoute);

        truckId += 1;
    }

    //std::cout<<"*************\n";


    truckId = 0;

    auto vetPackingType = {PackingType::Complete, PackingType::Complete};//, PackingType::NoSupport, PackingType::LoadingOnly};
    std::map<PackingType, std::string> mapPackingTypeToString;

    mapPackingTypeToString[PackingType::Complete] = "Complete";
    mapPackingTypeToString[PackingType::NoSupport] = "NoSupport";
    mapPackingTypeToString[PackingType::LoadingOnly] = "LoadingOnly";
//    mapPackingTypeToString[PackingType::NoSupportNoFragility] = "NoSupport";

    std::map<StatusOroloc3D, std::string> mapStatusOroloc3D_ToString;
    mapStatusOroloc3D_ToString[INFEASIBLE] = "INFEASIBLE";
    mapStatusOroloc3D_ToString[TIME_LIMIT] = "TIME_LIMIT";
    mapStatusOroloc3D_ToString[FEASIBLE]   = "FEASIBLE";

    std::printf("\n*****************\n\n");

    GRBEnv env;
    //exit(-1);
    int numCompleteFeasible = 0;
    for(auto& vet:vetVetRoutes)
    {

        //if(truckId == 2)
        {
        std::string output = std::format("{}; {}; ",  input.strInst , truckId);
        //std::cout<<"Truck "<<truckId<<": "<<vet<<"; \n";
        //VectorI vetItems;
        //copiaItensClientes(vet, vet.size(), vetItems, true);

        /*
        for(int i=0; i < vetItems.size(); ++i)
        {
            Item& item = instanciaG.vetItens[vetItems[i]];
            std::printf("%s\n", item.print().c_str());
        }
        */

        //PRINT_DEBUGG("", "");
        //exit(-1);

        Bin& bin = vetBin[truckId];
        //bin.vetItemId = vetVetItems[truckId];
        //bin.numItens  = (int)vetVetItems[truckId].size();

        std::printf("%d ", truckId);


        //semiTrailer.checkAxleWeights(bin);
        //truckId += 1;
        //continue;


        //PRINT_DEBUGG("", "");
        //exit(-1);

        //truckId += 1;
        //continue;

        //std::cout<<vetItems<<"\n";


        std::vector<Cuboid> vetCuboids;
        Collections::IdVector stopIds;

        convertVectorOfItensToVectorOfCuboids(bin.vetItemId, vetCuboids, bin.numItens);
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

        bool inputAxleWights = input.axleWights;

        int n = 0;
        for(PackingType type:vetPackingType)
        {

            double ompStart = omp_get_wtime();
            //PackingType::LoadingOnly
            if(n == 0)
                input.axleWights = true;
            else
                input.axleWights = false;

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

                if(n == 0)
                    numCompleteFeasible += 1;
            }

            output += std::format("{:.1f}; ", tempoCpu);

            if(doBreak)
                break;

            n += 1;
        }

        input.axleWights = inputAxleWights;
        if(n == 0)
            output += "NO_RUN; Inf";
        /*
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
        */

        //output += std::format("{:.1f}; {}; {}", tempoCpu, mapPackingTypeToString[lastType], mapStatusOroloc3D_ToString[statusOroc3D]);
        std::cout<<output<<"\n";
        appendToFile("../oroloc3D.csv", output);


        }


        truckId += 1;
//        break;



    }


    std::printf("numCompleteFeasible: %d\n", numCompleteFeasible);

    file.close();
    fileItems.close();

}


void TesteOroloc3D_NS::convertVectorOfItensToVectorOfCuboids(const VectorI& vetItens, std::vector<Cuboid>& vetCuboids,
                                                             int numItems)
{
    vetCuboids = std::vector<Cuboid>(numItems);

    int pos = numItems;

    for(int i=0; i < numItems; ++i)
    {
        Cuboid& cuboid = vetCuboids[i];
        Item& item = instanciaG.vetItens[vetItens[i]];

        if(i > 0)
        {
            Item& item_1 = instanciaG.vetItens[vetItens[i-1]];
            if(item.customer != item_1.customer)
                pos -= 1;


        }

        //std::cout<<pos<<" ";

        cuboid = Cuboid((size_t)i, (size_t)vetItens[i], item.vetDim[0], item.vetDim[1], item.vetDim[2], true,
                        Fragility::None, pos, item.weight);

    }

    //std::cout<<"\n";

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

void TesteOroloc3D_NS::testeOroloc3D_2()
{

    Solucao sol;
    readSolOroloc3D_2(sol);

    auto vetPackingType = {PackingType::Complete, PackingType::Complete};//, PackingType::NoSupport, PackingType::LoadingOnly};

    std::map<PackingType, std::string> mapPackingTypeToString;
    mapPackingTypeToString[PackingType::Complete] = "Complete";
    mapPackingTypeToString[PackingType::NoSupport] = "NoSupport";
    mapPackingTypeToString[PackingType::LoadingOnly] = "LoadingOnly";

    std::map<StatusOroloc3D, std::string> mapStatusOroloc3D_ToString;
    mapStatusOroloc3D_ToString[INFEASIBLE] = "INFEASIBLE";
    mapStatusOroloc3D_ToString[TIME_LIMIT] = "TIME_LIMIT";
    mapStatusOroloc3D_ToString[FEASIBLE]   = "FEASIBLE";

    int numCompleteFeasible = 0;

    std::cout<<"Veic: "<<instanciaG.vetDimVeiculo<<"\n";

    for(int veic=0; veic < sol.vetBin.size(); ++veic)
    {
        std::string output = std::format("{}; {}; ",  input.strInst , veic);


        Bin& bin = sol.vetBin[veic];
        Bin bin2;
        if(bin.vazio())
            continue;

        VectorI vetItems = bin.vetItemId;
        std::reverse(vetItems.begin(), vetItems.begin()+bin.numItens);

        bool feasibleSol = construtivoBinPacking(bin2, vetItems, bin.numItens, input.aphaBin, 10);
        std::printf("Extreme Point Heuristic: %d\n", feasibleSol);


        bool axleWeights = semiTrailer.checkAxleWeights(bin);

        std::printf("Veic %d:\n", veic);
        std::printf("AxleWeights: %d\n", (int)axleWeights);

        bool feasible = bin.verificaViabilidade();
        std::printf("Bin check: %d\n", feasible);

        std::vector<Cuboid> vetCuboids;
        Collections::IdVector stopIds;

        convertVectorOfItensToVectorOfCuboids(bin.vetItemId, vetCuboids, bin.numItens);
        for(int i=1; i < (sol.vetRota[veic].numPos-1); ++i)
            stopIds.push_back(sol.vetRota[veic].vetRota[i]);


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
        GRBEnv env;
        GRBModel model(env);
        model.set(GRB_IntParam_Threads, 4);
        model.set(GRB_IntParam_SolutionLimit, 1);
        model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_DoubleParam_TimeLimit, 60*10);

        //model.set(GRB_DoubleParam_Heuristics, 1);
        //model.set(GRB_IntParam_PumpPasses, 10000);
        //model.set(GRB_IntParam_ZeroObjNodes, 100000);

        Variables variables(model, bin.vetItemId, bin.numItens);
        //Bin bin;

        //bin.numItens = vetItems.size();
        //bin.vetItemId = vetItems;

        addBasicConstraints(model, variables, bin);
        std::printf("Truck %d: \n\n", veic);
        model.optimize();

        continue;
        */

        bool inputAxleWights = input.axleWights;

        int n = 0;
        for(PackingType type:vetPackingType)
        {

            double ompStart = omp_get_wtime();
            //PackingType::LoadingOnly
            if(n == 0)
                input.axleWights = true;
            else
                input.axleWights = false;

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

                if(n == 0)
                    numCompleteFeasible += 1;
            }

            output += std::format("{:.1f}; ", tempoCpu);

            if(doBreak)
                break;

            n += 1;
        }

        input.axleWights = inputAxleWights;
        if(n == 0)
            output += "NO_RUN; Inf";
        /*
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
        */

        //output += std::format("{:.1f}; {}; {}", tempoCpu, mapPackingTypeToString[lastType], mapStatusOroloc3D_ToString[statusOroc3D]);
        std::cout<<output<<"\n";
        appendToFile("../oroloc3D.csv", output);

    }

}

void TesteOroloc3D_NS::readSolOroloc3D_2(SolucaoNS::Solucao& sol)
{

    std::printf("FILE: %s\n", input.strSolOroloc3D_2.c_str());

    std::ifstream file(input.strSolOroloc3D_2);
    assertm(!file.is_open(), "Cant open the file: : "<<input.instOroloc3D_2);

    std::printf("Before reset\n");
    sol.reset();
    std::printf("After reset\n");

    std::string str;
    int numItems, numCust, numVeic, veicId;


    std::map<int, InstanceNS::Rotation> mapRotation_axisToRot;
    mapRotation_axisToRot[0]  = InstanceNS::Rot0;
    mapRotation_axisToRot[3]  = InstanceNS::Rot2;
    mapRotation_axisToRot[1]  = InstanceNS::Rot1;


    std::getline(file, str);    // Name:                   	S10_240T
    std::getline(file, str);    // Problem:                	3L-CVRP

    file>>str>>numVeic;
    file>>str>>sol.distTotal;

    std::printf("numVeic: %d; distTotal: %f\n", numVeic, sol.distTotal);

    for(int i=0; i < 6; ++i)
    {
        std::getline(file, str);
        std::printf("LINE: %s\n", str.c_str());
    }

    for(int veic=0; veic < numVeic; ++veic)
    {

        std::getline(file, str);  // Tour_Id:
        file>>str>>numCust;
        file>>str>>numItems;
        file>>str;               // Customer_Sequence:

        sol.vetRota[veic].vetRota[0] = 0;

        for(int i=0; i < numCust; ++i)
            file>>sol.vetRota[veic].vetRota[i+1];

        sol.vetRota[veic].vetRota[numCust+1] = 0;
        sol.vetRota[veic].numPos = numCust + 2;
        sol.vetRota[veic].computeDistance();

        std::printf("Rota: %s\n", sol.vetRota[veic].printRota().c_str());

        std::getline(file, str);
        std::getline(file, str);
        std::getline(file, str);  // CustId	Id	TypeId	Rotated	x	y	z	Length	Width	Height	mass	Fragilit

        std::printf("LINE: %s\n", str.c_str());

        if(str[0] == '-')
            continue;

        int custId, id, typeId, rot;
        Array<int, 3> pos;
        Bin& bin = sol.vetBin[veic];
        bin.numItens = numItems;

        for(int i=0; i < numItems; ++i)
        {
            file>>custId>>id>>typeId>>rot>>pos[0]>>pos[1]>>pos[2];
            bin.vetPosItem[i].set(pos[0], pos[1], pos[2]);
            bin.vetRotacao[i] = mapRotation_axisToRot[rot];
            bin.vetItemId[i]  = instanciaG.mapItem_IdItem[typeId];

            std::printf("Rot: %d; x(%d), y(%d), z(%d)\n", (int)bin.vetRotacao[i], pos[0], pos[1], pos[2]);

            std::getline(file, str);


        }

        std::getline(file, str);  // new line
        std::getline(file, str);  // new line
        std::getline(file, str);  //-----------------------------

        std::printf("LINE: %s\n", str.c_str());

    }

    file.close();


}
