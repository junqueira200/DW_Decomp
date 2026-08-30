/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: main.cpp
 * ****************************************
 * ****************************************/

#include "AuxT.h"
#include "BinPackingCP.h"
#include "Construtivo.h"
#include "ConstrutivoBin.h"
#include "IBM_CpOptimizer.h"
#include "Ig.h"
#include "InputOutput.h"
#include "Instancia.h"
#include "MILP.h"
#include "TesteOroloc3D.h"
#include "rand.h"
#include "c_api.h"
#include "DualFeasibleFunctions.h"

#include "AuxT.h"
#include "BCRoutingParams.h"
#include "LoadingChecker.h"
#include "ProblemParameters.h"
#include "ortools/base/version.h"

using namespace InstanceNS;
using namespace ConstrutivoBinNS;
using namespace SolucaoNS;
using namespace RandNs;
using namespace ParseInputNS;
using namespace ConstrutivoNS;
using namespace IgNs;
using namespace BinPackingCP_NS;
using namespace TesteOroloc3D_NS;
using namespace DualFeasibleFunctionsNS;

using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace MILP_NS;


void testMLIFO()
{
    std::printf("*****************************************\n");
    std::printf("**************TESTING MLIFO**************\n\n");
    PRINT_THROW();

    static Eigen::Matrix<int, -1, -1, Eigen::RowMajor>
        matSupportItems(instanciaG.numItens, instanciaG.numItens);
    matSupportItems.setConstant(0);

    int index0 = instanciaG.matCliItensIniFim(1, 0);
    int index1 = instanciaG.matCliItensIniFim(2, 0);
    const Item item0 = instanciaG.vetItens[index0];
    const Item item1 = instanciaG.vetItens[index1];

    double d0 = 73.25;
    double d1 = 55.75;
    double d2 = 61.0;

    instanciaG.vetItens[index0].set(d0, d1, d2, instanciaG.vetItens[index0].itemId);
    instanciaG.vetItens[index1].set(d0, d1, d2, instanciaG.vetItens[index1].itemId);

    Bin bin;
    bin.vetItemId[0] = index0;
    bin.vetItemId[1] = index1;

    bin.vetRotacao[0] = InstanceNS::Rot0;
    bin.vetRotacao[1] = InstanceNS::Rot0;
    bin.numItens = 2;

    Rota rota;
    rota.vetRota[1] = 1;
    rota.vetRota[2] = 2;
    rota.vetRota[3] = 0;
    rota.numPos = 4;

    bin.vetPosItem[0].set(326.75, 75.0, 0.0);
    bin.vetPosItem[1].set(400.0, 75.0, 0.0);
    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Pass\n");
    else
        std::printf("Not Pass\n");

    bin.vetPosItem[0].set(400.0, 75.0, 0.0);
    bin.vetPosItem[1].set(326.75, 75.0, 0.0);
    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Pass\n");
    else
        std::printf("Not Pass\n");


    bin.vetPosItem[0].set(400.0, 130.75, 0.0);
    bin.vetPosItem[1].set(400.0, 75.0, 0.0);
    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Pass\n");
    else
        std::printf("Not Pass\n");

    bin.vetPosItem[0].set(400.0, 75.0, 0.0);
    bin.vetPosItem[1].set(400.0, 130.75, 0.0);
    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Not Pass\n");
    else
        std::printf("Pass\n");


    bin.vetPosItem[0].set(400.0, 75.0, 61.0);
    bin.vetPosItem[1].set(400.0, 75.0, 0.0);
    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Pass\n");
    else
        std::printf("Not Pass\n");


    bin.vetPosItem[0].set(400.0, 75.0, 0.0);
    bin.vetPosItem[1].set(400.0, 75.0, 61.0);

    matSupportItems(index0, index1) = true;

    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Not Pass\n");
    else
        std::printf("Pass\n");

    bin.vetPosItem[0].set(400.0, 75.0, 0.0);
    bin.vetPosItem[1].set(400.0, 75.0, 61.0);

    matSupportItems(index0, index1) = false;

    if(SolucaoNS::checkUnloadingSequence(bin, rota, matSupportItems))
        std::printf("Pass\n");
    else
        std::printf("Not Pass\n");




    PRINT_THROW();
}

int main(int argc, const char *argv[])
{
    std::cout << "OR-Tools version: "
              //<< operations_research::OrToolsMajorVersion()
              << operations_research::OrToolsVersionString()
              << std::endl;

    // std::cout<<"main\n";
    // Item item(1, 2, 3, 1);

    // std::cout<<"Rot1: "<<getDim(item, Rot0)<<"\n";
    // return 0;

    ParseInputNS::parseInput(argc, argv);
    output.setup();
    std::cout << "INST: " << input.strInst << " SEED: " << RandNs::estado_ << " "
              << output.data << "";


    if(input.instOroloc3D_2)
        InstanceNS::readOroloc3D2(input.strInstCompleto);
    else if(input.instOroloc3D)
        InstanceNS::readOroloc3D(input.strInstCompleto);
    else if(input.inst2d)
        InstanceNS::read2dInstance(input.strInstCompleto);
    else
        InstanceNS::read3dInstance(input.strInstCompleto);

    //testMLIFO();

    for(int i=0; i < instanciaG.numItens; ++i)
    {
        //Item &item = instanciaG.vetItens[i];
        //for(int d=0; d < 3; ++d)
        //    maxDim = std::max(maxDim, item.vetDim[d]);


        if(i != instanciaG.vetItens[i].itemId)
        {
            std::printf("Error, i(%d) != itemId(%d)\n", i, instanciaG.vetItens[i].itemId);
            PRINT_THROW();
        }

        //std::printf("%d: %s\n", i, instanciaG.vetItens[i].print(InstanceNS::Rot0, true).c_str());
    }

    double vol = 0.0;
    double volNormal = 0.0;
    const double volVeich = instanciaG.vetDimVeiculo[0]*instanciaG.vetDimVeiculo[1]*
                            instanciaG.vetDimVeiculo[2];

    const double maxDim = std::cbrt(volVeich);

    double m0 = 3;
    double m1 = 4;
    double m2 = 1;

    const double volVeichNor = (instanciaG.vetDimVeiculo[0]/maxDim)*(instanciaG.vetDimVeiculo[1]/maxDim)*
                               (instanciaG.vetDimVeiculo[2]/maxDim);

    std::printf("Max dim: %.1f\n", maxDim);
    std::printf("volVeichNor: %.1f\n\n", volVeichNor);


    VectorI vetItems2;

    for(int i=0; i < instanciaG.numItens; ++i)
    {
        vetItems2.push_back(i);
        Item &item = instanciaG.vetItens[i];

        double dim0 = item.vetDim[0]/maxDim;
        double dim1 = item.vetDim[1]/maxDim;
        double dim2 = item.vetDim[2]/maxDim;

        vol += item.vetDim[0]*item.vetDim[1]*item.vetDim[2];
        volNormal += funcionU((1.0/8), dim0) * funcionU((1.0/7), dim1) *
                     funcionU((1.0/9), dim2);
        bool doBreak = false;

        if(vol > volVeich)
        {
            std::printf("VOL MAIOR\n");
            doBreak = true;
        }

        if(volNormal > 1.0)
        {
            std::printf("VOL NORMAL MAIOR\n\n");
            doBreak = true;
        }

        if(doBreak)
            break;
    }

    std::printf("volNormal: %f\n", volNormal);

    if(input.mlifo && !input.lifo)
    {
        std::printf("Error!, mlifo(1), and lifo(0)\n\n");
        PRINT_THROW();
    }


    /*
    int route[] = {0, 2, 12, 10, 0};

    input.instOroloc3D_2  		= false;
    input.axleWights      		= false;
    input.balancedLoading 		= false;
    input.compactness     		= false;
    input.lifo			  		= true;
    input.mlifo			  		= false;
    input.removeFromShortSide	= true;

    int result = testRoute(route, 5);
    std::printf("result: %d\n", result);

    PRINT_THROW();
    */

    startConstGlobalVaribles();
    std::printf("input.instOroloc3D_2: %d\n", input.instOroloc3D_2);

    if(!input.instOroloc3D_2)
    {
        setClassical3DPackingProblem();

        //bool result = testRoute(&vetItems2[0], vetItems2.size(), 0);
        //std::printf("Resultado: %d\n", (int)result);

        //exit(0);


        Solucao sol(instanciaG);
        Solucao best(instanciaG);

        if(metaheuristicaIg(best))
        {
            std::printf("Did find a feasible solution");
            std::printf("Dist; %.1f\nSol:\n\n%s\n\n", best.distTotal, best.printSol().c_str());
        }
        else
        {
            std::printf("Didnt find a solution\n\n");
        }

        return 0;

        best.distTotal = INF_Double;

        bool heuristicSol;
        bool updateBest = false;

        for(int i=0; i < 5000; ++i)
        {

            sol.reset();
            heuristicSol = construtivoVrp(sol, input.alphaVrp, input.aphaBin);
            if(heuristicSol && sol.distTotal < best.distTotal)
            {
                best.reset();
                best.copiaSolucao(sol);
                updateBest = true;
                std::printf("Dist; %.1f\n", best.distTotal);
            }
        }

        if(updateBest)
        {
            std::printf("Did find a feasible solution");
            std::printf("Dist; %.1f\nSol:\n\n%s\n\n", best.distTotal, best.printSol().c_str());


        }
        else
        {
            std::printf("Didnt find a solution\n\n");
        }

        return 0;
    }

    /*
    instanciaG.vetItens[0].set(2380.0, 1414.0, 934.0);   // 133
    Ponto p0(8120.0, 518.0, 1018.0);

    instanciaG.vetItens[1].set(3320.0, 1780.0, 610.0);   // 63
    Ponto p1(4800.0, 335.0, 1580.0);

    //instanciaG.vetItens[2].set(2031.0, 1270.0, 1170.0); // 75
    //Ponto p2(11469.0, 590.0, 0.0);
    InstanceNS::Rotation r = InstanceNS::Rot0;

    double area0 = getIntercetion(1, p1, r, Right, 0, p0, r, Left); // area 75 and 52
    double areaYZ = instanciaG.vetItens[1].vetDim[1]*instanciaG.vetItens[1].vetDim[2];

    std::printf("areaInterc(133 and 63): %.1f\nareaYZ: %.1f\n\n", area0, areaYZ);

    EXIT_PRINT();
    */

    /*
    bool toch = tochRightSideOfTruck(2, p2, r);
    if(toch)
        std::printf("75 tochs the right side\n");
    else
        std::printf("75 dont tochs the right side\n\n");
    */

    // EXIT_PRINT();

    testeOroloc3D_2();
    // IBM_CpOptimizerNS::testSCIP();
    return 0;

    GRBEnv   env;
    GRBModel model(env);
    model.set(GRB_IntParam_Threads, 4);
    model.set(GRB_IntParam_SolutionLimit, 1);

    VectorI vetItems;

    int numItems = generateRandomListOfItems(20, vetItems);

    std::cout << vetItems << "\n";
    for(int i = 0; i < numItems; ++i)
        std::cout << instanciaG.vetItens[vetItems[i]].print() << "\n";

    Variables variables(model, vetItems, numItems);
    Bin       bin;

    bin.numItens = numItems;
    bin.vetItemId = vetItems;

    addBasicConstraints(model, variables, bin);
    model.optimize();

    variables.vetPosX.setVetDoubleAttr_X(model, false);
    for(int i = 0; i < numItems; ++i)
        std::printf("posX[%d] = %.1f\n", i, variables.vetPosX.getX_value(i));

    std::printf("\n\n");
    variables.vetPosY.setVetDoubleAttr_X(model, false);
    for(int i = 0; i < numItems; ++i)
        std::printf("posY[%d] = %.1f\n", i, variables.vetPosY.getX_value(i));

    std::printf("\n\n");
    variables.vetPosZ.setVetDoubleAttr_X(model, false);
    for(int i = 0; i < numItems; ++i)
        std::printf("posZ[%d] = %.1f\n", i, variables.vetPosZ.getX_value(i));

    std::printf("\n\n");
    variables.vetDX.setVetDoubleAttr_X(model, false);
    for(int i = 0; i < numItems; ++i)
        std::printf("DX[%d] = %.1f\n", i, variables.vetDX.getX_value(i));

    variables.matRot.setVetDoubleAttr_X(model, false);
    for(int i = 0; i < numItems; ++i)
    {
        for(auto r : vetRot)
        {
            std::printf(
                "r[%i, %i] = %.0f\n", i, (int)r, variables.matRot.getX_value(i, (int)r));
        }
    }

    if(model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE)
    {
        std::cout << "\nFound a solution!\n";
    }

    // testeOroloc3D();
    return 0;
}
