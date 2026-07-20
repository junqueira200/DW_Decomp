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

#include "AuxT.h"
#include "BCRoutingParams.h"
#include "LoadingChecker.h"
#include "ProblemParameters.h"

using namespace InstanceNS;
using namespace ConstrutivoBinNS;
using namespace SolucaoNS;
using namespace RandNs;
using namespace ParseInputNS;
using namespace ConstrutivoNS;
using namespace IgNs;
using namespace BinPackingCP_NS;
using namespace TesteOroloc3D_NS;

using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;
using namespace MILP_NS;

int main(int argc, const char *argv[])
{
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


    for(int i=0; i < instanciaG.numItens; ++i)
    {
        if(i != instanciaG.vetItens[i].itemId)
        {
            std::printf("Error, i(%d) != itemId(%d)\n", i, instanciaG.vetItens[i].itemId);
            PRINT_THROW();
        }
    }

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

    if(!input.instOroloc3D_2)
    {
        setClassical3DPackingProblem();


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
