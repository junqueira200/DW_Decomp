#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"
#include <filesystem>
#include "VrpTW_DecompLabeling.h"
#include "LabelingAlgorithm.h"
#include "MemoryPool.h"
#include "BranchAndPrice.h"
#include "Alarm.h"
#include "Test.h"
#include <bits/stdc++.h>


//import teste;
// http://vrp.galgos.inf.puc-rio.br

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;
using namespace LabelingAlgorithmNS;
using namespace VrpTW_DecompLabelingNS;
using namespace BranchAndPriceNS;
using namespace SearchStrategyNS;
using namespace PrimalHeuristicNS;
using namespace BranchNS;
using namespace StatisticsNS;
using namespace TestNS;



int main(int argv, char **argc)
{

    try
    {
        std::feclearexcept(FE_OVERFLOW);
        std::feclearexcept(FE_UNDERFLOW);

        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";

        InstVRP_TW instVrpTw;
        std::string strFile(argc[1]);

        std::filesystem::path p(strFile);
        std::string fileName = p.filename();

        if(fileName[0] == 'A' || fileName[0] == 'P')
            leInstanciaAugerat(strFile, instVrpTw);
        else
            leInstanciaSalomon(strFile, instVrpTw);

        ptr_instVrpG = &instVrpTw;


        //RouteHash routeHash;
        //enumerateRoutes(instVrpTw, 5, routeHash);
        //return 0;

        //getSubInstancia(15, instVrpTw);

        VrpLabelingSubProb vrpLabelingSubProb(instVrpTw, instVrpTw.sumDist());

        if(argv == 3)
        {
            int option = atoi(argc[2]);
            if(option == 1)
            {
                vrpLabelingSubProb.setTypeLabelToBackward();
                std::cout<<"Seting typeLabel to backward\n";
            }
            else
                std::cout<<"Seting typeLabel to forward\n";
        }

        GRBEnv grbEnv;
        GRBModel model(grbEnv);
        criaMestre(instVrpTw, model);

        double distVarA = 100*std::max(somaDist(instVrpTw), 0.0);
        //VrpSubProb vrpSubProb(grbEnv, instVrpTw);

        DW_DecompNS::AuxData auxVectors;
        auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes * instVrpTw.numClientes));

        setAlarm(30.0*60); // 1.5 min timer

        std::cout << "Cria decompNode\n";
        DW_DecompNS::DW_DecompNode decompNode(grbEnv, model, distVarA, (DW_DecompNS::SubProb*)&vrpLabelingSubProb, 1, auxVectors);
        decompNode.rhsConv = instVrpTw.numVeic;
        //DepthFirst depthFirst;
        MinFuncObj minFuncObj;
        SimpleDiving simpleDiving;
        //StrongBranch branch;
        SimpleStrongBranch branch;

        StatisticsData statisticD;

        Eigen::VectorXd vetSol;

        CapacityCut capacityCut(instVrpTw, 5, 10, 0.00005);

        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";
        std::cout<<"Cost Var A.: "<<distVarA<<"\n\n";

        vetSol = branchAndPrice(decompNode,
                                auxVectors,
                                (SearchDataInter*)&minFuncObj,
                                (PrimalHeuristicInter*)&simpleDiving,
                                (BranchInter*)&branch,
                                (RobustCutGenerator*)&capacityCut,
                                statisticD);



/*        Eigen::MatrixXi matSol(instVrpTw.numVeic, instVrpTw.numClientes/2);
        matSol.setZero();

        int pos = 1;
        int i   = 0;
        int r   = 0;

        for(int j=1; j < instVrpTw.numClientes; ++j)
        {
            if(vetSol[VrpTW_DecompNS::getIndex(i, j, instVrpTw.numClientes)] >= 0.99)
            {
                matSol.coeffRef(r, pos) = j;
                r += 1;
            }

        }

        for(r=0; r < instVrpTw.numVeic; ++r)
        {
            pos = 2;
            i = matSol(r, 1);

            while(i != 0)
            {
                for(int j=0; j < instVrpTw.numClientes; ++j)
                {
                    if(vetSol[VrpTW_DecompNS::getIndex(i, j, instVrpTw.numClientes)] >= 0.99)
                    {
                        matSol.coeffRef(r, pos) = j;
                        pos += 1;
                        i = j;
                        break;
                    }
                }
            }

        }

        for(int k=0; k < vetSol.size(); ++k)
        {
            if(vetSol[k] >= 0.99)
            {
                int i = k/instVrpTw.numClientes;
                int j = k%instVrpTw.numClientes;

                std::cout<<"("<<i<<" "<<j<<"), ";
            }
        }

        std::cout<<"\n\n"<<matSol<<"\n";*/

        statisticD.inst = fileName;
        statisticD.numNodes = instVrpTw.numClientes;
        std::string extraHead = "Alg";
        std::string extraCont;

        if(vrpLabelingSubProb.typeLabel == LabelingAlgorithmNS::AlgForward)
            extraCont += "forward";
        else if(vrpLabelingSubProb.typeLabel == LabelingAlgorithmNS::AlgBackward)
            extraCont += "backward";
        else
            extraCont += "bidirectional";

        writeToFile(statisticD, "result.csv", extraHead, extraCont);

        //decompNode.columnGeneration(auxVectors);


        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";

        std::cout<<"DistMin: "<<minDistG<<"\nDistMax: "<<maxDistG<<"\n\n";

        if((bool)std::fetestexcept(FE_OVERFLOW))
            std::cout << "Overflow flag after: " << (bool)std::fetestexcept(FE_OVERFLOW) << std::endl;
        if((bool)std::fetestexcept(FE_UNDERFLOW))
            std::cout << "Underflow flag after: " << (bool)std::fetestexcept(FE_UNDERFLOW) << std::endl;

        /*
        std::cout<<"0 < -1("<<doubleLess(0, -1.0, FloatEp)<<")\n";
        std::cout<<"-1 == 0("<<doubleEqual(0, -1.0, FloatEp)<<")\n";
        std::cout<<"0 <= -1("<<doubleLessEqual(0, -1.0, FloatEp)<<")\n";
        */

    }
/*    catch(char const* str)
    {
        std::cout<<"catch(char* ):\n";
        std::printf("%s", str);
        std::cout<<"\n\n";
    }*/
    catch(GRBException &e)
    {
        std::cout<<"GRBException:\n"<<e.getMessage()<<"\n";
        std::cout<<"Code: "<<e.getErrorCode()<<"\n\n";
    }

    /* ****************************************************************************************************************
     * ****************************************************************************************************************
     */

    return 0;
}
