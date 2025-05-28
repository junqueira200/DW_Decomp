#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"
#include <filesystem>
#include "VrpTW_DecompLabeling.h"
#include "LabelingAlgorithm.h"
#include "MemoryPool.h"
#include "BranchAndPrice.h"
#include "Alarm.h"

//import teste;

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;
using namespace LabelingAlgorithmNS;
using namespace VrpTW_DecompLabelingNS;
using namespace BranchAndPriceNS;
using namespace SearchStrategyNS;
using namespace PrimalHeuristicNS;
using namespace BranchNS;
using namespace StatisticsNS;

int main(int argv, char **argc)
{
    /*
    std::bitset<10> bitset0 = 0;
    std::bitset<10> bitset1 = 0;

    bitset0[4] = true;
    bitset1[5] = true;

    std::bitset<10> result = bitset0&bitset1;
    std::cout<<result<<"\n";

    std::cout<<"==0: "<<((bitset0&bitset1) == 0)<<"\n";

    return 0;
    */

    try
    {
        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";

        InstVRP_TW instVrpTw;
        std::string strFile(argc[1]);

        std::filesystem::path p(strFile);
        std::string fileName = p.filename();

        if(fileName[0] == 'A' || fileName[0] == 'P')
            leInstanciaAugerat(strFile, instVrpTw);
        else
            leInstanciaSalomon(strFile, instVrpTw);

        getSubInstancia(16, instVrpTw);

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

        double distVarA = std::max(somaDist(instVrpTw), 0.0);
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


        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";
        std::cout<<"Cost Var A.: "<<distVarA<<"\n\n";

        vetSol = branchAndPrice(decompNode, auxVectors, (SearchDataInter*)&minFuncObj,
                                (PrimalHeuristicInter*)&simpleDiving, (BranchInter*)&branch, statisticD);



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
        writeToFile(statisticD, "result.csv");

        //decompNode.columnGeneration(auxVectors);


        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";

        std::cout<<"DistMin: "<<minDistG<<"\nDistMax: "<<maxDistG<<"\n\n";



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

    return 0;
}
