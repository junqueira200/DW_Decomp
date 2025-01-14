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

        VrpLabelingSubProb vrpLabelingSubProb(instVrpTw, 0.01*instVrpTw.sumDist());

        GRBEnv grbEnv;
        GRBModel model(grbEnv);
        criaMestre(instVrpTw, model);

        double distVarA = somaDist(instVrpTw);
        //VrpSubProb vrpSubProb(grbEnv, instVrpTw);

        DW_DecompNS::AuxData auxVectors;
        auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes * instVrpTw.numClientes));

        setAlarm(2.5*60*60); // 2.5H timer

        std::cout << "Cria decompNode\n";
        DW_DecompNS::DW_DecompNode decompNode(grbEnv, model, distVarA, (DW_DecompNS::SubProb*)&vrpLabelingSubProb, 1, auxVectors);
        //DepthFirst depthFirst;
        MinFuncObj minFuncObj;
        SimpleDiving simpleDiving;
        StrongBranch branch;

        StatisticsData statisticD;

        branchAndPrice(decompNode,
                       auxVectors,
                       (SearchDataInter*)&minFuncObj,
                       (PrimalHeuristicInter*)&simpleDiving,
                       (BranchInter*)&branch,
                       statisticD);


        statisticD.inst = fileName;
        writeToFile(statisticD, "result.csv");

        //decompNode.columnGeneration(auxVectors);


        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";




    }
    catch(char const* str)
    {
        std::cout<<"catch(char* ):\n";
        std::printf("%s", str);
        std::cout<<"\n\n";
    }
    catch(GRBException &e)
    {
        std::cout<<"GRBException:\n"<<e.getMessage()<<"\n";
        std::cout<<"Code: "<<e.getErrorCode()<<"\n\n";
    }

    return 0;
}