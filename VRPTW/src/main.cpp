#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"
#include <filesystem>

//import teste;

#include <boost/dynamic_bitset.hpp>

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;

#include "LabelingAlgorithm.h"
#include "MemoryPool.h"

typedef  std::bitset<5> BitSet;

using namespace LabelingAlgorithmNS;

int main(int argv, char **argc)
{
    Eigen::Vector<Step, 2> vetStepSize;

    vetStepSize[0].stepSize = 1;
    vetStepSize[0].start    = 99;
    vetStepSize[0].end      = 0;

    LabelingData labelingData(vetStepSize, 1, 17);

    return 0;

/*    BitSet b0(0);
    //std::vector<bool> vet(5);
    b0[4] = true;
//    b0[3] = true;

    BitSet b1(0);
    b1[3] = true;
    b1[4] = true;

    BitSet r(0);
    r = b1 & b0;
    bool equal = r == b0;

    std::cout<<"b0: "<<b0<<"\nb1: "<<b1<<"\n";
    std::cout<<"equal: "<<equal<<"\n";*/

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

        NgSet ngSet(instVrpTw.numClientes, NgSetSize);
        ngSet.setNgSets(instVrpTw.matDist);

        std::cout<<"\ncontain: "<<ngSet.contain(0, 4)<<"\n";

        return 0;

//        exit(-1);
        GRBEnv grbEnv;
        GRBModel model(grbEnv);
        //model.set(GRB_IntParam_Threads, 4);
        //model.set(GRB_DoubleParam_TimeLimit, 30.0);

        GRBModel modelComp(grbEnv);
        criaVRP_TW_CompleteModel(instVrpTw, modelComp);
        modelComp.optimize();

        return 0;

        criaMestre(instVrpTw, model);

        double distVarA = somaDist(instVrpTw);
        VrpSubProb vrpSubProb(grbEnv, instVrpTw);

        DW_DecompNS::AuxVectors auxVectors;
        auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes * instVrpTw.numClientes));

        DW_DecompNS::Info info;

        std::cout << "Cria decompNode\n";
        DW_DecompNS::DW_DecompNode decompNode(grbEnv, model, distVarA, (DW_DecompNS::SubProb *) &vrpSubProb, 1,
                                              auxVectors, info);


        decompNode.columnGeneration(auxVectors, info);

        std::cout << "..";
        std::cout << "Num de veic: " << instVrpTw.numVeic << "\n\n";

    }
    catch(char const* str)
    {
        std::cout<<"catch(char* ):\n";
        std::printf("%s", str);
        std::cout<<"\n\n";
    }

    return 0;
}