#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"
#include <filesystem>
#include "VrpTW_DecompLabeling.h"

//import teste;

#include <boost/dynamic_bitset.hpp>

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;

#include "LabelingAlgorithm.h"
#include "MemoryPool.h"

typedef  std::bitset<5> BitSet;

using namespace LabelingAlgorithmNS;
using namespace VrpTW_DecompLabelingNS;

int main(int argv, char **argc)
{
//    boost::dynamic_bitset<> dynamicBitset(8);
//    std::cout<<dynamicBitset<<"\n";
//    return 0;

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

        VrpLabelingSubProb vrpLabelingSubProb(instVrpTw, 0.01*instVrpTw.sumDist());

        GRBEnv grbEnv;
        GRBModel model(grbEnv);
        criaMestre(instVrpTw, model);

        double distVarA = somaDist(instVrpTw);
        //VrpSubProb vrpSubProb(grbEnv, instVrpTw);

        DW_DecompNS::AuxData auxVectors;
        auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes * instVrpTw.numClientes));

        DW_DecompNS::Info info;

        std::cout << "Cria decompNode\n";
        DW_DecompNS::DW_DecompNode decompNode(grbEnv, model, distVarA, (DW_DecompNS::SubProb *) &vrpLabelingSubProb, 1, auxVectors, info);

        int numSol;
        double maxDist = -10000;
        double redCost;

        /*
        forwardLabelingAlgorithm(2,
                                 instVrpTw.numClientes + 1,
                                 vrpLabelingSubProb.vetMatResCost,
                                 vrpLabelingSubProb.vetVetResBound,
                                 instVrpTw.numClientes,
                                 vrpLabelingSubProb.ngSet,
                                 vrpLabelingSubProb.labelingData,
                                 auxVectors.matColX_solSubProb,
                                 numSol,
                                 0.0,
                                 1,
                                 true,
                                 maxDist,
                                 redCost);
        */


        //std::cout<<"maxDist: "<<maxDist<<"\n";
        //std::cout<<"0.6*maxDist+1: "<<0.6*maxDist+1<<"\n\n";
        //vrpLabelingSubProb = VrpLabelingSubProb(instVrpTw, (0.6*maxDist+1));


        //return 0;

        decompNode.columnGeneration(auxVectors, info);

        std::cout << "..";
        std::cout << "Num de veic: " << instVrpTw.numVeic << "\n\n";

        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"Num Clientes: "<<instVrpTw.numClientes<<"\n";

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
    }

    return 0;
}