#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"
#include <filesystem>

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;

int main(int argv, char **argc)
{
    try
    {

        InstVRP_TW instVrpTw;
        std::string strFile(argc[1]);

        std::filesystem::path p(strFile);
        std::string fileName = p.filename();

        if(fileName[0] == 'A' || fileName[0] == 'P')
            leInstanciaAugerat(strFile, instVrpTw);
        else
            leInstanciaSalomon(strFile, instVrpTw);

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