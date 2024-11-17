#include <iostream>
#include "Instancia.h"
#include "VrpTW_Decomp.h"

using namespace InstanciaNS;
using namespace VrpTW_DecompNS;

int main(int argv, char **argc)
{
    InstVRP_TW instVrpTw;
    std::string strFile(argc[1]);

    leInstanciaAugerat(strFile, instVrpTw);

    GRBEnv grbEnv;
    GRBModel model(grbEnv);

    GRBModel modelComp(grbEnv);
    criaVRP_TW_CompleteModel(instVrpTw, modelComp);
    modelComp.optimize();

    return 0;

    criaMestre(instVrpTw, model);

    double distVarA = somaDist(instVrpTw);
    VrpSubProb vrpSubProb(grbEnv, instVrpTw);

    DW_DecompNS::AuxVectors auxVectors;
    auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes*instVrpTw.numClientes));

    DW_DecompNS::Info info;

    std::cout<<"Cria decompNode\n";
    DW_DecompNS::DW_DecompNode decompNode(grbEnv,
                                          model,
                                          distVarA,
                                          (DW_DecompNS::SubProb*)&vrpSubProb,
                                          1,
                                          auxVectors,
                                          info);


    decompNode.columnGeneration(auxVectors, info);

    std::cout<<"..";
    std::cout<<"Num de veic: "<<instVrpTw.numVeic<<"\n\n";

    return 0;
}