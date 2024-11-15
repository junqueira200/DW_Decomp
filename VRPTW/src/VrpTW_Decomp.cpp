//
// Created by igor on 15/11/24.
//
#include "VrpTW_Decomp.h"

void VrpTW_DecompNS::criaMestre(const InstanciaNS::InstVRP_TW &instVrpTw, GRBModel &model)
{
    const int NumClie = instVrpTw.numClientes;

    GRBVar *varX = model.addVars(NumClie*NumClie, GRB_BINARY);

    for(int i=0; i < NumClie; ++i)
    {
        for(int j=0; j < NumClie; ++j)
        {
            int index = getIndex(i, j, NumClie);
            varX[index].set(GRB_StringAttr_VarName, "X_"+std::to_string(i)+"_"+std::to_string(j));
            varX[index].set(GRB_DoubleAttr_Obj, instVrpTw.matDist(i, j));
        }
    }

    for(int i=0; i < NumClie; ++i)
    {
        GRBLinExpr linExpr;
        for(int j=0; j < NumClie; ++j)
        {
            if(i==j)
                continue;

            linExpr += 1*varX[getIndex(j, i, NumClie)];
        }

        model.addConstr(linExpr >= 1, "Clie_"+std::to_string(i));

    }

    model.update();
    model.write("vrp.lp");

}

int VrpTW_DecompNS::getIndex(int i, int j, int numClie)
{
    return i*numClie+j;
}


int64_t VrpTW_DecompNS::VrpSubProb::getNumberOfConvConstr()
{
    std::cout<<"VrpSubProb::getNumberOfConvConstr\n\n";
    return 0;
}


VrpTW_DecompNS::VrpSubProb::VrpSubProb(GRBEnv &e, InstanciaNS::InstVRP_TW &instVrpTw_)
{

    instVrpTw = &instVrpTw_;

}

void VrpTW_DecompNS::VrpSubProb::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{

}

int VrpTW_DecompNS::VrpSubProb::resolveSubProb(const Eigen::VectorXd &vetC,
                                               const Eigen::RowVectorXd &vetRowPi,
                                               GRBModel &mestre,
                                               Eigen::VectorXd &vetX,
                                               int itCG,
                                               bool &custoRedNeg,
                                               void *data,
                                               const int iniConv,
                                               int indSubProb,
                                               Eigen::VectorXd &vetCooefRestConv,
                                               const std::pair<int, int> &pairSubProb)
{

    return 0;
}