/* ****************************************
 * ****************************************
 *  Data:    21/02/26
 *  Arquivo: MILP.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L_SDVRP
 * ****************************************
 * ****************************************/

#include "MILP.h"
#include "Instancia.h"

using namespace MILP_NS;
using namespace InstanceNS;

VectorGRBVar::VectorGRBVar(GRBModel &model, int num_, const std::string &&name, char type)
{
    started = true;
    num     = num_;
    varType = type;
    vetVar  = model.addVars(num, type);

    for(int i=0; i < num; ++i)
        vetVar[num].set(GRB_StringAttr_VarName, name+"_("+ std::to_string(i)+")");
}

void VectorGRBVar::setUB(double ub)
{
    for(int i=0; i < num; ++i)
        vetVar[i].set(GRB_DoubleAttr_UB, ub);
}

void VectorGRBVar::setLB(double lb)
{
    for(int i=0; i < num; ++i)
        vetVar[i].set(GRB_DoubleAttr_UB, lb);
}

void VectorGRBVar::setUB_LB(double ub, double lb)
{
    for(int i=0; i < num; ++i)
    {
        vetVar[i].set(GRB_DoubleAttr_UB, lb);
        vetVar[i].set(GRB_DoubleAttr_UB, ub);
    }
}


void VectorGRBVar::setUB_LB(double ub, double lb, int i)
{
    vetVar[i].set(GRB_DoubleAttr_UB, lb);
    vetVar[i].set(GRB_DoubleAttr_UB, ub);
}

void VectorGRBVar::printVars()
{
    for(int i=0; i < num; ++i)
        std::cout<<vetVar[i].get(GRB_StringAttr_VarName)<<" ";

    std::cout<<"\n";
}

void VectorGRBVar::setVetDoubleAttr_X(GRBModel &model, bool Xn)
{
    delete []vetDoubleAttr_X;

    if(Xn)
        vetDoubleAttr_X = model.get(GRB_DoubleAttr_Xn, vetVar, num);
    else
        vetDoubleAttr_X = model.get(GRB_DoubleAttr_X, vetVar, num);
}

void VectorGRBVar::start(GRBModel &model, int num_, const std::string &&name, char type)
{
    if(started)
    {
        PRINT_DEBUGG("", "");
        std::cout<<"ERRO, vetor ja foi inicializada";
        exit(-1);
    }


    started = true;
    num     = num_;
    varType = type;
    vetVar  = model.addVars(num, type);

    for(int i=0; i < num; ++i)
        vetVar[i].set(GRB_StringAttr_VarName, name+"_("+ std::to_string(i)+")");

}

MatrixGRBVar::MatrixGRBVar(GRBModel &model, int numLin_, int numCol_, const std::string &&nome, char type, bool zeroI_EqualJ)
{
    inicializado = true;
    numLin = numLin_;
    numCol = numCol_;
    typeVar = type;
    vetVar = model.addVars(numLin*numCol, type);

    for(int i=0; i < numLin; ++i)
    {

        for(int j=0; j < numCol; ++j)
        {
            const int index = i*numCol+j;
            vetVar[index].set(GRB_StringAttr_VarName, nome+"_("+ std::to_string(i)+","+ std::to_string(j)+")");

            if(i == j && zeroI_EqualJ)
            {
                vetVar[index].set(GRB_DoubleAttr_LB, 0.0);
                vetVar[index].set(GRB_DoubleAttr_UB, 0.0);
            }
        }
    }

}

void MatrixGRBVar::setUB(const double ub)
{
    for(int i=0; i < numCol*numLin; ++i)
        vetVar[i].set(GRB_DoubleAttr_UB, ub);
}

void MatrixGRBVar::setLB(double lb)
{

    for(int i=0; i < numCol*numLin; ++i)
        vetVar[i].set(GRB_DoubleAttr_LB, lb);
}

void MatrixGRBVar::setUB_LB(double ub, double lb)
{

    for(int i=0; i < numCol*numLin; ++i)
    {
        vetVar[i].set(GRB_DoubleAttr_UB, ub);
        vetVar[i].set(GRB_DoubleAttr_LB, lb);
    }
}

void MatrixGRBVar::printVars()
{
    for(int i=0; i < numLin; ++i)
    {
        for(int j=0; j < numCol; ++j)
        {
            const int index = i * numCol + j;
            const auto nome = vetVar[index].get(GRB_StringAttr_VarName);
            std::cout<<nome<<" ";
        }
        std::cout<<"\n";
    }
}

void MatrixGRBVar::start(GRBModel &model, int numLin_, int numCol_, const std::string &&nome, char type, bool zeroI_EqualJ)
{
    if(inicializado)
    {
        PRINT_DEBUGG("", "");
        std::cout<<"ERRO, matrix ja foi inicializada";
        throw "ERROR";
    }

    inicializado = true;

    numLin = numLin_;
    numCol = numCol_;
    typeVar = type;
    vetVar = model.addVars(numLin*numCol, type);
    if(!zeroI_EqualJ)
    {
        std::printf("numLin: %d; numCol: %d\n", numLin, numCol);
        std::printf("numVars: %d\n", numLin*numCol);
    }
    for(int i=0; i < numLin; ++i)
    {

        for(int j=0; j < numCol; ++j)
        {
            const int index = i*numCol+j;
            vetVar[index].set(GRB_StringAttr_VarName, nome+"_("+ std::to_string(i)+","+ std::to_string(j)+")");

            if(i == j && zeroI_EqualJ)
            {
                vetVar[index].set(GRB_DoubleAttr_LB, 0.0);
                vetVar[index].set(GRB_DoubleAttr_UB, 0.0);
            }
        }
    }

}

void MatrixGRBVar::setVetDoubleAttr_X(GRBModel &model, bool X_n)
{
    delete []vetDoubleAttr_X;

    if(X_n)
        vetDoubleAttr_X  = model.get(GRB_DoubleAttr_Xn, vetVar, numCol * numLin);
    else
        vetDoubleAttr_X = model.get(GRB_DoubleAttr_X, vetVar, numCol * numLin);

}

Variables::Variables(GRBModel& model, VectorI& vetItems, int numItems)
{

    vetPosX.start(model, numItems, "X", GRB_CONTINUOUS);
    vetPosY.start(model, numItems, "Y", GRB_CONTINUOUS);
    vetPosZ.start(model, numItems, "Z", GRB_CONTINUOUS);

    vetPosX.setUB_LB(instanciaG.vetDimVeiculo[0], 0.0);
    vetPosY.setUB_LB(instanciaG.vetDimVeiculo[1], 0.0);
    vetPosZ.setUB_LB(instanciaG.vetDimVeiculo[2], 0.0);

    vetDX.start(model, numItems, "DX", GRB_CONTINUOUS);
    vetDY.start(model, numItems, "Dy", GRB_CONTINUOUS);
    vetDZ.start(model, numItems, "Dz", GRB_CONTINUOUS);

    for(int i=0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        double min = std::min(item.vetDim[0], item.vetDim[1]);
        double max = std::max(item.vetDim[0], item.vetDim[1]);

        vetDX(i).set(GRB_DoubleAttr_LB, min);
        vetDX(i).set(GRB_DoubleAttr_UB, max);

        vetDY(i).set(GRB_DoubleAttr_LB, min);
        vetDY(i).set(GRB_DoubleAttr_UB, max);

        vetDZ(i).set(GRB_DoubleAttr_LB, item.vetDim[2]);
        vetDZ(i).set(GRB_DoubleAttr_UB, item.vetDim[2]);
    }


    matX_pos.start(model, numItems, numItems, "X_pos", GRB_BINARY, true);
    matY_pos.start(model, numItems, numItems, "Y_pos", GRB_BINARY, true);
    matZ_pos.start(model, numItems, numItems, "Z_pos", GRB_BINARY, true);

    matX_neg.start(model, numItems, numItems, "X_neg", GRB_BINARY, true);
    matY_neg.start(model, numItems, numItems, "Y_neg", GRB_BINARY, true);
    matZ_neg.start(model, numItems, numItems, "Z_neg", GRB_BINARY, true);

    matRot.start(model, numItems, vetRot.size(), "R", GRB_BINARY, false);

}

void MILP_NS::addBasicConstraints(GRBModel& model, Variables& variables, SolucaoNS::Bin& bin)
{

    for(int i=0; i < bin.numItens; ++i)
    {

        model.addConstr(variables.vetPosX(i)+variables.vetDX(i), '<', instanciaG.vetDimVeiculo[0]);
        model.addConstr(variables.vetPosY(i)+variables.vetDY(i), '<', instanciaG.vetDimVeiculo[1]);
        model.addConstr(variables.vetPosZ(i)+variables.vetDZ(i), '<', instanciaG.vetDimVeiculo[2]);
    }

    double mX = instanciaG.vetDimVeiculo[0];
    double mY = instanciaG.vetDimVeiculo[1];

    //std::cout<<"Constraint in matRot\n";
    for(int i=0; i < bin.numItens; ++i)
    {
        Item& item = instanciaG.vetItens[bin.vetItemId[i]];
        GRBLinExpr sumR;



        for(Rotation r:vetRot)
        {
            //std::cout<<"\t"<<r<<"\n";

            double dx = item.getDimRotacionada(0, r);
            double dy = item.getDimRotacionada(1, r);
            double dz = item.getDimRotacionada(2, r);

            GRBVar& varR = variables.matRot(i, (int)r);

            model.addGenConstrIndicator(varR, 1, variables.vetDX(i) == dx);
            model.addGenConstrIndicator(varR, 1, variables.vetDY(i) == dy);
            model.addGenConstrIndicator(varR, 1, variables.vetDZ(i) == dz);

            sumR += varR;
        }

        std::string name = std::format("sum r item {}", i);
        model.addConstr(sumR, '=', 1, name);
    }


    for(int i=0; i < bin.numItens; ++i)
    {
        for(int j=i+1; j < bin.numItens; ++j)
        {
            GRBLinExpr linExp;

            model.addGenConstrIndicator(variables.matX_pos(i, j), 1,  variables.vetPosX(i) + variables.vetDX(i) <=
                                                                      variables.vetPosX(j));
            model.addGenConstrIndicator(variables.matX_neg(i, j), 1,  variables.vetPosX(j) + variables.vetDX(j) <=
                                                                      variables.vetPosX(i));

            model.addGenConstrIndicator(variables.matY_pos(i, j), 1,  variables.vetPosY(i) + variables.vetDY(i) <=
                                                                      variables.vetPosY(j));
            model.addGenConstrIndicator(variables.matY_neg(i, j), 1,  variables.vetPosY(j) + variables.vetDY(j) <=
                                                                      variables.vetPosY(i));


            model.addGenConstrIndicator(variables.matZ_pos(i, j), 1,  variables.vetPosZ(i) + variables.vetDZ(i) <=
                                                                      variables.vetPosZ(j));
            model.addGenConstrIndicator(variables.matZ_neg(i, j), 1,  variables.vetPosZ(j) + variables.vetDZ(j) <=
                                                                      variables.vetPosZ(i));

            linExp = variables.matX_neg(i, j) + variables.matX_pos(i, j) + variables.matY_neg(i, j) +
                     variables.matY_pos(i, j) + variables.matZ_neg(i, j) + variables.matZ_pos(i, j);
            model.addConstr(linExp, '=', 1);

        }
    }

    GRBVar max = model.addVar(0.0, GRB_INFINITY, 1, GRB_CONTINUOUS);
    for(int i=0; i < bin.numItens; ++i)
    {
        model.addConstr(variables.vetPosX(i) <= max);
        model.addConstr(variables.vetPosY(i) <= max);
        model.addConstr(variables.vetPosZ(i) <= max);
    }

    GRBLinExpr obj = max;

    model.setObjective(obj, GRB_MINIMIZE);

    model.update();
    model.write("model.lp");

}
