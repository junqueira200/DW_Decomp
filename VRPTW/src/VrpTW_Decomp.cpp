//
// Created by igor on 15/11/24.
//
#include "VrpTW_Decomp.h"

void VrpTW_DecompNS::criaMestre(const InstanceVRPTW_NS::InstanceVRPTW &instVrpTw, GRBModel &model)
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
    delete []varX;
}

int VrpTW_DecompNS::getIndex(int i, int j, int numClie)
{
    return i*numClie+j;
}


/*int64_t VrpTW_DecompNS::VrpSubProb::getNumberOfConvConstr()
{
    std::cout<<"VrpSubProb::getNumberOfConvConstr\n\n";
    return 0;
}*/

/*

VrpTW_DecompNS::VrpSubProb::VrpSubProb(GRBEnv &e, InstanciaNS::InstVRP_TW &instVrpTw_)
{

    instVrpTw = &instVrpTw_;
    subProb = std::make_unique<GRBModel>(e);

    buildSubProbModel();

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
                                               void* data,
                                               const int iniConv,
                                               int indSubProb,
                                               Eigen::VectorXd &vetCooefRestConv,
                                               const std::pair<int, int> &pairSubProb)
{

    static Eigen::VectorXd vetRedCost(vetC.size());
    static Eigen::VectorXi vetRota(instVrpTw->numClientes+1);

    vetRedCost.setZero();
    const int numClie = instVrpTw->numClientes;

    for(int i=0; i < numClie; ++i)
    {
        for(int j=0; j < numClie; ++j)
        {
            if(i == j)
                continue;

            int index = getIndex(i, j, numClie);
            vetRedCost[index] = vetC[index] - vetRowPi[j];
        }
    }

    int routeTam = 0;
    static int num = 0;

    double costHeur = geraSolHeuristica(*instVrpTw, vetRedCost, vetRota, routeTam);
    if(costHeur < -DW_DecompNS::TolObjSubProb)
    {
        std::cout<<"\tRoute with negative reduced cost was found\n";

        for(int i=0; i < (routeTam-1); ++i)
        {
            std::cout<<vetRota[i]<<" "<<vetRota[i+1]<<"\n";

            vetX[getIndex(vetRota[i], vetRota[i+1], instVrpTw->numClientes)] = 1;
        }

        custoRedNeg = true;

        num += 1;

        return 0;
    }

    num += 1;

//std::cout<<"Reduced costs: "<<vetRedCost.transpose()<<"\n\n";

    for(int i=0; i < numClie*numClie; ++i)
    {
        grbVarX[i].set(GRB_DoubleAttr_Obj, vetRedCost[i]);
    }

    static double time = 2.0;
    subProb->set(GRB_DoubleParam_TimeLimit, time);

    subProb->update();
    subProb->write("subProb.lp");
    subProb->optimize();

    while(subProb->get(GRB_DoubleAttr_ObjVal) >= -DW_DecompNS::TolObjSubProb
          && subProb->get(GRB_IntAttr_Status) != GRB_OPTIMAL)
    {
        time += 2.0;
        subProb->set(GRB_DoubleParam_TimeLimit, time);

        subProb->optimize();
    }


    if(subProb->get(GRB_DoubleAttr_ObjVal) < -DW_DecompNS::TolObjSubProb)
    {
        vetX.setZero();
        int i = 0;
        int j = 0;
        int dem = 0;

        std::cout<<"0 ";

        do
        {
            for(j=0; j < numClie; ++j)
            {
                if(j == i)
                    continue;
                int index = getIndex(i, j, numClie);
                if(grbVarX[index].get(GRB_DoubleAttr_X) >= 0.99)
                    break;
            }

            std::cout<<j<<" ";
            dem += instVrpTw->vetClieDem[j];

            int index = getIndex(i, j, numClie);
            vetX[index] = 1;
            i = j;

        }
        while(j != 0);

        std::cout<<"\n\ndem: "<<dem<<"\n";
        custoRedNeg = true;

    }
    else
        custoRedNeg = false;

    return 0;
}

VrpTW_DecompNS::VrpSubProb::~VrpSubProb()
{
    delete []grbVarX;
    delete []grbVarF;
    delete []grbVarT_Cheg;
    delete []grbVarT_Saida;
    delete []grbVarU;
}

void VrpTW_DecompNS::VrpSubProb::buildSubProbModel()
{

    subProb->set(GRB_IntParam_Threads, 4);
    subProb->set(GRB_DoubleParam_TimeLimit, 2.0);

    delete []grbVarX;
    delete []grbVarF;
    delete []grbVarT_Cheg;
    delete []grbVarT_Saida;
    delete []grbVarU;

    const int numClie = instVrpTw->numClientes;
    grbVarX       = subProb->addVars(numClie*numClie, GRB_BINARY);
    grbVarF       = subProb->addVars(numClie*numClie, GRB_CONTINUOUS);
    //grbVarT_Cheg  = subProb->addVars(numClie, GRB_CONTINUOUS);
    //grbVarT_Saida = subProb->addVars(numClie, GRB_CONTINUOUS);
    grbVarU       = subProb->addVars(numClie, GRB_CONTINUOUS);

    //grbVarT_Saida[0].set(GRB_DoubleAttr_UB, 0.0);
    //grbVarT_Saida[0].set(GRB_DoubleAttr_LB, 0.0);

    GRBModel &model = *subProb;




    // \sum_{j \in V, j \not = 0} x_{0,j} = numVeic
    GRBLinExpr linExpr;
    for(int j=1; j < numClie; ++j)
    {

        linExpr += grbVarX[getIndex(0, j, numClie)];
    }

    model.addConstr(linExpr, '=', 1, "Constr_2");

    // \sum_{j \in V, j \not = 0} x_{j,0} = numVeic
    linExpr = 0;
    for(int j=1; j < numClie; ++j)
    {

        linExpr += grbVarX[getIndex(j, 0, numClie)];
    }

    model.addConstr(linExpr, '=', 1, "Constr_3");



    for(int j=1; j < numClie; ++j)
    {

        GRBLinExpr grbLinExpr = 0;

        for(int i = 0; i < numClie; ++i)
        {
            if(i != j)
                grbLinExpr += grbVarX[getIndex(i, j, numClie)];
        }

        GRBLinExpr grbLinExpr1 = 0;

        for(int i = 0; i < numClie; ++i)
        {
            if(i != j)
                grbLinExpr += -grbVarX[getIndex(j, i, numClie)];
        }


        //if(modelo3Index)
        model.addConstr(grbLinExpr + grbLinExpr1 == 0, "Restricao_4_j_" + std::to_string(j));
    }

    const int capacidade = instVrpTw->capVeic;
    for(int i=1; i < numClie; ++i)
    {
        for(int j=1; j < numClie; ++j)
        {
            if(i != j)
            {
                GRBLinExpr linExpr = 0;

                // uj - ui + Q(1-xijk)>=  qj
                linExpr = -grbVarU[i] + grbVarU[j] + capacidade * (1-grbVarX[getIndex(i, j, numClie)]);
                double rhs = instVrpTw->vetClieDem[j];
                model.addConstr(linExpr >= rhs, "Restricao_5_ij_" + std::to_string(i)+"_"+std::to_string(j));
            }
        }
    }

    model.addConstr(grbVarU[0] <= 0, "Restricao_6_i_" + std::to_string(0));
    model.addConstr(grbVarU[0] >= 0, "Restricao_7_i_" + std::to_string(0));


    // TW

    /*

    for(int j=0; j < numClie; ++j)
    {
        grbVarT_Cheg[j].set(GRB_StringAttr_VarName, "t_cheg_"+std::to_string(j));
        grbVarT_Saida[j].set(GRB_StringAttr_VarName, "t_saida_"+std::to_string(j));

    }

    for(int j=0; j < numClie; ++j)
    {

        for(int i=0; i < numClie; ++i)
        {   GRBLinExpr linExpr;
            if(i == j)
                continue;

            linExpr += grbVarT_Cheg[j] - grbVarT_Saida[i] - instVrpTw->matDist(i, j) + instVrpTw->vetClieTime[0].dueTime*(1-grbVarX[getIndex(i, j, numClie)]);
            subProb->addConstr(linExpr, '>', 0, "Rest5_"+std::to_string(i)+"_"+std::to_string(j));
        }

    }

    for(int i=1; i < numClie; ++i)
    {
        GRBLinExpr linExpr;
        linExpr += grbVarT_Cheg[i] + instVrpTw->vetClieTime[i].servTime - grbVarT_Saida[i];
        subProb->addConstr(linExpr, '<', 0, "Rest6_"+std::to_string(i));

        linExpr = 0;
        linExpr += instVrpTw->vetClieTime[i].readyTime - grbVarT_Cheg[i];
        subProb->addConstr(linExpr, '<', 0, "Rest7_"+std::to_string(i));

        linExpr = 0;
        linExpr += -instVrpTw->vetClieTime[i].dueTime + grbVarT_Cheg[i];
        subProb->addConstr(linExpr, '<', 0, "Rest8_"+std::to_string(i));


    }


    subProb->update();
    subProb->write("subProb.lp");

}

void VrpTW_DecompNS::criaVRP_TW_CompleteModel(const InstanciaNS::InstVRP_TW &instVrpTw, GRBModel &model)
{

    const int numClie = instVrpTw.numClientes;
    GRBVar *grbVarX       = model.addVars(numClie*numClie, GRB_BINARY);
    GRBVar *grbVarF       = model.addVars(numClie*numClie, GRB_CONTINUOUS);
    GRBVar *grbVarT_Cheg  = model.addVars(numClie, GRB_CONTINUOUS);
    GRBVar *grbVarT_Saida = model.addVars(numClie, GRB_CONTINUOUS);
    GRBVar *grbVarU       = model.addVars(numClie, GRB_CONTINUOUS);


    for(int i=0; i < numClie; ++i)
    {
        for(int j=0; j < numClie; ++j)
        {
            grbVarX[getIndex(i, j, numClie)].set(GRB_StringAttr_VarName, "x_"+std::to_string(i) + "_"+std::to_string(j));
            grbVarX[getIndex(i, j, numClie)].set(GRB_DoubleAttr_Obj, instVrpTw.matDist(i, j));
        }

        grbVarU[i].set(GRB_StringAttr_VarName, "U_"+std::to_string(i));
    }


    // \sum_{j \in V, j \not = i} x_{i,j} = 1  \forall i \in V_+
    for(int i=1; i < numClie; ++i)
    {
        GRBLinExpr linExpr;
        for(int j=0; j < numClie; ++j)
        {
            if(i == j)
                continue;

            linExpr += grbVarX[getIndex(i, j, numClie)];
        }

        model.addConstr(linExpr, '=', 1, "Constr_1");
    }


    // \sum_{j \in V, j \not = 0} x_{0,j} = numVeic
    GRBLinExpr linExpr;
    for(int j=1; j < numClie; ++j)
    {

        linExpr += grbVarX[getIndex(0, j, numClie)];
    }

    model.addConstr(linExpr, '=', instVrpTw.numVeic, "Constr_2");

    // \sum_{j \in V, j \not = 0} x_{j,0} = numVeic
    linExpr = 0;
    for(int j=1; j < numClie; ++j)
    {

        linExpr += grbVarX[getIndex(j, 0, numClie)];
    }

    model.addConstr(linExpr, '=', instVrpTw.numVeic, "Constr_3");



    for(int j=1; j < numClie; ++j)
    {

        GRBLinExpr grbLinExpr = 0;

        for(int i = 0; i < numClie; ++i)
        {
            if(i != j)
                grbLinExpr += grbVarX[getIndex(i, j, numClie)];
        }

        GRBLinExpr grbLinExpr1 = 0;

        for(int i = 0; i < numClie; ++i)
        {
            if(i != j)
                grbLinExpr += -grbVarX[getIndex(j, i, numClie)];
        }


        //if(modelo3Index)
        model.addConstr(grbLinExpr + grbLinExpr1 == 0, "Restricao_4_j_" + std::to_string(j));
    }

    const int capacidade = instVrpTw.capVeic;
    for(int i=1; i < numClie; ++i)
    {
        for(int j=1; j < numClie; ++j)
        {
            if(i != j)
            {
                GRBLinExpr linExpr = 0;

                // uj - ui + Q(1-xijk)>=  qj
                linExpr = -grbVarU[i] + grbVarU[j] + capacidade * (1-grbVarX[getIndex(i, j, numClie)]);
                double rhs = instVrpTw.vetClieDem[j];
                model.addConstr(linExpr >= rhs, "Restricao_5_ij_" + std::to_string(i)+"_"+std::to_string(j));
            }
        }
    }

    model.addConstr(grbVarU[0] <= 0, "Restricao_6_i_" + std::to_string(0));
    model.addConstr(grbVarU[0] >= 0, "Restricao_7_i_" + std::to_string(0));




    model.set(GRB_IntParam_Threads, 4);
    //model.set(GRB_IntParam_MIPFocus, 2);
    model.update();
    model.write("vrp.lp");


    delete []grbVarX;
    delete []grbVarF;
    delete []grbVarT_Cheg;
    delete []grbVarT_Saida;
    delete []grbVarU;

}

double VrpTW_DecompNS::geraSolHeuristica(const InstanciaNS::InstVRP_TW &instVrpTw, Eigen::VectorXd &vetRedCost,
                                         Eigen::VectorXi &rota, int &routeTam)
{
std::cout<<"ini geraSolHeuristica"<<"\n";

    static Eigen::VectorXi vetClieVisi(instVrpTw.numClientes);
    vetClieVisi.setZero();
    rota.setZero();

    int clienteI = 0;

    int next = 1;
    int carga = 0;
    double cost = 0.0;
    double dist = 0.0;
    routeTam = 1;

    std::cout<<"\tRota criada: ";
    std::cout<<clienteI<<" ";

    do
    {
        int clienteJ = -1;
        double menor = std::numeric_limits<double>::infinity();

        for(int i=0; i < instVrpTw.numClientes; ++i)
        {
            if(vetClieVisi[i] != 0 || i == clienteI)
                continue;

            if((carga+instVrpTw.vetClieDem[i]) > instVrpTw.capVeic)
                continue;

            const int index = getIndex(clienteI, i, instVrpTw.numClientes);
            if(vetRedCost[index] < menor)
            {
                clienteJ = i;
                menor = vetRedCost[index];
            }
        }

        if(clienteJ == -1)
        {
            clienteJ = 0;
            cost += vetRedCost[getIndex(clienteI, clienteJ, instVrpTw.numClientes)];
        }
        else
            cost += menor;
        std::cout<<clienteJ<<" ";

        rota[next] = clienteJ;

        vetClieVisi[clienteJ] = 1;
        next += 1;
        routeTam += 1;
        clienteI = clienteJ;
        carga += carga+instVrpTw.vetClieDem[clienteJ];

    }while(clienteI != 0);


    std::cout<<"; custo: "<<cost<<"\n";

    return cost;
}
*/