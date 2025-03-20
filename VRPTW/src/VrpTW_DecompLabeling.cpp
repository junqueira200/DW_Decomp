/*  *****************************************************************
 *  *****************************************************************
 *  File:    VrpTW_DecompLabeling.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    20/11/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "VrpTW_DecompLabeling.h"

using namespace LabelingAlgorithmNS;


VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb(InstanciaNS::InstVRP_TW &instVrpTw_, double startDist)
{
    if(instVrpTw_.numClientes > NumMaxCust)
    {
        PRINT_DEBUG("", "ERRO; Const NumMaxCust("<<NumMaxCust<<") < NumCust("<<instVrpTw_.numClientes<<")\n\t");
        throw "ERRO";
    }

    instVrpTw = &instVrpTw_;

    //vetStepSize[0].stepSize = 400;
    vetStepSize[0].stepSize = 250.0; // 200
    vetStepSize[0].start    = (FloatType)-startDist;
    vetStepSize[0].end      = (FloatType)startDist;

    vetStepSize[1].stepSize = 30.0;  //30
    vetStepSize[1].start    = 0;
    vetStepSize[1].end      = (FloatType)instVrpTw->capVeic;


    labelingData  = LabelingAlgorithmNS::LabelingData(vetStepSize, 2, instVrpTw->numClientes+1);
    vetMatResCost = Vet3D_ResCost(instVrpTw->numClientes+1, instVrpTw->numClientes+1, 2);

    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        for(int j=0; j < (instVrpTw->numClientes+1); ++j)
        {
            if(i == j)
                continue;

            vetMatResCost.get(i, j, 1) = (FloatType)instVrpTw->vetClieDem[j];
        }
    }

    //vetVetResBound = LabelingAlgorithmNS::VetVetResBound(2);

    vetVetResBound = LabelingAlgorithmNS::MatBoundRes(instVrpTw->numClientes+1, 2);

    //vetVetResBound[0].resize(instVrpTw->numClientes+1);
    //vetVetResBound[1].resize(instVrpTw->numClientes+1);

    double distTotal = instVrpTw->sumDist();

    Bound bound0;
    bound0.lowerBound = -std::numeric_limits<FloatType>::infinity();
    bound0.upperBound =  std::numeric_limits<FloatType>::infinity();

    Bound bound1;
    bound1.lowerBound = 0;
    bound1.upperBound = (FloatType)instVrpTw->capVeic;

    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        vetVetResBound(i, 0) = bound0;
        if(i > 0 && i < instVrpTw->numClientes)
            bound1.lowerBound = (FloatType)instVrpTw->vetClieDem[i];
        else
            bound1.lowerBound = 0.0;

        vetVetResBound(i, 1) = bound1;
    }

    ngSet = NgSet(instVrpTw->numClientes+1, NgSetSize);
    ngSet.setNgSets(instVrpTw->matDist);
    ngSet.active = true;

}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{


    GRBLinExpr linExpr;
    GRBVar a = rmlp.addVar(0, GRB_INFINITY, custoVarA, GRB_CONTINUOUS);
    linExpr += -a;

    rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");


}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::resolveSubProb(const Eigen::VectorXd &vetC,
                                                               Eigen::RowVectorXd &vetRowPi,
                                                               GRBModel &mestre,
                                                               int itCG,
                                                               bool &custoRedNeg,
                                                               void *data,
                                                               const int iniConv,
                                                               int indSubProb,
                                                               Eigen::VectorXd &vetCooefRestConv,
                                                               const std::pair<int, int> &pairSubProb,
                                                               Eigen::MatrixXd &matColX,
                                                               int &numSol,
                                                               double& redCost,
                                                               double constPiValue,
                                                               const VectorI &vetVar0,
                                                               const VectorI &vetVar1,
                                                               DW_DecompNS::PhaseStatus phaseStatus)
{

    static DW_DecompNS::PhaseStatus phase = phaseStatus;

    //if(phase != phaseStatus && phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusColGen)
    {

        minDistG = std::numeric_limits<FloatType>::max();
        maxDistG = std::numeric_limits<FloatType>::min();
    }

    phase = phaseStatus;

//std::cout<<"constPiValue: "<<constPiValue<<"\n";
    static Eigen::VectorX<FloatType> vetRedCostFT(DW_DecompNS::NumMaxSolSubProb);

    //vetMatResCost.setVal(0.0);
    //constPiValue += -vetRowPi[0];

    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        for(int j=0; j < instVrpTw->numClientes; ++j)
        {
            if(i == j)
                continue;

            if(phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
                vetMatResCost.get(i, j, 0) =  - (FloatType)vetRowPi[j+1];
            else
                vetMatResCost.get(i, j, 0) = (FloatType)instVrpTw->matDist(i, j) - (FloatType)vetRowPi[j+1];
        }

        if(i != 0)
        {
            //vetMatResCost[0](i, 0) = instVrpTw->matDist(i, 0);
            if(phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
                vetMatResCost.get(i, instVrpTw->numClientes, 0) = - (FloatType)vetRowPi[1];
            else
                vetMatResCost.get(i, instVrpTw->numClientes, 0) = (FloatType)instVrpTw->matDist(i, 0) -
                                                                  (FloatType)vetRowPi[1];
        }
    }

    for(int varId:vetVar0)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        if(j == 0)
            j = instVrpTw->numClientes;

        vetMatResCost.get(i, j, 0) = std::numeric_limits<FloatType>::infinity();
        //std::cout<<varId<<"\n";
    }

    std::set<int> setI;
    int t = 0;

    /*
    for(int varId:vetVar1)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        //std::cout<<i<<" "<<j<<"\n\n";

        if(setI.contains(i))
        {
            vetMatResCost[0](i, j) = (FloatType)instVrpTw->matDist(i, j) - (FloatType)(vetRowPi[j+1] + vetRowPi[1+instVrpTw->numClientes+t]);
            if(j == 0)
                vetMatResCost[0](i, instVrpTw->numClientes) = vetMatResCost[0](i, j);
        }
        else
        {
            for(int ii=0; ii < instVrpTw->numClientes+1; ++ii)
            {
                if(ii == i)
                    continue;

                vetMatResCost[0](ii, j) = std::numeric_limits<FloatType>::infinity();
            }


            if(j == 0)
            {

                for(int ii=0; ii < instVrpTw->numClientes; ++ii)
                {
                    if(ii == i)
                        continue;

                    vetMatResCost[0](ii, instVrpTw->numClientes) = std::numeric_limits<FloatType>::infinity();
                }
            }

            //vetMatResCost[0](i, j) += -vetRowPi[1+instVrpTw->numClientes+t];
            vetMatResCost[0](i, j) = (FloatType)instVrpTw->matDist(i, j) - (FloatType)(vetRowPi[j+1] + vetRowPi[1+instVrpTw->numClientes+t]);
            if(j == 0)
                vetMatResCost[0](i, instVrpTw->numClientes) = vetMatResCost[0](i, j);

            setI.insert(i);
        }

        t += 1;
    }
    */

    //std::cout<<"Custo Reduzido: \n"<<vetMatResCost[0]<<"\n\n";
    //std::cout<<"Peso: "<<vetMatResCost[1]<<"\n\n";

    int it = 0;
    FloatType maxDist;

    ngSet.active = true;
    numSol = 0;

    exactLabelingG = false;

    for(int i=2; i <= 2; i += 1)
    {
        //std::cout<<"forwardLabelingAlgorithm: "<<i<<"\n\n";
        matColX.setZero();
        custoRedNeg = forwardLabelingAlgorithm(2,
                                               instVrpTw->numClientes + 1,
                                               vetMatResCost,
                                               vetVetResBound,
                                               instVrpTw->numClientes,
                                               ngSet,
                                               labelingData,
                                               matColX,
                                               numSol,
                                               (FloatType)constPiValue,
                                               i,
                                               true,
                                               maxDist,
                                               vetRedCostFT);

        it += 1;
        if(custoRedNeg)
            break;

    }


    if(!custoRedNeg)
    {
        matColX.setZero();

        //std::cout<<"EXACT LABELING\n";
        //exactLabelingG = true;
        custoRedNeg = forwardLabelingAlgorithm(2,
                                               instVrpTw->numClientes+1,
                                               vetMatResCost,
                                               vetVetResBound,
                                               instVrpTw->numClientes,
                                               ngSet,
                                               labelingData,
                                               matColX,
                                               numSol,
                                               (FloatType)constPiValue,
                                               -1,
                                               true,
                                               maxDist,
                                               vetRedCostFT);


        /*
        if(!custoRedNeg)
        {
            matColX.setZero();
            ngSet.active = false;
            custoRedNeg = forwardLabelingAlgorithm(2,
                                                   instVrpTw->numClientes+1,
                                                   vetMatResCost,
                                                   vetVetResBound,
                                                   instVrpTw->numClientes,
                                                   ngSet,
                                                   labelingData,
                                                   matColX,
                                                   numSol,
                                                   (FloatType)constPiValue,
                                                   -1,
                                                   true,
                                                   maxDist,
                                                   vetRedCostFT);

            if(custoRedNeg)
                std::cout<<"IMPROVE\n";
        }
        */

    }

    //redCost = (double)redCostFT;
    vetCooefRestConv[0] = 1;


/*    if(numSol == 0)
    {
        std::cout<<"CALL EXACT PRICING!\n\n";

        Eigen::Matrix<double, -1, 1, Eigen::ColMajor> matCol(instVrpTw->numClientes*instVrpTw->numClientes, 1);

        if(exactPricing(vetMatResCost, constPiValue, matCol, *instVrpTw, vetRedCostFT[0]))
        {

            custoRedNeg = forwardLabelingAlgorithm(2,
                                                   instVrpTw->numClientes+1,
                                                   vetMatResCost,
                                                   vetVetResBound,
                                                   instVrpTw->numClientes,
                                                   ngSet,
                                                   labelingData,
                                                   matColX,
                                                   numSol,
                                                   (FloatType)constPiValue,
                                                   -1,
                                                   true,
                                                   maxDist,
                                                   vetRedCostFT);

            if(custoRedNeg)
            {
                std::cout<<"WTF?\n";
                PRINT_DEBUG("", "");
                throw "ERROR";
            }
            else
            {   std::cout<<"\nERROR, labeling falhou na ultima execucao!";
                PRINT_DEBUG("", "");
                throw "ERROR";
            }


            for(int i=0; i < instVrpTw->numClientes*instVrpTw->numClientes; ++i)
                matColX(i, 0) = matCol[i];

            custoRedNeg = true;
            numSol = 1;



        }
    }*/

    // Check if solution have a negative reduced cost
    FloatType redCostTemp = 0.0;
    for(int j=0; j < numSol; ++j)
    {
        redCostTemp = constPiValue;
        std::cout<<j<<": ";
        for(int i=0; i < (instVrpTw->numClientes*instVrpTw->numClientes); ++i)
        {
            if(matColX(i,j) != 1.0)
                continue;

            int ii = i/instVrpTw->numClientes;
            int jj = i%instVrpTw->numClientes;

            redCostTemp += vetMatResCost(ii, jj, 0);//instVrpTw->matDist(ii, jj);
            std::cout<<"("<<ii<<","<<jj<<")["<<i<<"]; ";

        }
        std::cout<<"\n\n";
        if(redCostTemp >= -DW_DecompNS::TolObjSubProb || !doubleEqual(vetRedCostFT[j], redCostTemp, 1E-3))
        {
            std::cout<<"\nERROR, custo reduzido calculado: ("<<redCostTemp<<") \n";
            std::cout<<"dif: "<<(vetRedCostFT[j]-redCostTemp)<<"\n";
            std::cout<<"redCost: "<<vetRedCostFT[j]<<"\n\n";
            PRINT_DEBUG("", "");
            std::cout<<"phaseStatus: "<<(int)phaseStatus<<"\n\n";
            std::cout<<(redCostTemp >= -DW_DecompNS::TolObjSubProb)<<"\n";
            throw "ERROR";
        }

    }

    //std::cout<<"["<<minDistG<<", "<<maxDistG<<"]; numSol("<<numSol<<")\n";

    return 0;
}


bool VrpTW_DecompLabelingNS::exactPricing(const LabelingAlgorithmNS::Vet3D_ResCost&          vetMatResCost,
                                          const FloatType                                    startVal,
                                          Eigen::Matrix<double, -1, 1, Eigen::ColMajor>&     vetColSolX,
                                          InstanciaNS::InstVRP_TW&                           instVrpTw,
                                          double&                                            cost)
{

    static GRBEnv env;
    static GRBModel* ptrModel = nullptr;
    static GRBVar* vetVarX    = nullptr;
    static GRBVar* vetVarU    = nullptr;

    if(!ptrModel)
    {

        ptrModel = new GRBModel(env);
        vetVarX  = ptrModel->addVars(instVrpTw.numClientes*instVrpTw.numClientes, GRB_BINARY);
        vetVarU  = ptrModel->addVars(instVrpTw.numClientes, GRB_CONTINUOUS);



        GRBLinExpr linExpr;
        for(int j=1; j < instVrpTw.numClientes; ++j)
        {

            linExpr += vetVarX[getIndex(0, j, instVrpTw.numClientes)];
        }

        ptrModel->addConstr(linExpr, '=', 1, "Constr_0");

        // \sum_{j \in V, j \not = 0} x_{j,0} = numVeic
        linExpr = 0;
        for(int j=1; j < instVrpTw.numClientes; ++j)
        {

            linExpr += vetVarX[getIndex(j, 0, instVrpTw.numClientes)];
        }

        ptrModel->addConstr(linExpr, '=', 1, "Constr_1");



        for(int j=1; j < instVrpTw.numClientes; ++j)
        {

            GRBLinExpr grbLinExpr = 0;

            for(int i = 0; i < instVrpTw.numClientes; ++i)
            {
                if(i != j)
                    grbLinExpr += vetVarX[getIndex(i, j, instVrpTw.numClientes)];
            }

            GRBLinExpr grbLinExpr1 = 0;

            for(int i = 0; i < instVrpTw.numClientes; ++i)
            {
                if(i != j)
                    grbLinExpr += -vetVarX[getIndex(j, i, instVrpTw.numClientes)];
            }


            //if(modelo3Index)
            ptrModel->addConstr(grbLinExpr + grbLinExpr1 == 0, "Restricao_2_j_" + std::to_string(j));
        }

        const int capacidade = instVrpTw.capVeic;
        for(int i=1; i < instVrpTw.numClientes; ++i)
        {
            for(int j=1; j < instVrpTw.numClientes; ++j)
            {
                if(i != j)
                {
                    GRBLinExpr linExpr = 0;

                    // uj - ui + Q(1-xijk)>=  qj
                    linExpr = -vetVarU[i] + vetVarU[j] + capacidade * (1-vetVarX[getIndex(i, j, instVrpTw.numClientes)]);
                    double rhs = instVrpTw.vetClieDem[j];
                    ptrModel->addConstr(linExpr >= rhs, "Restricao_3_ij_" + std::to_string(i)+"_"+std::to_string(j));
                }
            }
        }

        ptrModel->addConstr(vetVarU[0] <= 0, "Restricao_6_i_" + std::to_string(0));
        ptrModel->addConstr(vetVarU[0] >= 0, "Restricao_7_i_" + std::to_string(0));

    }

    GRBLinExpr funcObj = startVal;

    for(int i=0; i < instVrpTw.numClientes; ++i)
    {
        for(int j=0; j < instVrpTw.numClientes; ++j)
        {
            if(i == j)
                continue;

            FloatType coeef = vetMatResCost(i, j, 0);
            funcObj += coeef*vetVarX[getIndex(i, j, instVrpTw.numClientes)];
        }
    }

    ptrModel->setObjective(funcObj, GRB_MINIMIZE);
    ptrModel->set(GRB_DoubleParam_BestObjStop, -0.00000009);
    ptrModel->set(GRB_IntParam_MIPFocus, 1);
    ptrModel->set(GRB_IntParam_Threads, 4);
    ptrModel->set(GRB_StringParam_NodefileDir, "nodeDir");

    ptrModel->update();
    ptrModel->optimize();


    vetColSolX.setZero();


    int i = 0;
    const int num = instVrpTw.numClientes*instVrpTw.numClientes;

    std::cout<<"0 ";
    vetRoteG.push_back(0);

    do
    {
        int j;
        for(j=0; j < num; ++j)
        {
            if(vetVarX[getIndex(i, j, instVrpTw.numClientes)].get(GRB_DoubleAttr_X) >= 0.99)
            {
                vetColSolX[getIndex(i, j, instVrpTw.numClientes)] = 1.0;
                std::cout<<j<<" ";
                vetRoteG.push_back(j);
                break;
            }
        }

        i = j;
    }
    while(i != 0);

    std::cout<<"\n";

    cost = ptrModel->get(GRB_DoubleAttr_ObjVal);

    if(cost < -DW_DecompNS::TolObjSubProb)
        return true;
    else
        return false;

}
