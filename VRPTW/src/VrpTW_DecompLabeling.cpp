//
// Created by igor on 24/11/24.
//
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
    const double sumDist = instVrpTw->sumDist();
    const int sumDem     = instVrpTw->sumDem();

    //vetStepSize[0].stepSize = 400;
    vetStepSize[0].stepSize = 50; // 10
    vetStepSize[0].start    = -startDist;
    vetStepSize[0].end      = startDist;

    //vetStepSize[1].stepSize = 5;
    vetStepSize[1].stepSize = 5;
    vetStepSize[1].start    = 0;
    vetStepSize[1].end      = instVrpTw->capVeic;

    labelingData = LabelingAlgorithmNS::LabelingData(vetStepSize, 2, instVrpTw->numClientes+1);


    vetMatResCost = LabelingAlgorithmNS::VetMatResCost(2);
    vetMatResCost[0] = instVrpTw->matDist;
    vetMatResCost[1].resize(instVrpTw->numClientes+1, instVrpTw->numClientes+1);
    vetMatResCost[1].setZero();


/*    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        for(int j=0; j < instVrpTw->numClientes+1; ++j)
        {

            if(i == j || (i==0 && j == instVrpTw->numClientes) || (j==0 && i == instVrpTw->numClientes))
                continue;

            vetMatResCost[0](i, j) += -sumDist;
        }
    }*/


    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        for(int j=0; j < (instVrpTw->numClientes+1); ++j)
        {
            if(i == j)
                continue;

            vetMatResCost[1](i, j) = instVrpTw->vetClieDem[j];
        }
    }


    vetVetResBound = LabelingAlgorithmNS::VetVetResBound(2);

    vetVetResBound[0].resize(instVrpTw->numClientes+1);
    vetVetResBound[1].resize(instVrpTw->numClientes+1);

    double distTotal = instVrpTw->sumDist();

    Bound bound0;
    bound0.lowerBound = -std::numeric_limits<double>::infinity();
    bound0.upperBound =  std::numeric_limits<double>::infinity();

    Bound bound1;
    bound1.lowerBound = 0;
    bound1.upperBound = instVrpTw->capVeic;

    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        vetVetResBound[0][i] = bound0;
        if(i > 0 && i < instVrpTw->numClientes)
            bound1.lowerBound = instVrpTw->vetClieDem[i];

        else
            bound1.lowerBound = 0.0;

        vetVetResBound[1][i] = bound1;
    }

    ngSet = NgSet(instVrpTw->numClientes+1, NgSetSize);
    ngSet.setNgSets(instVrpTw->matDist);
    ngSet.active = true;


    Eigen::VectorXd vetX(instVrpTw->numClientes*instVrpTw->numClientes);
    vetX.setZero();

    //forwardLabelingAlgorithm(2, instVrpTw->numClientes + 1, vetMatResCost, vetVetResBound, instVrpTw->numClientes,
    //                         ngSet, labelingData, vetX, 0, 0);

}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{

    //GRBLinExpr linExpr;
    //rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");

}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::resolveSubProb(const Eigen::VectorXd &vetC,
                                                               const Eigen::RowVectorXd &vetRowPi,
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
                                                               double &redCost,
                                                               double constPiValue,
                                                               const VectorI &vetDelVar)
{

    double pi0 = vetRowPi[0];
    vetMatResCost[0].setZero();

    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        for(int j=1; j < instVrpTw->numClientes; ++j)
        {
            if(i == j)
                continue;

            vetMatResCost[0](i, j) = instVrpTw->matDist(i, j) - vetRowPi[j-1];
        }

        if(i != 0)
        {
            vetMatResCost[0](i, 0) = instVrpTw->matDist(i, 0);
            vetMatResCost[0](i, instVrpTw->numClientes) = instVrpTw->matDist(i, 0);// - pi0;
        }
    }

    for(int varId:vetDelVar)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        if(j == 0)
            j = instVrpTw->numClientes;

        vetMatResCost[0](i, j) = std::numeric_limits<double>::infinity();
        //std::cout<<varId<<"\n";
    }

    //std::cout<<"Custo Reduzido: \n"<<vetMatResCost[0]<<"\n\n";
    //std::cout<<"Peso: "<<vetMatResCost[1]<<"\n\n";

    int it = 0;
    double maxDist;

    ngSet.active = true;

    for(int i=1; i <= 16; i += 5)
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
                                               constPiValue,
                                               i,
                                               true,
                                               maxDist,
                                               redCost);

        it += 1;
        if(custoRedNeg)
            break;

    }

    if(!custoRedNeg)
    {
        matColX.setZero();
        custoRedNeg = forwardLabelingAlgorithm(2,
                                               instVrpTw->numClientes+1,
                                               vetMatResCost,
                                               vetVetResBound,
                                               instVrpTw->numClientes,
                                               ngSet,
                                               labelingData,
                                               matColX,
                                               numSol,
                                               0.0,
                                               -1,
                                               true,
                                               maxDist,
                                               redCost);


        if(!custoRedNeg)
        {
            matColX.setZero();
            ngSet.active = false;
            //std::cout<<"Ultimo forwardLabelingAlgorithm\n\n";
            custoRedNeg = forwardLabelingAlgorithm(2,
                                                   instVrpTw->numClientes+1,
                                                   vetMatResCost,
                                                   vetVetResBound,
                                                   instVrpTw->numClientes,
                                                   ngSet,
                                                   labelingData,
                                                   matColX,
                                                   numSol,
                                                   0.0,
                                                   -1,
                                                   true,
                                                   maxDist,
                                                   redCost);
        }
    }


    return 0;
}
