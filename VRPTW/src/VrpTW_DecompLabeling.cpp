//
// Created by igor on 24/11/24.
//
#include "VrpTW_DecompLabeling.h"

using namespace LabelingAlgorithmNS;


VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb(InstanciaNS::InstVRP_TW &instVrpTw_)
{
    if(instVrpTw_.numClientes > NumMaxCust)
    {
        PRINT_DEBUG("", "ERRO; Const NumMaxCust("<<NumMaxCust<<") < NumCust("<<instVrpTw_.numClientes<<")\n\t");
        throw "ERRO";
    }

    instVrpTw = &instVrpTw_;
    const double sumDist = instVrpTw->sumDist();
    const int sumDem     = instVrpTw->sumDem();

    vetStepSize[0].stepSize = 100;
    vetStepSize[0].start    = -instVrpTw->sumDist();
    vetStepSize[0].end      = instVrpTw->sumDist();

    vetStepSize[1].stepSize = 5;
    vetStepSize[1].start    = 0;
    vetStepSize[1].end      = instVrpTw->capVeic;

    labelingData = LabelingAlgorithmNS::LabelingData(vetStepSize, 2, instVrpTw->numClientes+1);


    vetMatResCost = LabelingAlgorithmNS::VetMatResCost(2);
    vetMatResCost[0] = instVrpTw->matDist;
    vetMatResCost[1].resize(instVrpTw->numClientes+1, instVrpTw->numClientes+1);
    vetMatResCost[1].setZero();


    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        for(int j=0; j < instVrpTw->numClientes+1; ++j)
        {

            if(i == j || (i==0 && j == instVrpTw->numClientes) || (j==0 && i == instVrpTw->numClientes))
                continue;

            vetMatResCost[0](i, j) += -sumDist;
        }
    }


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
        vetVetResBound[1][i] = bound1;
    }

    ngSet = NgSet(instVrpTw->numClientes+1, NgSetSize);
    ngSet.setNgSets(instVrpTw->matDist);
    ngSet.active = false;


    Eigen::VectorXd vetX(instVrpTw->numClientes*instVrpTw->numClientes);
    vetX.setZero();
    forwardLabelingAlgorithm(2, instVrpTw->numClientes+1, vetMatResCost, vetVetResBound, instVrpTw->numClientes, ngSet, labelingData, vetX);

}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{

    GRBLinExpr linExpr;
    rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");

}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::resolveSubProb(const Eigen::VectorXd &vetC,
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

    double pi0 = vetRowPi[0];

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

//    std::cout<<"Custo Reduzido: \n"<<vetMatResCost[0]<<"\n\n";

    custoRedNeg = forwardLabelingAlgorithm(2,
                                           instVrpTw->numClientes+1,
                                           vetMatResCost,
                                           vetVetResBound,
                                           instVrpTw->numClientes,
                                           ngSet,
                                           labelingData,
                                           vetX);

    vetCooefRestConv[0] = 1;

    return 0;
}
