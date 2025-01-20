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
    vetStepSize[0].stepSize = 50; // 50
    vetStepSize[0].start    = (FloatType)-startDist;
    vetStepSize[0].end      = (FloatType)startDist;

    //vetStepSize[1].stepSize = 5;
    vetStepSize[1].stepSize = 5;
    vetStepSize[1].start    = 0;
    vetStepSize[1].end      = (FloatType)instVrpTw->capVeic;

    labelingData = LabelingAlgorithmNS::LabelingData(vetStepSize, 2, instVrpTw->numClientes+1);


    vetMatResCost = LabelingAlgorithmNS::VetMatResCost(2);
    vetMatResCost[0].resize(instVrpTw->matDist.rows(), instVrpTw->matDist.cols());
    vetMatResCost[1].resize(instVrpTw->numClientes+1, instVrpTw->numClientes+1);
    vetMatResCost[0].setZero();
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

            vetMatResCost[1](i, j) = (FloatType)instVrpTw->vetClieDem[j];
        }
    }


    vetVetResBound = LabelingAlgorithmNS::VetVetResBound(2);

    vetVetResBound[0].resize(instVrpTw->numClientes+1);
    vetVetResBound[1].resize(instVrpTw->numClientes+1);

    double distTotal = instVrpTw->sumDist();

    Bound bound0;
    bound0.lowerBound = -std::numeric_limits<FloatType>::infinity();
    bound0.upperBound =  std::numeric_limits<FloatType>::infinity();

    Bound bound1;
    bound1.lowerBound = 0;
    bound1.upperBound = (FloatType)instVrpTw->capVeic;

    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        vetVetResBound[0][i] = bound0;
        if(i > 0 && i < instVrpTw->numClientes)
            bound1.lowerBound = (FloatType)instVrpTw->vetClieDem[i];

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


    GRBLinExpr linExpr;
    //GRBVar a = rmlp.addVar(0, GRB_INFINITY, custoVarA, GRB_CONTINUOUS);
    //linExpr += -a;

    rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");


}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::resolveSubProb(const Eigen::VectorXd&     vetC,
                                                               const Eigen::RowVectorXd&  vetRowPi,
                                                               GRBModel&                  mestre,
                                                               int                        itCG,
                                                               bool&                      custoRedNeg,
                                                               void*                      data,
                                                               const int                  iniConv,
                                                               int                        indSubProb,
                                                               Eigen::VectorXd&           vetCooefRestConv,
                                                               const std::pair<int, int>& pairSubProb,
                                                               Eigen::MatrixXd&           matColX,
                                                               int&                       numSol,
                                                               double&                    redCost,
                                                               double                     constPiValue,
                                                               const VectorI&             vetDelVar)
{

//std::cout<<"constPiValue: "<<constPiValue<<"\n";
    static Eigen::VectorX<FloatType> vetRedCostFT(DW_DecompNS::NumMaxSolSubProb);

    vetMatResCost[0].setZero();
    //constPiValue += -vetRowPi[0];

    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        for(int j=0; j < instVrpTw->numClientes; ++j)
        {
            if(i == j)
                continue;

            vetMatResCost[0](i, j) = (FloatType)instVrpTw->matDist(i, j) - (FloatType)vetRowPi[j+1];
        }

        if(i != 0)
        {
            //vetMatResCost[0](i, 0) = instVrpTw->matDist(i, 0);
            vetMatResCost[0](i, instVrpTw->numClientes) = (FloatType)instVrpTw->matDist(i, 0) - (FloatType)vetRowPi[1];
        }
    }

    for(int varId:vetDelVar)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        if(j == 0)
            j = instVrpTw->numClientes;

        vetMatResCost[0](i, j) = std::numeric_limits<FloatType>::infinity();
        //std::cout<<varId<<"\n";
    }

    //std::cout<<"Custo Reduzido: \n"<<vetMatResCost[0]<<"\n\n";
    //std::cout<<"Peso: "<<vetMatResCost[1]<<"\n\n";

    int it = 0;
    FloatType maxDist;

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
                                                   (FloatType)constPiValue,
                                                   -1,
                                                   true,
                                                   maxDist,
                                                   vetRedCostFT);
        }
    }

    //redCost = (double)redCostFT;
    vetCooefRestConv[0] = 1;

    // Check if solution have a negative reduced cost
    FloatType redCostTemp = 0.0;
    for(int j=0; j < numSol; ++j)
    {
        redCostTemp = constPiValue;

        for(int i=0; i < (instVrpTw->numClientes*instVrpTw->numClientes); ++i)
        {
            if(doubleEqual(matColX(i,j), 0.0))
                continue;

            int ii = i/instVrpTw->numClientes;
            int jj = i%instVrpTw->numClientes;

            redCostTemp += (instVrpTw->matDist(ii, jj) - vetRowPi[jj+1]);
            //std::cout<<"("<<ii<<","<<jj<<"); ";

        }

        if(redCostTemp >= -DW_DecompNS::TolObjSubProb || !doubleEqual(redCostTemp, vetRedCostFT[j], 1E-4))
        {
            std::cout<<"\nERROR, custo reduzido calculado: ("<<redCostTemp<<") \n";
            std::cout<<"redCost: "<<vetRedCostFT[j]<<"\n\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

    }


    return 0;
}
