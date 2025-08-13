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
#include "LowerBound.h"

using namespace LabelingAlgorithmNS;
using namespace TestNS;
using namespace BranchAndPriceNS;

VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb(InstanceVRPTW_NS::InstanceVRPTW &instVrpTw_, double startDist,
                                                               LabelingAlgorithmNS::ArrayResources& vetMaxResouces)
{
    if(instVrpTw_.numClientes > NumMaxCust)
    {
        PRINT_DEBUG("", "ERRO; Const NumMaxCust("<<NumMaxCust<<") < NumCust("<<instVrpTw_.numClientes<<")\n\t");
        throw "ERRO";
    }

    instVrpTw = &instVrpTw_;
    double mult = 0.5;
    double numSteps = 1.0; // 2

    //vetStepSize[0].stepSize = 400;
    vetStepSize[0].stepSize = 5000000.0;//  10 // 5 //((2.0*mult)*startDist)/numSteps; // 1700
    vetStepSize[0].start    = -70;// -200 -50 // (FloatType)-mult*startDist;  // 1.0
    vetStepSize[0].end      = 400;//  200 50 // (FloatType) mult*startDist; // 1.0



    vetStepSize[1].stepSize = instVrpTw->capVeic/numSteps;  //50
    vetStepSize[1].start    = 0;
    vetStepSize[1].end      = (FloatType)instVrpTw->capVeic;


    labelingData  = LabelingAlgorithmNS::LabelingData(vetStepSize, 2, instVrpTw->numClientes+1, vetMaxResouces);
    labelingData.matDist = instVrpTw->matDist;
    vetMatResCostForward  = Vet3D_ResCost(instVrpTw->numClientes+1, instVrpTw->numClientes+1, 2);
    vetMatResCostForward.setVal(0.0);

    vetMatResCostBackward = Vet3D_ResCost(instVrpTw->numClientes+1, instVrpTw->numClientes+1, 2);
    vetMatResCostBackward.setVal(0.0);


    /*double val = vetStepSize[0].start;
    for(int i=0; i < labelingData.vetNumSteps(0); ++i)
    {
        //std::cout<<"["<<val<<";"<<val+vetStepSize[0].stepSize<<") ";
        val += vetStepSize[0].stepSize;
    }

    //std::cout<<"\n\n";

    val = vetStepSize[1].start;
    for(int i=0; i < labelingData.vetNumSteps(1); ++i)
    {
        //std::cout<<"["<<val<<";"<<val+vetStepSize[1].stepSize<<") ";
        val += vetStepSize[1].stepSize;
    }
    */
    //std::cout<<"\n";

    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        for(int j=0; j < (instVrpTw->numClientes+1); ++j)
        {
            if(i == j)
                continue;

            vetMatResCostForward.get(i, j, 1)  = (FloatType)instVrpTw->vetClieDem[j];
            vetMatResCostBackward.get(j, i, 1) = (FloatType)instVrpTw->vetClieDem[j];
        }

    }

    //std::cout<<"\n"<<instVrpTw->vetClieDem<<"\n\n";
    //PRINT_EXIT();

    //vetVetResBound = LabelingAlgorithmNS::VetVetResBound(2);

    vetVetResBound = LabelingAlgorithmNS::MatBoundRes(instVrpTw->numClientes+1, 2);

    //vetVetResBound[0].resize(instVrpTw->numClientes+1);
    //vetVetResBound[1].resize(instVrpTw->numClientes+1);


    Bound bound0;
    bound0.lowerBound = -std::numeric_limits<FloatType>::infinity();
    bound0.upperBound =  std::numeric_limits<FloatType>::infinity();

    Bound bound1;
    bound1.lowerBound = 0;
    bound1.upperBound = (FloatType)instVrpTw->capVeic;

    for(int i=0; i < instVrpTw->numClientes+1; ++i)
    {
        vetVetResBound(i, 0) = bound0;
        /*
        if(i > 0 && i < instVrpTw->numClientes)
            bound1.lowerBound = (FloatType)instVrpTw->vetClieDem[i];
        else
            bound1.lowerBound = 0.0;
        */
        vetVetResBound(i, 1) = bound1;
    }

    ngSet = NgSet(instVrpTw->numClientes+1, NgSetSize);
    ngSet.setNgSets(instVrpTw->matDist);
    ngSet.active = true;

    //enumerateRoutes(*instVrpTw, instVrpTw->numClientes-1, routeHash);


}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{


    GRBLinExpr linExpr;
    GRBVar a = rmlp.addVar(0, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
    //linExpr += -a;

    rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");


}

VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb()
{

    //enumerateRoutes(*instVrpTw, instVrpTw->numClientes-1, routeHash);

}

bool VrpTW_DecompLabelingNS::VrpLabelingSubProb::checkEnumeratedRoutesFinal(Eigen::RowVectorXd &vetRowPi)
{
    return false;

    cleanVetRouteG();
    std::cout<<"Start checkEnumeratedRoutes\n\n";
    int t = 0;
    bool ret = false;
    for(const Route& route:routeHash)
    {
        FloatType redCost = computeReducedCost(route, vetRowPi);
        if(t < 10)
            std::cout<<redCost<<"\n";
        t += 1;
        if(redCost < -DW_DecompNS::TolObjSubProb)
        {
            std::cout<<"Route("<<route.vetRoute<<") have a NEGATIVE("<<redCost<<") reduced cost!\n";
            if(redCost <= -1.0)
            {
                addToVetRoutesG(route.vetRoute);
                ret = true;
            }
        }
    }

    return ret;

}

bool  VrpTW_DecompLabelingNS::VrpLabelingSubProb::
      checkEnumeratedRoutesMid(Eigen::RowVectorXd &vetRowPi, Eigen::MatrixXd &matColX, int &numSol)
{

    //Remove all remaining labels
    labelingData.flushLabel();
    static LabelHeap labelHeap(1);

    for(int l=0; l < numSol; ++l)
    {
        Label* label = getLabel();

        // Seting label data
        label->vetRoute[0] = 0;
        label->tamRoute    = 1;
        label->typeLabel = LabelingAlgorithmNS::Forward;
        label->vetResources.setZero();
        label->vetResources[0] = -vetRowPi[0];
        label->typeLabel = Forward;
        int lastCust = 0;

        // Seting route
        do
        {
            int cliJ = -1;
            for(int i=0; i < instVrpTw->numClientes; ++i)
            {
                int index = getIndex(lastCust, i, instVrpTw->numClientes);
                if(matColX(index, l))
                {
                    cliJ = i;
                    break;
                }
            }

            if(cliJ == -1)
            {
                std::cout<<"Erro!; diden't find arc ("<<lastCust<<", _)\n";
                PRINT_EXIT();
            }

            label->vetResources[0] += instVrpTw->matDist(lastCust, cliJ) - vetRowPi[cliJ+1];
            label->vetResources[1] += instVrpTw->vetClieDem[cliJ];
            label->vetRoute[label->tamRoute] = cliJ;
            label->tamRoute += 1;
            lastCust = cliJ;
        }
        while(lastCust != 0);

        label->vetRoute[label->tamRoute-1] = instVrpTw->numClientes;
        label->cust = instVrpTw->numClientes;
        label->active = true;

        label->i = 0;
        label->j = 0;
        Bucket* bucket = labelingData.getBucket(label);
        bucket->addLabel(label);
    }

    Bucket& bucket = labelingData.vetMatBucketForward[instVrpTw->numClientes].mat(0, 0);
    if(bucket.sizeVetPtrLabel != numSol)
    {
        std::cout<<"Error, bucket size("<<bucket.sizeVetPtrLabel<<") is different from "<<numSol<<"\n\n";
        PRINT_EXIT();
    }

    bool addRoute = false;
    cleanVetRouteG();
    // Go through routes in hashRoutes
    for(const TestNS::Route& route:routeHash)
    {
        FloatType redCost = computeReducedCost(route, vetRowPi);
        if(redCost < -DW_DecompNS::TolObjSubProb)
        {
            Label* label = getLabel();
            label->cust = instVrpTw->numClientes;
            label->vetResources[0] = redCost;
            label->vetResources[1] = route.demand;
            label->active = true;
            label->typeLabel = Forward;

            writeNgSet(label, ngSet);
            int pos = -1;
            Bucket* bucket = dominanceIntraBucket(instVrpTw->numClientes, label, labelingData, &labelHeap, 2,
                                                  instVrpTw->numClientes, pos);

            if(bucket != nullptr)
            {
                //std::cout<<"Error!!\nRoute("<<route.vetRoute<<") have a negative reduced cost ("<<redCost<<
                //           ") and was inserted into the bucket!\n\n";

                //Vector<VectorI> vetRoute;
                //vetRoute.push_back(route.vetRoute);
                //startGlobalMemory(vetRoute);
                //return true;
                addToVetRoutesG(route.vetRoute);
                addRoute = true;

            }
        }
    }

    return addRoute;
}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::
    resolveSubProb(const Eigen::VectorXd &vetC, Eigen::RowVectorXd &vetRowPi, GRBModel &mestre, int itCG,
                   bool &custoRedNeg, void *data, const int iniConv, int indSubProb, Eigen::VectorXd &vetCooefRestConv,
                   const std::pair<int, int> &pairSubProb, Eigen::MatrixXd &matColX, int &numSol,
                   Eigen::Array<double, 1, DW_DecompNS::NumMaxSolSubProb>& vetRedCost, double constPiValue,
                   const VectorI &vetVar0, const VectorI &vetVar1, DW_DecompNS::PhaseStatus phaseStatus, bool exact)
{
    //std::cout<<"star value: "<<constPiValue<<"\n\n";
    //static Eigen::VectorXd vetDist(instVrpTw->numClientes+1);
    //LowerBoundNS::getDistLowerBound(vetMatResCostForward, vetDist, instVrpTw->numClientes, &labelingData);

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

            FloatType value0 = 0.0;
            FloatType value1 = 0.0;


            if(phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
            {
                value0 =  - (FloatType)vetRowPi[j+1];
                value1 =  - (FloatType)vetRowPi[i+1];
            }
            else
            {
                value0 = (FloatType)instVrpTw->matDist(i, j) - (FloatType)vetRowPi[j+1];
                value1 = (FloatType)instVrpTw->matDist(j, i) - (FloatType)vetRowPi[i+1];
            }

             vetMatResCostForward.get(i, j, 0)  = value0;
             vetMatResCostBackward.get(i, j, 0) = value0;
        }

        if(i != 0)
        {

            FloatType value = 0.0;


            //vetMatResCost[0](i, 0) = instVrpTw->matDist(i, 0);
            if(phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
                value = - (FloatType)vetRowPi[1];
            else
                value = (FloatType)instVrpTw->matDist(i, 0) - (FloatType)vetRowPi[1];

              vetMatResCostForward.get(i, instVrpTw->numClientes, 0)  = value;
              vetMatResCostBackward.get(i, instVrpTw->numClientes, 0) = value;
        }
    }

    for(int i=0; i < instVrpTw->numClientes; ++i)
    {
        if(!doubleEqual(vetMatResCostBackward(0, i, 0), vetMatResCostForward(0, i, 0)))
        {
            std::cout<<"backward("<<vetMatResCostBackward(0, i, 0)<<") != forward("<<vetMatResCostForward(0, i, 0)<<")\n";
            exit(-1);
        }
    }

    //std::cout<<"("<<instVrpTw->numClientes<<",1): "<<vetMatResCostBackward(instVrpTw->numClientes, 1, 0)<<"\n";
    //std::cout<<"(1,"<<instVrpTw->numClientes<<"): "<<vetMatResCostBackward(1, instVrpTw->numClientes, 0)<<"\n";


    for(int varId:vetVar0)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        if(j == 0)
            j = instVrpTw->numClientes;

        vetMatResCostForward.get(i, j, 0) = std::numeric_limits<FloatType>::infinity();
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

// TODO remove comments!
    /*
    if(!exact)
    {
        for(int i = 2; i <= 4; i += 2)
        {
            std::cout<<i<<": ";
            //std::cout<<"forwardLabelingAlgorithm: "<<i<<"\n\n";
            matColX.setZero();
            custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward,
                                                 vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes, ngSet,
                                                 labelingData, matColX, numSol, (FloatType) constPiValue, i, true,
                                                 maxDist, vetRedCostFT, exact, typeLabel, true, nullptr, vetSortRoute,
                                                 setSortRoute, 10);

            it += 1;
            if(custoRedNeg)
                break;

        }
    }
    */


    if(!custoRedNeg && !exact)
    {
        for(int i = 2; i <= 4; i += 2)
        {
            std::cout<<i<<": ";
            //std::cout<<"forwardLabelingAlgorithm: "<<i<<"\n\n";
            matColX.setZero();
            numSol = 0;
            custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward,
                                                 vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes, ngSet,
                                                 labelingData, matColX, numSol, (FloatType) constPiValue, i, true,
                                                 maxDist, vetRedCostFT, exact, typeLabel, true, nullptr, vetSortRoute,
                                                 setSortRoute, MaxInt, MaxRedCost);

            it += 1;
            if(custoRedNeg && numSol >= DW_DecompNS::NumMaxSolSubProb/2)
                break;

        }
    }


    /*
    if(!custoRedNeg && !exact)
    {
        for(int i = 2; i <= 8; i += 2)
        {
            std::cout<<i<<": ";
            //std::cout<<"forwardLabelingAlgorithm: "<<i<<"\n\n";
            matColX.setZero();
            custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward,
                                                 vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes, ngSet,
                                                 labelingData, matColX, numSol, (FloatType) constPiValue, i, true,
                                                 maxDist, vetRedCostFT, exact, typeLabel, true, nullptr, vetSortRoute,
                                                 setSortRoute, MaxInt, MaxFloatType);

            it += 1;
            if(custoRedNeg)
                break;

        }
    }
    */

    if(!custoRedNeg || numSol < DW_DecompNS::NumMaxSolSubProb/2)
    {
        matColX.setZero();
        numSol = 0;
        custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward, vetMatResCostBackward,
                                             vetVetResBound, instVrpTw->numClientes, ngSet, labelingData, matColX,
                                             numSol, (FloatType)constPiValue, -1, true, maxDist, vetRedCostFT, exact,
                                             typeLabel, true, nullptr,  vetSortRoute, setSortRoute, MaxInt, MaxFloatType);
    }



    if(!custoRedNeg)
    {
        if(checkEnumeratedRoutesFinal(vetRowPi))
        {
            PRINT_EXIT();
        }
    }
    /*
    else


        if(checkEnumeratedRoutesMid(vetRowPi, matColX, numSol))
        {
            numSol = 0;
            custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward, vetMatResCostBackward,
                                                 vetVetResBound, instVrpTw->numClientes, ngSet, labelingData, matColX,
                                                 numSol, (FloatType)constPiValue, -1, true, maxDist, vetRedCostFT, exact,
                                                 typeLabel, true);

            //std::println("**************************************\n");
            //PRINT_EXIT();
        }
    }
    */

    /*
    if(!custoRedNeg)
    {
        std::cout<<"Change Alg Type\n";
        changeTypeAlg(typeLabel);

        custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1), vetMatResCostForward, vetMatResCostBackward,
                                             vetVetResBound, instVrpTw->numClientes, ngSet, labelingData, matColX,
                                             numSol, (FloatType)constPiValue, -1, true, maxDist, vetRedCostFT, exact,
                                             typeLabel);

        if(custoRedNeg)
        {
            std::cout<<"ERROR, alterar o algoritimo fez com que fosse gerado uma coluna com custo red. neg.\n";
        }

        changeTypeAlg(typeLabel);

    }
    */

    //redCost = (double)redCostFT;
    vetCooefRestConv[0] = 1;


    std::cout<<"*************\n\n";

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
    /*
    FloatType redCostTemp = 0.0;
    for(int j=0; j < numSol; ++j)
    {
        redCostTemp = constPiValue;
        //std::cout<<j<<": ";
        for(int i=0; i < (instVrpTw->numClientes*instVrpTw->numClientes); ++i)
        {
            if(matColX(i,j) != 1.0)
                continue;

            int ii = i/instVrpTw->numClientes;
            int jj = i%instVrpTw->numClientes;

            redCostTemp += vetMatResCostForward(ii, jj, 0);//instVrpTw->matDist(ii, jj);
            //std::cout<<"("<<ii<<","<<jj<<")["<<i<<"]; ";

        }
        //std::cout<<"\n\n";
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

    */

    for(int i=0; i < numSol; ++i)
        vetRedCost[i] = (double)vetRedCostFT[i];

    //std::cout<<"["<<minDistG<<", "<<maxDistG<<"]; numSol("<<numSol<<")\n";

    return 0;
}


bool VrpTW_DecompLabelingNS::exactPricing(const LabelingAlgorithmNS::Vet3D_ResCost&          vetMatResCost,
                                          const FloatType                                    startVal,
                                          Eigen::Matrix<double, -1, 1, Eigen::ColMajor>&     vetColSolX,
                                          InstanceVRPTW_NS::InstanceVRPTW&                           instVrpTw,
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
    /*
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
    */
    std::cout<<"\n";

    cost = ptrModel->get(GRB_DoubleAttr_ObjVal);

    if(cost < -DW_DecompNS::TolObjSubProb)
        return true;
    else
        return false;

}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToForward()
{
    typeLabel = LabelingAlgorithmNS::AlgForward;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToBackward()
{
    typeLabel = LabelingAlgorithmNS::AlgBackward;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToBidirectional()
{
    typeLabel = LabelingAlgorithmNS::AlgBidirectional;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::
     convertRouteIntoLabel(const TestNS::Route& route, LabelingAlgorithmNS::Label* label)
{
    //label->vetResources.setZero();
    label->bitSetNg = 0;
    for(int i=0; i < (int)route.vetRoute.size(); ++i)
    {
        int cliI = route.vetRoute[i];
        label->vetRoute[i] = cliI;

        if((i+1) < (int)route.vetRoute.size())
        {
            int cliJ = route.vetRoute[i+1];

            for(int r=0; r < 2; ++r)
            {
                label->vetResources[r] += vetMatResCostForward(cliI, cliJ, r);
            }

            if(ngSet.contain(cliJ, cliI))
                label->bitSetNg[cliJ] = 1;
        }
    }

    label->tamRoute = (int)route.vetRoute.size();
}


VrpTW_DecompLabelingNS::CapacityCut::CapacityCut(InstanceVRPTW_NS::InstanceVRPTW &instVrpTw, int dim_, int maxNoOfCuts_, double eps)
{
    capacity 	= instVrpTw.capVeic;
    dim      	= dim_;
    maxNoOfCuts = maxNoOfCuts_;
    numCust  	= instVrpTw.numClientes - 1; // Don't count with the deposit
    edgeSize 	= ((numCust+1)*(numCust))/2;
    epsForIntegrality = eps;

    std::printf("EdgeSize: %d\n\n", edgeSize);

    edgeHead = new int[edgeSize];
    edgeTail = new int[edgeSize];
    edgeX    = new double[edgeSize];
    demand   = new int[numCust+1];

    //std::cout<<"Demand: ";
    for(int i=0; i <= numCust; ++i)
    {
        demand[i] = instVrpTw.vetClieDem[i];
        //std::printf("%d ", demand[i]);
    }

    //std::println();

    CMGR_CreateCMgr(&cutsCMP, dim);
    CMGR_CreateCMgr(&oldCutsCMP, dim);

    bool createMap = true;
    //if(createMap)
    mapArcToIndex = new std::map<std::pair<int,int>, int>();

    int next = 0;
    for(int i=0; i < numCust+1; ++i)
    {
        for(int j=(i+1); j < numCust+1; ++j)
        {
            if(next >= edgeSize)
            {
                std::printf("Error, next(%i) >= edgeSize(%i)\n", next, edgeSize);
                PRINT_EXIT();
            }

            if(createMap)
            {	//std::printf("(%d; %d): %i\n", i, j, next);
                (*mapArcToIndex)[{i,j}] = next;
            }

            edgeHead[next] = i;
            edgeTail[next] = j;
            //std::printf("Add: (%d;%d)\n", i, j);

            next += 1;
        }
    }

    //std::printf("\n\n");

    if(next != edgeSize)
        std::printf("next(%d) != edgeSize(%d)\n\n", next, edgeSize);
    if(createMap)
        list = Eigen::VectorXi(numCust+1);
}

VrpTW_DecompLabelingNS::CapacityCut::~CapacityCut()
{
    delete []edgeHead;
    delete []edgeTail;
    delete []edgeX;

    CMGR_FreeMemCMgr(&cutsCMP);
    CMGR_FreeMemCMgr(&oldCutsCMP);
}
int VrpTW_DecompLabelingNS::CapacityCut::operator()(DW_DecompNS::DW_DecompNode& decompNode)
{
    std::printf("*************************************************************************\n");
    std::printf("***************************CAPACITY_CUT**********************************\n\n");

    // Convert the arc solution for the edge one
    convertArcSolution(decompNode.vetSolX);

    // Separate the capacity cuts
    CAPSEP_SeparateCapCuts(numCust, demand, capacity, edgeSize, edgeTail, edgeHead, edgeX, oldCutsCMP, maxNoOfCuts,
                           epsForIntegrality, &integerAndFeasible, &maxViolation, cutsCMP);

    list.setConstant(-1);

    // Runs through the cuts
    for(int i=0; i < cutsCMP->Size; ++i)
    {

        std::printf("Add cut %d: \n\n\t", i);

        int listSize = 0;
        // Recover the customers
        for(int j=0; j < cutsCMP->CPL[i]->IntListSize; ++j)
        {
            list[listSize] = cutsCMP->CPL[i]->IntList[j];
            std::printf("%d ", list[listSize]);
            listSize += 1;
        }

        double rhs = cutsCMP->CPL[i]->RHS;
        std::printf("< %f", rhs);
        RobustCut robustCut(decompNode.vetSolX.size());
        robustCut.sense = '<';
        robustCut.rhs = rhs;

        // Add the arcs of the clique of nodes of the list
        createInducedSubGraphArcs(listSize, robustCut);
        addMasterCut(robustCut, decompNode, i, false);
    }

    // Move cuts in cutsCMP for oldCutsCMP
    for(int i=0; i < cutsCMP->Size; ++i)
        CMGR_MoveCnstr(cutsCMP, oldCutsCMP, i, 0);

    int size = cutsCMP->Size;
    cutsCMP->Size = 0;

    return size;
};

void VrpTW_DecompLabelingNS::CapacityCut::convertArcSolution(const Eigen::VectorXd& vetX)
{

    int numNodes = numCust+1;

    for(int i=0; i < edgeSize; ++i)
        edgeX[i] = 0.0;

    int index;

    for(int i=0; i < numNodes; ++i)
    {
        for(int j=0; j < numNodes; ++j)
        {
        //reverseIndex(t, numNodes, i, j);
            if(i == j)
                continue;

            index = (*mapArcToIndex).at(getEdge(i,j));
            edgeX[index] += vetX[getIndex(i, j, numNodes)];

            if(edgeX[index] > 1.0)
            {
                if(j > i)
                    std::cout<<"Edge: ("<<i<<", "<<j<<") have a value greter then 1\n";
                else
                    std::cout<<"Edge: ("<<j<<", "<<i<<") have a value greter then 1\n";
            }    
        }
    }

//    for(int i=0; i < edgeSize; ++i)
//        edgeX[i] = std::min(edgeX[i], 1.0);


    std::printf("edgeX: ");
    for(int i=0; i < edgeSize; ++i)
        std::printf("%.1f ", edgeX[i]);

    std::cout<<"\n\n";

}

std::pair<int,int> VrpTW_DecompLabelingNS::CapacityCut::getEdge(int i, int j)
{

    if(i < j)
        return {i,j};
    else
        return {j,i};
}

void  VrpTW_DecompLabelingNS::CapacityCut::createInducedSubGraphArcs(int listSize, BranchAndPriceNS::RobustCut& cut)
{

    int numNodes = numCust + 1;

    for(int i=0; i < listSize; ++i)
    {
        for(int j=0; j < listSize; ++j)
        {
            if(i == j)
                continue;

            int index0 = getIndex(i, j, numNodes);
            int index1 = getIndex(j, i, numNodes);

            cut.vetX.coeffRef(index0) = 1.0;
            cut.vetX.coeffRef(index1) = 1.0;
        }
    }
}









