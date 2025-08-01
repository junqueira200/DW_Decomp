#include <iostream>
#include "InstanceVRPTW.h"
#include "VrpTW_Decomp.h"
#include <filesystem>
#include "VrpTW_DecompLabeling.h"
#include "Cvrp_DecompLabeling.h"
#include "LabelingAlgorithm.h"
#include "MemoryPool.h"
#include "BranchAndPrice.h"
#include "Alarm.h"
#include "Test.h"
#include <bits/stdc++.h>
#include "Instancia.h"
#include "InputOutput.h"
#include "Ig.h"

//import teste;
// http://vrp.galgos.inf.puc-rio.br

using namespace InstanceVRPTW_NS;
using namespace VrpTW_DecompNS;
using namespace LabelingAlgorithmNS;
using namespace VrpTW_DecompLabelingNS;
using namespace BranchAndPriceNS;
using namespace SearchStrategyNS;
using namespace PrimalHeuristicNS;
using namespace BranchNS;
using namespace StatisticsNS;
using namespace TestNS;
using namespace InstanciaNS;
using namespace ParseInputNS;
using namespace IgNs;

void convertInstance(const InstanceVRPTW& instanceVrptw, Instancia& instancia);
void computeMeanMaxMin();

int main(int argc, const char **argv)
{


    try
    {
        std::feclearexcept(FE_OVERFLOW);
        std::feclearexcept(FE_UNDERFLOW);

        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";

        ParseInputNS::parseInput(argc, argv);
        output.setup();
        std::cout << "INST: " << input.strInst << " SEMENTE: " << RandNs::estado_ << " " << output.data << "";


        InstanceVRPTW instVrpTw;
        std::string strFile(input.strInstCompleto);

        std::filesystem::path p(strFile);
        std::string fileName = p.filename();

        if(fileName[0] == 'A' || fileName[0] == 'P')
            leInstanciaAugerat(strFile, instVrpTw);
        else
            leInstanciaSalomon(strFile, instVrpTw);

        ptr_instVrpG = &instVrpTw;

        instanciaG = Instancia(instVrpTw.numClientes, instVrpTw.numClientes-1, instVrpTw.numVeic);
        convertInstance(instVrpTw, instanciaG);

        //RouteHash routeHash;
        //enumerateRoutes(instVrpTw, 5, routeHash);
        //return 0;

        //getSubInstancia(15, instVrpTw);
        ArrayResources vetMaxResources;
        vetMaxResources[1] = instVrpTw.capVeic;

        VrpLabelingSubProb vrpLabelingSubProb(instVrpTw, instVrpTw.sumDist(), vetMaxResources);
        //Cvrp_DecompLabelingNS::CvrpLabelingSubProb vrpLabelingSubProb(instVrpTw, instVrpTw.sumDist());


        /*
        ArrayResources vetMaxResouces;
        vetMaxResouces[1] = instVrpTw.capVeic;
        Label label;
        label.typeLabel = Backward;
        label.active = true;
        label.cust   = 5;
        label.vetResources[0] = 0.0;
        label.vetResources[1] = 149.0;

        Index index = vrpLabelingSubProb.labelingData.getListOfIndexForMerge(label, vetMaxResouces);

        std::printf("Start: %d; %d\n", index.start(0), index.start(1));
        std::printf("End: %d; %d\n", index.end(0), index.end(1));

        return 0;
        */



        int option = input.labelingType;
        if(option == 1)
        {
            vrpLabelingSubProb.setTypeLabelToBackward();
            std::cout<<"Seting typeLabel to backward\n";
        }
        else if(option == 0)
            std::cout<<"Seting typeLabel to forward\n";
        else
        {
            vrpLabelingSubProb.setTypeLabelToBidirectional();
            std::cout<<"Seting typeLabel to bidirectional\n";
        }


        GRBEnv grbEnv;
        GRBModel model(grbEnv);
        criaMestre(instVrpTw, model);
        //vrpLabelingSubProb.createMaster(model);

        double distVarA = 100*std::max(somaDist(instVrpTw), 0.0);
        //VrpSubProb vrpSubProb(grbEnv, instVrpTw);

        DW_DecompNS::AuxData auxVectors;
        auxVectors.vetPairSubProb.push_back(std::make_pair(0, instVrpTw.numClientes * instVrpTw.numClientes));
        //auxVectors.vetPairSubProb.push_back(std::make_pair(0, ((instVrpTw.numClientes * (instVrpTw.numClientes-1))/2)));

        for(std::pair<int,int>& pair:auxVectors.vetPairSubProb)
        {
            std::printf("%d; %d\n", pair.first, pair.second);
        }

        setAlarm(30.0*60); // 30 min timer

        vetValueOfReducedCostsG.reserve(99999999);

        std::cout << "Cria decompNode\n";
        DW_DecompNS::DW_DecompNode decompNode(grbEnv, model, distVarA, (DW_DecompNS::SubProb*)&vrpLabelingSubProb, 1,
                                              auxVectors);

        std::cout<<"Num Vars Model: "<<decompNode.uRmlp->get(GRB_IntAttr_NumVars);
        std::cout<<"\nNum Vars vetVarLambdaCol:"<<decompNode.vetVarLambdaCol.size()<<"\n";

        input.comprimentoAlturaIguais1 = true;
        input.numItIG = 2;
        int numRotas = 0;
        for(int i=0; i < 50; ++i)
        {
            SolucaoNS::Solucao best(instanciaG);
            metaheuristicaIg(best);

            //std::cout<<"IG: "<<best.distTotal<<"\n";

            for(int j=0; j < (int)best.vetRota.size(); ++j)
            {
                Eigen::VectorXd vetX(decompNode.info.numVarMaster);
                vetX.setZero();
                SolucaoNS::Rota& rota = best.vetRota[j];


                for(int t=0; t < (int)rota.numPos-1; ++t)
                {

                    int cliI = rota.vetRota[t];
                    int cliJ = rota.vetRota[t+1];
                    //std::cout<<"\t"<<cliI<<", "<<cliJ<<"\n";
                    int index = VrpTW_DecompLabelingNS::getIndex(cliI, cliJ, instVrpTw.numClientes);
                    vetX[index] = 1;
                }
                auxVectors.vetColConvCooef[0] = 1;
                decompNode.addColumnX(rota.distTotal, numRotas, auxVectors, vetX);
                numRotas += 1;

                //std::cout<<"\n\n";
            }
        }

        decompNode.uRmlp->update();

        std::cout<<"Num Vars Model: "<<decompNode.uRmlp->get(GRB_IntAttr_NumVars);
        std::cout<<"\nNum Vars vetVarLambdaCol:"<<decompNode.vetVarLambdaCol.size()<<"\n";
        std::cout<<"NumRotas: "<<numRotas<<"\n\n";



        decompNode.rhsConv = instVrpTw.numVeic;
        //DepthFirst depthFirst;
        MinFuncObj minFuncObj;
        SimpleDiving simpleDiving;
        //StrongBranch branch;
        SimpleStrongBranch branch;

        StatisticsData statisticD;

        Eigen::VectorXd vetSol;

        CapacityCut capacityCut(instVrpTw, 5000, 5000, 0.00001);

        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";
        std::cout<<"Cost Var A.: "<<distVarA<<"\n\n";

        vetSol = branchAndPrice(decompNode,
                                auxVectors,
                                (SearchDataInter*)&minFuncObj,
                                (PrimalHeuristicInter*)&simpleDiving,
                                (BranchInter*)&branch,
                                nullptr,//(RobustCutGenerator*)&capacityCut,
                                statisticD);

        computeMeanMaxMin();

/*        Eigen::MatrixXi matSol(instVrpTw.numVeic, instVrpTw.numClientes/2);
        matSol.setZero();

        int pos = 1;
        int i   = 0;
        int r   = 0;

        for(int j=1; j < instVrpTw.numClientes; ++j)
        {
            if(vetSol[VrpTW_DecompNS::getIndex(i, j, instVrpTw.numClientes)] >= 0.99)
            {
                matSol.coeffRef(r, pos) = j;
                r += 1;
            }

        }

        for(r=0; r < instVrpTw.numVeic; ++r)
        {
            pos = 2;
            i = matSol(r, 1);

            while(i != 0)
            {
                for(int j=0; j < instVrpTw.numClientes; ++j)
                {
                    if(vetSol[VrpTW_DecompNS::getIndex(i, j, instVrpTw.numClientes)] >= 0.99)
                    {
                        matSol.coeffRef(r, pos) = j;
                        pos += 1;
                        i = j;
                        break;
                    }
                }
            }

        }

        for(int k=0; k < vetSol.size(); ++k)
        {
            if(vetSol[k] >= 0.99)
            {
                int i = k/instVrpTw.numClientes;
                int j = k%instVrpTw.numClientes;

                std::cout<<"("<<i<<" "<<j<<"), ";
            }
        }

        std::cout<<"\n\n"<<matSol<<"\n";*/

        statisticD.inst = fileName;
        statisticD.numNodes = instVrpTw.numClientes;
        std::string extraHead = "Alg";
        std::string extraCont;

        if(vrpLabelingSubProb.typeLabel == LabelingAlgorithmNS::AlgForward)
            extraCont += "forward";
        else if(vrpLabelingSubProb.typeLabel == LabelingAlgorithmNS::AlgBackward)
            extraCont += "backward";
        else
            extraCont += "bidirectional";

        writeToFile(statisticD, "result.csv", extraHead, extraCont);

        //decompNode.columnGeneration(auxVectors);


        std::cout<<"sizeof(Label): "<<sizeof(Label)<<"\n";
        std::cout<<"vetNumSteps r0: "<<vrpLabelingSubProb.labelingData.vetNumSteps[0]<<"\n";
        std::cout<<"start r0: "<<vrpLabelingSubProb.labelingData.vetStepSize[0].start<<"\n";
        std::cout<<"vetNumSteps r1: "<<vrpLabelingSubProb.labelingData.vetNumSteps[1]<<"\n";
        std::cout<<"numClien: "<<instVrpTw.numClientes<<"\n";

        std::cout<<"DistMin: "<<minDistG<<"\nDistMax: "<<maxDistG<<"\n\n";

        if((bool)std::fetestexcept(FE_OVERFLOW))
            std::cout << "Overflow flag after: " << (bool)std::fetestexcept(FE_OVERFLOW) << std::endl;
        if((bool)std::fetestexcept(FE_UNDERFLOW))
            std::cout << "Underflow flag after: " << (bool)std::fetestexcept(FE_UNDERFLOW) << std::endl;

        /*
        std::cout<<"0 < -1("<<doubleLess(0, -1.0, FloatEp)<<")\n";
        std::cout<<"-1 == 0("<<doubleEqual(0, -1.0, FloatEp)<<")\n";
        std::cout<<"0 <= -1("<<doubleLessEqual(0, -1.0, FloatEp)<<")\n";
        */

    }
/*    catch(char const* str)
    {
        std::cout<<"catch(char* ):\n";
        std::printf("%s", str);
        std::cout<<"\n\n";
    }*/
    catch(GRBException &e)
    {
        std::cout<<"GRBException:\n"<<e.getMessage()<<"\n";
        std::cout<<"Code: "<<e.getErrorCode()<<"\n\n";
    }

    /* ****************************************************************************************************************
     * ****************************************************************************************************************
     */

    return 0;
}

void convertInstance(const InstanceVRPTW& instanceVrptw, Instancia& instancia)
{

    instancia.vetItens = Vector<Item>(instancia.numItens);

    for(int i=0; i < instanceVrptw.numClientes; ++i)
    {
        for(int j=0; j < instanceVrptw.numClientes; ++j)
            instancia.matDist.get(i, j) = instanceVrptw.matDist(i, j);

        instancia.vetDemandaCliente[i] = instanceVrptw.vetClieDem[i];
        if(i != 0)
        {
            instancia.matCliItensIniFim.get(i, 0) = i-1;
            instancia.matCliItensIniFim.get(i, 1) = i-1;

            instancia.vetItens[i-1].set(1.0, 1.0, 1.0);
            instancia.vetItens[i-1].peso = instanceVrptw.vetClieDem[i];
            instancia.vetPesoItens[i-1] = instanceVrptw.vetClieDem[i];
            instancia.vetItens[i-1].volume = 1.0;
        }
    }

    instancia.packing = false;
    instancia.maxNumItensPorClie = 1;

    instancia.vetDimVeiculo[0] = 5000;
    instancia.vetDimVeiculo[1] = 5000;
    instancia.veicCap = instanceVrptw.capVeic;

}

void computeMeanMaxMin()
{
    double min, max, mean, median;

    min = MaxFloatType;
    max = MinFloatType;
    mean = 0.0;

    for(double val:vetValueOfReducedCostsG)
    {
        min = std::min(min, val);
        max = std::max(max, val);
        mean += val;
    }

    std::sort(vetValueOfReducedCostsG.begin(), vetValueOfReducedCostsG.end());

    mean = mean/(FloatType)vetValueOfReducedCostsG.size();

    size_t size = vetValueOfReducedCostsG.size();
    if(size % 2 == 0)
    {
        median = (vetValueOfReducedCostsG[size/2] + vetValueOfReducedCostsG[(size/2)+1])/2.0;
    }
    else
        median = vetValueOfReducedCostsG[size/2];

    std::printf("Min: %.2f\nMax: %.2f\nMean: %.2f\nMedian: %.2f\n\n", min, max, mean, median);
}









