/*  *****************************************************************
 *  *****************************************************************
 *  File:    BranchAndPrice.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    04/01/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#include <fstream>
#include <filesystem>
#include "BranchAndPrice.h"
#include "PrimalHeuristic.h"
#include "Alarm.h"

using namespace DW_DecompNS;
using namespace SearchStrategyNS;
using namespace PrimalHeuristicNS;
using namespace BranchNS;
using namespace StatisticsNS;

int BranchAndPriceNS::getMostFractionVariable(const Eigen::VectorXd &vetSolX)
{
    int id = 0;
    double val = std::abs(std::abs(vetSolX[0]-int(vetSolX[0])) - 0.5);

    for(int i=1; i < int(vetSolX.size()); ++i)
    {
        double valTemp = std::abs(std::abs(vetSolX[i]-int(vetSolX[i])) - 0.5);
        if(valTemp < val)
        {
            val = valTemp;
            id = i;
        }

    }

    return id;
}

void BranchAndPriceNS::addMasterCut(const Cut &cut, DW_DecompNS::DW_DecompNode &decompNode, int num)
{
    if(cut.rhs != 0 && cut.sense != '<')
    {

        // Update info
        decompNode.info.numConstrsMaster += 1;
        decompNode.info.numVarRmlpPi += 1;

        // Add cut to matA
        auto &matA = decompNode.matA;
        matA.conservativeResize(matA.rows() + 1, matA.cols());
        matA.row(matA.rows() - 1) = cut.vetX;

        // Update columns coef
        GRBLinExpr linExpr;
        //GRBVar var = decompNode.uRmlp->addVar(0, GRB_INFINITY, decompNode.info.costA_Var, GRB_CONTINUOUS, "a_"+std::to_string(num));
        //double mult = 1.0;
        //if(cut.sense == '<')
        //    mult = -1.0;

        //linExpr += mult*var;

        //decompNode.uRmlp->update();
        //decompNode.uRmlp->optimize();
        GRBVar *varRmlp = decompNode.uRmlp->getVars();


        int numVars = decompNode.uRmlp->get(GRB_IntAttr_NumVars);

        if(numVars != int(decompNode.vetVarLambdaCol.size()))
        {
            std::cout << "Num Vars is wrong;\n\t Model: " << numVars << "; vetVarLamdaCol: "
                      << decompNode.vetVarLambdaCol.size() << "\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }


        // Runs through columns and compute their coefficients
        for(int i = 0; i < int(decompNode.vetVarLambdaCol.size()); ++i)
        {
            auto &col = *decompNode.vetVarLambdaCol[i];
            double coef = (cut.vetX.dot(col));

            if(coef != 0.0)
                linExpr += coef * varRmlp[i];

            //exit(-1);
        }

        std::cout << "add cut\n";
        // Add cut to rmlp
        decompNode.uRmlp->addConstr(linExpr, cut.sense, cut.rhs, "masterCut_" + std::to_string(num));
        delete[]varRmlp;
    }
    else
    {
        int varId = -1;
        for(int i=0; i < (int)cut.vetX.size(); ++i)
        {
            if(cut.vetX.coeff(i) != 0.0)
            {
                varId = i;
                break;
            }
        }

        if(varId == -1)
        {
            std::cout<<"ERROR, cut.vetX is zero!\nvetX: "<<cut.vetX<<"\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

        std::cout<<"X_"<<varId<<" <= 0\n";

        VectorI vetRmColumns;
        vetRmColumns.reserve(50);

        // Removes all columns with X_<varId> variable
        for(int i=0; i < (int)decompNode.vetVarLambdaCol.size(); ++i)
        {
            Eigen::Matrix<double, -1, 1> &vetX = (*decompNode.vetVarLambdaCol[i]);
            if(vetX[varId] != 0.0)
                vetRmColumns.push_back(i);
        }

        std::reverse(vetRmColumns.begin(), vetRmColumns.end());
        GRBVar *varRmlp = decompNode.uRmlp->getVars();

        for(int rmId:vetRmColumns)
        {
            decompNode.setVarLamdaCol.erase(SolXHash((*decompNode.vetVarLambdaCol[rmId])));
            decompNode.vetVarLambdaCol.erase(decompNode.vetVarLambdaCol.begin() + rmId);
            decompNode.uRmlp->remove(varRmlp[rmId]);
        }

        decompNode.vetDelVar.push_back(varId);

        delete []varRmlp;
    }
}

Eigen::VectorXd BranchAndPriceNS::branchAndPrice(DW_DecompNS::DW_DecompNode &cRootNode,
                                                 DW_DecompNS::AuxData &auxVectors,
                                                 SearchDataInter* searchD,
                                                 PrimalHeuristicInter* ptrPrimalH,
                                                 BranchInter* branch,
                                                 StatisticsData &statisticD)
{

   clock_t start = clock();

    //PrimalHeuristicInter* ptrPrimalH = new SimpleDiving;

    //std::list<DW_DecompNode*> listDecomNode;
    //listDecomNode.emplace_back(new DW_DecompNode(cRootNode));


    if(cRootNode.uRmlp->get(GRB_IntAttr_ModelSense) != GRB_MINIMIZE)
    {
        delete ptrPrimalH;
        std::cout<<"ERROR, model should be a min problem\n";
        PRINT_DEBUG("", "");
        throw "ERROR IN MODEL";
    }

    double lowerBound = -std::numeric_limits<double>::infinity();
    double upperBound =  std::numeric_limits<double>::infinity();
    auto vetX_best    = cRootNode.vetSolX;
    vetX_best.setZero();
    double gap = std::numeric_limits<double>::infinity();

    int numIt = -1;

    DW_DecompNode *rootNode = new DW_DecompNode(cRootNode);
    int status = rootNode->columnGeneration(auxVectors);

    //rootNode->uRmlp->write("root.lp");

    /*
    auto rootCp = new DW_DecompNode(*rootNode);

    rootCp->uRmlp->write("rootCP.lp");

    rootCp->columnGeneration(auxVectors);*/

    clock_t end = clock();
    statisticD.rootTime = double(end-start)/CLOCKS_PER_SEC;
    statisticD.rootLB   = rootNode->funcObj;

    if(status != StatusSubProb_Otimo)
    {
        delete ptrPrimalH;
        std::cout<<"Root node can't be soved\n";
        delete rootNode;
        vetX_best.setZero();
        return vetX_best;
    }

    delete rootNode;
    return vetX_best;

    searchD->insert(rootNode);
    rootNode = nullptr;

    int it = -1;

    //PRINT_DEBUG("", "");
    //throw "NAO EH ERRO";

    while(!searchD->empty() && gap > gapLimit && !alarm_stopG)
    {

        it += 1;


        //lowerBound = computeLowerBaound(listDecomNode);
        lowerBound = searchD->getMin();

        //DW_DecompNode* ptrDecomNode = listDecomNode.back();
        DW_DecompNode* ptrDecomNode = searchD->pop();
        std::cout<<"Diving Heuristic\n";
        DW_DecompNode* primalNode   = (*ptrPrimalH)(ptrDecomNode, auxVectors, upperBound);
        if(primalNode)
        {
            if(primalNode->funcObj < upperBound)
            {
                upperBound = primalNode->funcObj;
                vetX_best  = primalNode->vetSolX;
            }

            delete primalNode;
            primalNode = nullptr;
        }

        gap = computeGap(lowerBound, upperBound);
        end = clock();

        std::cout<<"it("<<it<<") \t LB("<<lowerBound<<") \t UB("<<upperBound<<") \t gap("<<gap<<"%) \t T("<<double(end-start)/CLOCKS_PER_SEC<<" s)\n";

        std::cout<<"Processando NO: "<<ptrDecomNode<<"\n\n";
        //listDecomNode.pop_back();

        // TODO Colocar antes da chamada da heuristica
        if(ptrDecomNode->funcObj > upperBound)
        {
            delete ptrDecomNode;
            continue;
        }

        if(isInteger(ptrDecomNode->vetSolX))
        {
            if(ptrDecomNode->funcObj < upperBound)
            {
                upperBound = ptrDecomNode->funcObj;
                vetX_best  = ptrDecomNode->vetSolX;
            }

            delete ptrDecomNode;
            continue;

        }

        if(alarm_stopG)
        {
            if(!searchD->empty())
                lowerBound = searchD->getMin();

            delete ptrDecomNode;
            break;
        }

        int varId = (*branch)(ptrDecomNode, auxVectors);
        double varValue = ptrDecomNode->vetSolX[varId];

        Cut cut;
        cut.vetX.resize(ptrDecomNode->vetSolX.size());
        cut.vetX.coeffRef(varId) = 1;

        DW_DecompNode* ptrNodeGreater = new DW_DecompNode(*ptrDecomNode);
        DW_DecompNode* ptrNodeSmaller = ptrDecomNode;
        ptrDecomNode = nullptr;

        cut.sense = '>';
        cut.rhs   = std::ceil(varValue);
        std::cout<<"x_"<<varId<<" >= "<<cut.rhs<<"\n";
        addMasterCut(cut, *ptrNodeGreater, numIt);

        cut.sense = '<';
        cut.rhs   = std::floor(varValue);
        std::cout<<"x_"<<varId<<" <= "<<cut.rhs<<"\n";
        addMasterCut(cut,*ptrNodeSmaller, numIt);

        auxVectors.updateSizes(*ptrNodeGreater);

        int statusSubProb = ptrNodeSmaller->columnGeneration(auxVectors);
        if(statusSubProb == StatusSubProb_Otimo)
            searchD->insert(ptrNodeSmaller);
        else
            delete ptrNodeSmaller;

        statusSubProb = ptrNodeGreater->columnGeneration(auxVectors);
        if(statusSubProb == StatusSubProb_Otimo)
            searchD->insert(ptrNodeGreater);
        else
            delete ptrNodeGreater;

        delete ptrDecomNode;

    }

    lowerBound            = std::min(searchD->getMin(), lowerBound);
    gap                   = computeGap(lowerBound, upperBound);
    end                   = clock();
    statisticD.totalTime  = double(end-start)/CLOCKS_PER_SEC;
    statisticD.upperBound = upperBound;
    statisticD.lowerBound = lowerBound;
    statisticD.gap        = gap;
    statisticD.timeLimit  = alarm_stopG;
    statisticD.numIt      = numIt;


    std::cout<<"it("<<it<<") \t LB("<<lowerBound<<") \t UB("<<upperBound<<") \t gap("<<gap<<"%)\n";

    return vetX_best;

}

bool BranchAndPriceNS::isInteger(const Eigen::VectorXd &vet)
{
    for(const double &val:vet)
    {
        if(!isInteger(val))
            return false;
    }

    return true;
}

double BranchAndPriceNS::computeGap(double lb, double ub)
{
    return ((ub-lb)/lb)*100.0;
}

void BranchAndPriceNS::writeToFile(StatisticsNS::StatisticsData& statisticsD, const std::string& fileStr)
{

    if(fileStr.empty())
        return;

    statisticsD.inst = std::filesystem::path(statisticsD.inst);

    std::filesystem::path path(statisticsD.inst);
    path.replace_extension();
    statisticsD.inst = path.string();

    std::ofstream file;
    if(std::filesystem::exists(fileStr))
    {
        file.open(fileStr, std::ios_base::app);

        if(!file.is_open())
        {
            std::cout<<"It is not possible to open: "<<fileStr<<"\n";
            throw "ERROR FILE";
        }
    }
    else
    {
        file.open(fileStr, std::ios_base::out);

        if(!file.is_open())
        {
            std::cout<<"It is not possible to open: "<<fileStr<<"\n";
            throw "ERROR FILE";
        }

        file<<"#"<<statisticsD.date<<"\n";
        file<<"inst; rootLB; rootTime; LB; UB; gap; totalTime\n";
    }


    file<<statisticsD.inst<<"; "<<std::format("{:.2f}; {:.2f}; {:.2f}; {:.2f}; {:.2f}; {:.2f}\n", statisticsD.rootLB,
                                                                                                       statisticsD.rootTime,
                                                                                                       statisticsD.lowerBound,
                                                                                                       statisticsD.upperBound,
                                                                                                       statisticsD.gap,
                                                                                                       statisticsD.totalTime);

    file.close();


}
