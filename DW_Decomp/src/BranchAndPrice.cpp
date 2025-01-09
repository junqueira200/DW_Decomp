//
// Created by igor on 04/01/25.
//

#include "BranchAndPrice.h"
#include "PrimalHeuristic.h"

using namespace DW_DecompNS;
using namespace SearchStrategyNS;
using namespace PrimalHeuristicNS;

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

    // Update info
    decompNode.info.numConstrsMaster += 1;

    // Add cut to matA
    auto &matA = decompNode.matA;
    matA.conservativeResize(matA.rows()+1, matA.cols());
    matA.row(matA.rows()-1) = cut.vetX;

    // Update columns coef
    GRBLinExpr linExpr;

    //decompNode.uRmlp->update();
    //decompNode.uRmlp->optimize();
    GRBVar *varRmlp = decompNode.uRmlp->getVars();


    int numVars = decompNode.uRmlp->get(GRB_IntAttr_NumVars);

    if(numVars != int(decompNode.vetVarLambdaCol.size()))
    {
        std::cout<<"Num Vars is wrong;\n\t Model: "<<numVars<<"; vetVarLamdaCol: "<<decompNode.vetVarLambdaCol.size()<<"\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }


    // Runs through columns and compute their coefficients
    for(int i=0; i < int(decompNode.vetVarLambdaCol.size()); ++i)
    {
        auto &col    = *decompNode.vetVarLambdaCol[i];
        double coef = (cut.vetX.dot(col));

        if(coef != 0.0)
            linExpr += coef * varRmlp[i];

        //exit(-1);
    }

    std::cout<<"add cut\n";
    // Add cut to rmlp
    decompNode.uRmlp->addConstr(linExpr, cut.sense, cut.rhs, "masterCut_"+std::to_string(num));
    delete []varRmlp;
}

void BranchAndPriceNS::branchAndPrice(DW_DecompNS::DW_DecompNode &cRootNode,
                                      DW_DecompNS::AuxData &auxVectors,
                                      SearchDataInter* searchD)
{

    PrimalHeuristicInter* ptrPrimalH = new SimpleDiving;

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
    double gap = std::numeric_limits<double>::infinity();

    int numIt = -1;

    DW_DecompNode *rootNode = new DW_DecompNode(cRootNode);
    int status = rootNode->columnGeneration(auxVectors);
    if(status != StatusSubProb_Otimo)
    {
        delete ptrPrimalH;
        std::cout<<"Root node can't be soved\n";
        delete rootNode;
        return;
    }

    searchD->insert(rootNode);
    rootNode = nullptr;

    int it = -1;

    while(!searchD->empty() && gap > gapLimit)
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

        std::cout<<"it("<<it<<") \t LB("<<lowerBound<<") \t UB("<<upperBound<<") \t gap("<<gap<<"%)\n";

        std::cout<<"Processando NO: "<<ptrDecomNode<<"\nNumVars: "<<ptrDecomNode->uRmlp->get(GRB_IntAttr_NumVars)<<"\n\n";
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

        int varId = getMostFractionVariable(ptrDecomNode->vetSolX);
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

    delete ptrPrimalH;

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

