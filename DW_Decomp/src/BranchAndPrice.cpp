//
// Created by igor on 04/01/25.
//

#include "BranchAndPrice.h"
using namespace DW_DecompNS;

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
    GRBVar *varRmlp = decompNode.uRmlp->getVars();

    // Runs through columns and compute their coefficients
    for(int i=0; i < int(decompNode.vetVarLambdaCol.size()); ++i)
    {
        auto &col    = *decompNode.vetVarLambdaCol[i];
        double coef = (cut.vetX*col)[0];
        if(coef != 0.0)
            linExpr += coef*varRmlp[i];
    }

    // Add cut to rmlp
    decompNode.uRmlp->addConstr(linExpr, cut.sense, cut.rhs, "masterCut_"+std::to_string(num));

}

void BranchAndPriceNS::branchAndPrice(const DW_DecompNS::DW_DecompNode &rootNode, DW_DecompNS::AuxData &auxVectors)
{

    std::list<DW_DecompNode*> listDecomNode;
    listDecomNode.emplace_back(new DW_DecompNode(rootNode));

    if(rootNode.uRmlp->get(GRB_IntAttr_ModelSense) != GRB_MINIMIZE)
    {
        std::cout<<"ERROR, model should be a min problem\n";
        PRINT_DEBUG("", "");
        throw "ERROR IN MODEL";
    }

    double lowerBound = -std::numeric_limits<double>::infinity();
    double upperBound =  std::numeric_limits<double>::infinity();
    auto vetX_best    = rootNode.vetSolX;

    int numIt = -1;

    while(!listDecomNode.empty())
    {
        numIt += 1;

        DW_DecompNode* ptrDecomNode = listDecomNode.back();
        listDecomNode.pop_back();


        int statusSubProb = ptrDecomNode->columnGeneration(auxVectors);
        if(statusSubProb != StatusSubProb_Otimo || ptrDecomNode->funcObj > upperBound)
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

        if(numIt == 1)
        {
            std::cout<<"NAO EH ERROR\n";
            PRINT_DEBUG("", "");
            throw "NAO_EH_ERROR";
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
        addMasterCut(cut, *ptrNodeGreater, numIt);

        cut.sense = '<';
        cut.rhs   = std::floor(varValue);
        addMasterCut(cut,*ptrNodeSmaller, numIt);

        listDecomNode.push_back(ptrNodeSmaller);
        listDecomNode.push_back(ptrNodeGreater);

    }

}

bool BranchAndPriceNS::isInteger(const Eigen::VectorXd &vet)
{
    for(const double &val:vet)
    {
        double ceil     = std::ceil(val);
        double floor    = std::floor(val);
        double difCeil  = std::abs(ceil-val);
        double difFloor = std::abs(val-floor);
        double smaller  = difCeil<difFloor? difCeil:difFloor;

        if(smaller > IntFeasTol)
            return false;
    }

    return true;
}

