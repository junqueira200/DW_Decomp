//
// Created by igor on 04/01/25.
//

#include "BranchAndPrice.h"

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