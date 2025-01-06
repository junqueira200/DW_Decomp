//
// Created by igor on 04/01/25.
//

#ifndef DW_BRANCHANDPRICE_H
#define DW_BRANCHANDPRICE_H

#include "DW_Decomp.h"


namespace BranchAndPriceNS
{
    constexpr double IntFeasTol = 1E-5;

    struct Cut
    {
        Eigen::SparseVector<double> vetX;
        char sense;
        double rhs;
    };

    int getMostFractionVariable(const Eigen::VectorXd &vetSolX);
    bool isInteger(const Eigen::VectorXd &vet);
    void addMasterCut(const Cut &cut, DW_DecompNS::DW_DecompNode &decompNode, int num);
    void branchAndPrice(const DW_DecompNS::DW_DecompNode &cRootNode, DW_DecompNS::AuxData &auxVectors);
    double computeGap(double lb, double ub);
    double computeLowerBaound(std::list<DW_DecompNS::DW_DecompNode*> &list);

}

#endif //DW_BRANCHANDPRICE_H
