//
// Created by igor on 04/01/25.
//

#ifndef DW_BRANCHANDPRICE_H
#define DW_BRANCHANDPRICE_H

#include "DW_Decomp.h"


namespace BranchAndPriceNS
{

    struct Cut
    {
        Eigen::SparseVector<double> vetX;
        char sense;
        double rhs;
    };

    int getMostFractionVariable(const Eigen::VectorXd &vetSolX);
    void addMasterCut(const Cut &cut, DW_DecompNS::DW_DecompNode &decompNode, int num);

}

#endif //DW_BRANCHANDPRICE_H
