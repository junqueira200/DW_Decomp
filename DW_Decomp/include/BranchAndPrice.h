//
// Created by igor on 04/01/25.
//

#ifndef DW_BRANCHANDPRICE_H
#define DW_BRANCHANDPRICE_H

#include "DW_Decomp.h"
#include "SearchStrategy.h"


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
    void branchAndPrice(DW_DecompNS::DW_DecompNode &cRootNode,
                        DW_DecompNS::AuxData &auxVectors,
                        SearchStrategyNS::SearchDataInter* searchD);
    double computeGap(double lb, double ub);

    inline __attribute__((always_inline))
    bool isInteger(double val)
    {

        double ceil     = std::ceil(val);
        double floor    = std::floor(val);
        double difCeil  = std::abs(ceil-val);
        double difFloor = std::abs(val-floor);
        double smaller  = difCeil<difFloor? difCeil:difFloor;

        return (smaller <= IntFeasTol);
    }

}

#endif //DW_BRANCHANDPRICE_H
