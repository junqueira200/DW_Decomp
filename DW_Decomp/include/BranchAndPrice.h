/*  *****************************************************************
 *  *****************************************************************
 *  File:    BranchAndPrice.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    04/01/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_BRANCHANDPRICE_H
#define DW_BRANCHANDPRICE_H

#include "DW_Decomp.h"
#include "SearchStrategy.h"
#include "PrimalHeuristic.h"
#include "Branch.h"
#include "Statistics.h"

namespace BranchAndPriceNS
{
    constexpr double IntFeasTol = 1E-5;

    struct RobustCut
    {
        Eigen::SparseVector<double> vetX;
        char sense  = '<';
        double rhs  = 0.0;

        RobustCut()=default;
        RobustCut(int sizeVetX){vetX = Eigen::SparseVector<double>(sizeVetX);};
    };

    class RobustCutGenerator
    {
    public:
        ~RobustCutGenerator()=default;
        virtual int operator()(DW_DecompNS::DW_DecompNode& decompNode)=0;
    };

    int getMostFractionVariable(const Eigen::VectorXd &vetSolX);
    bool isInteger(const Eigen::VectorXd &vet);
    void addMasterCut(const RobustCut &cut, DW_DecompNS::DW_DecompNode &decompNode, int num,
                      bool isBranching);

    Eigen::VectorXd branchAndPrice(DW_DecompNS::DW_DecompNode &cRootNode,
                                   DW_DecompNS::AuxData &auxVectors,
                                   SearchStrategyNS::SearchDataInter* searchD,
                                   PrimalHeuristicNS::PrimalHeuristicInter* ptrPrimalH,
                                   BranchNS::BranchInter* branch,
                                   RobustCutGenerator* robustCutGenerator,
                                   StatisticsNS::StatisticsData& statisticD);

    double computeGap(double lb, double ub);
    void writeToFile(StatisticsNS::StatisticsData& statisticsD, const std::string& fileStr, std::string& extraHead,
                     std::string& extraCont);

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
