/*  *****************************************************************
 *  *****************************************************************
 *  File:    Branch.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    10/01/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_BRANCH_H
#define DW_BRANCH_H

#include "DW_Decomp.h"

namespace BranchNS
{
    class BranchInter
    {
    public:

        ~BranchInter()=default;
        virtual int operator()(const DW_DecompNS::DW_DecompNode* const node, DW_DecompNS::AuxData &auxVet)=0;
    };

    class StrongBranch : public BranchInter
    {
    public:
        ~StrongBranch()=default;
        int operator()(const DW_DecompNS::DW_DecompNode* const node, DW_DecompNS::AuxData &auxVet) override;
    };

    class SimpleStrongBranch : public BranchInter
    {
    public:
        ~SimpleStrongBranch()=default;
        int operator()(const DW_DecompNS::DW_DecompNode* const node, DW_DecompNS::AuxData &auxVet) override;
    };

    struct BranchVar
    {
        int varId    = -1;
        double minLB = -std::numeric_limits<double>::infinity();
    };
}

#endif //DW_BRANCH_H
