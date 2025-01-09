//
// Created by igor on 09/01/25.
//

#ifndef DW_PRIMALHEURISTIC_H
#define DW_PRIMALHEURISTIC_H

#include "DW_Decomp.h"

namespace PrimalHeuristicNS
{

    class PrimalHeuristicInter
    {
    public:

        virtual ~PrimalHeuristicInter()=default;
        virtual DW_DecompNS::DW_DecompNode* operator()(DW_DecompNS::DW_DecompNode* node,
                                                        DW_DecompNS::AuxData &auxVet,
                                                        double upperBound)=0;

    };

    class SimpleDiving : public PrimalHeuristicInter
    {
    public:

        DW_DecompNS::DW_DecompNode* operator()(DW_DecompNS::DW_DecompNode* node,
                                               DW_DecompNS::AuxData &auxVet,
                                               double upperBound) override;

        SimpleDiving()=default;
        ~SimpleDiving() override=default;

    };

}

#endif //DW_PRIMALHEURISTIC_H
