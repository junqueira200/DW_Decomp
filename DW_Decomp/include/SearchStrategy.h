//
// Created by igor on 08/01/25.
//
#include <list>
#include "DW_Decomp.h"

#ifndef DW_SEARCHSTRATEGY_H
#define DW_SEARCHSTRATEGY_H

namespace SearchStrategyNS
{
    class SearchDataInter
    {
    public:

        virtual bool empty()=0;
        virtual DW_DecompNS::DW_DecompNode* pop()= 0;
        virtual void insert(DW_DecompNS::DW_DecompNode*)=0;
        virtual int size()=0;

        virtual ~SearchDataInter()=default;

    };
};

#endif //DW_SEARCHSTRATEGY_H
