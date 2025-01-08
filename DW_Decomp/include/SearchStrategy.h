//
// Created by igor on 08/01/25.
//

#ifndef DW_SEARCHSTRATEGY_H
#define DW_SEARCHSTRATEGY_H

#include <list>
#include "DW_Decomp.h"

namespace SearchStrategyNS
{
    class SearchDataInter
    {
    public:

        virtual bool empty()=0;
        virtual DW_DecompNS::DW_DecompNode* pop()= 0;
        virtual void insert(DW_DecompNS::DW_DecompNode* node)=0;
        virtual int size()=0;
        virtual double getMin()=0;
        virtual double getMax()=0;

        virtual ~SearchDataInter()=default;

    };

    class DepthFirst : public SearchDataInter
    {
    public:

        std::list<DW_DecompNS::DW_DecompNode*> nodeList;

        bool empty() override{return nodeList.empty();};
        DW_DecompNS::DW_DecompNode* pop() override;
        void insert(DW_DecompNS::DW_DecompNode* node) override{nodeList.push_back(node);};
        int size() override{return int(nodeList.size());};
        double getMin() override;
        double getMax() override;
        ~DepthFirst() override=default;

    };

};

#endif //DW_SEARCHSTRATEGY_H
