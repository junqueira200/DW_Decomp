/*  *****************************************************************
 *  *****************************************************************
 *  File:    SearchStrategy.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    08/01/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_SEARCHSTRATEGY_H
#define DW_SEARCHSTRATEGY_H

#include <list>
#include <set>
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
        ~DepthFirst() override;

    };

    class BreadthFirst : public SearchDataInter
    {
    public:

        std::list<DW_DecompNS::DW_DecompNode*> nodeList;

        bool empty() override{return nodeList.empty();};
        DW_DecompNS::DW_DecompNode* pop() override;
        void insert(DW_DecompNS::DW_DecompNode* node) override{nodeList.push_back(node);};
        int size() override{return int(nodeList.size());};
        double getMin() override;
        double getMax() override;
        ~BreadthFirst() override;
    };

    class DecompNodeCompMin
    {
    public:

        bool operator()(DW_DecompNS::DW_DecompNode *node0, DW_DecompNS::DW_DecompNode *node1) const
        {return node0->funcObj < node1->funcObj;}

    };

    class MinFuncObj : public SearchDataInter
    {
    public:

        std::multiset<DW_DecompNS::DW_DecompNode*, DecompNodeCompMin> setDecompNode;


        bool empty() override{return setDecompNode.empty();};
        DW_DecompNS::DW_DecompNode* pop() override;
        void insert(DW_DecompNS::DW_DecompNode* node) override{setDecompNode.insert(node);};
        int size() override{return int(setDecompNode.size());};
        double getMin() override;
        double getMax() override;
        ~MinFuncObj() override;
    };

};

#endif //DW_SEARCHSTRATEGY_H
