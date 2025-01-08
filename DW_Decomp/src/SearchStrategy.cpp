//
// Created by igor on 08/01/25.

#include "SearchStrategy.h"

using namespace DW_DecompNS;


DW_DecompNS::DW_DecompNode *SearchStrategyNS::DepthFirst::pop()
{
    if(nodeList.empty())
        return nullptr;

    DW_DecompNode* node = nodeList.back();
    nodeList.pop_back();
    return node;
}


double SearchStrategyNS::DepthFirst::getMin()
{
    double min = std::numeric_limits<double>::infinity();

    for(DW_DecompNode* node:nodeList)
    {
        if(node->funcObj < min)
            min = node->funcObj;
    }

    return min;
}

double SearchStrategyNS::DepthFirst::getMax()
{

    double max = -std::numeric_limits<double>::infinity();

    for(DW_DecompNode* node:nodeList)
    {
        if(node->funcObj > max)
            max = node->funcObj;
    }

    return max;
}
