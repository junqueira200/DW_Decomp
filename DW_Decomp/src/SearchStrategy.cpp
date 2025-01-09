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

SearchStrategyNS::DepthFirst::~DepthFirst()
{
    for(DW_DecompNode* node:nodeList)
        delete node;
}

DW_DecompNS::DW_DecompNode *SearchStrategyNS::BreadthFirst::pop()
{
    if(nodeList.empty())
        return nullptr;

    DW_DecompNode* node = nodeList.front();
    nodeList.pop_front();

    return node;
}

double SearchStrategyNS::BreadthFirst::getMin()
{
    double min = std::numeric_limits<double>::infinity();

    for(DW_DecompNode* node:nodeList)
    {
        if(node->funcObj < min)
            min = node->funcObj;
    }

    return min;
}

double SearchStrategyNS::BreadthFirst::getMax()
{

    double max = -std::numeric_limits<double>::infinity();

    for(DW_DecompNode* node:nodeList)
    {
        if(node->funcObj > max)
            max = node->funcObj;
    }

    return max;
}

SearchStrategyNS::BreadthFirst::~BreadthFirst()
{

    for(DW_DecompNode* node:nodeList)
        delete node;
}

DW_DecompNS::DW_DecompNode *SearchStrategyNS::MinFuncObj::pop()
{
   if(setDecompNode.empty())
       return nullptr;

    DW_DecompNode* node = (*setDecompNode.begin());
    setDecompNode.erase(node);
    return node;
}

double SearchStrategyNS::MinFuncObj::getMin()
{
    double min = std::numeric_limits<double>::infinity();

    if(!setDecompNode.empty())
        min = (*setDecompNode.begin())->funcObj;

    return min;
}

double SearchStrategyNS::MinFuncObj::getMax()
{

    double max = -std::numeric_limits<double>::infinity();

    if(!setDecompNode.empty())
    {
        auto itEnd = setDecompNode.end();
        --itEnd;
        max = (*itEnd)->funcObj;
    }

    return max;
}

SearchStrategyNS::MinFuncObj::~MinFuncObj()
{

    for(DW_DecompNode* node:setDecompNode)
        delete node;
}