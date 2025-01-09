//
// Created by igor on 09/01/25.
//
#include "PrimalHeuristic.h"
#include "SearchStrategy.h"
#include "BranchAndPrice.h"

using namespace DW_DecompNS;
using namespace SearchStrategyNS;
using namespace BranchAndPriceNS;

DW_DecompNode* PrimalHeuristicNS::SimpleDiving::operator()(DW_DecompNode *node,
                                                           AuxData &auxVet,
                                                           double upperBound)
{
    DepthFirst dfs;
    DW_DecompNode* root = new DW_DecompNode(*node);

    int status = root->columnGeneration(auxVet);
    if(status != StatusSubProb_Otimo)
    {
        delete root;
        return nullptr;
    }

    dfs.insert(root);
    root = nullptr;
    int numIt = -1;

    while(true)
    {
        numIt += 1;

        DW_DecompNode* nodeTemp = dfs.pop();

        if(nodeTemp->funcObj > upperBound)
        {
            delete nodeTemp;
            return nullptr;
        }

        int id = -1;
        double val = std::numeric_limits<double>::infinity();

        // Chose a variable to branch
        for(int i=0; i < int(nodeTemp->vetSolX.size()); ++i)
        {
            double dif = (std::ceil(nodeTemp->vetSolX[i])-nodeTemp->vetSolX[i]);

            if(dif < val && !isInteger(nodeTemp->vetSolX[i]))
            {
                val = dif;
                id = i;
            }
        }

        // Check is solution is integer
        if(id == -1)
            return nodeTemp;


        Cut cut;
        cut.vetX.resize(nodeTemp->vetSolX.size());
        cut.vetX.coeffRef(id) = 1;
        cut.rhs = std::ceil(nodeTemp->vetSolX[id]);
        cut.sense = '>';

        addMasterCut(cut, *nodeTemp, numIt);
        auxVet.updateSizes(*nodeTemp);

        status = nodeTemp->columnGeneration(auxVet);
        if(status != StatusSubProb_Otimo)
        {
            delete nodeTemp;
            return nullptr;
        }

        dfs.insert(nodeTemp);
    }

}
