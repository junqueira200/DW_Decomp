/*  *****************************************************************
 *  *****************************************************************
 *  File:    Branch.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    10/01/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "Branch.h"
#include "BranchAndPrice.h"

using namespace DW_DecompNS;
using namespace BranchAndPriceNS;

int BranchNS::StrongBranch::operator()(const DW_DecompNS::DW_DecompNode *const node, DW_DecompNS::AuxData &auxVet)
{
    std::cout<<"StrongBranch\n";

    Eigen::VectorXd vetSol = node->vetSolX;
    BranchVar best;

    const int numCandC = std::min((int)vetSol.size(), NumCandidatesBranch);

    for(int i=0; i < numCandC; ++i)
    {
        BranchVar bVar;

        int id = getMostFractionVariable(vetSol);
        if(id == -1)
        {
            std::cout<<"Solution should not be integer!\n";
            throw "ERROR";
        }

        // Branch and compute the lowerBound
        Cut cut;
        cut.vetX.resize(vetSol.size());
        cut.vetX.setZero();
        cut.vetX.coeffRef(id) = 1;

        cut.sense = '>';
        cut.rhs   = std::ceil(vetSol[id]);
        DW_DecompNode* nodeAux = new DW_DecompNode(*node);
        addMasterCut(cut, *nodeAux, -1);
        auxVet.updateSizes(*nodeAux);

        int status = nodeAux->columnGeneration(auxVet);
        if(status == StatusSubProb_Otimo)
        {
            bVar.minLB = nodeAux->funcObj;
            bVar.varId = id;
        }

        delete nodeAux;
        nodeAux = new DW_DecompNode(*node);

        cut.sense = '<';
        cut.rhs   = std::floor(vetSol[id]);

        addMasterCut(cut, *nodeAux, -1);
        status = nodeAux->columnGeneration(auxVet);
        if(status == StatusSubProb_Otimo)
        {
            bVar.minLB = std::min(bVar.minLB, nodeAux->funcObj);
            bVar.varId = id;
        }

        if(bVar.minLB > best.minLB)
            best = bVar;

        vetSol[id] = 0.0;

        delete nodeAux;
    }

    std::cout<<"Var("<<best.varId<<") LB("<<best.minLB<<")\n\n";

    return best.varId;

}

int BranchNS::SimpleStrongBranch::operator()(const DW_DecompNS::DW_DecompNode *const node, AuxData &auxVet)
{

    return getMostFractionVariable(node->vetSolX);

}
