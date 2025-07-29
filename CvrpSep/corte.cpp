#include <iostream>

extern "C" {

#include "src/basegrph.h"
#include <stdio.h>
#include "src/cnstrmgr.h"
#include "src/capsep.h"
}

int main()
{
    int dim = 10;
    int maxNoOfCuts = 10;
    CnstrMgrPointer cutsCMP;
    CnstrMgrPointer oldCutsCMP;

    CMGR_CreateCMgr(&cutsCMP, dim);
    CMGR_CreateCMgr(&oldCutsCMP, dim);

    double EpsForIntegrality = 0.000001;
    double MaxViolation;
    char IntegerAndFeasible;

    int n_customers = 4;
    int noOfEdeges = 10;
    int demand[] = {0, 8, 18, 1, 13};
    int capacity = 20;

    int edge_head[] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 3};
    int edge_tail[] = {1, 2, 3, 4, 2, 3, 4, 3, 4, 4};

    double edge_x[] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5};

    CAPSEP_SeparateCapCuts(n_customers, demand, capacity, noOfEdeges, edge_tail, edge_head, edge_x, oldCutsCMP,
                           maxNoOfCuts, EpsForIntegrality, &IntegerAndFeasible, &MaxViolation, cutsCMP);

    std::cout<<"Size: "<<cutsCMP->Size<<"\n";
    std::cout<<"IntegerAndFeasible: "<<(int)IntegerAndFeasible<<"\n";

    int i, j, Listsize;
    double rhs;
    int List[] = {-1, -1, -1, -1, -1};


    for(i=0; i < cutsCMP->Size; ++i)
    {
        Listsize = 0;
        std::cout<<"IntListSize: "<<cutsCMP->CPL[i]->IntListSize<<"\nList: ";
        for(j=1; j <= cutsCMP->CPL[i]->IntListSize; ++j)
        {
            List[++Listsize] = cutsCMP->CPL[i]->IntList[j];
        }

        for(int j=1; j <= Listsize; ++j)
            std::cout<<List[j]<<" ";

        rhs = cutsCMP->CPL[i]->RHS;
        std::cout<<"\nRhs: "<<rhs<<"\n\n";

    }

    CMGR_FreeMemCMgr(&cutsCMP);
    CMGR_FreeMemCMgr(&oldCutsCMP);

    return 0;
}
