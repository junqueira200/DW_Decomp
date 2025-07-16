#include <stdio.h>
#include "src/basegrph.h"
#include "src/cnstrmgr.h"
#include "src/capsep.h"


int main()
{
    int dim = 5;
    int maxNoOfCuts = 10;
    CnstrMgrPointer cutsCMP;
    CnstrMgrPointer oldCutsCMP;

    CMGR_CreateCMgr(&cutsCMP, dim);
    CMGR_CreateCMgr(&oldCutsCMP, dim);

    double EpsForIntegrality = 0.000001;
    double MaxViolation;
    char IntegerAndFeasible;

    int n_customers = 4;
    int noOfEdeges = 6;
    int demand[] = {0, 8, 18, 1, 13};
    int capacity = 20;
    int edge_tail[] = {1, 1, 1, 2, 2, 3};

    int edge_head[] = {2, 3, 4, 3, 4, 4};

    double edge_x[] = {0.5, 0.5, 0.0, 0.5, 0.5, 0.5};

    CAPSEP_SeparateCapCuts(n_customers, demand, capacity, noOfEdeges, edge_tail, edge_head, edge_x, oldCutsCMP,
                           maxNoOfCuts, EpsForIntegrality, &IntegerAndFeasible, &MaxViolation, cutsCMP);

}
