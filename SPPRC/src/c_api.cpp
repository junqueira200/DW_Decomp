#include "c_api.h"

using namespace InstanceVRPTW_NS;
using namespace LabelingAlgorithmNS;
using namespace VrpTW_DecompLabelingNS;

extern "C"
void setInstance(int n)
{
    PRINT_DEBUG("", "");
    std::printf("n: %d\n", n);

    instanceVrptw = InstanceVRPTW(n);
}

extern "C"
void setDistance(int i, int j, int dist)
{
    //std::printf("\t%d %d: %d\n", i, j, dist);
    instanceVrptw.matDist(i, j) = dist;
}

extern "C"
void setDemand(int i, int dem)
{
    //std::printf("%d: %d\n", i, dem);
    instanceVrptw.vetClieDem[i] = dem;
}

extern "C"
int solveSpprc(double*	reducedCosts,
               int* 	vetRoutesOutput,
               int* 	sizeVetRoutesOutput,
               double* 	reducedCostRoutes,
               int		exact,
               int*     isComplete)
{
    static ArrayReducedCost arrayReducedCost;
    //PRINT_DEBUG("", "");
    int n = instanceVrptw.numClientes;

    //for(int i=0; i < n*n; ++i)
    //    std::printf("%.1f ", reducedCosts[i]);

    //std::printf("\n\n");

    bool complete = false;
    int numSol =
        vrpLabelingSubProb->resolveSubProbColuna(reducedCosts,
                                                 vetRoutesOutput,
                                                 *sizeVetRoutesOutput,
                                                 arrayReducedCost,
                                                 (bool)exact,
                                                 complete);
    *isComplete = complete ? 1 : 0;

    for(int i=0; i < numSol; ++i)
        reducedCostRoutes[i] = arrayReducedCost[i];

    return numSol;
}

extern "C"
double fixArcsSpprc(double* reducedCosts, double gap, double* arcMinOut)
{
    return vrpLabelingSubProb->fixArcsColuna(reducedCosts, gap, arcMinOut);
}

extern "C"
void setVrpLabelingSubProb()
{

    ArrayResources vetMaxResources;
    vetMaxResources[1] = instanceVrptw.capVeic;

    vrpLabelingSubProb = new
        VrpLabelingSubProb(instanceVrptw, instanceVrptw.sumDist(), vetMaxResources);

    //std::printf("END OF setVrpLabelingSubProb\n");
}

extern "C"
void setCapacityVehicle(int cap)
{
    instanceVrptw.capVeic = cap;
}

extern "C"
void setNumberOfVehicles(int k)
{
    instanceVrptw.numVeic = k;
}

extern "C"
int getSizeVetRoutesOutput()
{
    return NumMaxSolSubProb*NumMaxRoute;
}

extern "C"
int getMaxNumberOfRoutes()
{
    return NumMaxSolSubProb;
}
