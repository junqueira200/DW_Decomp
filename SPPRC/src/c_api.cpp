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

extern "C"
int setRank1Cuts(int numCuts, const int* custs, const double* duals, const int* ids,
                 const int* memStart, const int* memCusts)
{
    LabelingData& ld = vrpLabelingSubProb->labelingDataForward;
    const int n = instanceVrptw.numClientes;

    int numNeg = 0;
    for(int c = 0; c < numCuts; ++c)
        numNeg += (duals[c] < -1e-9);

    if(numNeg > MaxActiveCuts)
    {
        ld.sizeVetRank1Cut = 0;
        return -1;
    }

    int k = 0;
    for(int c = 0; c < numCuts; ++c)
    {
        const double pen = -duals[c];           // sign fixed here: penalty >= 0
        if(pen <= 1e-9)
            continue;

        Rank1Cut& r = ld.vetRank1Cut[k];
        r.S.reset();
        r.memory.reset();                       // limited memory (vertices)
        for(int p = memStart[c]; p < memStart[c + 1]; ++p)
        {
            const int v = memCusts[p];
            if(v <= 0 || v >= n)
            {
                std::printf("setRank1Cuts: memory customer %d out of range (n = %d)\n", v, n);
                PRINT_THROW();
            }
            r.memory.set(v);
        }
        for(int i = 0; i < 3; ++i)
        {
            const int cust = custs[3*c + i];
            if(cust <= 0 || cust >= n)
            {
                std::printf("setRank1Cuts: customer %d out of range (n = %d)\n", cust, n);
                PRINT_THROW();
            }
            r.vetS[i] = cust;
            r.S.set(cust);
            r.memory.set(cust);                 // S ⊆ memory
        }
        r.sizeVetS = 3;
        r.penalty  = pen;
        r.colunaId = ids[c];
        ++k;
    }

    ld.sizeVetRank1Cut = k;                     // updateRank1CutsData() runs in the pricing
    return k;
}

extern "C"
int getMaxActiveCuts()
{
    return MaxActiveCuts;
}
