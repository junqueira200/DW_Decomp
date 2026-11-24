#ifndef C_API_H
#define C_API_H

#include "VrpTW_DecompLabeling.h"

InstanceVRPTW_NS::InstanceVRPTW instanceVrptw;
inline VrpTW_DecompLabelingNS::VrpLabelingSubProb* vrpLabelingSubProb = nullptr;

extern "C"
{

    void setInstance(int n);
    void setDistance(int i, int j, int dist);
    void setDemand(int i, int dem);
    void setCapacityVehicle(int cap);
    void setNumberOfVehicle(int k);
    int solveSpprc(double* 	reducedCosts,
                   int* 	vetRoutesOutput,
                   int* 	sizeVetRoutesOutput,
                   double* 	reducedCostRoutes,
                   int		exact,
                   int*     isComplete);      // out: 1 = exact search finished

    // Reduced-cost arc fixing. reducedCosts: n*n, same layout as solveSpprc.
    // arcMinOut[i*n+j] = lower bound on the reduced cost of any route using (i,j)
    // (Inf if no route with value < gap uses it). Returns the min route reduced cost.
    double fixArcsSpprc(double* reducedCosts,
                        double  gap,
                        double* arcMinOut);

    void setVrpLabelingSubProb();
    int getSizeVetRoutesOutput();
    int getMaxNumberOfRoutes();

    // Rank-1 cuts (subset-row, |S| = 3, multiplier 1/2) for the next pricing calls.
    // custs: 3*numCuts customers (0 = depot, 1..n-1), duals: sigma <= 0 of each cut,
    // ids: any id (kept in Rank1Cut::colunaId). Cuts with sigma >= -1e-9 are skipped.
    // Returns the number of cuts loaded, or -1 if more than MaxActiveCuts have sigma < 0
    // (then none is loaded: the exact pricing would not be exact).
    int setRank1Cuts(int numCuts, const int* custs, const double* duals, const int* ids);
    int getMaxActiveCuts();
}

#endif // C_API_H
