/*  *****************************************************************
 *  *****************************************************************
 *  File:    VrpTW_DecompLabeling.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    20/11/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "VrpTW_DecompLabeling.h"
#include "sefe_array.h"
#include "LowerBound.h"
#include <set>
#include <algorithm>
#include "LowerBound.h"

using namespace LabelingAlgorithmNS;
using namespace TestNS;
using namespace LowerBoundNS;
//using namespace BranchAndPriceNS;

VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb(
    InstanceVRPTW_NS::InstanceVRPTW     &instVrpTw_,
    double                               startDist,
    LabelingAlgorithmNS::ArrayResources &vetMaxResouces)
{
    PRINT_DEBUG("", "");

    // Nodes are 0..numClientes (numClientes = destination, a copy of the depot),
    // so the ng bitset needs numClientes + 1 positions.
    if(instVrpTw_.numClientes + 1 > NumMaxCust)
    {
        PRINT_DEBUG("",
                    "ERRO; Const NumMaxCust(" << NumMaxCust << ") < NumCust + 1("
                                              << instVrpTw_.numClientes + 1 << ")\n\t");
        throw "ERRO";
    }

    instVrpTw = &instVrpTw_;
    double numSteps = 10.0; // 20 // 10

    // vetStepSize[0].stepSize = 400;
    vetStepSize[0].stepSize = 500.0; // 200
    vetStepSize[0].start = -1000;       // -500
    vetStepSize[0].end = 1000;          // 500

    std::printf("instVrpTw->capVeic: %d\n", instVrpTw->capVeic);

    vetStepSize[1].stepSize = instVrpTw->capVeic / numSteps; // 50
    vetStepSize[1].start = 0;
    vetStepSize[1].end = (FloatType)instVrpTw->capVeic+1;

    labelingDataForward = LabelingAlgorithmNS::LabelingData(
        vetStepSize, 2, instVrpTw->numClientes + 1, vetMaxResouces);
    labelingDataForward.matDist = instVrpTw->matDist;

    labelingDataBackward = LabelingAlgorithmNS::LabelingData(
        vetStepSize, 2, instVrpTw->numClientes + 1, vetMaxResouces);
    labelingDataBackward.matDist = instVrpTw->matDist;

    vetMatResCostForward = Vet3D_ResCost(instVrpTw->numClientes + 1, instVrpTw->numClientes + 1, 2);
    vetMatResCostForward.setVal(0.0);

    vetMatResCostBackward =
        Vet3D_ResCost(instVrpTw->numClientes + 1, instVrpTw->numClientes + 1, 2);
    vetMatResCostBackward.setVal(0.0);

    double val = vetStepSize[0].start;
    std::printf("Resorce 0: %d\n", labelingDataForward.vetNumSteps(0));
    std::printf("\nResorce 1: %d\n", labelingDataForward.vetNumSteps(0));
    std::printf("start: %f\n", vetStepSize[0].start);

    for(int i=0; i < labelingDataForward.vetNumSteps(0); ++i)
    {
        std::cout<<"["<<val<<";"<<val+vetStepSize[0].stepSize<<") ";
        val += vetStepSize[0].stepSize;
    }

    std::cout<<"\n\n";

    std::printf("\nResorce 1: %d\n", labelingDataForward.vetNumSteps(1));
    std::printf("start: %f\n", vetStepSize[1].start);

    val = vetStepSize[1].start;
    for(int i=0; i < labelingDataForward.vetNumSteps(1); ++i)
    {
        std::cout<<"["<<val<<";"<<val+vetStepSize[1].stepSize<<") ";
        val += vetStepSize[1].stepSize;
    }

    // std::cout<<"\n";

    const int dest = instVrpTw->numClientes;   // copy of the depot

    for(int i = 0; i < instVrpTw->numClientes + 1; ++i)
    {
        for(int j = 0; j < (instVrpTw->numClientes + 1); ++j)
        {
            if(i == j)
                continue;

            // vetClieDem has only numClientes entries (0..numClientes-1).
            // The destination node (dest == numClientes) is a copy of the depot,
            // so it has no demand of its own.
            const FloatType demJ = (j == dest) ? (FloatType)0.0
                                               : (FloatType)instVrpTw->vetClieDem[j];

            vetMatResCostForward.get(i, j, 1)  = demJ;
            vetMatResCostBackward.get(j, i, 1) = demJ;
        }
    }

    // Same size as vetMatResCostForward; filled on every heuristic call.
    vetMatResCostForwardHeur =
        Vet3D_ResCost(instVrpTw->numClientes + 1, instVrpTw->numClientes + 1, 2);

    // Completion bound of the exact pricing and of the arc fixing.
    matDpColuna = EigenMatrixRowD(instVrpTw->numClientes + 1, instVrpTw->capVeic + 1);
    matDpHeur   = EigenMatrixRowD(instVrpTw->numClientes + 1, instVrpTw->capVeic + 1);

    // Bidirectional labeling (backward run). Resource 1 is copied from the forward
    // graph (demand of the head) and never changes; resource 0 is set on every exact
    // call by updateVetMatResCostBackward().
    vetMatResCostTransposed =
        Vet3D_ResCost(instVrpTw->numClientes + 1, instVrpTw->numClientes + 1, 2);
    vetMatResCostTransposed.copy(vetMatResCostForward);
    matDpBackward = EigenMatrixRowD(instVrpTw->numClientes + 1, instVrpTw->capVeic + 1);

    // std::cout<<"\n"<<instVrpTw->vetClieDem<<"\n\n";
    // PRINT_EXIT();

    // vetVetResBound = LabelingAlgorithmNS::VetVetResBound(2);

    vetVetResBound = LabelingAlgorithmNS::MatBoundRes(instVrpTw->numClientes + 1, 2);

    // vetVetResBound[0].resize(instVrpTw->numClientes+1);
    // vetVetResBound[1].resize(instVrpTw->numClientes+1);

    Bound bound0;
    bound0.lowerBound = -std::numeric_limits<FloatType>::infinity();
    bound0.upperBound = std::numeric_limits<FloatType>::infinity();

    Bound bound1;
    bound1.lowerBound = 0;
    bound1.upperBound = (FloatType)instVrpTw->capVeic;

    for(int i = 0; i < instVrpTw->numClientes + 1; ++i)
    {
        vetVetResBound(i, 0) = bound0;
        /*
        if(i > 0 && i < instVrpTw->numClientes)
            bound1.lowerBound = (FloatType)instVrpTw->vetClieDem[i];
        else
            bound1.lowerBound = 0.0;
        */
        vetVetResBound(i, 1) = bound1;
    }

    ngSet = NgSet(instVrpTw->numClientes + 1, NgSetSize);
    ngSet.setNgSets(instVrpTw->matDist);
    ngSet.active = true;

    // enumerateRoutes(*instVrpTw, instVrpTw->numClientes-1, routeHash);
}
/*
void VrpTW_DecompLabelingNS::VrpLabelingSubProb::iniConvConstr(GRBModel    &rmlp,
                                                               void        *data,
                                                               const double custoVarA)
{

    GRBLinExpr linExpr;
    GRBVar     a = rmlp.addVar(0, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
    // linExpr += -a;

    rmlp.addConstr(linExpr, '<', instVrpTw->numVeic, "convConstr");
}
*/

VrpTW_DecompLabelingNS::VrpLabelingSubProb::VrpLabelingSubProb()
{

    // enumerateRoutes(*instVrpTw, instVrpTw->numClientes-1, routeHash);
}

bool VrpTW_DecompLabelingNS::VrpLabelingSubProb::checkEnumeratedRoutesFinal(
    Eigen::RowVectorXd &vetRowPi)
{
    return false;

    cleanVetRouteG();
    std::cout << "Start checkEnumeratedRoutes\n\n";
    int  t = 0;
    bool ret = false;
    for(const Route &route : routeHash)
    {
        FloatType redCost = computeReducedCost(route, vetRowPi);
        if(t < 10)
            std::cout << redCost << "\n";
        t += 1;
        if(redCost < -TolObjSubProb)
        {
            std::cout << "Route(" << route.vetRoute << ") have a NEGATIVE(" << redCost
                      << ") reduced cost!\n";
            if(redCost <= -1.0)
            {
                addToVetRoutesG(route.vetRoute);
                ret = true;
            }
        }
    }

    return ret;
}

// Reduced costs of the transposed graph (resource 0 of vetMatResCostTransposed), from
// vetMatResCostForward. Nodes: 0 = source (depot), 1..n-1 = customers, n = destination
// (depot copy). A path 0 -> v -> ... -> n here is a route suffix read backwards.
void VrpTW_DecompLabelingNS::VrpLabelingSubProb::updateVetMatResCostBackward()
{
    const int n = instVrpTw->numClientes;          // dest = n

    // Customer arcs: (i, j) here is (j, i) in the forward graph.
    for(int i = 1; i < n; ++i)
    {
        for(int j = 1; j < n; ++j)
        {
            if(i == j)
                continue;

            vetMatResCostTransposed.get(i, j, 0) = vetMatResCostForward(j, i, 0);
        }
    }

    for(int i = 1; i < n; ++i)
    {
        // First arc of the backward path: forward arc (i, depot) = (i, dest).
        vetMatResCostTransposed.get(0, i, 0) = vetMatResCostForward(i, n, 0);
        // Last arc of the backward path: forward arc (depot, i) = (0, i).
        vetMatResCostTransposed.get(i, n, 0) = vetMatResCostForward(0, i, 0);
        // Arcs into the source are never used.
        vetMatResCostTransposed.get(i, 0, 0) = InfFloatType;
    }

    // Empty route (source -> destination): same arc in both directions.
    vetMatResCostTransposed.get(0, n, 0) = vetMatResCostForward(0, n, 0);
}

// Merge of the bidirectional labeling.
//
// Forward label F at u: path 0 -> ... -> u, load(F) includes d_u. Backward label B at v:
// path v -> ... -> depot (built on the transposed graph), load(B) includes d_v. Through
// the arc (u, v) they give the route 0 ... u v ... 0 with
//     rc = rc(F) + c(u, v) + rc(B),   load = load(F) + load(B) <= Q.
// ng: the memories must not share a customer (this also excludes v in F and u in B).
//
// Completeness (every elementary route is found or dominated): for a route, let u be the
// last customer whose prefix load is <= splitLoad (u = 0, the start label, if the first
// customer alone is above it) and v the next one. The prefix up to u was extended (load
// <= splitLoad), so F exists (or a label that dominates it). The suffix from v exists in
// the backward run: its parent has load = total - (prefix load up to v) < Q - splitLoad
// = splitLoad (split at Q/2), so the parent was extended. Routes whose load stays <=
// splitLoad up to the last customer are found by the forward run itself.
// Dominance keeps this true: a dominating label has lower or equal rc and load and an
// ng memory contained in the dominated one, so the same join stays feasible.
//
// So only the forward labels with load <= splitLoad are joined (they are the closed
// ones; the others were kept but not extended), with every backward label.
int VrpTW_DecompLabelingNS::VrpLabelingSubProb::mergeBidirectional(
    FloatType                        splitLoad,
    int                              numSolForward,
    const Eigen::VectorX<FloatType> &vetRedCostForward,
    int                             *vetRoutesOut,
    int                             &sizeVetRoutesOut,
    ArrayReducedCost                &vetRedCost)
{
    const int    n   = instVrpTw->numClientes;     // nodes 0..n-1, dest = n
    const int    Q   = instVrpTw->capVeic;
    const int    K   = NumMaxSolSubProb;
    const double inf = std::numeric_limits<double>::infinity();

    LabelingData &lF = labelingDataForward;
    LabelingData &lB = labelingDataBackward;
    const int numI = lF.vetNumSteps(0);
    const int numJ = lF.vetNumSteps(1);

    // 1. Backward labels per customer v, sorted by reduced cost, and
    //    bestB(v, r) = best reduced cost of a backward label at v with load <= r.
    struct BackLabel
    {
        double       rc;
        int          load;
        const Label *label;
    };
    static std::vector<std::vector<BackLabel>> vetBack;
    static EigenMatrixRowD                     bestB;
    vetBack.resize(n);
    bestB.resize(n, Q + 1);
    bestB.setConstant(inf);

    for(int v = 1; v < n; ++v)
    {
        std::vector<BackLabel> &list = vetBack[v];
        list.clear();

        for(int i = 0; i < numI; ++i)
        {
            for(int j = 0; j < numJ; ++j)
            {
                // Open buckets (labels above the split, never extended) and closed ones.
                for(const Bucket *bucket : {&lB.vetMatBucketForward[v].mat(i, j),
                                            &lB.vetMatBucketForwardClosed[v].mat(i, j)})
                {
                    for(int t = 0; t < bucket->sizeVetPtrLabel; ++t)
                    {
                        const Label *label = bucket->vetPtrLabel[t];
                        if(!label || !label->active)
                            continue;

                        const int load = (int)std::lround(label->vetResources[1]);
                        if(load > Q)
                            continue;

                        const double rc = label->vetResources[0];
                        list.push_back({rc, load, label});
                        bestB(v, load) = std::min(bestB(v, load), rc);
                    }
                }
            }
        }

        std::sort(list.begin(), list.end(),
                  [](const BackLabel &a, const BackLabel &b) { return a.rc < b.rc; });

        for(int r = 1; r <= Q; ++r)
            bestB(v, r) = std::min(bestB(v, r), bestB(v, r - 1));
    }

    // 2. G(u, R) = min over v of c(u, v) + bestB(v, R): best completion of a forward
    //    label at u with remaining capacity R (skips whole forward labels).
    static EigenMatrixRowD G;
    G.resize(n, Q + 1);
    G.setConstant(inf);
    for(int u = 0; u < n; ++u)
    {
        for(int v = 1; v < n; ++v)
        {
            if(v == u)
                continue;
            const double c = vetMatResCostForward(u, v, 0);
            if(c == InfFloatType)
                continue;
            for(int r = 0; r <= Q; ++r)
                G(u, r) = std::min(G(u, r), c + bestB(v, r));
        }
    }

    // 3. Candidate routes: the forward run's own routes, then the joins. Kept: the K
    //    best distinct routes; once there are K of them, the K-th reduced cost is the
    //    threshold a new join must beat.
    struct Candidate
    {
        double           rc;
        std::vector<int> route;    // 0 c1 ... ck 0
    };
    static std::vector<Candidate> cands;
    cands.clear();
    double threshold = -TolObjSubProb;

    auto compact = [&]()
    {
        std::sort(cands.begin(), cands.end(),
                  [](const Candidate &a, const Candidate &b) { return a.rc < b.rc; });
        std::vector<Candidate>       kept;
        std::set<std::vector<int>>   seen;
        for(Candidate &c : cands)
        {
            if((int)kept.size() == K)
                break;
            if(seen.insert(c.route).second)
                kept.push_back(std::move(c));
        }
        cands.swap(kept);
        if((int)cands.size() == K)
            threshold = std::min(threshold, cands.back().rc);
    };

    // Forward routes (already complete): numSolForward routes, each 0 ... 0.
    {
        int pos = 0;
        for(int k = 0; k < numSolForward && pos < sizeVetRoutesOut; ++k)
        {
            int end = pos + 1;
            while(end < sizeVetRoutesOut && vetRoutesOut[end] != 0)
                ++end;
            cands.push_back({(double)vetRedCostForward[k],
                             std::vector<int>(vetRoutesOut + pos, vetRoutesOut + end + 1)});
            pos = end + 1;
        }
    }

    std::vector<int> route;
    auto buildRoute = [&](const Label *F, const Label *B)
    {
        route.clear();
        for(const Label *l = F; l; l = l->ptrLabel)      // u ... 0
            route.push_back(l->cust);
        std::reverse(route.begin(), route.end());        // 0 ... u
        for(const Label *l = B; l; l = l->ptrLabel)      // v ... 0
            route.push_back(l->cust);
    };

    for(int u = 0; u < n; ++u)
    {
        for(int i = 0; i < numI; ++i)
        {
            for(int j = 0; j < numJ; ++j)
            {
                const Bucket &bucket = lF.vetMatBucketForwardClosed[u].mat(i, j);
                for(int t = 0; t < bucket.sizeVetPtrLabel; ++t)
                {
                    const Label *F = bucket.vetPtrLabel[t];
                    if(!F || !F->active || F->vetResources[1] > splitLoad)
                        continue;

                    const int R = Q - (int)std::lround(F->vetResources[1]);
                    if(R < 0)
                        continue;

                    const double rcF = F->vetResources[0];
                    if(rcF + G(u, R) >= threshold)
                        continue;

                    for(int v = 1; v < n; ++v)
                    {
                        if(v == u)
                            continue;
                        const double c = vetMatResCostForward(u, v, 0);
                        if(c == InfFloatType)
                            continue;

                        const double base = rcF + c;
                        if(base + bestB(v, R) >= threshold)
                            continue;

                        for(const BackLabel &b : vetBack[v])
                        {
                            if(base + b.rc >= threshold)
                                break;                // sorted by reduced cost
                            if(b.load > R)
                                continue;

                            auto common = F->bitSetNg & b.label->bitSetNg;
                            common.reset(0);
                            common.reset(n);
                            if(common.any())
                                continue;

                            buildRoute(F, b.label);
                            cands.push_back({base + b.rc, route});
                            if((int)cands.size() >= 4 * K)
                                compact();
                        }
                    }
                }
            }
        }
    }

    compact();

    // 4. Output.
    int next = 0;
    int num  = 0;
    for(const Candidate &c : cands)
    {
        if(next + (int)c.route.size() > NumMaxSolSubProb * NumMaxRoute)
            break;
        for(int node : c.route)
            vetRoutesOut[next++] = node;
        vetRedCost[num] = c.rc;
        ++num;
    }
    sizeVetRoutesOut = next;

    if(PrintIterInfoG)
        std::printf("merge: %d forward routes, %d routes returned, best rc %.6f\n",
                    numSolForward, num, num > 0 ? cands[0].rc : 0.0);

    return num;
}

/*
bool VrpTW_DecompLabelingNS::VrpLabelingSubProb::checkEnumeratedRoutesMid(
    Eigen::RowVectorXd &vetRowPi,
    Eigen::MatrixXd    &matColX,
    int                &numSol)
{

    // Remove all remaining labels
    labelingData.flushLabel();
    static LabelHeap labelHeap(1);

    for(int l = 0; l < numSol; ++l)
    {
        Label *label = getLabel();

        // Seting label data
        label->vetRoute[0] = 0;
        label->tamRoute = 1;
        label->typeLabel = LabelingAlgorithmNS::Forward;
        // label->vetResources.setZero();
        setVetResources0(*label);
        label->vetResources[0] = -vetRowPi[0];
        label->typeLabel = Forward;
        int lastCust = 0;

        // Seting route
        do
        {
            int cliJ = -1;
            for(int i = 0; i < instVrpTw->numClientes; ++i)
            {
                int index = getIndex(lastCust, i, instVrpTw->numClientes);
                if(matColX(index, l))
                {
                    cliJ = i;
                    break;
                }
            }

            if(cliJ == -1)
            {
                std::cout << "Erro!; diden't find arc (" << lastCust << ", _)\n";
                PRINT_EXIT();
            }

            label->vetResources[0] += instVrpTw->matDist(lastCust, cliJ) - vetRowPi[cliJ + 1];
            label->vetResources[1] += instVrpTw->vetClieDem[cliJ];
            label->vetRoute[label->tamRoute] = cliJ;
            label->tamRoute += 1;
            lastCust = cliJ;
        } while(lastCust != 0);

        label->vetRoute[label->tamRoute - 1] = instVrpTw->numClientes;
        label->cust = instVrpTw->numClientes;
        label->active = true;

        label->i = 0;
        label->j = 0;
        Bucket *bucket = labelingData.getBucket(label);
        bucket->addLabel(label);
    }

    Bucket &bucket = labelingData.vetMatBucketForward[instVrpTw->numClientes].mat(0, 0);
    if(bucket.sizeVetPtrLabel != numSol)
    {
        std::cout << "Error, bucket size(" << bucket.sizeVetPtrLabel << ") is different from "
                  << numSol << "\n\n";
        PRINT_EXIT();
    }

    bool addRoute = false;
    cleanVetRouteG();
    // Go through routes in hashRoutes
    for(const TestNS::Route &route : routeHash)
    {
        FloatType redCost = computeReducedCost(route, vetRowPi);
        if(redCost < -TolObjSubProb)
        {
            Label *label = getLabel();
            label->cust = instVrpTw->numClientes;
            label->vetResources[0] = redCost;
            label->vetResources[1] = route.demand;
            label->active = true;
            label->typeLabel = Forward;

            writeNgSet(label, ngSet);
            int     pos = -1;
            Bucket *bucket = dominanceIntraBucketSlow(instVrpTw->numClientes,
                                                      label,
                                                      labelingData,
                                                      &labelHeap,
                                                      2,
                                                      instVrpTw->numClientes,
                                                      pos);

            if(bucket != nullptr)
            {
                // std::cout<<"Error!!\nRoute("<<route.vetRoute<<") have a negative
                // reduced cost
                // ("<<redCost<<
                //            ") and was inserted into the bucket!\n\n";

                // Vector<VectorI> vetRoute;
                // vetRoute.push_back(route.vetRoute);
                // startGlobalMemory(vetRoute);
                // return true;
                addToVetRoutesG(route.vetRoute);
                addRoute = true;
            }
        }
    }

    return addRoute;
}
*/

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::resolveSubProb(
    const Eigen::VectorXd        &vetC,
    Eigen::RowVectorXd           &vetRowPi,
    //GRBModel                                               &mestre,
    int                           itCG,
    bool                         &custoRedNeg,
    void                         *data,
    const int                     iniConv,
    int                           indSubProb,
    Eigen::VectorXd              &vetCooefRestConv,
    const std::pair<int, int>    &pairSubProb,
    Eigen::MatrixXd              &matColX,
    int                          &numSol,
    ArrayReducedCost			 &vetRedCost,
    double                        constPiValue,
    const VectorI                &vetVar0,
    const VectorI                &vetVar1,
    //DW_DecompNS::PhaseStatus                                phaseStatus,
    bool                          exact)
{
    // std::cout<<"star value: "<<constPiValue<<"\n\n";
    static Eigen::VectorXd vetDist(instVrpTw->numClientes + 1);
    Eigen::VectorXd       *ptrVetDist = nullptr; //&vetDist;
    vetDist.setConstant(-MaxFloatType);

    // if(LowerBoundNS::getDistLowerBound(vetMatResCostForward, vetDist,
    // instVrpTw->numClientes, &labelingData))
    //     ptrVetDist = &vetDist;

    //static DW_DecompNS::PhaseStatus phase = phaseStatus;

    // if(phase != phaseStatus && phaseStatus ==
    // DW_DecompNS::PhaseStatus::PhaseStatusColGen)
    {

        minDistG = std::numeric_limits<FloatType>::max();
        maxDistG = std::numeric_limits<FloatType>::min();
    }

    //phase = phaseStatus;

    // std::cout<<"constPiValue: "<<constPiValue<<"\n";
    static Eigen::VectorX<FloatType> vetRedCostFT(NumMaxSolSubProb);

    // vetMatResCost.setVal(0.0);
    // constPiValue += -vetRowPi[0];

    for(int i = 0; i < instVrpTw->numClientes; ++i)
    {
        for(int j = 0; j < instVrpTw->numClientes; ++j)
        {
            if(i == j)
                continue;

            FloatType value0 = 0.0;

            FloatType dual = 0.0;
            if(j != 0)
                dual = (FloatType)vetRowPi[j];

            value0 = -dual;

            //if(phaseStatus != DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
            {
                value0 += (FloatType)instVrpTw->matDist(i, j); // - (FloatType)vetRowPi[j];
            }

            vetMatResCostForward.get(i, j, 0) = value0;
            vetMatResCostBackward.get(i, j, 0) = value0;
        }

        if(i != 0)
        {

            FloatType value = 0.0;

            // vetMatResCost[0](i, 0) = instVrpTw->matDist(i, 0);
            //if(phaseStatus == DW_DecompNS::PhaseStatus::PhaseStatusTwoPhase)
            //    value = 0.0;
            //else
            value = (FloatType)instVrpTw->matDist(i, 0);

            vetMatResCostForward.get(i, instVrpTw->numClientes, 0) = value;
            vetMatResCostBackward.get(i, instVrpTw->numClientes, 0) = value;
        }
    }

    /*
    if(LowerBoundNS::getDistLowerBound(vetMatResCostForward, vetDist,
    instVrpTw->numClientes, &labelingData))
    {
        std::printf("Seting Lower Bound\n");
        ptrVetDist = &vetDist;
    }
    */

    for(int i = 0; i < instVrpTw->numClientes; ++i)
    {
        if(!doubleEqual(vetMatResCostBackward(0, i, 0), vetMatResCostForward(0, i, 0)))
        {
            std::cout << "backward(" << vetMatResCostBackward(0, i, 0) << ") != forward("
                      << vetMatResCostForward(0, i, 0) << ")\n";
            exit(-1);
        }
    }

    // std::cout<<"("<<instVrpTw->numClientes<<",1):
    // "<<vetMatResCostBackward(instVrpTw->numClientes, 1, 0)<<"\n";
    // std::cout<<"(1,"<<instVrpTw->numClientes<<"): "<<vetMatResCostBackward(1,
    // instVrpTw->numClientes, 0)<<"\n";

    for(int varId : vetVar0)
    {
        int i = varId / instVrpTw->numClientes;
        int j = varId % instVrpTw->numClientes;

        if(j == 0)
            j = instVrpTw->numClientes;

        vetMatResCostForward.get(i, j, 0) = std::numeric_limits<FloatType>::infinity();
        // std::cout<<varId<<"\n";
    }

    std::set<int> setI;
    int           t = 0;

    /*
    for(int varId:vetVar1)
    {
        int i = varId/instVrpTw->numClientes;
        int j = varId%instVrpTw->numClientes;

        //std::cout<<i<<" "<<j<<"\n\n";

        if(setI.contains(i))
        {
            vetMatResCost[0](i, j) = (FloatType)instVrpTw->matDist(i, j) -
    (FloatType)(vetRowPi[j+1] + vetRowPi[1+instVrpTw->numClientes+t]); if(j == 0)
                vetMatResCost[0](i, instVrpTw->numClientes) = vetMatResCost[0](i, j);
        }
        else
        {
            for(int ii=0; ii < instVrpTw->numClientes+1; ++ii)
            {
                if(ii == i)
                    continue;

                vetMatResCost[0](ii, j) = std::numeric_limits<FloatType>::infinity();
            }


            if(j == 0)
            {

                for(int ii=0; ii < instVrpTw->numClientes; ++ii)
                {
                    if(ii == i)
                        continue;

                    vetMatResCost[0](ii, instVrpTw->numClientes) =
    std::numeric_limits<FloatType>::infinity();
                }
            }

            //vetMatResCost[0](i, j) += -vetRowPi[1+instVrpTw->numClientes+t];
            vetMatResCost[0](i, j) = (FloatType)instVrpTw->matDist(i, j) -
    (FloatType)(vetRowPi[j+1] + vetRowPi[1+instVrpTw->numClientes+t]); if(j == 0)
                vetMatResCost[0](i, instVrpTw->numClientes) = vetMatResCost[0](i, j);

            setI.insert(i);
        }

        t += 1;
    }
    */

    // std::cout<<"Custo Reduzido: \n"<<vetMatResCost[0]<<"\n\n";
    // std::cout<<"Peso: "<<vetMatResCost[1]<<"\n\n";

    int       it = 0;
    FloatType maxDist;

    ngSet.active = true;
    numSol = 0;

    exactLabelingG = false;

    // LowerBoundNS::getDistLowerBound(vetMatResCostForward, vetDist,
    // instVrpTw->numClientes, &labelingData);

    /*
    matColX.setZero();
    vetRedCostFT.setZero();
    bidirectionalAlgorithm(1, (instVrpTw->numClientes+1), vetMatResCostForward,
                           vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes,
    ngSet, labelingData, matColX, numSol, (FloatType) constPiValue, -1, true, maxDist,
    vetRedCostFT, true, typeLabel, true, ptrVetDist);

    std::cout<<vetRedCostFT.transpose()<<"\n";
    PRINT_EXIT();
    */

    G_numRecDominance = 2;
    // static Array<double, 3> arrayPercentage

    // TODO remove comments!
    if(!exact)
    {
        for(int i = 1; i <= 2; i += 1)
        {
            if(PrintIterInfoG)
                std::cout << i << ": ";
            // std::cout<<"forwardLabelingAlgorithm: "<<i<<"\n\n";
            matColX.setZero();
            /*
            custoRedNeg = bidirectionalAlgorithm(2,
                                                 (instVrpTw->numClientes + 1),
                                                 vetMatResCostForward,
                                                 vetMatResCostBackward,
                                                 vetVetResBound,
                                                 instVrpTw->numClientes,
                                                 ngSet,
                                                 labelingData,
                                                 matColX,
                                                 numSol,
                                                 (FloatType)constPiValue,
                                                 i,
                                                 true,
                                                 maxDist,
                                                 vetRedCostFT,
                                                 exact,
                                                 typeLabel,
                                                 true,
                                                 ptrVetDist,
                                                 1.0);
            */

            it += 1;
            if(custoRedNeg && numSol >= NumMaxSolSubProb) // 20
                break;
        }
    }

    if(!custoRedNeg || numSol < NumMaxSolSubProb)
    {
        matColX.setZero();
        /*
        custoRedNeg = bidirectionalAlgorithm(2,
                                             (instVrpTw->numClientes + 1),
                                             vetMatResCostForward,
                                             vetMatResCostBackward,
                                             vetVetResBound,
                                             instVrpTw->numClientes,
                                             ngSet,
                                             labelingData,
                                             matColX,
                                             numSol,
                                             (FloatType)constPiValue,
                                             -1,
                                             true,
                                             maxDist,
                                             vetRedCostFT,
                                             exact,
                                             typeLabel,
                                             true,
                                             ptrVetDist,
                                             .6);
        */
    }

    if(!custoRedNeg)
    {
        G_numRecDominance = 2;

        /*
        matColX.setZero();
        custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1),
        vetMatResCostForward, vetMatResCostBackward, vetVetResBound,
        instVrpTw->numClientes, ngSet, labelingData, matColX, numSol,
        (FloatType)constPiValue, -1, true, maxDist, vetRedCostFT, exact, typeLabel, true,
        ptrVetDist, .61);
        */

        matColX.setZero();
        /*
        custoRedNeg = bidirectionalAlgorithm(2,
                                             (instVrpTw->numClientes + 1),
                                             vetMatResCostForward,
                                             vetMatResCostBackward,
                                             vetVetResBound,
                                             instVrpTw->numClientes,
                                             ngSet,
                                             labelingData,
                                             matColX,
                                             numSol,
                                             (FloatType)constPiValue,
                                             -1,
                                             true,
                                             maxDist,
                                             vetRedCostFT,
                                             exact,
                                             typeLabel,
                                             true,
                                             ptrVetDist,
                                             1.0);
        */
    }

    if(!custoRedNeg)
    {
        if(checkEnumeratedRoutesFinal(vetRowPi))
        {
            PRINT_EXIT();
        }
    }
    /*
    else
    {

        if(checkEnumeratedRoutesMid(vetRowPi, matColX, numSol))
        {
            numSol = 0;
            custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1),
    vetMatResCostForward, vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes,
    ngSet, labelingData, matColX, numSol, (FloatType)constPiValue, -1, true, maxDist,
    vetRedCostFT, exact, typeLabel, true);

            //std::println("**************************************\n");
            //PRINT_EXIT();
        }
    }
    */

    /*
    if(!custoRedNeg)
    {
        std::cout<<"Change Alg Type\n";
        changeTypeAlg(typeLabel);

        custoRedNeg = bidirectionalAlgorithm(2, (instVrpTw->numClientes+1),
    vetMatResCostForward, vetMatResCostBackward, vetVetResBound, instVrpTw->numClientes,
    ngSet, labelingData, matColX, numSol, (FloatType)constPiValue, -1, true, maxDist,
    vetRedCostFT, exact, typeLabel);

        if(custoRedNeg)
        {
            std::cout<<"ERROR, alterar o algoritimo fez com que fosse gerado uma coluna
    com custo red. neg.\n";
        }

        changeTypeAlg(typeLabel);

    }
    */

    // redCost = (double)redCostFT;
    vetCooefRestConv[0] = 1;

    if(PrintIterInfoG)
        std::cout << "*************\n\n";

    /*    if(numSol == 0)
        {
            std::cout<<"CALL EXACT PRICING!\n\n";

            Eigen::Matrix<double, -1, 1, Eigen::ColMajor>
       matCol(instVrpTw->numClientes*instVrpTw->numClientes, 1);

            if(exactPricing(vetMatResCost, constPiValue, matCol, *instVrpTw,
       vetRedCostFT[0]))
            {

                custoRedNeg = forwardLabelingAlgorithm(2,
                                                       instVrpTw->numClientes+1,
                                                       vetMatResCost,
                                                       vetVetResBound,
                                                       instVrpTw->numClientes,
                                                       ngSet,
                                                       labelingData,
                                                       matColX,
                                                       numSol,
                                                       (FloatType)constPiValue,
                                                       -1,
                                                       true,
                                                       maxDist,
                                                       vetRedCostFT);

                if(custoRedNeg)
                {
                    std::cout<<"WTF?\n";
                    PRINT_DEBUG("", "");
                    throw "ERROR";
                }
                else
                {   std::cout<<"\nERROR, labeling falhou na ultima execucao!";
                    PRINT_DEBUG("", "");
                    throw "ERROR";
                }


                for(int i=0; i < instVrpTw->numClientes*instVrpTw->numClientes; ++i)
                    matColX(i, 0) = matCol[i];

                custoRedNeg = true;
                numSol = 1;



            }
        }*/

    // Check if solution have a negative reduced cost

    /*
    FloatType redCostTemp = 0.0;
    for(int j=0; j < numSol; ++j)
    {
        redCostTemp = constPiValue;
        //std::cout<<j<<": ";
        for(int i=0; i < (instVrpTw->numClientes*instVrpTw->numClientes); ++i)
        {
            if(matColX(i,j) != 1.0)
                continue;

            int ii = i/instVrpTw->numClientes;
            int jj = i%instVrpTw->numClientes;

            redCostTemp += vetMatResCostForward(ii, jj, 0);//instVrpTw->matDist(ii, jj);
            //std::cout<<"("<<ii<<","<<jj<<")["<<i<<"]; ";

        }
        //std::cout<<"\n\n";
        if(redCostTemp >= -DW_DecompNS::TolObjSubProb || !doubleEqual(vetRedCostFT[j],
    redCostTemp, 1E-3))
        {
            std::cout<<"\nERROR, custo reduzido calculado: ("<<redCostTemp<<") \n";
            std::cout<<"dif: "<<(vetRedCostFT[j]-redCostTemp)<<"\n";
            std::cout<<"redCost: "<<vetRedCostFT[j]<<"\n\n";
            PRINT_DEBUG("", "");
            std::cout<<"phaseStatus: "<<(int)phaseStatus<<"\n\n";
            std::cout<<(redCostTemp >= -DW_DecompNS::TolObjSubProb)<<"\n";
            throw "ERROR";
        }

    }
    */

    for(int i = 0; i < numSol; ++i)
        vetRedCost[i] = (double)vetRedCostFT[i];

    // std::cout<<"["<<minDistG<<", "<<maxDistG<<"]; numSol("<<numSol<<")\n";

    return 0;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::
    loadReducedCostsColuna(const double *vetReducedCostInput)
{
    const int n = instVrpTw->numClientes;          // nodes 0..n-1, dest = n

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            int k = getIndex(i, j, n);
            vetMatResCostForward.get(i, j, 0) = vetReducedCostInput[k];
        }
    }

    // Arc (i, 0) of Coluna is the arc (i, dest) of the labeling.
    for(int i = 1; i < n; ++i)
    {
        int k = getIndex(i, 0, n);
        vetMatResCostForward.get(i, n, 0) = vetReducedCostInput[k];
    }
}

/// Reduced-cost arc fixing. Exact forward labeling where a label is pruned only when
/// rc(L) + completion bound >= gap (instead of >= 0). For each arc (i,j) it records the
/// minimum of rc(L) + c̄(i,j) + matDp(j, Q - load(L) - d_j) over the labels L at i.
/// Any route through (i,j) has reduced cost >= vetArcMinOut[i*n+j], so the arc can be
/// removed if vetArcMinOut[i*n+j] > gap (decided by the caller).
double VrpTW_DecompLabelingNS::VrpLabelingSubProb::
    fixArcsColuna(const double *vetReducedCostInput, double gap, double *vetArcMinOut)
{
    const int n = instVrpTw->numClientes;          // nodes 0..n-1, dest = n

    loadReducedCostsColuna(vetReducedCostInput);

    int capacityGcd = 1;
    getDistLowerBound(vetMatResCostForward, n, vetVetResBound, matDpColuna, capacityGcd);
    matDpColunaValid = true;
    matDpColunaGcd   = capacityGcd;

    static ArcFixingData fixing;
    fixing.gap = (FloatType)gap;
    fixing.matArcMin.setConstant(n + 1, n + 1, InfFloatType);
    fixing.minRouteRc = InfFloatType;

    static Eigen::VectorX<FloatType> vetRedCostFT(NumMaxSolSubProb);
    static int vetRoutesDummy[NumMaxSolSubProb * NumMaxRoute];   // nothing is written
    int        sizeRoutesDummy = 0;
    int        numSol = 0;
    FloatType  maxDist;

    // Forward only: the per-arc minimum is recorded on forward extensions.
    bidirectionalAlgorithm(2,
                           n + 1,
                           vetMatResCostForward,
                           vetMatResCostBackward,
                           vetVetResBound,
                           n,
                           ngSet,
                           labelingDataForward,
                           numSol,
                           MaxInt,                  // no label limit per bucket
                           true,
                           maxDist,
                           vetRedCostFT,
                           true,                    // exact: never stops early
                           LabelingAlgorithmNS::AlgForward,
                           true,
                           &matDpColuna,
                           1.0,
                           vetRoutesDummy,
                           sizeRoutesDummy,
                           instVrpTw->capVeic,
                           capacityGcd,
                           nullptr,
                           &fixing);

    // Back to Coluna's indexing: arc (i, dest) is (i, 0).
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            const int jj = (j == 0) ? n : j;
            vetArcMinOut[getIndex(i, j, n)] =
                (i == j) ? InfDouble : (double)fixing.matArcMin(i, jj);
        }
    }

    fflush(stdout);
    std::cout.flush();

    return (double)fixing.minRouteRc;
}

int VrpTW_DecompLabelingNS::VrpLabelingSubProb::
    resolveSubProbColuna(double 		  *vetReducedCostInput,
                         int              *vetRoutesOut,
                         int              &sizeVetRoutesOut,
                         ArrayReducedCost &vetRedCost,
                         bool			   exact,
                         bool             &complete)
{
    // A heuristic result is never complete; the exact call sets it below.
    complete = !exact;

    //std::print("exact: {}\n", exact);
    static Eigen::VectorX<FloatType> vetRedCostFT(NumMaxSolSubProb);
    EigenMatrixRowD       *ptrMatDp = &matDpColuna;

    loadReducedCostsColuna(vetReducedCostInput);

    // The completion bound is computed only for the exact labeling: it is
    // expensive (O(n^2 * Q)) and the heuristic does not need it.
    int capacityGcd = 1;
    if(exact && ExactUseCompletionBound)
    {
        getDistLowerBound(vetMatResCostForward, instVrpTw->numClientes, vetVetResBound,
                          matDpColuna, capacityGcd);
        matDpColunaValid = true;
        matDpColunaGcd   = capacityGcd;
    }
    else
        ptrMatDp = nullptr;

    int numSol;
    FloatType maxDist;

    if(!exact)
    {
        // Reduced graph for the heuristic: every customer keeps only its
        // HeurNumBestArcs cheapest arcs to other customers. The copy keeps
        // resource 1 (demands), arcs out of the depot and arcs to the destination.
        const int n = instVrpTw->numClientes;          // nodes 0..n-1, dest = n
        vetMatResCostForwardHeur.copy(vetMatResCostForward);

        if(HeurNumBestArcs > 0 && HeurNumBestArcs < (n - 2))
        {
            static std::vector<std::pair<FloatType, int>> vetCand;

            for(int i = 1; i < n; ++i)
            {
                vetCand.clear();
                for(int j = 1; j < n; ++j)
                {
                    if(j == i)
                        continue;

                    const FloatType c = vetMatResCostForward(i, j, 0);
                    if(c != InfFloatType)               // skip arcs removed by branching
                        vetCand.emplace_back(c, j);
                }

                if((int)vetCand.size() <= HeurNumBestArcs)
                    continue;

                // The HeurNumBestArcs cheapest arcs end up in [0, HeurNumBestArcs).
                std::nth_element(vetCand.begin(),
                                 vetCand.begin() + HeurNumBestArcs,
                                 vetCand.end());

                for(size_t t = HeurNumBestArcs; t < vetCand.size(); ++t)
                    vetMatResCostForwardHeur.get(i, vetCand[t].second, 0) = InfFloatType;
            }
        }

        // Without a reduced graph, the "reduced" stages already run on the full
        // graph, so the full-graph stages would only repeat work: skip them.
        const bool reducedGraphActive = (HeurNumBestArcs > 0 && HeurNumBestArcs < (n - 2));

        int stageId = 0;
        const Vet3D_ResCost *boundGraph = nullptr;   // graph of matDpHeur (this call)
        int                  boundGcd   = 1;         // its capacity divisor
        for(const HeurStage &stage : HeurStages)
        {
            ++stageId;

            if(!stage.reducedGraph && !reducedGraphActive)
                continue;

            const Vet3D_ResCost &graph =
                stage.reducedGraph ? vetMatResCostForwardHeur : vetMatResCostForward;

            // Completion bound of the heuristic (see HeurCompletionBound).
            ptrMatDp    = nullptr;
            capacityGcd = 1;
            if(HeurCompletionBound == 1)
            {
                if(boundGraph != &graph)
                {
                    getDistLowerBound(graph, instVrpTw->numClientes, vetVetResBound,
                                      matDpHeur, boundGcd);
                    boundGraph = &graph;
                }
                capacityGcd = boundGcd;
                ptrMatDp    = &matDpHeur;
            }
            else if(HeurCompletionBound == 2 && matDpColunaValid)
            {
                ptrMatDp    = &matDpColuna;
                capacityGcd = matDpColunaGcd;
            }

            bool reducedCostNeg =
            bidirectionalAlgorithm(2,
                                   (instVrpTw->numClientes + 1),
                                   graph,
                                   vetMatResCostBackward,
                                   vetVetResBound,
                                   instVrpTw->numClientes,
                                   ngSet,
                                   labelingDataForward,
                                   //matColX,
                                   numSol,
                                   //(FloatType)constPiValue,
                                   stage.labelsPerBucket,
                                   true,
                                   maxDist,
                                   vetRedCostFT,
                                   exact,
                                   typeLabel,
                                   true,
                                   ptrMatDp,
                                   1.0,
                                   vetRoutesOut,
                                   sizeVetRoutesOut,
                                   instVrpTw->capVeic,
                                   capacityGcd,
                                   nullptr,
                                   nullptr,
                                   stage.useNg);

            if(reducedCostNeg)
            {
                if(PrintIterInfoG)
                    std::printf("heuristic stage %d found %d columns\n", stageId, numSol);

                fflush(stdout);
                std::cout.flush();
                return numSol;
            }
        }

        sizeVetRoutesOut = 0;
        fflush(stdout);
        std::cout.flush();
        return 0;
    }

    // Bidirectional labeling (exact pricing): forward run up to the split load, backward
    // run (forward labeling on the transposed graph) up to the split load, then the merge
    // of the two writes the routes that go back to Coluna.
    const bool      bidirMerge = BidirExactBackwardRun && !BidirBackwardFullCheck;
    const FloatType splitLoad  = (FloatType)(BidirSplitFraction * instVrpTw->capVeic);
    const FloatType forwardStopLoad = bidirMerge ? splitLoad : (FloatType)InfFloatType;

    // Forward run.
    bidirectionalAlgorithm(2,
                           (instVrpTw->numClientes + 1),
                           vetMatResCostForward,
                           vetMatResCostBackward,
                           vetVetResBound,
                           instVrpTw->numClientes,
                           ngSet,
                           labelingDataForward,
                           //matColX,
                           numSol,
                           //(FloatType)constPiValue,
                           MaxInt,
                           true,
                           maxDist,
                           vetRedCostFT,
                           exact,
                           typeLabel,
                           true,
                           ptrMatDp,
                           1.0,
                           vetRoutesOut,
                           sizeVetRoutesOut,
                           instVrpTw->capVeic,      // full Q: only used by the matDp pruning
                           capacityGcd,
                           &complete,
                           nullptr,
                           true,
                           true,
                           // Bidirectional: stops extending at the split load (labels
                           // above it are kept, not extended); the merge below completes
                           // the routes. Otherwise (or in the check mode): full Q.
                           forwardStopLoad);

    if(BidirExactBackwardRun)
    {
        // Backward run (forward labeling on the transposed graph). Its own routes are
        // not used (scratch variables): the merge uses its labels, which stay alive in
        // labelingDataBackward (and the forward ones in labelingDataForward, since the
        // label pool is not reset).
        static Eigen::VectorX<FloatType> vetRedCostBackward(NumMaxSolSubProb);
        static std::vector<int>          vetRoutesBackward(NumMaxSolSubProb * NumMaxRoute);
        int       numSolBackward     = 0;
        int       sizeRoutesBackward = 0;
        bool      completeBackward   = false;
        FloatType maxDistBackward;

        updateVetMatResCostBackward();

        // Completion bound of the transposed graph (the forward matDp is not valid
        // here unless the reduced costs are symmetric). getDistLowerBound keeps the
        // demands and the capacity in static variables from its first call; they are
        // the same for both graphs.
        int capacityGcdBackward = 1;
        if(ExactUseCompletionBound)
            getDistLowerBound(vetMatResCostTransposed, instVrpTw->numClientes, vetVetResBound,
                              matDpBackward, capacityGcdBackward);

        // Labels with a load above this are kept but not extended.
        const FloatType backwardStopLoad =
            BidirBackwardFullCheck ? (FloatType)InfFloatType : splitLoad;

        bidirectionalAlgorithm(2,
                               (instVrpTw->numClientes + 1),
                               vetMatResCostTransposed,
                               vetMatResCostTransposed,
                               vetVetResBound,
                               instVrpTw->numClientes,
                               ngSet,
                               labelingDataBackward,
                               numSolBackward,           // scratch
                               MaxInt,
                               true,
                               maxDistBackward,          // scratch
                               vetRedCostBackward,       // scratch
                               exact,
                               typeLabel,
                               true,
                               ExactUseCompletionBound ? &matDpBackward : nullptr,
                               1.0,
                               vetRoutesBackward.data(), // scratch
                               sizeRoutesBackward,       // scratch
                               instVrpTw->capVeic,       // full Q: only used by the matDp pruning
                               capacityGcdBackward,
                               &completeBackward,        // scratch
                               nullptr,
                               true,
                               false,                    // keep the forward labels
                               backwardStopLoad);

        if(bidirMerge)
        {
            numSol   = mergeBidirectional(splitLoad, numSol, vetRedCostFT,
                                          vetRoutesOut, sizeVetRoutesOut, vetRedCost);
            complete = complete && completeBackward;
        }

        if(BidirBackwardFullCheck)
        {
            // Best route of each run (routes are written one after the other, each
            // as 0 c1 ... ck 0) and whether it visits a customer twice. ng-routes
            // with a cycle can differ between the two directions, so the runs must
            // agree only when the better route is elementary.
            auto best = [](const int *vet, int size, const Eigen::VectorX<FloatType> &rc,
                           int num, double &bestRc, bool &elementary)
            {
                bestRc = 0.0;
                elementary = true;
                if(num <= 0)
                    return;

                int bestK = 0;
                for(int k = 1; k < num; ++k)
                    if(rc[k] < rc[bestK])
                        bestK = k;
                bestRc = rc[bestK];

                int k = -1, pos = 0;
                std::vector<int> route;
                while(pos < size)
                {
                    // Route starts at a 0; it ends at the next 0.
                    int end = pos + 1;
                    while(end < size && vet[end] != 0)
                        ++end;
                    ++k;
                    if(k == bestK)
                    {
                        route.assign(vet + pos + 1, vet + end);
                        break;
                    }
                    pos = end + 1;
                }
                std::sort(route.begin(), route.end());
                elementary = std::adjacent_find(route.begin(), route.end()) == route.end();
            };

            double fwdRc, bwdRc;
            bool   fwdElem, bwdElem;
            best(vetRoutesOut, sizeVetRoutesOut, vetRedCostFT, numSol, fwdRc, fwdElem);
            best(vetRoutesBackward.data(), sizeRoutesBackward, vetRedCostBackward,
                 numSolBackward, bwdRc, bwdElem);

            const double tol = 1e-6;
            // A run is wrong if the other one found a better elementary route.
            const bool bug = (fwdElem && numSol > 0 && bwdRc > fwdRc + tol) ||
                             (bwdElem && numSolBackward > 0 && fwdRc > bwdRc + tol);
            const char *status = bug ? "MISMATCH"
                               : (std::abs(fwdRc - bwdRc) > tol ? "differ (ng cycle)" : "ok");

            std::printf("BIDIR check: %s; forward %d routes best %.6f%s, backward %d routes best %.6f%s\n",
                        status, numSol, fwdRc, fwdElem ? "" : " (cycle)",
                        numSolBackward, bwdRc, bwdElem ? "" : " (cycle)");
        }
    }

    fflush(stdout);
    std::cout.flush();

    return numSol;

}


void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToForward()
{
    typeLabel = LabelingAlgorithmNS::AlgForward;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToBackward()
{
    typeLabel = LabelingAlgorithmNS::AlgBackward;
}

void VrpTW_DecompLabelingNS::VrpLabelingSubProb::setTypeLabelToBidirectional()
{
    typeLabel = LabelingAlgorithmNS::AlgBidirectional;
}

/*
void VrpTW_DecompLabelingNS::VrpLabelingSubProb::convertRouteIntoLabel(
    const TestNS::Route        &route,
    LabelingAlgorithmNS::Label *label)
{
    // label->vetResources.setZero();
    label->bitSetNg = 0;
    for(int i = 0; i < (int)route.vetRoute.size(); ++i)
    {
        int cliI = route.vetRoute[i];
        label->vetRoute[i] = cliI;

        if((i + 1) < (int)route.vetRoute.size())
        {
            int cliJ = route.vetRoute[i + 1];

            for(int r = 0; r < 2; ++r)
            {
                label->vetResources[r] += vetMatResCostForward(cliI, cliJ, r);
            }

            if(ngSet.contain(cliJ, cliI))
                label->bitSetNg[cliJ] = 1;
        }
    }

    label->tamRoute = (int)route.vetRoute.size();
}
*/

/*
VrpTW_DecompLabelingNS::CapacityCut::CapacityCut(InstanceVRPTW_NS::InstanceVRPTW &instVrpTw,
                                                 int                              dim_,
                                                 int                              maxNoOfCuts_,
                                                 double                           eps)
{
    capacity = instVrpTw.capVeic;
    dim = dim_;
    maxNoOfCuts = maxNoOfCuts_;
    numCust = instVrpTw.numClientes - 1; // Don't count with the deposit
    edgeSize = ((numCust + 1) * (numCust)) / 2;
    epsForIntegrality = eps;

    std::printf("EdgeSize: %d\n\n", edgeSize);

    edgeHead = new int[edgeSize];
    edgeTail = new int[edgeSize];
    edgeX = new double[edgeSize];
    demand = new int[numCust + 1];

    // std::cout<<"Demand: ";
    for(int i = 0; i <= numCust; ++i)
    {
        demand[i] = instVrpTw.vetClieDem[i];
        // std::printf("%d ", demand[i]);
    }

    // std::println();

    CMGR_CreateCMgr(&cutsCMP, dim);
    CMGR_CreateCMgr(&oldCutsCMP, dim);

    bool createMap = true;
    // if(createMap)
    mapArcToIndex = new std::map<std::pair<int, int>, int>();

    int next = 0;
    for(int i = 0; i < numCust + 1; ++i)
    {
        for(int j = (i + 1); j < numCust + 1; ++j)
        {
            if(next >= edgeSize)
            {
                std::printf("Error, next(%i) >= edgeSize(%i)\n", next, edgeSize);
                PRINT_EXIT();
            }

            if(createMap)
            { // std::printf("(%d; %d): %i\n", i, j, next);
                (*mapArcToIndex)[{i, j}] = next;
            }

            edgeHead[next] = i;
            edgeTail[next] = j;
            // std::printf("Add: (%d;%d)\n", i, j);

            next += 1;
        }
    }

    // std::printf("\n\n");

    if(next != edgeSize)
        std::printf("next(%d) != edgeSize(%d)\n\n", next, edgeSize);
    if(createMap)
        list = Eigen::VectorXi(numCust + 1);
}

VrpTW_DecompLabelingNS::CapacityCut::~CapacityCut()
{
    delete[] edgeHead;
    delete[] edgeTail;
    delete[] edgeX;

    CMGR_FreeMemCMgr(&cutsCMP);
    CMGR_FreeMemCMgr(&oldCutsCMP);
}
int VrpTW_DecompLabelingNS::CapacityCut::operator()(DW_DecompNS::DW_DecompNode &decompNode)
{
    std::printf("*************************************************************************\n");
    std::printf("***************************CAPACITY_CUT**********************************\n\n");

    // Convert the arc solution for the edge one
    convertArcSolution(decompNode.vetSolX);

    // Separate the capacity cuts
    CAPSEP_SeparateCapCuts(numCust,
                           demand,
                           capacity,
                           edgeSize,
                           edgeTail,
                           edgeHead,
                           edgeX,
                           oldCutsCMP,
                           maxNoOfCuts,
                           epsForIntegrality,
                           &integerAndFeasible,
                           &maxViolation,
                           cutsCMP);

    list.setConstant(-1);

    // Runs through the cuts
    for(int i = 0; i < cutsCMP->Size; ++i)
    {

        std::printf("Add cut %d: \n\n\t", i);

        int listSize = 0;
        // Recover the customers
        for(int j = 0; j < cutsCMP->CPL[i]->IntListSize; ++j)
        {
            list[listSize] = cutsCMP->CPL[i]->IntList[j];
            std::printf("%d ", list[listSize]);
            listSize += 1;
        }

        double rhs = cutsCMP->CPL[i]->RHS;
        std::printf("< %f", rhs);
        RobustCut robustCut(decompNode.vetSolX.size());
        robustCut.sense = '<';
        robustCut.rhs = rhs;

        // Add the arcs of the clique of nodes of the list
        createInducedSubGraphArcs(listSize, robustCut);
        addMasterCut(robustCut, decompNode, i, false);
    }

    // Move cuts in cutsCMP for oldCutsCMP
    for(int i = 0; i < cutsCMP->Size; ++i)
        CMGR_MoveCnstr(cutsCMP, oldCutsCMP, i, 0);

    int size = cutsCMP->Size;
    cutsCMP->Size = 0;

    return size;
};

void VrpTW_DecompLabelingNS::CapacityCut::convertArcSolution(const Eigen::VectorXd &vetX)
{

    int numNodes = numCust + 1;

    for(int i = 0; i < edgeSize; ++i)
        edgeX[i] = 0.0;

    int index;

    for(int i = 0; i < numNodes; ++i)
    {
        for(int j = 0; j < numNodes; ++j)
        {
            // reverseIndex(t, numNodes, i, j);
            if(i == j)
                continue;

            index = (*mapArcToIndex).at(getEdge(i, j));
            edgeX[index] += vetX[getIndex(i, j, numNodes)];

            if(edgeX[index] > 1.0)
            {
                if(j > i)
                    std::cout << "Edge: (" << i << ", " << j << ") have a value greter then 1\n";
                else
                    std::cout << "Edge: (" << j << ", " << i << ") have a value greter then 1\n";
            }
        }
    }

    //    for(int i=0; i < edgeSize; ++i)
    //        edgeX[i] = std::min(edgeX[i], 1.0);

    std::printf("edgeX: ");
    for(int i = 0; i < edgeSize; ++i)
        std::printf("%.1f ", edgeX[i]);

    std::cout << "\n\n";
}

std::pair<int, int> VrpTW_DecompLabelingNS::CapacityCut::getEdge(int i, int j)
{

    if(i < j)
        return {i, j};
    else
        return {j, i};
}

void VrpTW_DecompLabelingNS::CapacityCut::createInducedSubGraphArcs(
    int                          listSize,
    BranchAndPriceNS::RobustCut &cut)
{

    int numNodes = numCust + 1;

    for(int i = 0; i < listSize; ++i)
    {
        for(int j = 0; j < listSize; ++j)
        {
            if(i == j)
                continue;

            int index0 = getIndex(i, j, numNodes);
            int index1 = getIndex(j, i, numNodes);

            cut.vetX.coeffRef(index0) = 1.0;
            cut.vetX.coeffRef(index1) = 1.0;
        }
    }
}
*/
