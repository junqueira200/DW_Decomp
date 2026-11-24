#ifndef LABELINGCONSTANTS_H
#define LABELINGCONSTANTS_H

#include "safe_3D_matrix.h"
#include <Eigen/Core>
#include <cfenv>
#include <iostream>
#include <vector>

typedef double FloatType;

#define HEAP_KEY vetResources[0]
#define KEY_TYPE FloatType
#define KEY_EXTREME -std::numeric_limits<KEY_TYPE>::max()

constexpr bool MinHeap = true;

#define INLINE inline __attribute__((always_inline))

namespace LabelingAlgorithmNS
{
inline FloatType minDistG = std::numeric_limits<FloatType>::max();
inline FloatType maxDistG = std::numeric_limits<FloatType>::min();
inline bool      exactLabelingG = false;
inline bool      PrintG = false;
// Prints per-call information of the pricing (number of iterations, routes, ...)
inline bool      PrintIterInfoG = false;

constexpr int       NumMaxResources 	= 2;
constexpr int       NumMaxRoute 		= 100;
constexpr int       NumMaxCust 			= 128; // 100; 128 uses the same 2 words of bitset<>
constexpr int       MaxActiveCuts 		= 64;
constexpr int       NgSetSize 			= 8;    // 8
constexpr int       NumBuckets 			= 10;
constexpr int       vetPtrLabelSize 	= 10;
constexpr bool      NullFlush 			= false;
constexpr bool      Print 				= false;
constexpr int       numMaxLabelG 		= 40; // 2000
constexpr bool      DominaIterBuckets 	= false;
constexpr FloatType FloatEp 			= 1E-6;
constexpr bool      TrackingRoutes 		= false;
constexpr bool      SortLabels 			= false;
// (ClosedBucketWindow removed 2026-11-15: the closed buckets are now scanned in all
// lower buckets, with the best reduced cost per row (vetMatClosedCBest) to stop early.)
// TEMP test: false = extended labels are not kept in closed buckets and never checked.
constexpr bool      UseClosedBuckets 	= true;

// NOT USED since 2026-11-17 (removed from funcStop: the heuristic now stops only at
// HeurMaxSolSubProb columns). Kept for reference.
// Early stop of the heuristic labeling (exact labeling never stops early).
// Each rule {minIt, minCols}: stop when at least minIt labels were processed and
// at least minCols columns were found.
struct HeurStopRule
{
    int minIt;
    int minCols;
};

constexpr HeurStopRule HeurStopRules[] = {
    {500, 20},
    {1000, 15},
    {10000, 10},
    {100000, 5},
    {1000000, 1},
};

// Heuristic labeling only: each customer keeps only its HeurNumBestArcs outgoing
// arcs (to other customers) with the lowest reduced cost. Arcs to the destination
// and arcs out of the depot are always kept. 0 = use the full graph.
// (see also HeurStages below)
constexpr int HeurNumBestArcs = 10;

// Heuristic pricing cascade: the stages run in this order until one of them finds
// columns with negative reduced cost. Exact pricing runs only if all stages fail.
//  - reducedGraph:    true  = graph with only the HeurNumBestArcs best arcs per customer
//                     false = full graph
//  - labelsPerBucket: maximum number of labels kept per bucket (NumMaxLabePerBucket)
//  - useNg:           true  = full dominance test (ng-sets included), as in the exact
//                             labeling
//                     false = dominance ignores the ng-sets (reduced cost and load only):
//                             more labels are pruned and each test is cheaper; routes stay
//                             ng-feasible (the ng-set is still checked on extension)
struct HeurStage
{
    bool reducedGraph;
    int  labelsPerBucket;
    bool useNg;
};

constexpr HeurStage HeurStages[] = {
    {true,  2, false},  // reduced graph, cheap, no ng in dominance
    {true,  4, false},  // reduced graph, more labels, no ng in dominance
    {true,  4, true},   // reduced graph, full dominance: catches what the stages above miss
    // {false, 4, true},   // full graph: catches what the reduced graph misses (off)
};

// Exact labeling: stop as soon as at least ExactStopMinCols columns with negative
// reduced cost were found. The search is then incomplete, so its best reduced cost is
// not a valid dual bound (the caller is told through `complete`).
// 0 = never stop early (always explore everything, as before).
constexpr int ExactStopMinCols = 0;

// Exact pricing with bidirectional labeling: forward run and backward run (forward
// labeling on the transposed graph) both stop at BidirSplitFraction * Q, then
// mergeBidirectional() joins them. Exact (same best route as the forward-only
// labeling on all recorded inputs), but OFF: each half costs 70-90% of a full run, so
// it is slower than ExactDestKeepAll (E-n76-k7 4-min tree, 2026-11-16: bound 672.85
// vs 672.92, 47 vs 53 nodes). false = forward-only exact labeling up to Q.
constexpr bool BidirExactBackwardRun = false;
// Split load of the bidirectional labeling, as a fraction of the capacity: labels
// above it are kept (for the merge) but not extended.
constexpr double BidirSplitFraction = 0.5;
// Debug (with BidirExactBackwardRun): both runs go up to Q, no merge; prints one
// "BIDIR check" line per exact call comparing their best routes.
constexpr bool BidirBackwardFullCheck = false;

// Exact (forward-only) labeling keeps every negative route at the destination (no
// dominance there), so up to NumMaxSolSubProb routes go back per exact call.
constexpr bool ExactDestKeepAll = true;
// Set by bidirectionalAlgorithm on each call (exact && ExactDestKeepAll).
inline bool destNoDominanceG = false;

// Exact pricing: prune labels with the completion bound (matDp). false = no completion
// bound in the exact labeling (forward and backward runs); test only.
constexpr bool ExactUseCompletionBound = true;

// Completion bound in the heuristic labeling (it only has to find some negative
// columns, so an invalid prune is allowed: exact pricing catches what it misses).
//   0 = none;
//   1 = computed on every heuristic call, on the graph the stage searches (reduced
//       graph or full graph): a valid bound for that graph;
//   2 = the bound of the last exact call (or arc fixing run), reused as it is: free,
//       but computed from older duals.
// 2026-11-16: mode 1 made the heuristic ~3x slower (the bound costs more than a call);
// mode 2: heuristic -35% on the E-n76-k7 4-min tree (19.4 vs 30.0 s), 55 vs 52 nodes,
// bound 673.17 vs 672.86.
constexpr int HeurCompletionBound = 2;

constexpr double MaxFloatType = std::numeric_limits<FloatType>::max();
constexpr double MinFloatType = std::numeric_limits<FloatType>::min();
constexpr double InfFloatType = std::numeric_limits<FloatType>::infinity();

constexpr double InfDouble    = std::numeric_limits<double>::infinity();

constexpr int MaxInt = std::numeric_limits<int>::max();
constexpr int MinInt = std::numeric_limits<int>::min();
constexpr int InfInt = std::numeric_limits<int>::infinity();

constexpr int64_t MaxInt64 = std::numeric_limits<int64_t>::max();
constexpr int64_t MinInt64 = std::numeric_limits<int64_t>::min();
constexpr int64_t InfInt64 = std::numeric_limits<int64_t>::infinity();

struct Bound
{
    FloatType lowerBound = -std::numeric_limits<FloatType>::infinity();
    FloatType upperBound = std::numeric_limits<FloatType>::infinity();
};

/// Access (i,j,r) where (i,j) is an arc and r a resource
typedef Vector3D<FloatType, false> Vet3D_ResCost;
/// Access (i,r) where i is a customer and r a resource
typedef Eigen::Matrix<Bound, -1, -1, Eigen::RowMajor> MatBoundRes;
typedef Eigen::Array<FloatType, 1, NumMaxResources>   ArrayResources;
std::ostream &operator<<(std::ostream &out, const Bound &bound);

} // namespace LabelingAlgorithmNS
#endif // LABELINGCONSTANTS_H
