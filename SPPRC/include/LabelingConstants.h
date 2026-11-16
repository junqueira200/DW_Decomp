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

constexpr int       NumMaxResources = 2;
constexpr int       NumMaxRoute = 100;
constexpr int       NumMaxCust = 128; // 100; 128 uses the same 2 words of bitset<>
constexpr int       NgSetSize = 8;    // 8
constexpr int       NumBuckets = 10;
constexpr int       vetPtrLabelSize = 10;
constexpr bool      NullFlush = false;
constexpr bool      Print = false;
constexpr int       numMaxLabelG = 40; // 2000
constexpr bool      DominaIterBuckets = false;
constexpr FloatType FloatEp = 1E-6;
constexpr bool      TrackingRoutes = false;
constexpr bool      SortLabels = false;
// (ClosedBucketWindow removed 2026-11-15: the closed buckets are now scanned in all
// lower buckets, with the best reduced cost per row (vetMatClosedCBest) to stop early.)
// TEMP test: false = extended labels are not kept in closed buckets and never checked.
constexpr bool      UseClosedBuckets = true;

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

// Exact pricing: also run the labeling on the transposed graph (bidirectional
// labeling, work in progress: no merge yet, its routes are not returned).
constexpr bool BidirExactBackwardRun = true;
// Fraction of the capacity where the backward run stops extending labels (labels
// above it are kept, not extended). The forward run will stop there too once the
// merge exists; until then it must explore up to Q (its routes go to Coluna).
constexpr double BidirSplitFraction = 0.5;
// Debug: run the backward labeling up to Q (no split) and compare its best route
// with the forward one on every exact call (prints one "BIDIR check" line per call).
constexpr bool BidirBackwardFullCheck = false;

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
