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
constexpr int       NumMaxCust = 100; // 100
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
// How many buckets back (in each resource) the closed buckets are scanned when
// checking if a new label is dominated. 0 = own bucket only.
constexpr int       ClosedBucketWindow = 1;

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
struct HeurStage
{
    bool reducedGraph;
    int  labelsPerBucket;
};

constexpr HeurStage HeurStages[] = {
    {true,  2},     // reduced graph, cheap
    {true,  4},     // reduced graph, more labels
    // {false, 4},  // full graph: catches what the reduced graph misses (cascade, off)
};

// Exact labeling: stop as soon as at least ExactStopMinCols columns with negative
// reduced cost were found. The search is then incomplete, so its best reduced cost is
// not a valid dual bound (the caller is told through `complete`).
// 0 = never stop early (always explore everything, as before).
constexpr int ExactStopMinCols = 0;

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
