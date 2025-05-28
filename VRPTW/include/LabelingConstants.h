#ifndef LABELINGCONSTANTS_H
#define LABELINGCONSTANTS_H

#include <vector>
#include <iostream>
#include <Eigen/Eigen>
#include "safe_3D_matrix.h"

typedef double FloatType;

namespace LabelingAlgorithmNS
{
    inline std::vector<int> vetRoteG;
    inline FloatType minDistG = std::numeric_limits<FloatType>::max();
    inline FloatType maxDistG = std::numeric_limits<FloatType>::min();
    inline bool exactLabelingG = false;

    constexpr int   NumMaxResources   = 2;
    constexpr int   NumMaxRoute       = 256;
    constexpr int   NumMaxCust        = 25; //100
    constexpr int   NgSetSize         = 5;
    constexpr int   NumBuckets        = 10;
    constexpr int   vetPtrLabelSize   = 10;
    constexpr bool  NullFlush         = true;
    constexpr bool  Print             = false;
    constexpr int   numMaxLabelG      = 2000; // 2000
    constexpr bool  DominaIterBuckets = true;
    constexpr FloatType FloatEp       = 1E-4;

    struct Bound
    {
        FloatType lowerBound = -std::numeric_limits<FloatType>::infinity();
        FloatType upperBound =  std::numeric_limits<FloatType>::infinity();
    };


    /// Access (i,j,r) where (i,j) is an arc and r a resource
    typedef Vector3D<FloatType, false> Vet3D_ResCost;
    /// Access (i,r) where i is a customer and r a resource
    typedef Eigen::Matrix<Bound, -1, -1, Eigen::RowMajor> MatBoundRes;
    typedef Eigen::Array<FloatType, 1, NumMaxResources> ArrayResources;
    std::ostream& operator<< (std::ostream& out, const Bound &bound);

}
#endif // LABELINGCONSTANTS_H
