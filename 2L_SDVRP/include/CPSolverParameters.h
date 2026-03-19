#pragma once

namespace ContainerLoading
{
namespace Algorithms
{
struct CPSolverParams
{
    int Threads = 4;
    int Seed = 0;
    bool LogFlag = false;
    bool Presolve = true;

    bool EnableCumulativeDimensions = false;
    bool EnableNoOverlap2DFloor = false;
};

}
}
