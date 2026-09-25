#ifndef DUALFEASIBLEFUNCTIONS_H
#define DUALFEASIBLEFUNCTIONS_H

#include "safe_vector.h"

namespace DualFeasibleFunctionsNS
{

    double funcionU(double ep, double x);
    bool check(const VectorI& vetItems, int sizeVetItems=-1);
}


#endif // DUALFEASIBLEFUNCTIONS_H
