#include "DualFeasibleFunctions.h"
#include "AuxT.h"
#include <iostream>

using namespace DualFeasibleFunctionsNS;


double DualFeasibleFunctionsNS::funcionU(double ep, double x)
{
    if(ep > 0.5 || ep < 0.0)
    {
        std::printf("Erro, ep shuld be bething [0;0.5]\n");
        PRINT_THROW()
    }

    if(x < ep)
        return 0;

    else if(x >= ep && x <= (1-ep))
        return x;

    else// if(x > (1-ep))
        return 1;

    /*
    else
    {

        std::printf("Erro, it shudend be here; x: %f; ep: %f\n\n", x, ep);
        PRINT_THROW();
        return -1;
    }
    */

}
