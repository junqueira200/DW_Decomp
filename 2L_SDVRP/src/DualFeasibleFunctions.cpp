#include "DualFeasibleFunctions.h"
#include "AuxT.h"
#include "c_api.h"

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

bool DualFeasibleFunctionsNS::check(const VectorI &vetItems, int sizeVetItems)
{
    static const Array<double, 5> arrayEp{1.0/3, 1.0/4, 1.0/5, 1.0/6, 1.0/7};

    if(sizeVetItems <= 0)
        sizeVetItems = vetItems.size();

    for(double ep0:arrayEp)
    {
        for(double ep1:arrayEp)
        {
            setDualFeasibleFunction(ep0, ep1);
            double vol = 0.0;

            for(int i=0; i < sizeVetItems; ++i)
                vol += InstanceNS::instanciaG.vetItens[vetItems[i]].volNormDual;

            if(doubleGreater(vol, 1.0))
            {
                std::printf("Dual function detects a infeasible packing! vol: %.2f\n", vol);
                doBreakTestRoute = true;
                return false;
            }
        }
    }

    return true;

}
