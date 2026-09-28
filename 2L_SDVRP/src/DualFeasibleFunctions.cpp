#include "DualFeasibleFunctions.h"
#include "AuxT.h"
#include "c_api.h"
#include "Instancia.h"

#include <iostream>

using namespace DualFeasibleFunctionsNS;
using namespace InstanceNS;

double DualFeasibleFunctionsNS::funcionU(double ep, double x)
{
    if(ep > 0.5 || ep < 0.0)
    {
        std::printf("Erro, ep shuld be bething [0;0.5]\n");
        PRINT_THROW()
    }

    constexpr double tol = 1E-9;

    if(x < ep - tol)
        return 0.0;

    else if(x > (1-ep) - tol)
        return 1.0;

    else// if(x > (1-ep))
        return x;

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
    static const Array<double, 5> arrayEp{1.0/3.0, 1.0/4.0, 1.0/5.0, 1.0/6.0, 1.0/7.0};

    if(sizeVetItems <= 0)
        sizeVetItems = vetItems.size();

    for(double ep0:arrayEp)
    {
        for(double ep1:arrayEp)
        {
            for(double ep2:arrayEp)
            {
                //setDualFeasibleFunction(ep0, ep1);
                double vol = 0.0;

                for(int i=0; i < sizeVetItems; ++i)
                {
                    Item& item = instanciaG.vetItens[vetItems[i]];
                    double d0 = funcionU(ep0, item.vetDimNor[0]);
                    double d1 = funcionU(ep1, item.vetDimNor[1]);
                    double d2 = funcionU(ep2, item.vetDimNor[2]);
                    vol += d0*d1*d2;
                }

                if(doubleGreater(vol, 1.0, 1E-5))
                {
                    //std::printf("Dual function detects a infeasible packing! vol: %.2f\n", vol);
                    //std::printf("ep0(%f), ep1(%f), ep2(%f)\n", ep0, ep1, ep2);


                    for(int i=0; i < sizeVetItems; ++i)
                    {
                        //std::printf("Item %d: ", vetItems[i]);

                        Item& item = instanciaG.vetItens[vetItems[i]];
                        double d0 = funcionU(ep0, item.vetDimNor[0]);
                        double d1 = funcionU(ep1, item.vetDimNor[1]);
                        double d2 = funcionU(ep2, item.vetDimNor[2]);

                        //std::printf("n0(%f), n1(%f), n2(%f)\n", item.vetDimNor[0],
                        //                                        item.vetDimNor[1],
                        //                                        item.vetDimNor[2]);
                        //std::printf("d0(%f), d1(%f), d2(%f); vol(%f)\n\n", d0, d1, d2,
                        //                                                 d0*d1*d2);

                        vol += d0*d1*d2;
                    }


                    doBreakTestRoute = true;
                    return false;
                }
            }
        }
    }

    return true;

}
