#include "AxleWeights.h"
#include "Instancia.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin& bin)
{
    int numItens = bin.numItens;
    static VectorD vetR(instanciaG.numItens);
    vetR.setAll(distanceCargoSpaceTrailerAxle);

    for(int i=0; i < numItens; ++i)
    {
        vetR[i] = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                  instanciaG.vetItens[bin.vetItens[i]].getDimRotacionada(0, bin.vetRotacao[i])/2.0;
    }

    double fK  = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;
    double sumF = 0.0;

    for(int i=0; i < numItens; ++i)
        fK += instanciaG.vetItens[bin.vetItens[i]].weight * vetR[i];

    fFA = (1.0/(double)wheelBase) * fK * distanceKingpinRearAxle;

    fRA = fK - fFA;

    //if(fK <= )

    return false;
}
