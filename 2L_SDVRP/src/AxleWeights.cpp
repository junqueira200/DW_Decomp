#include "AxleWeights.h"
#include "Instancia.h"
#include "Constantes.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin& bin)
{
    double fK  = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;

    for(int i=0; i < bin.numItens; ++i)
    {
        double f = instanciaG.vetItens[bin.vetItens[i]].weightForce;
        fTA += f;
        double r = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                   instanciaG.vetItens[bin.vetItens[i]].getDimRotacionada(0, bin.vetRotacao[i])/2.0;

        fK += f*r;
    }


    fK *= distanceKingpinTrailerAxle;
    fFA = (1.0/(double)wheelBase) * fK * distanceKingpinRearAxle;
    fRA = fK - fFA;
    fTA = fTA - fK;

    if(fFA > maxMassFrontAxle*Gravity || fRA > maxMassRearAxle*Gravity || fTA > maxMassTrailerAxle*Gravity)
        return false;

    return true;
}

double AxleWeightsNS::SemiTrailer::computeMaxFK(SolucaoNS::Bin& bin)
{
    double maxFK = 0.0;
    for(int i=0; i < bin.numItens; ++i)
    {
        Item& item = instanciaG.vetItens[bin.vetItens[i]];

        double f = item.weightForce;
        double r =  (double)distanceCargoSpaceTrailerAxle - std::min(item.vetDim[0], item.vetDim[1]);
        maxFK += f*r;
    }

    return distanceKingpinTrailerAxle*maxFK;
}
