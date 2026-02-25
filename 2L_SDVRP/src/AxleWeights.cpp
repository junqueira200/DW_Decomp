#include "AxleWeights.h"
#include "Instancia.h"
#include "Constantes.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin& bin) const
{
    double fK  = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;
    //std::printf("numItems: %d\n", bin.numItens);
    for(int i=0; i < bin.numItens; ++i)
    {
        double f = instanciaG.vetItens[bin.vetItemId[i]].weightForce;
        fTA += f;
        double r = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                   instanciaG.vetItens[bin.vetItemId[i]].getDimRotacionada(0, bin.vetRotacao[i])/2.0;
        //r = std::abs(r);

        //std::printf("F: %.1f; px: %.1f; R: %.1f\n", f, bin.vetPosItem[i].vetDim[0], r);

        fK += f*r;
    }


    fK *= 1.0/distanceKingpinTrailerAxle;
    fFA = (1.0/(double)wheelBase) * fK * distanceKingpinRearAxle;
    fRA = fK - fFA;
    fTA = fTA - fK;

    //std::printf("fK: %.1f; fFA: %.1f; fRA: %.1f; FTA: %.1f\n", fK, fFA, fRA, fTA);

    if(fFA > maxMassFrontAxle*Gravity || fRA > maxMassRearAxle*Gravity || fTA > maxMassTrailerAxle*Gravity ||
       fK < 0.0 || fFA < 0.0 || fRA < 0.0 || fTA < 0.0)
    {
        std::printf("\tINFEASIBLE\n");
        return false;  // physically infeasible loading
    }

    std::printf("\tFEASIBLE\n");

    return true;
}

double AxleWeightsNS::SemiTrailer::computeMaxFK(SolucaoNS::Bin& bin) const
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
