#include "AxleWeights.h"
#include "Instancia.h"
#include "Constantes.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin& bin) const
{
    //std::printf("\n");
    double fK  = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;

    double sumF = 0.0;
    double sumM = 0.0;
    //std::printf("numItems: %d\n", bin.numItens);
    for(int i=0; i < bin.numItens; ++i)
    {

        double f = instanciaG.vetItens[bin.vetItemId[i]].weightForce;
        sumF += f;
        double r = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                   instanciaG.vetItens[bin.vetItemId[i]].getDimRotacionada(0, bin.vetRotacao[i])/2.0;
        //r = std::abs(r);

        //std::printf("W: %0.f; F: %.2f; Px: %.1f; Dx/2: %.0f; r: %.2f; M: %.2f\n",  instanciaG.vetItens[bin.vetItemId[i]].weight, f, bin.vetPosItem[i].vetDim[0], instanciaG.vetItens[bin.vetItemId[i]].getDimRotacionada(0, bin.vetRotacao[i])/2.0, r, f*r);

        //std::printf("F: %.1f; px: %.1f; R: %.1f\n", f, bin.vetPosItem[i].vetDim[0], r);

        sumM += f*r;
    }

    //std::printf("\nsumF: %.2f; sumM: %.2f\n", sumF, sumM);

    fK  = (1.0/distanceKingpinTrailerAxle)*(sumM + massTrailer*Gravity*distanceMassTrailerTrailerAxle);
    fFA = (1.0/(double)wheelBase) * (fK * distanceKingpinRearAxle + massTractor*Gravity*distanceMassTractorRearAxle);
    fRA = fK + massTractor*Gravity - fFA;
    fTA = sumF + massTrailer*Gravity - fK;

    //std::printf("fK: %.1f; fFA: %.1f; fRA: %.1f; FTA: %.1f\n", fK, fFA, fRA, fTA);

    if(fFA > maxMassFrontAxle*Gravity || fRA > maxMassRearAxle*Gravity || fTA > maxMassTrailerAxle*Gravity ||
       fK < 0.0 || fFA < 0.0 || fRA < 0.0 || fTA < 0.0)
    {
        std::printf("\tINFEASIBLE\t");

        if(fFA > maxMassFrontAxle*Gravity)
            std::printf("fFa >; ");
        if(fK < 0.0)
            std::printf("fk < 0; ");
        if(fFA < 0.0)
            std::printf("fFA < 0; ");
        if(fRA > maxMassRearAxle*Gravity)
            std::printf("fRA >; ");
        if(fRA < 0.0)
            std::printf("fRA < 0; ");
        if(fTA > maxMassTrailerAxle*Gravity)
            std::printf("fTA >; ");
        if(fTA < 0.0)
            std::printf("fTA < 0; ");

        std::printf("\n");
        return false;
    }

    //std::printf("\tFEASIBLE: Fk: %.1f; fFA: %.1f \n", fK, fFA);

    return true;
}

double AxleWeightsNS::SemiTrailer::computeMaxFK(SolucaoNS::Bin& bin) const
{
    double maxFK = massTrailer*10*distanceMassTrailerTrailerAxle;
    for(int i=0; i < bin.numItens; ++i)
    {
        Item& item = instanciaG.vetItens[bin.vetItens[i]];

        double f = item.weight*10;
        double r =  (double)distanceCargoSpaceTrailerAxle;// - std::min(item.vetDim[0], item.vetDim[1]);
        maxFK += f*r;
    }

    return (1.0/distanceKingpinTrailerAxle)*maxFK;
}
