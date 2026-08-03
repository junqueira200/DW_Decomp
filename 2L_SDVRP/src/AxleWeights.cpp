#include "AxleWeights.h"
#include "Constants.h"
#include "Instancia.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;
/*
 * fFA: 6343672.8 fRA: 7727029.6(Limit: 11845575.0); fTA: 23665424.6(Limit: 23661720.0)
 *      8240400.0
 *
 * fRA:   7727029.6
 *Limit: 11845575.0
 *
 * fTA:   23662139.0(Limit: 23661720.0)
 * Limit: 23661720.0
 */

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin &bin, bool print) const
{
    GravityCm = GravityCmConst;

    double fK = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;

    double sumF = 0.0;
    double sumM = 0.0;

    for(int i = 0; i < bin.numItens; ++i)
    {

        double f = instanciaG.vetItens[bin.vetItemId[i]].weight * GravityCm;
        sumF += f;
        double r = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                   instanciaG.vetItens[bin.vetItemId[i]].getDimRotacionada(
                       0, bin.vetRotacao[i]) /
                       2.0;

        sumM += f * r;
    }

    fK = (1.0 / (double)distanceKingpinTrailerAxle) *
         (sumM + (double)massTrailer * GravityCm * distanceMassTrailerTrailerAxle);
    fFA = (1.0 / (double)wheelBase) *
          (fK * (double)distanceKingpinRearAxle +
           (double)massTractor * GravityCm * distanceMassTractorRearAxle);
    fRA = fK + (double)massTractor * GravityCm - fFA;
    fTA = sumF + (double)massTrailer * GravityCm - fK;

    // 1.01
    if(fFA > (double)maxMassFrontAxle* GravityCm ||
       fRA > (double)maxMassRearAxle * GravityCm ||
       fTA > (double)maxMassTrailerAxle * GravityCm) // ||
    {
        if(print)
        {
            std::printf("GravityMM: %.1f\n", GravityCm);
            std::printf("FROM AxleWeights: fK: %.1f; fFA: %.1f(Limit: %.1f); fRA: %.1f(Limit: %.1f); fTA: %.1f(Limit: %.1f)\n",
                        fK, fFA, (double)maxMassFrontAxle * GravityCm,
                        fRA, (double)maxMassRearAxle * GravityCm,
                        fTA, (double)maxMassTrailerAxle * GravityCm);
        }
        return false;
    }

    return true;
}

double AxleWeightsNS::SemiTrailer::computeMaxFK(SolucaoNS::Bin &bin) const
{
    double maxFK = massTrailer * 10 * distanceMassTrailerTrailerAxle;
    for(int i = 0; i < bin.numItens; ++i)
    {
        Item &item = instanciaG.vetItens[bin.vetItens[i]];

        double f = item.weight * 10;
        double r = (double)
            distanceCargoSpaceTrailerAxle; // - std::min(item.vetDim[0], item.vetDim[1]);
        maxFK += f * r;
    }

    return (1.0 / (double)distanceKingpinTrailerAxle) * maxFK;
}
