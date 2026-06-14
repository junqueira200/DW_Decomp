#include "AxleWeights.h"
#include "Constants.h"
#include "Instancia.h"

using namespace AxleWeightsNS;
using namespace InstanceNS;

bool AxleWeightsNS::SemiTrailer::checkAxleWeights(SolucaoNS::Bin &bin, bool print) const
{
    GravityMM = GravityMM_const;

    double fK = 0.0;
    double fFA = 0.0;
    double fRA = 0.0;
    double fTA = 0.0;

    double sumF = 0.0;
    double sumM = 0.0;

    for(int i = 0; i < bin.numItens; ++i)
    {

        double f = instanciaG.vetItens[bin.vetItemId[i]].weight * GravityMM;
        sumF += f;
        double r = (double)distanceCargoSpaceTrailerAxle - bin.vetPosItem[i].vetDim[0] -
                   instanciaG.vetItens[bin.vetItemId[i]].getDimRotacionada(
                       0, bin.vetRotacao[i]) /
                       2.0;

        sumM += f * r;
    }

    fK = (1.0 / (double)distanceKingpinTrailerAxle) *
         (sumM + (double)massTrailer * GravityMM * distanceMassTrailerTrailerAxle);
    fFA = (1.0 / (double)wheelBase) *
          (fK * (double)distanceKingpinRearAxle +
           (double)massTractor * GravityMM * distanceMassTractorRearAxle);
    fRA = fK + (double)massTractor * GravityMM - fFA;
    fTA = sumF + (double)massTrailer * GravityMM - fK;

    if(fFA > (double)maxMassFrontAxle*1.01* GravityMM ||
       fRA > (double)maxMassRearAxle*1.01 * GravityMM ||
       fTA > (double)maxMassTrailerAxle*1.01 * GravityMM) // ||
    {
        if(print)
        {
            std::printf("GravityMM: %d\n", GravityMM);
            std::printf("FROM AxleWeights: fK: %.1f; fFA: %.1f(Limit: %.1f); fRA: %.1f(Limit: %.1f); fTA: %.1f(Limit: %.1f)\n",
                        fK, fFA, (double)maxMassFrontAxle * GravityMM,
                        fRA, (double)maxMassRearAxle * GravityMM,
                        fTA, (double)maxMassTrailerAxle * GravityMM);
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
