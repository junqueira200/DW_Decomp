// #ifndef AXLEWEIGHTS_H
// #define AXLEWEIGHTS_H
#pragma once

#include "Instancia.h"
#include "Solucao.h"

namespace AxleWeightsNS
{

class AxleData
{
  public:
    AxleData() {};
    virtual ~AxleData() {};
    virtual bool checkAxleWeights(SolucaoNS::Bin &bin, bool print=false) const = 0;
};

class SemiTrailer : public AxleData
{
  public:
    SemiTrailer() {}
    ~SemiTrailer() {}
    bool   checkAxleWeights(SolucaoNS::Bin &bin, bool print=false) const override;
    double computeMaxFK(SolucaoNS::Bin &bin) const;

    int wheelBase = 379; // WB in cm
    int maxMassFrontAxle = 8400;
    int maxMassRearAxle = 12075;
    int maxMassTrailerAxle = 24120;
    int distanceKingpinRearAxle =
        56; // The distance between the kingpin and the rear axle (l_K|RA) in cm
    int distanceKingpinTrailerAxle =
        816; // The distance between the kingpin and the resultant trailer axle (l_K|TA)cm
    int distanceCargoSpaceTrailerAxle =
        955; // The distance between the cargo area to the resultant trailer axle (l_TA)cm
    int distanceMassTractorRearAxle = 257;    // ? in cm
    int distanceMassTrailerTrailerAxle = 163; // ? in cm
    int massTractor = 8197;
    int massTrailer = 6472;
};

inline const SemiTrailer semiTrailer;

} // namespace AxleWeightsNS

// #endif // AXLEWEIGHTS_H
