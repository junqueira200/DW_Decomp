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
    virtual bool checkAxleWeights(SolucaoNS::Bin &bin,
                                  bool            print = false,
                                  double*         fk    = nullptr,
                                  double*         fFa   = nullptr,
                                  double*         fRa   = nullptr,
                                  double*	      fTa   = nullptr) const = 0;
};

class SemiTrailer : public AxleData
{
  public:
    SemiTrailer() {}
    ~SemiTrailer() {}
    bool   checkAxleWeights(SolucaoNS::Bin &bin,
                            bool            print    = false,
                            double*         ptrFk    = nullptr,
                            double*         ptrFFa   = nullptr,
                            double*         ptrFRa   = nullptr,
                            double*			ptrFTa   = nullptr) const override;

    double computeMaxFK(SolucaoNS::Bin &bin) const;

    int wheelBase = 3790; // WB in mm
    int maxMassFrontAxle = 8400;
    int maxMassRearAxle = 12075;
    int maxMassTrailerAxle = 24120;
    int distanceKingpinRearAxle =
        560; // The distance between the kingpin and the rear axle (l_K|RA) in mm
    int distanceKingpinTrailerAxle =
        8160; // The distance between the kingpin and the resultant trailer axle (l_K|TA)mm
    int distanceCargoSpaceTrailerAxle =
        9550; // The distance between the cargo area to the resultant trailer axle (l_TA)mm
    int distanceMassTractorRearAxle = 2570;    // ? in mm
    int distanceMassTrailerTrailerAxle = 1630; // ? in mm
    int massTractor = 8197;
    int massTrailer = 6472;
};

inline const SemiTrailer semiTrailer;

} // namespace AxleWeightsNS

// #endif // AXLEWEIGHTS_H
