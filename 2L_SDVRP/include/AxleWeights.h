// #ifndef AXLEWEIGHTS_H
// #define AXLEWEIGHTS_H
#pragma once

#include "Instancia.h"
#include "Solucao.h"
#include "InputOutput.h"

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
    SemiTrailer()
    {
        if(ParseInputNS::input.useCm)
            return;

        wheelBase                      *= 10;
        distanceKingpinRearAxle        *= 10;
        distanceKingpinTrailerAxle     *= 10;
        distanceCargoSpaceTrailerAxle  *= 10;
        distanceMassTractorRearAxle    *= 10;
        distanceMassTrailerTrailerAxle *= 10;



    }
    ~SemiTrailer() {}
    bool   checkAxleWeights(SolucaoNS::Bin &bin,
                            bool            print    = false,
                            double*         ptrFk    = nullptr,
                            double*         ptrFFa   = nullptr,
                            double*         ptrFRa   = nullptr,
                            double*			ptrFTa   = nullptr) const override;

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
