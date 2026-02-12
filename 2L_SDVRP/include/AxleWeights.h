//#ifndef AXLEWEIGHTS_H
//#define AXLEWEIGHTS_H
#pragma once

#include "Instancia.h"
#include "Solucao.h"

namespace AxleWeightsNS
{

    class AxleData
    {
    public:
        AxleData();
        virtual ~AxleData(){};
        virtual bool checkAxleWeights(SolucaoNS::Bin& bin)=0;
    };

    class SemiTrailer : public AxleData
    {
    public:

        SemiTrailer(){}
        ~SemiTrailer(){}
        bool checkAxleWeights(SolucaoNS::Bin& bin) override;

        int wheelBase                      = 36;        // WB
        int maxMassFrontAxle               = 10000;
        int maxMassRearAxle                = 11500;
        int maxMassTrailerAxle             = 24000;
        int distanceKingpinRearAxle        = 6;         // The distance between the kingpin and the rear axle (l_K|RA)
        int distanceKingpinTrailerAxle     = 76;		// The distance between the kingpin and the resultant trailer axle (l_K|TA)
        int distanceCargoSpaceTrailerAxle  = 92;		// The distance between the cargo area to the resultant trailer axle (l_TA)
        int distanceMassTractorRearAxle    = 25;        // ?
        int distanceMassTrailerTrailerAxle = 16;        // ?
        int massTractor				       = 7300;
        int massTrailer				       = 6750;

    };

}

//#endif // AXLEWEIGHTS_H
