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
        AxleData(){};
        virtual ~AxleData(){};
        virtual bool checkAxleWeights(SolucaoNS::Bin& bin)const=0;
    };

    class SemiTrailer : public AxleData
    {
    public:

        SemiTrailer(){}
        ~SemiTrailer(){}
        bool checkAxleWeights(SolucaoNS::Bin& bin) const override;
        double computeMaxFK(SolucaoNS::Bin& bin) const;

        int wheelBase                      = 3790;        // WB
        int maxMassFrontAxle               = 8400;
        int maxMassRearAxle                = 12075;
        int maxMassTrailerAxle             = 24120;
        int distanceKingpinRearAxle        = 560;         // The distance between the kingpin and the rear axle (l_K|RA)
        int distanceKingpinTrailerAxle     = 8160;		// The distance between the kingpin and the resultant trailer axle (l_K|TA)
        int distanceCargoSpaceTrailerAxle  = 9550;		// The distance between the cargo area to the resultant trailer axle (l_TA)
        int distanceMassTractorRearAxle    = 2568;        // ?
        int distanceMassTrailerTrailerAxle = 1632;        // ?
        int massTractor				       = 8197;
        int massTrailer				       = 6472;

    };

    inline const SemiTrailer semiTrailer;

}

//#endif // AXLEWEIGHTS_H
