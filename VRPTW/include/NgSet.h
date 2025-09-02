#ifndef NGSET_H
#define NGSET_H

#include "Aux.h"
#include "LabelingConstants.h"
#include <bitset>


typedef std::vector<std::bitset<LabelingAlgorithmNS::NumMaxCust>> VetBitSet;

namespace LabelingAlgorithmNS
{

    class NgSet
    {
    public:

        //Eigen::VectorX<std::bitset<NumMaxCust>>
        VetBitSet vetNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);

        inline __attribute__((always_inline))
        bool contain(int i, int j) const
        {
            return false;
            //if(!active)
            //    return true;
            return (bool)vetNgSet[i][j];
        };

        void setNgSets(const EigenMatrixRowD &matDist);
    };

    struct Node
    {
        int       i;
        FloatType dist;

        bool operator<(const Node &node) const
        {
            return dist < node.dist;
        }
    };


}

#endif // NGSET_H
