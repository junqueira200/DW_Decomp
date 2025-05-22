#ifndef NGSET_H
#define NGSET_H

#include "Aux.h"
#include "LabelingConstants.h"

namespace LabelingAlgorithmNS
{

    class NgSet
    {
    public:

        EigenMatrixRowI matNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);
        [[nodiscard]] bool contain(int i, int j) const;
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
