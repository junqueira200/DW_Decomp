//
// Created by igor on 20/11/24.
//

#ifndef DW_LABELINGALGORITHM_H
#define DW_LABELINGALGORITHM_H

#include <Eigen/Eigen>
#include "Grafo.h"
#include <bitset>

namespace LabelingAlgorithmNS
{

    constexpr int NumMaxResources = 1;
    constexpr int NumMaxRoute     = 300;
    constexpr int NumMaxCust      = 150;
    constexpr int NgSetSize       = 5;

    struct Bound
    {
        double lowerBound = -std::numeric_limits<double>::infinity();
        double upperBound = std::numeric_limits<double>::infinity();
    };

    class NgSet
    {
    public:

        Eigen::MatrixXi matNgSet;
        int numCust   = NumMaxCust;
        int ngSetSize = NgSetSize;
        bool active   = true;

        NgSet();
        NgSet(int numCust, int ngSetSize);
        bool contain(int i, int j) const;
        void setNgSets(const Eigen::MatrixXd &matDist);
    };

    struct Node
    {
        int i;
        double dist;

        bool operator<(const Node &node) const
        {
            return dist < node.dist;
        }
    };

    class Label
    {
    public:

        //int numResources = 1;
        Eigen::Array<double, 1, NumMaxResources> vetResources;
        Eigen::VectorXi vetRoute;
        int tamRoute = 0;
        std::bitset<NumMaxCust> bitSet;

        Label() = default;

    };


    class Bucket
    {
        // Bound: [lower;upper)
        Eigen::Array<Bound, 1, NumMaxResources> vetBound;
        Eigen::VectorX<Label*> vetPtrLabel;
        int sizeVetPtrLabel;
    };

    class MatBucket
    {

    };

    // Fist access the resource, and then the cost of an arc (i,j)
    typedef Eigen::Array<Eigen::Array<double, NumMaxCust, NumMaxCust>, 1, NumMaxResources> VetMatResCost;

    // Fist access the resource, and then the bound of a customer
    typedef Eigen::Array<Eigen::Array<Bound, 1, NumMaxCust>, 1, NumMaxResources> VetVetResBound;

    void forwardLabelingAlgorithm(const int numRes,
                                  const int numCust,
                                  const VetMatResCost& vetMatResCost,
                                  const VetVetResBound& vetVetBound,
                                  int dest,
                                  const NgSet &ngSet);

    bool checkDominance(const Label& l0, const Label& l1, int numResources);

    bool extendLabel(const Label &label,
                     Label &newLabel,
                     const VetMatResCost& vetMatResCost,
                     const VetVetResBound& vetVetBound,
                     int custI,
                     int custJ,
                     const NgSet &ngSet,
                     int numResources);

}
#endif //DW_LABELINGALGORITHM_H
