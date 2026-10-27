/*  *****************************************************************
 *  *****************************************************************
 *  File:    LowerBound.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    24/02/25
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "LowerBound.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>

using namespace LabelingAlgorithmNS;
using namespace boost;

// Define edge property: weight (double)
typedef property<edge_weight_t, double> EdgeWeightProperty;

// Define a directed graph with edge weights
typedef adjacency_list<vecS,
                       vecS,
                       directedS,
                       no_property, // vertex properties
                       EdgeWeightProperty>
                                                 graph_t;
typedef graph_traits<graph_t>::vertex_descriptor Vertex;
typedef graph_traits<graph_t>::edge_descriptor   edge_descriptor;
typedef std::pair<int, int>                      Edge;

/** ******************************************************************************************
 *  ******************************************************************************************
 *
 *  @param vetMatResCost
 *  @param vetDist          Returns the LB distance from the i-th customer to the last one
 *
 *  ******************************************************************************************
 *  ******************************************************************************************
 */


bool LowerBoundNS::getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost &vetMatResCost,
                                     int                                       dest,
                                     const LabelingAlgorithmNS::MatBoundRes   &matBoundRes,
                                     EigenMatrixRowD 						  &matDp)
{
    static int                 numNodes = matDp.rows();
    static int 				   lastNode = numNodes - 1;
    static bool allocateMat = false;
    static int capacity_Q = -InfInt;
    static Vector<int> demands(numNodes);

    demands[0] = 0;

    assert(matDp.cols() >= capacity_Q + 1);

    if(!allocateMat)
    {
        for(int i=1; i < (numNodes-1); ++i)
            capacity_Q = std::max(capacity_Q, (int)matBoundRes(i, 1).upperBound);

        for(int v=1; v < (numNodes-1); ++v)
            demands[v] = (int)vetMatResCost(0, v, 1);

        allocateMat = true;
    }


    matDp.setConstant(InfDouble);
    //vetDist.setConstant(-InfDouble);

    for(int i=0; i <= capacity_Q; ++i)
        matDp(dest, i) = 0.0;

    for(int r=1; r <= capacity_Q; ++r)
    {
        for(int u=0; u < lastNode; ++u)
        {
            if(u == dest)
                continue;

            matDp(u, r) = std::min(matDp(u, r), vetMatResCost(u, dest, 0));

            for(int v=1; v < lastNode; ++v)
            {

                if (v == u || v == dest)
                    continue;

                if(demands[v] <= r)
                {
                    int rem = r - demands[v];

                    double arcCost = vetMatResCost(u, v, 0);
                    double remCost = matDp(v, rem);

                    if (arcCost < InfDouble / 2.0 && remCost < InfDouble / 2.0)
                        matDp(u, r) = std::min(matDp(u, r), arcCost+remCost);
                }
            }

            if(r > 0)
                matDp(u, r) = std::min(matDp(u, r), matDp(u, r-1));

        }
    }

    return true;
}

void LowerBoundNS::copyDistMat(const LabelingAlgorithmNS::Vet3D_ResCost &vetMatResCost,
                               Eigen::MatrixXd                          &distMat)
{

    const int sizeDimI = vetMatResCost.getNumDimI();
    const int sizeDimJ = vetMatResCost.getNumDimJ();

    FloatType arcMin = InfFloatType;
    FloatType temp;

    for(int j = 0; j < sizeDimJ; ++j)
    {
        for(int i = 0; i < sizeDimI; ++i)
        {
            temp = vetMatResCost(i, j, 0);
            distMat(i, j) = temp;

            if(temp != InfFloatType && temp < arcMin)
                arcMin = temp;
        }
    }

    // std::cout<<distMat<<"\n";

    // if(arcMin >= 0.0)
    //     return;

    arcMin = -arcMin;

    for(int j = 0; j < sizeDimJ; ++j)
    {
        for(int i = 0; i < sizeDimI; ++i)
        {
            if(i == j)
                continue;

            distMat(i, j) += arcMin;
        }
    }

    // std::cout<<"\n\n"<<distMat<<"\n";
}
