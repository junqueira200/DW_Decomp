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

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>

using namespace LabelingAlgorithmNS;
using namespace boost;

// Define edge property: weight (double)
typedef property<edge_weight_t, double> EdgeWeightProperty;

// Define a directed graph with edge weights
typedef adjacency_list<vecS, vecS, directedS,
                       no_property,      // vertex properties
                       EdgeWeightProperty> graph_t;
typedef graph_traits < graph_t >::vertex_descriptor Vertex;
typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
typedef std::pair<int, int> Edge;


/** ******************************************************************************************
 *  ******************************************************************************************
 *
 *  @param vetMatResCost
 *  @param vetDist          Returns the LB distance from the i-th customer to the last one
 *
 *  ******************************************************************************************
 *  ******************************************************************************************
 */
bool LowerBoundNS::getDistLowerBound(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                                     Eigen::VectorXd&                          vetDist,
                                     int                                       dest,
                                     LabelingAlgorithmNS::LabelingData*        lDataPtr)

{
    //return false;
    //return false;
    //std::cout<<"ini\n";

    Eigen::MatrixXd distMat(vetDist.size(), vetDist.size());
    copyDistMat(vetMatResCost, distMat);
    static int numNodes = vetDist.size();
    static std::vector<Vertex> predecessor(numNodes);
    static std::vector<double> distance(numNodes, std::numeric_limits<double>::infinity());

    //std::cout<<distMat<<"\n";

    // Create the graph
    graph_t g(numNodes+1);

    /*
    // Reverse graph: edges j -> i
    for(int i=0; i < numNodes; ++i)
    {
        for(int j=0; j < numNodes; ++j)
        {
            if(i == j) continue;
            if(distMat(i, j) == std::numeric_limits<double>::infinity())
                continue;

            add_edge(i, j, EdgeWeightProperty(distMat(i,j)), g);
        }
    }


//    dijkstra_shortest_paths(g, dest, predecessor_map(make_iterator_property_map(predecessor.begin(), get(vertex_index, g)))
//            .distance_map(make_iterator_property_map(distance.begin(), get(vertex_index, g))));


    bool ok = bellman_ford_shortest_paths(g, numNodes, weight_map(get(edge_weight, g)).distance_map(make_iterator_property_map(
                                          distance.begin(), get(vertex_index, g))));


    if(!ok)
    {
        vetDist.setConstant(-MaxFloatType);
        return false;
    }


    // Copy result
    for(int i = 0; i < numNodes; ++i)
    {
        vetDist[i] = distance[i];
    }
    */

    for(int i=0; i < numNodes+1; ++i)
    {
        for(int j=0; j < numNodes+1; ++j)
        {
            if(i == j)
                continue;
            //std::cout<<i<<", "<<j<<"\n";
           add_edge(i, j, EdgeWeightProperty(distMat(i,j)), g);
        }
    }

    //std::cout<<"Create graph!\n";

    for(int i=0; i < numNodes; ++i)
    {
        //std::cout<<"Cust("<<i<<"\n";
        Vertex source = i;
        dijkstra_shortest_paths(g, source, predecessor_map(boost::make_iterator_property_map(
                                               predecessor.begin(), get(boost::vertex_index, g))));
        //std::cout<<"dijkstra\n";
        Vertex node = dest;
        VectorI vetRoute;
        vetRoute.reserve(5);
        vetRoute.push_back(dest);
        //std::cout<<"node: "<<dest<<"\n";
        while(node != source)
        {
            node = predecessor[node];
            //std::cout<<"node: "<<node<<"\n";
            vetRoute.push_back(node);
        }

        //std::cout<<"\n";

        std::reverse(vetRoute.begin(), vetRoute.end());
        //std::cout<<vetRoute<<"\n\n";
        double dist = 0.0;
\
        for(int i=0; i < (vetRoute.size()-1); ++i)
        {
            dist += vetMatResCost(vetRoute[i], vetRoute[i+1], 0);
        }

        vetDist[i] = dist;
    }



    std::cout<<"Dist: "<<vetDist.transpose()<<"\n";
    //PRINT_EXIT();

    return true;
}

void LowerBoundNS::copyDistMat(const LabelingAlgorithmNS::Vet3D_ResCost& vetMatResCost,
                               Eigen::MatrixXd& distMat)
{

    const int sizeDimI = vetMatResCost.getNumDimI();
    const int sizeDimJ = vetMatResCost.getNumDimJ();

    FloatType arcMin = InfFloatType;
    FloatType temp;

    for(int j=0; j < sizeDimJ; ++j)
    {
        for(int i=0; i < sizeDimI; ++i)
        {
            temp = vetMatResCost(i, j, 0);
            distMat(i, j) = temp;

            if(temp != InfFloatType && temp < arcMin)
                arcMin = temp;
        }
    }

    //std::cout<<distMat<<"\n";

    //if(arcMin >= 0.0)
    //    return;

    arcMin = -arcMin;

    for(int j=0; j < sizeDimJ; ++j)
    {
        for(int i=0; i < sizeDimI; ++i)
        {
            if(i == j)
                continue;

            distMat(i, j) += arcMin;
        }
    }

    //std::cout<<"\n\n"<<distMat<<"\n";

}
