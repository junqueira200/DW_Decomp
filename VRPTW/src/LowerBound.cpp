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
    std::cout<<"ini\n";

    Eigen::MatrixXd distMat(vetDist.size(), vetDist.size());
    copyDistMat(vetMatResCost, distMat);
    static int numNodes = vetDist.size();
    static std::vector<Vertex> predecessor(numNodes+1);

    //std::cout<<distMat<<"\n";

    // Create the graph
    graph_t g(numNodes+1);

    for(int i=0; i < numNodes; ++i)
    {
        for(int j=0; j < numNodes; ++j)
        {
            if(i == j)
                continue;
            //std::cout<<i<<", "<<j<<"\n";
            if(distMat(i, j) < -1E-5)
            {
                std::cout<<"dist("<<i<<", "<<j<<") : "<<distMat(i, j)<<"\n\n";
                std::cout<<distMat<<"\n";
                PRINT_EXIT();
            }
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

    for(double& val:vetDist)
    {
        if(val > 0.0)
            val = 0.0;
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
