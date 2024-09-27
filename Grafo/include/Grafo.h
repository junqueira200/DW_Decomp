//
// Created by igor on 25/09/24.
//

#ifndef MNFP_GRAFO_H
#define MNFP_GRAFO_H

#include "GAux.h"
#include <iostream>
#include <cinttypes>
#include <map>
#include <Eigen/Eigen>
#include <ranges>
#include <exception>
#include <limits>

namespace GraphNS
{

    class ArcDontExist : public std::exception
    {
    public:

        int64_t i,j;
        std::string msg;
        ArcDontExist(int64_t ii, int64_t jj):i(ii), j(jj)
        {
            msg = "arc ("+ std::to_string(i) + ", " + std::to_string(j) + ") don't exist!";
        }

        virtual const char* what() const throw()
        {
            return msg.c_str();
        }
    };

    class OutOfRange : public std::exception
    {
    public:

        std::string msg;
        OutOfRange()
        {
            msg = "Out of range";
        }


        virtual const char* what() const throw()
        {
            return msg.c_str();
        }
    };

    struct Arc
    {
        int64_t src    = -1;
        int64_t dest   = -1;
        std::reference_wrapper<double> weight;
    };

    template<typename T>
    requires std::is_arithmetic_v<T>
    struct Graph
    {

    //private:

        int64_t numVertices    = 0;
        //int64_t numArcs      = 0;
        std::vector<std::map<int64_t , T>> arcs;       // Arcs of vertice i
        std::vector<Arc> vetArcs;                      // All arcs of G

    //public:

        explicit Graph(int64_t numV):numVertices(numV)
        {
            arcs = std::vector<std::map<int64_t ,T>>(numVertices);

        }

        Graph()=default;

        // Remove all arcs from the graph
        void reset(int64_t numV)
        {
            numVertices = numV;
            arcs = std::vector<std::map<int64_t ,T>>(numVertices);
        }

        void loadVetArcs()
        {
            vetArcs = std::vector<Arc>();

            int pos = 0;
            for(int64_t i=0; i < numVertices; ++i)
            {
                for(auto &it: arcs[i])
                {
                    vetArcs.push_back(Arc(i, it.first, it.second));
                    std::cout<<i<<" "<<it.first<<" "<<(it.second)<<"\n";
                }
            }

        }

        std::string printVetArcs()
        {
            std::cout<<"printVetArcs\n";

            std::string str;
            for(auto &it:vetArcs)
            {
                str += "(" + std::to_string(it.src) + ", " + std::to_string(it.dest) + "): " +
                        std::to_string(it.weight.get());
            }

            return str;
        }

        void addVertice()
        {
            numVertices += 1;
            arcs.push_back(std::map<int64_t ,T>());
        }

        void addArc(int64_t i, int64_t j, T val)
        {
            if(i >= numVertices || i < 0 || j >= numVertices || j < 0)
                throw OutOfRange();

            (arcs[i])[j] = val;
        }


        // Get the arcs that starts with i, of the form (i,_); Returns a pair of begin and end
        auto getArcs(int64_t i)
        {

            if(i >= numVertices || i < 0)
                throw OutOfRange();

            return std::make_pair(arcs[i].begin(), arcs[i].end());
        }

        bool arcExist(int64_t i, int64_t j)
        {

            if(i >= numVertices || i < 0 || j >= numVertices || j < 0)
                throw OutOfRange();

            return arcs[i].count(j) == 1;
        }

        // Get the value storaged in (i,j)
        T getArc(int64_t i, int64_t j)
        {
            if(i >= numVertices || i < 0 || j >= numVertices || j < 0)
                throw OutOfRange();

            if(!arcExist(i, j))
                throw ArcDontExist(i,j);

            return (arcs[i])[j];
        }


        // Get the arcs that starts with i, of the form (i,_); Returns range
        auto getArcsRange(int64_t i)
        {

            if(i >= numVertices || i < 0)
                throw OutOfRange();

            return std::ranges::subrange(arcs[i].begin(), arcs[i].end());
        }

        // Get the number of arcs that starts with i, of the form (i,_)
        int64_t getNumArcsI(int64_t i)
        {
            if(i >= numVertices || i < 0)
                throw OutOfRange();

            return arcs[i].size();
        }

        int64_t calculateNumArcs()
        {
            int64_t numArcs = 0;
            for(auto &it:arcs)
                numArcs += it.size();

            return numArcs;
        }

        int64_t getNumVertices(){return numVertices;}
    };

    template<typename T>
    void bellmanFord(Graph<T> &graph,
                     int64_t src,
                     Eigen::VectorXd &vetDist,
                     Eigen::VectorX<int64_t> &vetPredecessor,
                     int64_t &verticeNegCycle)
    {
        std::cout<<"**************************************\n";
        std::cout<<"*************BELLMAN FORD*************\n\n";


        T Inf = std::numeric_limits<T>::infinity();

        if(!std::numeric_limits<T>::has_infinity)
            Inf = std::numeric_limits<T>::max();


        vetDist.setConstant(Inf);
        vetPredecessor.setConstant(-1);
        vetDist[src] = T(0);


        // Relax all edges |V| - 1 times.
        for(int i=0; i < graph.numVertices-1; ++i)
        {
            for(int j=0; j < graph.numVertices; ++j)
            {
                for(auto &it: graph.getArcsRange(j))
                {
                    int64_t u = j;
                    int64_t v = it.first;
                    double weight = it.second;
                    double sum = vetDist[u] + weight;

                    if(vetDist[u] != Inf && sum < vetDist[v])
                    {
                        vetDist[v] = sum;
                        vetPredecessor[v] = u;
                    }
                }
            }
        }

        bool negCycle = false;
        verticeNegCycle = -1;

        // Check for negative cycles
        for(int j=0; j < graph.numVertices; ++j)
        {
            for(auto &it: graph.getArcsRange(j))
            {
                int64_t u = j;
                int64_t v = it.first;
                double weight = it.second;
                double sum = vetDist[u] + weight;

                if(vetDist[u] != Inf && sum < vetDist[v])
                {
                    negCycle = true;
                    verticeNegCycle = u;
                    break;

                }
            }

            if(negCycle)
                break;
        }



        std::cout<<"negCycle: "<<negCycle<<"\n";
        std::cout<<"id: "<<verticeNegCycle<<"\n";
        std::cout<<vetDist.transpose()<<"\n";
        std::cout<<vetPredecessor.transpose()<<"\n";
    }

    template<typename T>
    void copyGraph(const Graph<T> &src, Graph<T> &dest)
    {
        dest.numVertices = src.numVertices;
        dest.reset(src.numVertices);
        dest.arcs = src.arcs;
    }

    template<typename T>
    std::string printGraph(const Graph<T> &graph)
    {
        std::string str;
        auto p = std::make_pair(0,0);

        for(int64_t i=0; i < graph.numVertices; ++i)
        {
            for(const auto &it:graph.arcs[i])
                str += "("+std::to_string(i)+", " + std::to_string(it.first)+"): "+std::to_string(it.second)+"\n";
        }

        return str;
    }
}

#endif //MNFP_GRAFO_H
