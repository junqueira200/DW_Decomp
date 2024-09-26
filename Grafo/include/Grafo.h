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



    template<typename T>
    requires std::is_arithmetic_v<T>
    struct Graph
    {

    //private:

        const int64_t numVertices;
        int64_t numArcs      = 0;
        Eigen::VectorX<std::map<int, T>> arcs;

    //public:

        explicit Graph(int64_t numV):numVertices(numV)
        {
            arcs = Eigen::VectorX<std::map<int,T>>(numVertices);

        }

        void addArc(int64_t i, int64_t j, T val)
        {
            if(i >= numVertices || i < 0 || j >= numVertices || j < 0)
                throw OutOfRange();

            (arcs[i])[j] = val;
        }


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

        T getArc(int64_t i, int64_t j)
        {
            if(i >= numVertices || i < 0 || j >= numVertices || j < 0)
                throw OutOfRange();

            if(!arcExist(i, j))
                throw ArcDontExist(i,j);

            return (arcs[i])[j];
        }

        auto getArcsRange(int64_t i)
        {

            if(i >= numVertices || i < 0)
                throw OutOfRange();

            return std::ranges::subrange(arcs[i].begin(), arcs[i].end());
        }

        int64_t getNumArcs(int64_t i)
        {
            if(i >= numVertices || i < 0)
                throw OutOfRange();

            return arcs[i].size();
        }
    };

    template<typename T>
    void bellmanFord(Graph<T> &graph, int src)
    {


        T Inf = std::numeric_limits<T>::infinity();

        if(!std::numeric_limits<T>::has_infinity)
            Inf = std::numeric_limits<T>::max();

        int64_t numArcs = graph.getNumArcs(src);
        Eigen::VectorX<T> vetDist(graph.numVertices);
        vetDist.setConstant(Inf);
        vetDist[src] = T(0);

        for(int i=0; i < graph.numVertices; ++i)
        {

            for(auto &it: graph.getArcsRange(i))
            {
                int j = it.first;
                T distJ = it.second;
                T sum = vetDist[i]+distJ;

                if(vetDist[i] != Inf && sum < vetDist[j])
                    vetDist[j] = sum;
            }

        }

        std::cout<<vetDist.transpose()<<"\n";
    }


}

#endif //MNFP_GRAFO_H
