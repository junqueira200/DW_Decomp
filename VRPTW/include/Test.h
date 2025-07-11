/*  ********************************************************************************************************************
 *  ********************************************************************************************************************
 *  File:    Test.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    06/06/25
 *
 *  ********************************************************************************************************************
 *  ********************************************************************************************************************/


#ifndef TEST_H
#define TEST_H

#include "safe_vector.h"
#include "Instancia.h"
#include <bits/stdc++.h>
#include "LabelingConstants.h"

namespace TestNS
{

    struct Route
    {
        VectorI vetRoute;
        double dist = 0.0;
        int demand  = 0;
        std::size_t valHash = 0;

        Route(int n){vetRoute.resize(n); vetRoute.setAll(0);}
        Route(){};
        void resize(int n){vetRoute.resize(n); vetRoute.setAll(0);}
        double getDistance(){return dist;}
        bool operator == (const Route& route) const;

    };

    struct Hash
    {
        std::size_t operator()(const Route& route) const{return route.valHash;}
    };

    typedef std::unordered_set<Route, Hash> RouteHash;

    void combinationUtil(int ind, int r, VectorI &data, Vector<VectorI> &result, VectorI &arr);
    Vector<VectorI> findCombination(VectorI &arr, int r);
    void enumerateRoutes(InstanciaNS::InstVRP_TW& instVrp, int numMax, RouteHash& routeHash);
    void computeDistance(Route& route);
    void computeHash(Route& route);
    int computeDemand(VectorI& route);
    FloatType computeReducedCost(const Route& route, const Eigen::RowVectorXd& vetRowPi);

}

#endif // TEST_H
