/* ****************************************
 * ****************************************
 *  Data:    22/07/26
 *  File:    SetPartition.cpp
 *  Author:  Igor de Andrade Junqueira
 *  Project: 2L-SDVRP
 * ****************************************
 * ****************************************/

#include "SetPartition.h"
#include "gurobi_c++.h"
#include "Instancia.h"
#include "c_api.h"

using namespace SetPartitionNS;
using namespace SolucaoNS;
using namespace InstanceNS;


void SetPartitionNS::setPatition(Solucao &startSol, Solucao &outPutSol)
{
    std::printf("Num rotas: %ld\n", routeData.routeSetFeasible.size());
    std::printf("Start Sol dist: %.2f\n", startSol.distTotal);

    //startSol.hashVal = hashRoute(startSol);
    for(int i=0; i < startSol.vetRota.size(); ++i)
        startSol.vetRota[i].hashVal = hashRoute(startSol.vetRota[i]);

    try {


    GRBEnv env;
    GRBModel model(env);
    //model.set(GRB_DoubleParam_MIPGap, 0.0); // Set 5% gap
    model.set(GRB_IntParam_Threads, 1);
    GRBVar* varX = model.addVars(routeData.routeSetFeasible.size(),GRB_BINARY);
    Vector<GRBLinExpr> vetLinExpr;


    for(int i=0; i < (instanciaG.numClientes-1); ++i)
        vetLinExpr.push_back(GRBLinExpr());

    Vector<int8_t> vetSelectedRoutes(startSol.vetRota.size());
    vetSelectedRoutes.setAll((int8_t)0);

    int64_t routeId = 0;
    GRBLinExpr obj;
    for(const Rota& route:routeData.routeSetFeasible)
    {
        size_t hashVal = hashRoute(route);

        for(int i=1; i < (route.numPos-1); ++i)
        {
            int cust = route.vetRota[i];
            vetLinExpr[cust-1] += varX[routeId];

        }

        int64_t routeIdSol = -1;
        for(int i=0; i < startSol.vetRota.size(); ++i)
        {
            if(startSol.vetRota[i].hashVal != hashVal)
                continue;

            if(startSol.vetRota[i].numPos != route.numPos)
                continue;

            bool equal = true;
            for(int j=1; j < (route.numPos-1); ++j)
            {
                if(route.vetRota[j] != startSol.vetRota[i].vetRota[j])
                {
                    equal = false;
                    break;
                }
            }

            if(equal)
            {
                routeIdSol = i;
                break;
            }
        }

        if(routeIdSol >= 0)
        {
            if(vetSelectedRoutes[routeIdSol] == (int8_t)1)
            {
                std::printf("Error, route: %s, was already selected\n",
                            startSol.vetRota[routeIdSol].printRota(false).c_str());
                PRINT_THROW();
            }

            vetSelectedRoutes[routeIdSol] = (int8_t)1;
            varX[routeId].set(GRB_DoubleAttr_Start, 1.0);
        }
        else
            varX[routeId].set(GRB_DoubleAttr_Start, 0.0);

        obj += route.distTotal*varX[routeId];
        routeId += 1;
    }

    for(int i=0; i < (instanciaG.numClientes-1); ++i)
        model.addConstr(vetLinExpr[i], '=', 1, "Cust_"+std::to_string(i+1));


    model.setObjective(obj, GRB_MINIMIZE);
    model.update();
    model.optimize();

    delete []varX;

    } catch (GRBException& e)
    {
        std::cout<<e.getMessage()<<"\n"<<e.getErrorCode()<<"\n\n\n";
        PRINT_THROW();
    }
}
