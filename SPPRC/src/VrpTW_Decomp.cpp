//
// Created by igor on 15/11/24.
//
#include "VrpTW_Decomp.h"
// #include "rand.h"
#include "InstanceVRPTW.h"
#include <iostream>

using namespace LabelingAlgorithmNS;

double VrpTW_DecompNS::geraSolHeuristica(const InstanceVRPTW_NS::InstanceVRPTW &instVrpTw,
                                          Eigen::VectorXd                       &vetRedCost,
                                          Eigen::VectorXi                       &rota,
                                          int                                   &routeTam)
{
std::cout<<"ini geraSolHeuristica"<<"\n";

    static Eigen::VectorXi vetClieVisi(instVrpTw.numClientes);
    vetClieVisi.setZero();
    rota.setZero();

    int clienteI = 0;

    int next = 1;
    int carga = 0;
    double cost = 0.0;
    double dist = 0.0;
    routeTam = 1;

    std::cout<<"\tRota criada: ";
    std::cout<<clienteI<<" ";

    do
    {
        int clienteJ = -1;
        double menor = std::numeric_limits<double>::infinity();

        for(int i=0; i < instVrpTw.numClientes; ++i)
        {
            if(vetClieVisi[i] != 0 || i == clienteI)
                continue;

            if((carga+instVrpTw.vetClieDem[i]) > instVrpTw.capVeic)
                continue;

            const int index = getIndex(clienteI, i, instVrpTw.numClientes);
            if(vetRedCost[index] < menor)
            {
                clienteJ = i;
                menor = vetRedCost[index];
            }
        }

        if(clienteJ == -1)
        {
            clienteJ = 0;
            cost += vetRedCost[getIndex(clienteI, clienteJ, instVrpTw.numClientes)];
        }
        else
            cost += menor;
        std::cout<<clienteJ<<" ";

        rota[next] = clienteJ;

        vetClieVisi[clienteJ] = 1;
        next += 1;
        routeTam += 1;
        clienteI = clienteJ;
        carga += carga+instVrpTw.vetClieDem[clienteJ];

    }while(clienteI != 0);


    std::cout<<"; custo: "<<cost<<"\n";

    return cost;
}

