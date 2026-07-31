/* ****************************************
 * ****************************************
 *  Data:    11/12/24
 *  File:    IG.cpp
 *  Author:  Igor de Andrade Junqueira
 *  Project: 2L-SDVRP
 * ****************************************
 * ****************************************/

#include "Ig.h"
#include "BuscaLocal.h"
#include "Construtivo.h"
#include "InputOutput.h"
#include "Instancia.h"
#include "rand.h"
#include "c_api.h"
#include "SetPartition.h"

using namespace SolucaoNS;
using namespace ConstrutivoNS;
using namespace InstanceNS;
using namespace ParseInputNS;
using namespace RandNs;
using namespace BuscaLocalNS;
using namespace SetPartitionNS;

bool IgNs::metaheuristicaIg(SolucaoNS::Solucao &best)
{

    std::string strError;
    Solucao     sol(instanciaG);
    int         ultimaA = 0;

    for(int i = 0; i < 500; ++i)
    {
        sol.reset();
        if(construtivoVrp(sol, input.alphaVrp, input.aphaBin))
        {
            if(!sol.verificaSol(strError))
            {
                std::cout << "ERROR\n" << strError << "\n";
                throw "ERROR";
            }

            best.copiaSolucao(sol);
            break;
        }
    }

    if(!best.verificaSol(strError))
    {
        std::cout << strError << "\n\n";
        return false;
    }

    std::printf("Dist; %.1f\n\n", best.distTotal);
    static int numRm = 0.2*instanciaG.numVeiculos+1;
    std::printf("Number of veich to remove: %d\n\n", numRm);

    for(int i = 0; i < input.numItIG; ++i)
    {
        if(i%100 == 0 && i > 0)
            std::printf("IT: %d; Number of Routes: %ld\n", i, routeData.routeSetFeasible.size());

        for(int k = 0; k < numRm; ++k)
        {
            int r = getRandInt(0, instanciaG.numVeiculos - 1);
            while(sol.vetRota[r].numPos == 2)
                r = (r + 1) % instanciaG.numVeiculos;

            rmRota(sol, r);
        }

        if(!construtivoVrp(sol, input.alphaVrp, input.aphaBin))
            sol.copiaSolucao(best);
        else
        {
            // std::cout<<"Viavel!\n";
            rvnd(sol);

            if(doubleLess(sol.distTotal, best.distTotal))
            {
                best.copiaSolucao(sol);
                ultimaA = i;

                std::printf("Dist: %.1f; i: %d\n", best.distTotal, i);
            }

            else if(sol.distTotal > best.distTotal)
            {
                double gap = ((sol.distTotal - best.distTotal) / best.distTotal);
                if(gap > input.gapIgReset)
                    sol.copiaSolucao(best);
            }

            // std::cout<<"\t"<<sol.distTotal<<"\n\n";
        }

    }

    for(int i=0; i < best.vetRota.size(); ++i)
    {
        Rota& rota = best.vetRota[i];
        if(rota.numPos <= 2)
            continue;

        if(testRoute(&rota.vetRota[0], rota.numPos, 1) == 0)
        {
            std::printf("Error, packing for route: %s not found\n",
                        rota.printRota(false).c_str());
            PRINT_THROW();
        }
    }

    Solucao setPartitionSol(instanciaG);
    SetPartitionNS::setPatition(best, setPartitionSol);

    // std::cout<<"MELHOR SOL: "<<best.distTotal<<"\n\n";
    return sol.verificaSol(strError);
}
