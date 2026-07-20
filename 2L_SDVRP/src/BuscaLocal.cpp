/* ****************************************
 * ****************************************
 *  Data:    17/12/24
 *  Arquivo: BuscaLocal.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#include "BuscaLocal.h"
#include "Construtivo.h"
#include "ConstrutivoBin.h"
#include "InputOutput.h"
#include "rand.h"
#include "c_api.h"

using namespace SolucaoNS;
using namespace RandNs;
using namespace ConstrutivoNS;
using namespace ConstrutivoBinNS;
using namespace ParseInputNS;

bool BuscaLocalNS::intraRotaShift(SolucaoNS::Solucao &sol, int rotaR, int pos0, int pos1)
{
    // std::cout<<"MvIntraRotaShifit pos0("<<pos0<<"); pos1("<<pos1<<");
    // rotaR("<<rotaR<<")\n";

    // Verifica se posicoes sao viaveis

    if((pos0 + 1) == pos1 || pos0 == pos1 || pos1 <= 0)
        return false;

    Rota &rota = sol.vetRota[rotaR];
    static Rota rotaTemp;

    // std::cout<<"\t"<<rota.printRota()<<"\n";

    if(pos1 >= (rota.numPos - 1))
        return false;

    if(pos0 >= (rota.numPos - 1))
        return false;

    const int clientePos1 = rota.vetRota[pos1];

    // Calcula nova distancia

    double novaDist = rota.distTotal;
    novaDist -= instanciaG.matDist(rota.vetRota[pos1 - 1], rota.vetRota[pos1]);
    novaDist -= instanciaG.matDist(rota.vetRota[pos1], rota.vetRota[pos1 + 1]);
    novaDist -= instanciaG.matDist(rota.vetRota[pos0], rota.vetRota[pos0 + 1]);
    novaDist += instanciaG.matDist(rota.vetRota[pos0], clientePos1);
    novaDist += instanciaG.matDist(clientePos1, rota.vetRota[pos0 + 1]);
    novaDist += instanciaG.matDist(rota.vetRota[pos1-1], rota.vetRota[pos1 + 1]);

    // std::cout<<"\tNova Dist: "<<novaDist<<"\n";

    if(doubleLess(novaDist, rota.distTotal))
    {
        /*
        std::cout<<"MvIntraRotaShifit pos0("<<pos0<<"); pos1("<<pos1<<");
        rotaR("<<rotaR<<")\n"; std::cout<<"\t"<<rota.printRota()<<"\n"; std::cout<<"\tNova
        Dist: "<<novaDist<<"\n"; std::cout<<"\tMenor\n";
        */

        // sol.distTotal += -rota.distTotal + novaDist;
        // rota.distTotal = novaDist;
        rotaTemp.reset();
        copiaRota(rota, rotaTemp);

        // Shifit rota
        if(pos0 < pos1)
        {
            // std::cout<<"\tpos0<pos1\n";
            for(int i = pos1; i > pos0; --i)
                rotaTemp.vetRota[i] = rotaTemp.vetRota[i-1];

            rotaTemp.vetRota[pos0+1] = clientePos1;
        }
        else
        {
            // std::cout<<"\tpos1<pos0\n";
            for(int i = pos1; i < pos0; ++i)
            {
                rotaTemp.vetRota[i] = rotaTemp.vetRota[i+1];
            }

            rotaTemp.vetRota[pos0] = clientePos1;
        }

        if(testRoute(&rotaTemp.vetRota[0], rotaTemp.numPos, 1))
        {

            sol.distTotal += -rota.distTotal + novaDist;

            copiaRota(rotaTemp, rota);
            rota.distTotal = novaDist;

            return true;
        }



        // std::cout<<"\t"<<rota.printRota()<<"\n\n";

    }

    return false;
}

bool BuscaLocalNS::mvIntraRotaSwap(SolucaoNS::Solucao &sol, int rotaR, int pos0, int pos1)
{

    if(pos0 >= pos1)
        return false;

    Rota &rota = sol.vetRota[rotaR];

    static Rota rotaTemp;

    if(pos0 <= 0 || pos1 <= 0)
        return false;

    if(pos0 >= (rota.numPos - 1) || pos1 >= (rota.numPos - 1))
        return false;

    double novaDist = rota.distTotal;
    int    cliente0 = rota.vetRota[pos0];
    int    cliente1 = rota.vetRota[pos1];

    const bool diff1 = abs(pos0 - pos1) == 1;

    novaDist -= instanciaG.matDist(rota.vetRota[pos0 - 1], cliente0);

    if(!diff1)
    {
        novaDist -= instanciaG.matDist(cliente0, rota.vetRota[pos0 + 1]);
        novaDist -= instanciaG.matDist(rota.vetRota[pos1 - 1], cliente1);
    }

    novaDist -= instanciaG.matDist(cliente1, rota.vetRota[pos1 + 1]);

    novaDist += instanciaG.matDist(rota.vetRota[pos0 - 1], cliente1);

    if(!diff1)
    {
        novaDist += instanciaG.matDist(cliente1, rota.vetRota[pos0 + 1]);
        novaDist += instanciaG.matDist(rota.vetRota[pos1 - 1], cliente0);
    }

    novaDist += instanciaG.matDist(cliente0, rota.vetRota[pos1 + 1]);

    if(doubleLess(novaDist, rota.distTotal))
    {

        rotaTemp.reset();
        copiaRota(rota, rotaTemp);
        /*
        std::cout<<"MvIntraRotaSwap pos0("<<pos0<<"); pos1("<<pos1<<");
        rotaR("<<rotaR<<")\n"; std::cout<<"\tDiff1: "<<diff1<<"\n";
        std::cout<<"\t"<<rota.printRota()<<"\n";
        std::cout<<"\tNova Dist: "<<novaDist<<"\n";
        std::cout<<"\tMenor\n";
        */

        std::swap(rotaTemp.vetRota[pos0], rotaTemp.vetRota[pos1]);

        if(testRoute(&rotaTemp.vetRota[0], rotaTemp.numPos, 1))
        {
            std::swap(rota.vetRota[pos0], rota.vetRota[pos1]);
            sol.distTotal -= rota.distTotal;
            rota.distTotal = novaDist;
            sol.distTotal += novaDist;

            return true;
        }

        // std::cout<<"\t"<<rota.printRota()<<"\n\n";

    }

    return false;
}

bool BuscaLocalNS::mvIntraRota2opt(SolucaoNS::Solucao &sol, int rotaR, int pos0, int pos1)
{

    if(pos0 <= 0 || pos1 <= 0 || pos0 > pos1)
        return false;

    Rota &rota = sol.vetRota[rotaR];
    static Rota rotaTemp;

    if(pos0 >= (rota.numPos - 1) || pos1 >= (rota.numPos - 1))
        return false;

    double novaDist = rota.distTotal;

    int clie0, clie1;
    clie0 = rota.vetRota[pos0];
    clie1 = rota.vetRota[pos1];

    novaDist -= instanciaG.matDist(rota.vetRota[pos0 - 1], clie0);
    novaDist -= instanciaG.matDist(clie1, rota.vetRota[pos1 + 1]);

    novaDist += instanciaG.matDist(clie0, rota.vetRota[pos1 + 1]);
    novaDist += instanciaG.matDist(rota.vetRota[pos0 - 1], clie1);

    if(doubleLess(novaDist, rota.distTotal))
    {

        /*
        std::cout<<"mvIntraRota2opt\n";
        std::cout<<"\tpos0("<<pos0<<"); pos1("<<pos1<<")\n";
        std::cout<<"\t"<<rota.printRota()<<"\n";
         */
        rotaTemp.reset();
        copiaRota(rota, rotaTemp);

        std::reverse((rotaTemp.vetRota.begin() + pos0),
                     (rotaTemp.vetRota.begin() + pos1 + 1));

        if(testRoute(&rotaTemp.vetRota[0], rotaTemp.numPos, 1))
        {
            copiaRota(rotaTemp, rota);

            sol.distTotal -= rota.distTotal;
            sol.distTotal += novaDist;
            rota.distTotal = novaDist;
        }

        // std::cout<<"\t"<<rota.printRota()<<"\n\n";

        return true;
    }

    return false;
}

bool BuscaLocalNS::mvInterRotasShift(SolucaoNS::Solucao &sol,
                                     double              alphaBin,
                                     int                 idRota0,
                                     int                 pos0,
                                     int                 idRota1,
                                     int                 pos1)
{
    // std::cout<<"rota0: "<<idRota0<<"; rota1: "<<idRota1<<"\n";
    // std::cout<<"pos0: "<<pos0<<"; pos1: "<<pos1<<"\n\n";
    Rota &rota0 = sol.vetRota[idRota0];
    Rota &rota1 = sol.vetRota[idRota1];

    static Rota rotaAux0;
    static Rota rotaAux1;

    // std::cout<<"Rota0: "<<rota0.printRota()<<"\n";
    // std::cout<<"Rota1: "<<rota1.printRota()<<"\n\n";

    // std::cout<<"BIN1: \n"<<sol.vetBin[idRota1]<<"\n\n";

    int clie = sol.vetRota[idRota1].vetRota[pos1];

    if(pos0 <= 0 || pos0 >= (rota0.numPos - 1))
        return false;

    if(pos1 <= 0 || pos1 >= (rota1.numPos - 1))
        return false;

    double novaDistRota0 = rota0.distTotal;
    double novaDistRota1 = rota1.distTotal;

    // Remove clie de rota1
    novaDistRota1 += -instanciaG.matDist(rota1.vetRota[pos1 - 1], clie);
    novaDistRota1 += -instanciaG.matDist(clie, rota1.vetRota[pos1 + 1]);
    novaDistRota1 += instanciaG.matDist(rota1.vetRota[pos1 - 1], rota1.vetRota[pos1 + 1]);

    // Add clie na rota0
    novaDistRota0 += -instanciaG.matDist(rota0.vetRota[pos0], rota0.vetRota[pos0 + 1]);
    novaDistRota0 += instanciaG.matDist(rota0.vetRota[pos0], clie);
    novaDistRota0 += instanciaG.matDist(clie, rota0.vetRota[pos0 + 1]);

    if(!doubleLess((novaDistRota0 + novaDistRota1), (rota0.distTotal + rota1.distTotal)))
        return false;


    copiaRota(rota0, rotaAux0);
    copiaRota(rota1, rotaAux1);

    // Atualiza Rotas
    // Rota0

    for(int i = (rotaAux0.numPos - 1); i > pos0; --i)
        rotaAux0.vetRota[i + 1] = rotaAux0.vetRota[i];

    rotaAux0.vetRota[pos0 + 1] = clie;
    rotaAux0.numPos += 1;

    // Rota1
    for(int i = pos1; i < (rota1.numPos - 1); ++i)
        rotaAux1.vetRota[i] = rotaAux1.vetRota[i + 1];

    rotaAux1.numPos -= 1;

    std::swap(rotaAux0.vetDemClie[clie], rotaAux1.vetDemClie[clie]);

    bool feasible = testRoute(&rotaAux0.vetRota[0], rotaAux0.numPos, 1);
    if(!feasible)
        return false;

    feasible = testRoute(&rotaAux1.vetRota[0], rotaAux1.numPos, 1);
    if(!feasible)
        return false;

    copiaRota(rotaAux0, rota0);
    copiaRota(rotaAux1, rota1);

    // Atualiza distancias
    sol.distTotal += -(rota0.distTotal + rota1.distTotal);
    sol.distTotal += (novaDistRota0 + novaDistRota1);

    rota0.distTotal = novaDistRota0;
    rota1.distTotal = novaDistRota1;

    // std::cout<<"Rota0: "<<rota0.printRota()<<"\n";
    // std::cout<<"Rota1: "<<rota1.printRota()<<"\n\n";

    // std::cout<<"BIN1: \n"<<sol.vetBin[idRota1]<<"\n\n";

    return true;
}

bool BuscaLocalNS::mvInterRotasSwap(SolucaoNS::Solucao &sol,
                                    double              alphaBin,
                                    int                 idRota0,
                                    int                 pos0,
                                    int                 idRota1,
                                    int                 pos1)
{

    //static Vector<Bin> bin1(1);
    //static Vector<Bin> bin0(1);
    //static VectorI     vetItens0(instanciaG.maxNumItensPorClie);
    //static VectorI     vetItens1(instanciaG.maxNumItensPorClie);

    Rota &rota0 = sol.vetRota[idRota0];
    Rota &rota1 = sol.vetRota[idRota1];

    static Rota rotaAux0;
    static Rota rotaAux1;

    if(pos0 <= 0 || pos1 <= 0)
        return false;

    if(pos0 >= (rota0.numPos - 1) || pos1 >= (rota1.numPos - 1))
        return false;

    const int cliente0 = rota0.vetRota[pos0];
    const int cliente1 = rota1.vetRota[pos1];

    double novaDist0 = rota0.distTotal;
    double novaDist1 = rota1.distTotal;

    // Calcula novas distancias
    // Rm cliente0
    novaDist0 += -instanciaG.matDist(rota0.vetRota[pos0 - 1], cliente0);
    novaDist0 += -instanciaG.matDist(cliente0, rota0.vetRota[pos0 + 1]);

    // Rm cliente1
    novaDist1 += -instanciaG.matDist(rota1.vetRota[pos1 - 1], cliente1);
    novaDist1 += -instanciaG.matDist(cliente1, rota1.vetRota[pos1 + 1]);

    // Add cliente1 a rota0
    novaDist0 += instanciaG.matDist(rota0.vetRota[pos0 - 1], cliente1);
    novaDist0 += instanciaG.matDist(cliente1, rota0.vetRota[pos0 + 1]);

    // Add cliente0 a rota1
    novaDist1 += instanciaG.matDist(rota1.vetRota[pos1 - 1], cliente0);
    novaDist1 += instanciaG.matDist(cliente0, rota1.vetRota[pos1 + 1]);

    if(!doubleLess((novaDist0 + novaDist1), (rota0.distTotal + rota1.distTotal)))
        return false;

    copiaRota(rota0, rotaAux0);
    copiaRota(rota1, rotaAux1);

    std::swap(rotaAux0.vetRota[pos0], rotaAux1.vetRota[pos1]);

    // demandas
    rotaAux1.vetDemClie[cliente0] = instanciaG.vetDemandaCliente[cliente0];
    rotaAux1.vetDemClie[cliente1] = 0.0;

    rotaAux0.vetDemClie[cliente1] = instanciaG.vetDemandaCliente[cliente1];
    rotaAux0.vetDemClie[cliente0] = 0.0;

    bool feasible = testRoute(&rotaAux0.vetRota[0], rotaAux0.numPos, 1);
    if(!feasible)
        return false;

    feasible = testRoute(&rotaAux1.vetRota[0], rotaAux1.numPos, 1);
    if(!feasible)
        return false;

    copiaRota(rotaAux0, rota0);
    copiaRota(rotaAux1, rota1);

    sol.distTotal += -(rota0.distTotal + rota1.distTotal);
    sol.distTotal += (novaDist0 + novaDist1);

    rota0.distTotal = novaDist0;
    rota1.distTotal = novaDist1;

    return true;
}

bool BuscaLocalNS::mvInterRotasCross(SolucaoNS::Solucao &sol,
                                     double              alphaBin,
                                     int                 idRota0,
                                     int                 pos0,
                                     int                 idRota1,
                                     int                 pos1)
{

    // std::cout<<"********************************************************\n\n";

    static VectorI vetClientesRota0(TamRota);
    static VectorI vetClientesRota1(TamRota);
    //static VectorI vetItensRota0(NumItensPorBin);
    //static VectorI vetItensRota1(NumItensPorBin);

    static Rota rotaAux0;
    static Rota rotaAux1;

    //static Vector<Bin> bin0(1);
    //static Vector<Bin> bin1(1);

    int tamVetC_Rota0 = 0;
    int tamVetC_Rota1 = 0;
    int numItensRota0 = 0;
    int numItensRota1 = 0;

    double distSubRota0  = 0.0;
    double distSubRota1  = 0.0;
    double novaDistRota0 = 0.0;
    double novaDistRota1 = 0.0;

    if(pos0 < 0 || pos1 < 0)
        return false;

    Rota &rota0 = sol.vetRota[idRota0];
    Rota &rota1 = sol.vetRota[idRota1];

    if((pos0 >= (rota0.numPos - 1)) || (pos1 >= (rota1.numPos - 1)))
        return false;

    if(rota0.numPos == 0 && rota1.numPos == 0)
        return false;

    // std::cout<<"Rota0: "<<rota0.printRota()<<"\npos0: "<<pos0<<"\n\n";
    // std::cout<<"Rota1: "<<rota1.printRota()<<"\npos1: "<<pos1<<"\n\n";

    auto getClientes = [](Rota &rota, VectorI &vet, int &num, int pos)
    {
        num = 0;
        for(int i = pos; i < (rota.numPos - 1); ++i)
        {
            vet[num] = rota.vetRota[i];
            num += 1;
        }
    };

    // Pega todos os clientes da rota0 apos pos0 e seus itens
    getClientes(rota0, vetClientesRota0, tamVetC_Rota0, pos0 + 1);
    //numItensRota0 = copiaItensClientes(vetClientesRota0, tamVetC_Rota0, vetItensRota0);
    vetClientesRota0[tamVetC_Rota0] = 0;
    tamVetC_Rota0 += 1;
    distSubRota0 = calculaDistancia(vetClientesRota0, tamVetC_Rota0);

    // std::cout<<"Clientes Rota0: "<<printVet(vetClientesRota0, tamVetC_Rota0)<<"\nItens:
    // "<<printVet(vetItensRota0, numItensRota0)<<"\n";
    //     std::cout<<"\n\n";

    // Pega todos os clientes da rota1 apos pos1 e seus itens
    getClientes(rota1, vetClientesRota1, tamVetC_Rota1, pos1 + 1);
    //numItensRota1 = copiaItensClientes(vetClientesRota1, tamVetC_Rota1, vetItensRota1);
    vetClientesRota1[tamVetC_Rota1] = 0;
    tamVetC_Rota1 += 1;
    distSubRota1 = calculaDistancia(vetClientesRota1, tamVetC_Rota1);

    // std::cout<<"Clientes Rota1: "<<printVet(vetClientesRota1, tamVetC_Rota1)<<"\nItens:
    // "<<printVet(vetItensRota1, numItensRota1)<<"\n\n";

    /*
    std::cout<<"Itens BIN0: "<<printVet(rota0.binPtr->vetItemId,
    rota0.binPtr->numItens)<<"\n"; std::cout<<"Itens BIN1:
    "<<printVet(rota1.binPtr->vetItemId, rota1.binPtr->numItens)<<"\n"; std::cout<<"\n\n";
    */
    novaDistRota0 = rota0.distTotal;
    novaDistRota1 = rota1.distTotal;

    // Remove a distancia dos clientes da subRota0
    novaDistRota0 += -distSubRota0;
    novaDistRota0 += -instanciaG.matDist(rota0.vetRota[pos0], rota0.vetRota[pos0 + 1]);

    // Remove a distancia dos clientes da subRota1
    novaDistRota1 += -distSubRota1;
    novaDistRota1 += -instanciaG.matDist(rota1.vetRota[pos1], rota1.vetRota[pos1 + 1]);

    // Add a distancia dos clientes da subRota1
    novaDistRota0 += distSubRota1;
    if(tamVetC_Rota1 > 0)
    {
        novaDistRota0 += instanciaG.matDist(rota0.vetRota[pos0], vetClientesRota1[0]);
        // novaDistRota0 += instanciaG.matDist(vetClientesRota1[tamVetC_Rota1 - 1], 0);
    }

    // Add a distancia dos clientes da subRota0
    novaDistRota1 += distSubRota0;
    if(tamVetC_Rota0 > 0)
    {
        novaDistRota1 += instanciaG.matDist(rota1.vetRota[pos1], vetClientesRota0[0]);
        // novaDistRota1 += instanciaG.matDist(vetClientesRota0[tamVetC_Rota0 - 1], 0);
    }

    if(!doubleLess((novaDistRota0 + novaDistRota1), (rota0.distTotal + rota1.distTotal)))
        return false;


    copiaRota(rota0, rotaAux0);
    copiaRota(rota1, rotaAux1);

    // Copia subRota1 para rota0
    int pos = pos0 + 1;
    for(int i = 0; i < tamVetC_Rota1; ++i)
    {
        int cliente = vetClientesRota1[i];
        rotaAux0.vetRota[pos] = cliente;
        pos += 1;
        rotaAux0.vetDemClie[cliente] = instanciaG.vetDemandaCliente[cliente];
    }

    rotaAux0.numPos = pos;
    // rota0.vetRota[rota0.numPos-1] = 0;

    // Remove demanda da subRota0 da rota0
    for(int i = 0; i < tamVetC_Rota0; ++i)
        rotaAux0.vetDemClie[vetClientesRota0[i]] = 0;

    // Copia subRota1 para rota0
    pos = pos1 + 1;
    for(int i = 0; i < tamVetC_Rota0; ++i)
    {
        int cliente = vetClientesRota0[i];
        rotaAux1.vetRota[pos] = cliente;
        pos += 1;
        rotaAux1.vetDemClie[cliente] = instanciaG.vetDemandaCliente[cliente];
    }

    rotaAux1.numPos = pos;
    // rota1.vetRota[rota1.numPos-1] = 0;

    // Remove demanda da subRota0 da rota0
    for(int i = 0; i < tamVetC_Rota1; ++i)
        rotaAux1.vetDemClie[vetClientesRota1[i]] = 0;

    bool feasible = testRoute(&rotaAux0.vetRota[0], rotaAux0.numPos, 1);
    if(!feasible)
        return false;

    feasible = testRoute(&rotaAux1.vetRota[0], rotaAux1.numPos, 1);
    if(!feasible)
        return false;

    copiaRota(rotaAux0, rota0);
    copiaRota(rotaAux1, rota1);

    // Corrige distancias
    sol.distTotal += -(rota0.distTotal + rota1.distTotal);
    sol.distTotal += novaDistRota0 + novaDistRota1;

    rota0.distTotal = novaDistRota0;
    rota1.distTotal = novaDistRota1;

    /*
    std::cout<<"*Itens BIN0: "<<printVet(rota0.binPtr->vetItemId,
    rota0.binPtr->numItens)<<"\n"; std::cout<<"*Itens BIN1:
    "<<printVet(rota1.binPtr->vetItemId, rota1.binPtr->numItens)<<"\n"; std::cout<<"Rota0:
    "<<rota0.printRota()<<"\n\n"; std::cout<<"Rota1: "<<rota1.printRota()<<"\n\n";
    */

    return true;
}

void BuscaLocalNS::rvnd(SolucaoNS::Solucao &sol)
{
    //NumMv
    static Array<MV, NumMv> arrayMv;
    // arrayMv[0] = MvInterRotasCross;



    for(int i = 0; i < NumMv; ++i)
    {
        auto checkFuncm = [&]() -> bool
        {
            int choice = arrayMv[i];
            for(int j = 0; j < i; ++j)
            {
                if(choice == arrayMv[j])
                    return false;
            }

            return true;
        };

        do
        {
            arrayMv[i] = static_cast<MV>(getRandInt(0, (NumMv-1)));
        }
        while(!checkFuncm());
    }


    bool improve = true;
    int  k = 0;


    while(k < NumMv)
    {

        switch(arrayMv[k])
        {
        case MvIntraRotaShifit:
            improve = buscaIntraRota(sol, intraRotaShift);
            break;

        case MvIntraRotaSwap:
            improve = buscaIntraRota(sol, mvIntraRotaSwap);
            break;

        case MvIntraRota2opt:
            improve = buscaIntraRota(sol, mvIntraRota2opt);
            break;

        case MvInterRotasShifit:
            improve = buscaInterRotas(
                sol, MvInterRotasShifit, input.aphaBin, mvInterRotasShift);
            break;

        case MvInterRotasSwap:
            improve =
                buscaInterRotas(sol, MvInterRotasSwap, input.aphaBin, mvInterRotasSwap);
            break;

        case MvInterRotasCross:
            improve =
                buscaInterRotas(sol, MvInterRotasCross, input.aphaBin, mvInterRotasCross);
            break;
        }

        if(improve)
            k = 0;
        else
            k += 1;
    }
}
