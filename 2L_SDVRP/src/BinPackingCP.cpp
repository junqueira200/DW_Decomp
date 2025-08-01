/* ****************************************
 * ****************************************
 *  Data:    10/02/25
 *  Arquivo: BinPackingCP.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#include "BinPackingCP.h"
#include "Instancia.h"
#include "rand.h"
#include "ConstrutivoBin.h"
#include "InputOutput.h"


using namespace SolucaoNS;
using namespace InstanciaNS;
using namespace RandNs;
using namespace BinPackingCP_NS;
using namespace ConstrutivoBinNS;
using namespace ParseInputNS;

bool BinPackingCP_NS::cpSatBinPacking(SolucaoNS::Bin &binResult, VectorI &vetItens, int tamVet)
{
    if(!input.cpSat)
    {   std::cout<<"cp-sat igual a false\n";
        return false;
    }

    binResult.reset();

    CpModelBuilder model;
    SatParameters params;
    params.set_num_workers(1);
    params.set_stop_after_first_solution(true);

    if(input.cpSatTime > 0.0)
        params.set_max_time_in_seconds(input.cpSatTime);

    static Vector<IntVar> vetX(instanciaG.numItens);
    static Vector<IntVar> vetY(instanciaG.numItens);

    static Vector<IntVar> vetDx(instanciaG.numItens);
    static Vector<IntVar> vetDy(instanciaG.numItens);

    static Vector<IntVar> vetR(instanciaG.numItens);


    for(int i=0; i < tamVet; ++i)
    {
        const Item item = instanciaG.vetItens[vetItens[i]];
        Domain domainX(0, ((int)instanciaG.vetDimVeiculo[0]));
        Domain domainY(0, ((int)instanciaG.vetDimVeiculo[1]));

        IntVar x = model.NewIntVar(domainX).WithName("x_" + std::to_string(i));
        vetX[i] = x;


        IntVar y = model.NewIntVar(domainY).WithName("y_" + std::to_string(i));
        vetY[i] = y;

        int max = (int)std::max(item.vetDim[0], item.vetDim[1]);
        int min = (int)std::min(item.vetDim[0], item.vetDim[1]);

        Domain domainDxDy(min, max);

        IntVar dx = model.NewIntVar(domainDxDy).WithName("dx_"+std::to_string(i));
        vetDx[i] = dx;


        IntVar dy = model.NewIntVar(domainDxDy).WithName("dy_"+std::to_string(i));
        vetDy[i] = dy;

        Domain domainRot(0, 1);
        IntVar r = model.NewIntVar(domainRot).WithName("r_"+std::to_string(i));
        vetR[i] = r;

        model.AddLessOrEqual(x+dx, (int)instanciaG.vetDimVeiculo[0]);
        model.AddLessOrEqual(y+dy, (int)instanciaG.vetDimVeiculo[1]);

        model.AddElement(r, {(int)item.vetDim[0], (int)item.vetDim[1]}, dx);
        model.AddElement(r, {(int)item.vetDim[1], (int)item.vetDim[0]}, dy);

        if((int)item.vetDim[0] == (int)item.vetDim[1])
            model.AddEquality(r, 0);

    }

    for(int i=0; i < tamVet-1; ++i)
    {
        for(int j=i+1; j < tamVet; ++j)
        {
            BoolVar no_overlap_x1 = model.NewBoolVar();
            BoolVar no_overlap_x2 = model.NewBoolVar();
            BoolVar no_overlap_y1 = model.NewBoolVar();
            BoolVar no_overlap_y2 = model.NewBoolVar();

            model.AddLessOrEqual(vetX[i] + vetDx[i], vetX[j]).OnlyEnforceIf(no_overlap_x1);
            model.AddLessOrEqual(vetX[j] + vetDx[j], vetX[i]).OnlyEnforceIf(no_overlap_x2);


            model.AddLessOrEqual(vetY[i] + vetDy[i], vetY[j]).OnlyEnforceIf(no_overlap_y1);
            model.AddLessOrEqual(vetY[j] + vetDy[j], vetY[i]).OnlyEnforceIf(no_overlap_y2);

            model.AddBoolOr({no_overlap_x1, no_overlap_x2, no_overlap_y1, no_overlap_y2});
        }
    }

    const CpSolverResponse response = SolveWithParameters(model.Build(), params);

    if(response.status() == FEASIBLE || response.status() == OPTIMAL)
    {
        for(int i = 0; i < tamVet; ++i)
        {
            Item item = instanciaG.vetItens[vetItens[i]];

            int x  = SolutionIntegerValue(response, vetX[i]);
            int y  = SolutionIntegerValue(response, vetY[i]);
            int r  = SolutionIntegerValue(response, vetR[i]);

            int dx = SolutionIntegerValue(response, vetDx[i]);
            int dy = SolutionIntegerValue(response, vetDy[i]);

            //bin.addItem(vetItens[i], )
            binResult.vetItemId[i] = vetItens[i];
            binResult.vetItens[vetItens[i]] = 1;
            binResult.vetPosItem[i] = Ponto(x, y, 0);
            binResult.vetRotacao[i] = Rotacao(r);
            binResult.volumeOcupado += instanciaG.vetItens[vetItens[i]].volume;
            binResult.demandaTotal  += instanciaG.vetItens[vetItens[i]].peso;
            binResult.numItens += 1;

            //std::cout<<x<<" "<<y<<"; "<<dx<<" "<<dy<<"\n";
            //std::cout<<instanciaG.vetItens[vetItens[i]].vetDim<<"\n\n";
        }

        if(!binResult.verificaViabilidade())
        {

            std::cout<<"CP-SAT Encontrou Solucao Viavel!\n";
            std::cout << "ERROR, Bin NAO eh Viavel!\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

        criaEPs(binResult);
        return true;
    }

    return false;

}

void BinPackingCP_NS::testaCpSatBinPacking(int numItens)
{
    VectorI vetItens(numItens);
    vetItens.setAll(0);

    Vector<int8_t> vetItensSelecionados(instanciaG.numItens);
    vetItensSelecionados.setAll((int8_t)0);

    double volumeOcupado = 0.0;
    double volumeVeiculo = 1.0;
    double demanda       = 0.0;

    for(int i=0; i < instanciaG.numDim; ++i)
        volumeVeiculo *= instanciaG.vetDimVeiculo[i];


    for(int t=0; t < numItens; ++t)
    {
        int itemId = getRandInt(0, instanciaG.numItens-1);
        const int itemIdIni = itemId;
        while(vetItensSelecionados[itemId] == (int8_t)1 ||
             (volumeOcupado+instanciaG.vetItens[itemId].volume) > volumeVeiculo ||
             demanda + instanciaG.vetItens[itemId].peso > instanciaG.veicCap)
        {
            itemId = (itemId+1)%instanciaG.numItens;

            if(itemId == itemIdIni)
            {
                itemId = -1;
                break;
            }
        }

        if(itemId == -1)
        {
            numItens = t;
            break;
        }

        vetItensSelecionados[itemId] = (int8_t)1;
        volumeOcupado += instanciaG.vetItens[itemId].volume;
        demanda       += instanciaG.vetItens[itemId].peso;
        vetItens[t] = itemId;

    }

    std::cout<<"\nNum itens: "<<numItens<<"\n\n";

    Bin binCpSat;
    if(cpSatBinPacking(binCpSat, vetItens, numItens))
        std::cout<<"CP-SAT Encontrou Solucao Viavel!\n";

    Bin binHeu;
    binHeu.reset();
    //binHeu.addEp(PontoZero);
    //binHeu.numEps = 1;

    std::sort(vetItens.begin(), vetItens.begin()+numItens);
    std::cout<<vetItens<<"\n\n";

    if(construtivoBinPacking(binHeu, vetItens, numItens, input.aphaBin, 2))
    {
        std::cout << "Construtivo Encontrou Solucao Viavel!\n";
        if(binHeu.verificaViabilidade())
            std::cout<<"Bin Viavel!\n";
        else
            std::cout<<"Bin NAO eh Viavel\n";
    }
    else
        std::cout<<"Construtivo NAO Encontrou Solucao Viavel!\n";

}

void BinPackingCP_NS::criaEPs(SolucaoNS::Bin &bin)
{

    for(int i=0; i < bin.numItens; ++i)
    {
        Ponto p    = bin.vetPosItem[i];
        Item& item = instanciaG.vetItens[bin.vetItemId[i]];

        Array<double, 3> vetD;

        for(int d=0; d < instanciaG.numDim; ++d)
            vetD[d] = item.getDimRotacionada(d, bin.vetRotacao[i]);

        for(int d=0; d < instanciaG.numDim; ++d)
        {
            Ponto ep(p);
            ep.vetDim[d] += vetD[d];
            bool colisao = false;

            for(int t=0; t < bin.numItens; ++t)
            {
                if(epColideItem(ep, bin.vetPosItem[t], bin.vetItemId[t]))
                {
                    colisao = true;
                    break;
                }
            }

            if(!colisao)
                bin.addEp(ep);
        }

    }

}