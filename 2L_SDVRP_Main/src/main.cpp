/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: main.cpp
 * ****************************************
 * ****************************************/

#include "InputOutput.h"
#include "Instancia.h"
#include "ConstrutivoBin.h"
#include "Construtivo.h"
#include "AuxT.h"
#include "rand.h"
#include "Ig.h"
#include "BinPackingCP.h"
#include "TesteOroloc3D.h"


#include "ProblemParameters.h"
#include "BCRoutingParams.h"
#include "LoadingChecker.h"

using namespace InstanceNS;
using namespace ConstrutivoBinNS;
using namespace SolucaoNS;
using namespace RandNs;
using namespace ParseInputNS;
using namespace ConstrutivoNS;
using namespace IgNs;
using namespace BinPackingCP_NS;
using namespace TesteOroloc3D_NS;

using namespace ContainerLoading;
using namespace VehicleRouting;
using namespace VehicleRouting::Algorithms;


Array<double, 3> getDim(Item& item, InstanceNS::Rotation r)
{
    Array<double,3> array;
    array[0] = item.getDimRotacionada(0, r);
    array[1] = item.getDimRotacionada(1, r);
    array[2] = item.getDimRotacionada(2, r);

    return array;
}



int main(int argc, const char* argv[])
{

    //Item item(1, 2, 3, 1);

    //std::cout<<"Rot1: "<<getDim(item, Rot0)<<"\n";
    //return 0;

    ParseInputNS::parseInput(argc, argv);
    output.setup();
    std::cout << "INST: " << input.strInst << " SEMENTE: " << RandNs::estado_ << " " << output.data << "";

    if(input.instOroloc3D)
        InstanceNS::readOroloc3D(input.strInstCompleto);
    else if(input.inst2d)
        InstanceNS::read2dInstance(input.strInstCompleto);
    else
        InstanceNS::read3dInstance(input.strInstCompleto);

    testeOroloc3D();
    return 0;

    //std::cout<<"clie: "<<instanciaG.numClientes<<"\t"<<instanciaG.numItens<<"\n";
    //std::cout<<instanciaG.matDist<<"\n\n";
    //escreveFileNum(input.file.fileNum);

    //return 0;



    try
    {
        clock_t start = clock();
        //double ompStart = omp_get_wtime();

        //Solucao best(instanciaG);
        //metaheuristicaIg(best);

        InputParameters inputParam;
        inputParam.ContainerLoading.LoadingProblem.Variant = LoadingProblemParams::VariantType::AllConstraints;
        inputParam.SetLoadingFlags();

        LoadingChecker loadingChecker(inputParam.ContainerLoading);
        Container container((int)instanciaG.vetDimVeiculo[0], (int)instanciaG.vetDimVeiculo[1],
                            (int)instanciaG.vetDimVeiculo[2], (int)instanciaG.maxPayload);

        int numItens = 5;

        VectorI vetItens(numItens);
        VectorI vetCustomers(numItens);
        vetItens.setAll(0);
        Collections::IdVector stopIds(numItens);

        for(int i=1; i <= numItens; ++i)
        {
            vetCustomers[i-1] = i;
            stopIds[i-1] = i;
        }

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
                 (volumeOcupado+instanciaG.vetItens[itemId].volume) > volumeVeiculo/2.0 ||
                  demanda + instanciaG.vetItens[itemId].weight > instanciaG.maxPayload)
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
            demanda       += instanciaG.vetItens[itemId].weight;
            vetItens[t] = itemId;

        }

        std::vector<Cuboid> vetCuboids;
        convertVectorOfItensToVectorOfCuboids(vetItens, vetCuboids);


        auto status = loadingChecker.ConstraintProgrammingSolver(PackingType::Complete, container, stopIds, vetCuboids,
                                                                 std::numeric_limits<double>::infinity());
        if(status != LoadingStatus::FeasOpt)
        {
            std::printf("Error in ConstraintProgrammingSolver\n");
            if(status == LoadingStatus::Infeasible)
                std::printf("INFEASIBLE\n");
        }
        else
        {
            std::printf("Loading was successful!\n");
        }

        int heuristica, exato, inviavel, total;
        heuristica = exato = inviavel = total = 0;
        /*
        for(int i=0; i < 600; ++i)
        {
            if((i%200) == 0)
                std::cout<<"IT: "<<i<<"\n";

            Resultado resul =  testaCpSatBinPacking(50);
            total += 1;
            switch (resul)
            {
            case HEURISTICA:
                heuristica += 1;
                break;

            case EXATO:
                exato += 1;
                break;

            case INVIAVEL:
                inviavel += 1;
                break;

            }
        }

        */
        clock_t end = clock();
        //double ompEnd = omp_get_wtime();

        output.tempoCpu = double(end-start)/CLOCKS_PER_SEC;
        //output.tempoRelogio = ompEnd - ompStart;

        //std::cout<<"Num Clientes: "<<instanciaG.numClientes<<"\n";
        //std::cout<<"Num Veic: "<<instanciaG.numVeiculos<<"\n\n";

        //escreveSaidas(best, output.tempoCpu);
        //int tamVet = std::min(instanciaG.numItens, 15);
        //testaCpSatBinPacking(tamVet);

        //std::cout<<"\nTotal: \t\t"<<total<<"\nHEURISTIC: \t"<<heuristica<<"\nCP-SAT: \t"<<exato<<"\nUNFEASIBLE: \t"<<inviavel<<"\n\n";

        //std::cout<<"TOTAL TIME : "<<output.tempoCpu<<" s\n\n";
        std::printf("TOTAL TIME: \t%.2f s\n\n", output.tempoCpu);

    }
    catch (char const* exeption)
    {
        std::cout<<"\nCAT EXEPTION: "<<exeption<<"\n";
        return -1;
    }

    return 0;
}
