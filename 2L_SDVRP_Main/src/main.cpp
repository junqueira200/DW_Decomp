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

using namespace InstanceNS;
using namespace ConstrutivoBinNS;
using namespace SolucaoNS;
using namespace RandNs;
using namespace ParseInputNS;
using namespace ConstrutivoNS;
using namespace IgNs;
using namespace BinPackingCP_NS;

Array<double, 3> getDim(Item& item, Rotation r)
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

    if(input.inst2d)
        InstanceNS::read2dInstance(input.strInstCompleto);
    else
        InstanceNS::read3dInstance(input.strInstCompleto);

    //std::cout<<"clie: "<<instanciaG.numClientes<<"\t"<<instanciaG.numItens<<"\n";
    //std::cout<<instanciaG.matDist<<"\n\n";
    //escreveFileNum(input.file.fileNum);

    //return 0;

    std::cout<<"\n";

    try
    {
        clock_t start = clock();
        //double ompStart = omp_get_wtime();

        //Solucao best(instanciaG);
        //metaheuristicaIg(best);


        int heuristica, exato, inviavel, total;
        heuristica = exato = inviavel = total = 0;

        for(int i=0; i < 600; ++i)
        {
            if((i%200) == 0)
                std::cout<<"IT: "<<i<<"\n";

            Resultado resul =  testaCpSatBinPacking(6);
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


        clock_t end = clock();
        //double ompEnd = omp_get_wtime();

        output.tempoCpu = double(end-start)/CLOCKS_PER_SEC;
        //output.tempoRelogio = ompEnd - ompStart;

        //std::cout<<"Num Clientes: "<<instanciaG.numClientes<<"\n";
        //std::cout<<"Num Veic: "<<instanciaG.numVeiculos<<"\n\n";

        //escreveSaidas(best, output.tempoCpu);
        //int tamVet = std::min(instanciaG.numItens, 15);
        //testaCpSatBinPacking(tamVet);

        std::cout<<"\nTotal: \t\t"<<total<<"\nHEURISTIC: \t"<<heuristica<<"\nCP-SAT: \t"<<exato<<"\nUNFEASIBLE: \t"<<inviavel<<"\n\n";

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
