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

using namespace InstanciaNS;
using namespace ConstrutivoBinNS;
using namespace SolucaoNS;
using namespace RandNs;
using namespace ParseInputNS;
using namespace ConstrutivoNS;
using namespace IgNs;
using namespace BinPackingCP_NS;


int main(int argc, const char* argv[])
{

    ParseInputNS::parseInput(argc, argv);
    output.setup();
    std::cout << "INST: " << input.strInst << " SEMENTE: " << RandNs::estado_ << " " << output.data << "";

    InstanciaNS::leInstancia(input.strInstCompleto);
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

        int heuristica, exato, inviavel, total;
        heuristica = exato = inviavel = total = 0;

        for(int i=0; i < 100; ++i)
        {
            if((i%10) == 0)
                std::cout<<i<<"\n";
            Resultado resul =  testaCpSatBinPacking(9999);
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

        std::cout<<"Total: \t\t"<<total<<"\nHeuristica: \t"<<heuristica<<"\nExato: \t\t"<<exato<<
                  "\nInviavel: \t"<<inviavel<<"\n\n";

        std::cout<<"Tempo Total: "<<output.tempoCpu<<"\n\n";

    }
    catch (char const* exeption)
    {
        std::cout<<"\nCAT EXEPTION: "<<exeption<<"\n";
        return -1;
    }

    return 0;
}
