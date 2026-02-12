/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: Instancia.cpp
 * ****************************************
 * ****************************************/

#include "Instancia.h"
#include "AuxT.h"
#include <cuchar>
#include <fstream>
#include "InputOutput.h"

using namespace InstanceNS;
using namespace ParseInputNS;

InstanceNS::Instance::Instance(int numClientes_, int numItens_, int numVeiculos_):
                                                                   matDist(numClientes_, numClientes_),
                                                                   matTempo(numClientes, numClientes),
                                                                   //vetItens(numItens_),
                                                                   vetPesoItens(numItens_),
                                                                   //vetItemLargura(numItens_),
                                                                   vetDemandaCliente(numClientes_),
                                                                   vetItemCliente(numItens_),
                                                                   matCliItensIniFim(0, 0),
                                                                   vetNumItensPorCli(numClientes_),
                                                                   vetMinDimItens(numItens_),
                                                                   vetTw(numClientes_)

{

    numClientes = numClientes_;
    numItens    = numItens_;
    numVeiculos = numVeiculos_;

    if(numClientes < 0 || numItens < (numClientes-1))
    {
        assertm(true, "numClientes("<<numClientes<<") < 0 || numItens("<<numItens<<") < numClientes-1("<<numClientes<<")");
    }

    matCliItensIniFim = Matrix<int>(numClientes, 2);
    matCliItensIniFim.setVal(-1);

}

InstanceNS::Instance::Instance():
                                    matDist(0, 0),
                                    matTempo(0,0),
                                    vetItens(0),
                                    vetDemandaCliente(0),
                                    vetItemCliente(0, 0),
                                    matCliItensIniFim(0, 0),
                                    vetPesoItens(0)

{

}


void InstanceNS::read2dInstance(const std::string &strFile)
{

    int numClientes, numVeiculos, numItens;

    std::ifstream file(strFile);
    assertm(!file.is_open(), "Nao foi possivel abrir o arquivo: "<<strFile);

    //std::cout<<"File("<<strFile<<") aberto\n";

    std::string lineLixo;
    int num = 2;
    if(input.splitInstancia)
        num = 1;
    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numClientes;
    numClientes += 1;

    for(int i=0; i < 1; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numVeiculos;

    for(int i=0; i < 1; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numItens;
    //numItens += 1;

//std::cout<<"numClientes: "<<numClientes<<"; numItens: "<<numItens<<" numVeic: "<<numVeiculos<<"\n";


    if(!input.splitInstancia)
    {
        for(int i = 0; i < 2; ++i)
        {
            getline(file, lineLixo);
            //std::cout << "line: " << lineLixo << "\n";
        }
    }

    instanciaG = Instance(numClientes, numItens, numVeiculos);
    instanciaG.matTempo.setVal(1.0);
    instanciaG.numDim = 2;
    instanciaG.nome = strFile;
    instanciaG.numRotation = 2;

    file >> instanciaG.maxPayload >> instanciaG.vetDimVeiculo[0] >> instanciaG.vetDimVeiculo[1];

//std::cout << "veicCap: " << instanciaG.veicCap << "; veicComprimento: " << instanciaG.vetDimVeiculo[0] << "; veicLargura: " << instanciaG.vetDimVeiculo[1] << "\n\n";

    num = 2;
    if(input.splitInstancia)
        num = 1;

    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    Matrix<double> matPontos(numClientes, 2);
    matPontos.setVal(-1.0);

    for(int i=0; i < numClientes; ++i)
    {
        double x, y, cli;
        file >> cli >> x >> y >> instanciaG.vetDemandaCliente[i];
        matPontos.get(i, 0) = x;
        matPontos.get(i, 1) = y;

//std::cout<<cli<<" "<<x<<" "<<y<<" "<<instanciaG.vetDemandaCliente[i]<<"\n";

    }

    for(int i=0; i < numClientes; ++i)
    {
        instanciaG.matDist.get(i, i) = 0.0;

        for(int j=i+1; j < numClientes; ++j)
        {
            double dist = std::sqrt(std::pow(matPontos(i, 0) - matPontos(j, 0), 2) +
                                      std::pow(matPontos(i, 1) - matPontos(j, 1), 2));

            instanciaG.matDist.get(i, j) = dist;
            instanciaG.matDist.get(j, i) = dist;
        }
    }

    //std::cout << instanciaG.matDist << "\n";

    num = 2;
    if(input.splitInstancia)
        num = 1;
    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    Vector<Vector<double>> vetVetItens(numClientes);
    int maxNumItensPorCli = 0;
    int nextItem = 0;


    //instanciaG.vetItemAltura.reserve(instanciaG.numItens);
    //instanciaG.vetItemLargura.reserve(instanciaG.numItens);
    //instanciaG.vetItemArea.reserve(instanciaG.numItens);

    instanciaG.vetPesoItens.reserve(instanciaG.numItens);
    instanciaG.vetItens.reserve(instanciaG.numItens);

    num = numClientes;
    if(input.splitInstancia)
        num -= 1;
    for(int i=0; i < num; ++i)
    {
        int node, numItensPorClie;

        file >> node >> numItensPorClie;
//std::cout<<node<<" "<<numItensPorClie<<"\n";

        maxNumItensPorCli = std::max(maxNumItensPorCli, numItensPorClie);

        if(node != 0)
        {
            instanciaG.matCliItensIniFim.get(node, 0) = nextItem;
            //instancia.matCliItensIniFim.get(node, 1) = nextItem+(numItensPorClie-1);
        }

        for(int item = 0; item < numItensPorClie; ++item)
        {
            double altura, comprimento;
            file >> altura >> comprimento;

            if(input.comprimentoAlturaIguais1)
            {
                altura      = 1.0;
                comprimento = 1.0;
            }

            vetVetItens[i].push_back(altura);
            vetVetItens[i].push_back(comprimento);
            instanciaG.vetItemCliente[nextItem] = node;

            instanciaG.matCliItensIniFim.get(node, 1) = nextItem;

            nextItem += 1;
        }

        //std::cout << vetVetItens[i] << "\n";

        for(int j = 0; j < numItensPorClie*2; j+=2)
        {
            /*
            instanciaG.vetItemAltura.push_back(vetVetItens[i][j]);
            instanciaG.vetItemLargura.push_back(vetVetItens[i][j + 1]);
            instanciaG.vetItemArea.push_back(vetVetItens[i][j] * vetVetItens[i][j + 1]);
            */

            instanciaG.vetItens.push_back(Item(vetVetItens[i][j], vetVetItens[i][j+1], 0.0, instanciaG.vetDemandaCliente[node]/numItensPorClie));

        }
    }

    instanciaG.maxNumItensPorClie = maxNumItensPorCli;
    instanciaG.atualizaVetMinDimItens();

    file.close();

}

void InstanceNS::readOroloc3D(const std::string &strFile)
{
    //std::printf("FILE: %s\n", strFile.c_str());


    int numClientes = 0, numVeiculos = 0, numItens = 0, numArcs = 0;
    double maxPayload = 0.0;

    std::ifstream file(strFile);
    assertm(!file.is_open(), "Cant open the file: : "<<strFile);

    std::string trash;
    std::getline(file, trash);	// #Date: 10/02/26 16:22:09
    std::getline(file, trash);  // #Instance Name: S20_288T


    file>>numClientes;
    std::getline(file, trash); // # Customers + deposit
    std::getline(file, trash); // #Trailer size (length, width, height) in millimeters
    //std::printf("Num of Customrs: %d\n", numClientes);

    Array<double, 3> dimTruck;// = instanciaG.vetDimVeiculo;
    file>>dimTruck[0] >> dimTruck[1] >> dimTruck[2];

    //std::cout<<"DimTruck: "<<dimTruck<<"\n";

    file>>numVeiculos;
    std::getline(file, trash);  // #Number of Trucks

    //std::cout<<dimTruck<<"\n";
    //std::cout<<"numTrucks: "<<numVeiculos<<"\n";

    file>>maxPayload;
    instanciaG.maxPayload = maxPayload;
    //std::printf("MaxPayload: %.1f\n", maxPayload);
    std::getline(file, trash);  // #Maximum truck payload

    file>>numItens;
    std::getline(file, trash);  // #Number of Itms

    //std::printf("numItems: %d\n", numItens);

    instanciaG = Instance(numClientes, numItens, numVeiculos);
    instanciaG.matTempo.setVal(1.0);
    instanciaG.numDim = 3;
    instanciaG.nome = strFile;
    instanciaG.maxPayload = maxPayload;
    instanciaG.vetDimVeiculo = dimTruck;

    std::getline(file, trash);  // #Distances - Simetric

    file>>numArcs;
    std::getline(file, trash);  //  #Number of arcs

    //std::printf("numArcs: %d\n", numArcs);

    //std::printf("Trash: %s\n", trash.c_str());

    for(int i=0; i < numArcs; ++i)
    {
        int custI = 0, custJ = 0;
        double dist = 0;

        file>>custI>>custJ>>dist;
        //std::printf("%d %d %.1f\n", custI, custJ, dist);
        instanciaG.matDist.get(custI, custJ) = dist;
        instanciaG.matDist.get(custJ, custI) = dist;
    }

    //std::cout<<instanciaG.matDist;
    std::getline(file, trash);  //  #Itens
    std::getline(file, trash);  //  ??
    //std::printf("Trash: %s\n", trash.c_str());

    instanciaG.vetItens.resize(numItens);
    int nextItem = 0;
    int maxNumItensPorCli = -1;

    for(int i=0; i < numClientes-1; ++i)
    {
        int cust=0, custNumItems=0;
        file>>cust>>custNumItems;
        maxNumItensPorCli = std::max(maxNumItensPorCli, custNumItems);

        //std::cout<<cust<<" "<<custNumItems<<"\n";
        instanciaG.matCliItensIniFim.get(cust, 0) = nextItem;

        //std::cout<<cust<<"\n";

        for(int j=0; j < custNumItems; ++j)
        {
            Item& item = instanciaG.vetItens[nextItem];
            file>>item.vetDim[0]>>item.vetDim[1]>>item.vetDim[2]>>item.weight;
            item.volume = item.vetDim[0] * item.vetDim[2] * item.vetDim[2];
            item.fragility = false;
            item.customer = cust;

            //std::cout<<item.print()<<"; ";

            instanciaG.matCliItensIniFim.get(cust, 1) = nextItem;
            instanciaG.vetItemCliente[nextItem] = cust;

            nextItem += 1;


        }
        //std::printf("\n");
        //break;
    }

    instanciaG.vetOrderId.resize(numItens);
    std::getline(file, trash); //
    std::getline(file, trash); // #Order Id

    for(int i=0; i < numItens; ++i)
    {
        file>>instanciaG.vetOrderId[i];
    }

    //std::cout<<instanciaG.vetOrderId<<"\n";

    for(int i=0; i < numItens; ++i)
    {
        int orderId = instanciaG.vetOrderId[i];
        if(instanciaG.mapOrderIdItem.count(orderId) == 0)
        {	VectorI vet;
            vet.push_back(i);
            instanciaG.mapOrderIdItem[orderId] = vet;
        }
        else
        {
            VectorI& vet = instanciaG.mapOrderIdItem[orderId];
            vet.push_back(i);
        }

        instanciaG.mapOrderIdCust[orderId] = instanciaG.vetItens[i].customer;
    }

    /*
    for(auto it:instanciaG.mapOrderIdItem)
    {
        std::cout<<it.first<<": "<<it.second<<"\n";
    }

    std::cout<<"OrderId  customer\n";
    for(auto it:instanciaG.mapOrderIdCust)
        std::cout<<it.first<<" "<<it.second<<"\n";
    */

    file.close();


    instanciaG.maxNumItensPorClie = maxNumItensPorCli;
    instanciaG.atualizaVetMinDimItens();

    int itemId = 0;
    /*
    for(Item& item:instanciaG.vetItens)
    {
        std::cout<<itemId<<": "<<item.customer<<"\n";
        itemId += 1;
    }
    */

    //PRINT_DEBUGG("", "")
    //exit(-1);
}

void InstanceNS::read3dInstance(const std::string &strFile)
{

    int numClientes, numVeiculos, numItens;

    std::ifstream file(strFile);
    assertm(!file.is_open(), "Cant open the file: : "<<strFile);

    //std::cout<<"File("<<strFile<<") aberto\n";

    std::string lineLixo;
    int num = 2;

    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numClientes;
    numClientes += 1;

    for(int i=0; i < 1; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numVeiculos;

    for(int i=0; i < 1; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    file>>numItens;
    //numItens += 1;

    std::cout<<"numClientes: "<<numClientes<<"; numItens: "<<numItens<<" numVeic: "<<numVeiculos<<"\n";


    {
        for(int i = 0; i < 2; ++i)
        {
            getline(file, lineLixo);
            //std::cout << "line: " << lineLixo << "\n";
        }
    }

    instanciaG = Instance(numClientes, numItens, numVeiculos);
    instanciaG.matTempo.setVal(1.0);
    instanciaG.numDim = 3;
    instanciaG.nome = strFile;

    file >> instanciaG.maxPayload >> instanciaG.vetDimVeiculo[2] >> instanciaG.vetDimVeiculo[1]
         >> instanciaG.vetDimVeiculo[0];

    std::cout << "veicCap: " << instanciaG.maxPayload << "; veicComprimento: " << instanciaG.vetDimVeiculo<<"\n\n";

    num = 2;

    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    Matrix<double> matPontos(numClientes, 2);
    matPontos.setVal(-1.0);

    for(int i=0; i < numClientes; ++i)
    {
        double x, y, cli;
        file >> cli >> x >> y >> instanciaG.vetDemandaCliente[i];
        matPontos.get(i, 0) = x;
        matPontos.get(i, 1) = y;

        //std::cout<<cli<<" "<<x<<" "<<y<<" "<<instanciaG.vetDemandaCliente[i]<<"\n";

    }

    for(int i=0; i < numClientes; ++i)
    {
        instanciaG.matDist.get(i, i) = 0.0;

        for(int j=i+1; j < numClientes; ++j)
        {
            double dist = std::sqrt(std::pow(matPontos(i, 0) - matPontos(j, 0), 2) +
                                    std::pow(matPontos(i, 1) - matPontos(j, 1), 2));

            instanciaG.matDist.get(i, j) = dist;
            instanciaG.matDist.get(j, i) = dist;
        }
    }

    //std::cout << instanciaG.matDist << "\n";

    num = 2;
    for(int i=0; i < num; ++i)
    {
        getline(file, lineLixo);
        //std::cout << "line: " << lineLixo << "\n";
    }

    //Vector<Vector<double>> vetVetItens(numClientes);
    int maxNumItensPorCli = 0;
    int nextItem = 0;


    //instanciaG.vetItemAltura.reserve(instanciaG.numItens);
    //instanciaG.vetItemLargura.reserve(instanciaG.numItens);
    //instanciaG.vetItemArea.reserve(instanciaG.numItens);

    instanciaG.vetPesoItens.reserve(instanciaG.numItens);
    instanciaG.vetItens.reserve(instanciaG.numItens);

    num = numClientes;
    if(input.splitInstancia)
        num -= 1;
    for(int i=0; i < num; ++i)
    {
        int node, numItensPorClie;

        file >> node >> numItensPorClie;
        //std::cout<<node<<" "<<numItensPorClie<<"\n";

        maxNumItensPorCli = std::max(maxNumItensPorCli, numItensPorClie);

        if(node != 0)
        {
            instanciaG.matCliItensIniFim.get(node, 0) = nextItem;
            //instancia.matCliItensIniFim.get(node, 1) = nextItem+(numItensPorClie-1);
        }

        double wight = instanciaG.vetDemandaCliente[node]/numItensPorClie;

        for(int item = 0; item < numItensPorClie; ++item)
        {
            int fragility;
            double altura, largura, comprimento;
            file >> altura >> largura >>comprimento >> fragility;
            //std::cout<<"\t"<<altura<<", "<<largura<<", "<<comprimento<<"\n";

            if(input.comprimentoAlturaIguais1)
            {
                altura      = 1.0;
                comprimento = 1.0;
            }

            /*
            vetVetItens[i].push_back(altura);
            vetVetItens[i].push_back(largura);
            vetVetItens[i].push_back(comprimento);
            */

            instanciaG.vetItemCliente[nextItem] = node;

            instanciaG.matCliItensIniFim.get(node, 1) = nextItem;
            instanciaG.vetItens.push_back(Item(largura, comprimento, altura, wight));
            instanciaG.vetItens[instanciaG.vetItens.size()-1].fragility = fragility;
            instanciaG.vetItens[instanciaG.vetItens.size()-1].customer = node;
            nextItem += 1;


        }

        //std::cout << vetVetItens[i] << "\n";
    }

    instanciaG.maxNumItensPorClie = maxNumItensPorCli;
    instanciaG.atualizaVetMinDimItens();

    file.close();

}

void InstanceNS::Instance::atualizaVetMinDimItens()
{
    for(int i=0; i < numItens; ++i)
    {
        auto tuple = getMinArray(vetItens[i].vetDim, numDim);
        vetMinDimItens[i] = get<1>(tuple);
    }

    //std::cout<<"vetMinDimItens: "<<vetMinDimItens<<"\n\n";
}

std::string InstanceNS::Item::print(bool printVol)
{
    std::string str = "(";
    for(int i=0; i < instanciaG.numDim; ++i)
    {
        str += std::format("{:.1f}", vetDim[i]);
        if(i < (instanciaG.numDim-1))
            str += ",";
    }

    if(printVol)
        str += "; " + std::to_string(volume);

    str += "; " + std::format("{:.1f}", weight);

    str += ")";
    return str;
}


/*std::string InstanceNS::printItem(int itemId)
{
    return std::to_string(instanciaG.vetItemLargura[itemId])+","+std::to_string(instanciaG.vetItemAltura[itemId]);
}*/

InstanceNS::Item::Item(double x, double y, double z, double peso_)
{
    vetDim[0] = x;
    vetDim[1] = y;
    vetDim[2] = z;
    weight = peso_;
    volume = 1.0;

    for(int d=0; d < 3; ++d)
    {
        volume *= vetDim[d];
        if((d+1) == instanciaG.numDim)
            break;
    }
}

void InstanceNS::Item::set(double x, double y, double z)
{
    vetDim[0] = x;
    vetDim[1] = y;
    vetDim[2] = z;

    volume = 1.0;

    for(int d=0; d < instanciaG.numDim; ++d)
    {
        volume *= vetDim[d];
    }
}

double InstanceNS::Item::getDimRotacionada(int d, Rotation r)
{
    assertm(static_cast<int>(r) >= instanciaG.numRotation, "Error r("<<r<<")");

    static const int perm[6][3] = {
        {0,1,2}, // LWH
        {1,0,2}, // WLH
        {0,2,1}, // LHW
        {1,2,0}, // WHL
        {2,0,1}, // HLW
        {2,1,0}  // HWL
    };

    return vetDim[perm[static_cast<int>(r)][d]];
}

int InstanceNS::copiaItensCliente(int cliente, VectorI &vetItens)
{
    int numItens = 0;

    for(int i=instanciaG.matCliItensIniFim(cliente, 0); i <= instanciaG.matCliItensIniFim(cliente, 1); ++i)
    {
        vetItens[numItens] = i;
        numItens += 1;
    }

    return numItens;
}


int InstanceNS::copiaItensClientes(VectorI& vetClientes, int tam, VectorI& vetItens, bool push)
{
    int tamVetItens = 0;
    for(int i=0; i < tam; ++i)
    {
        int cliente = vetClientes[i];
        const int ini = instanciaG.matCliItensIniFim(cliente, 0);
        const int fim = instanciaG.matCliItensIniFim(cliente, 1);

        //std::printf("Ini: %d, Fim: %d\n", ini, fim);

        for(int j=ini; j <= fim; ++j)
        {
            if(!push)
                vetItens[tamVetItens] = j;
            else
                vetItens.push_back(j);
            tamVetItens += 1;
        }
    }

    return tamVetItens;
}

double InstanceNS::calculaDistancia(VectorI& vet, int tam)
{
    double dist = 0.0;
    if(tam <= 1)
        return 0.0;

    for(int i=0; i < (tam-1); ++i)
        dist += instanciaG.matDist(vet[i], vet[i+1]);

    return dist;
}
