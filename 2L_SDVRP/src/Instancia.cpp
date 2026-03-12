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
#include "rand.h"
#include <filesystem>

using namespace InstanceNS;
using namespace ParseInputNS;
using namespace RandNs;

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
            file>>item.vetDim[0]>>item.vetDim[1]>>item.vetDim[2]>>item.weight>>item.oroloc3D_item_id;
            item.volume = item.vetDim[0] * item.vetDim[2] * item.vetDim[2];
            item.fragility = false;
            item.customer = cust;
            item.weightForce = item.weight*Gravity;

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

    int item_id = -1;


    for(int i=0; i < numItens; ++i)
    {
        file>>instanciaG.vetOrderId[i]>>item_id;
        instanciaG.mapItem_IdItem[item_id] = i;
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
    trash = "";
    std::getline(file, trash);
    trash = "";
    std::getline(file, trash);

    std::printf("Trash: %s\n", trash.c_str());

    int customer_id, customer;
    for(int i=0; i < instanciaG.numClientes; ++i)
    {
        file>>customer_id>>customer;

        if(instanciaG.mapCustomerToCustomer_id.contains(customer) ||
           instanciaG.mapCustomer_idToCustomer.contains(customer_id))
        {
            std::printf("Error while creating mapCustomerToCustomer_id and mapCustomer_idToCustomer to %d and %d\n",
                         customer, customer_id);
            PRINT_DEBUGG("", "");
            throw "ERROR";

            exit(EXIT_FAILURE);
        }

        instanciaG.mapCustomerToCustomer_id[customer] = customer_id;
        instanciaG.mapCustomer_idToCustomer[customer_id] = customer;
    }

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
    weightForce = weight*Gravity;
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


int InstanceNS::generateRandomListOfItems(int numItens, VectorI& vetItems)
{
    Vector<int8_t> vetItensSelecionados(instanciaG.numItens);
    vetItensSelecionados.setAll((int8_t)0);

    vetItems = VectorI();

    double volumeOcupado = 0.0;
    double volumeVeiculo = 1.0;
    double demanda       = 0.0;

    volumeVeiculo = instanciaG.vetDimVeiculo[0]*instanciaG.vetDimVeiculo[1]*instanciaG.vetDimVeiculo[2];

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

        vetItems.push_back(itemId);

    }

    return numItens;

}


void InstanceNS::readOroloc3D2(const std::string &strFile)
{
    std::printf("FILE: %s\n", strFile.c_str());


    int numClientes = 0, numVeiculos = 0, numItens = 0, numArcs = 0;
    double maxPayload = 0.0;
    Array<double, 3> veicDim;
    std::string trash;

    std::ifstream file(strFile);
    assertm(!file.is_open(), "Cant open the file: : "<<strFile);

    std::getline(file, trash);
    file>>trash>>numClientes;
    numClientes += 1;

    file>>trash>>numItens;
    std::getline(file, trash);
    std::getline(file, trash);
    file>>trash>>numVeiculos;

    std::getline(file, trash);
    std::getline(file, trash);
    std::getline(file, trash);
    std::getline(file, trash);

    file>>trash>>maxPayload;
    file>>trash>>veicDim[0];
    file>>trash>>veicDim[1];
    file>>trash>>veicDim[2];

    //std::printf("maxPayload: %f\n", maxPayload);
    //std::cout<<"Veic dim: "<<veicDim<<"\n";



    for(int i=0; i < 13; ++i)
    {
        std::getline(file, trash);
        //std::printf("Trash: %s\n", trash.c_str());
    }

    std::getline(file, trash);   // CUSTOMERS
    std::getline(file, trash);   // i	x	y	Demand	ReadyTime	DueDate	ServiceTime	DemandedMass	DemandedVolume

    /*
    std::printf("numClientes: %d\n", numClientes);
    std::printf("numItems: %d\n", numItens);
    std::printf("numVeic: %d\n", numVeiculos);
    std::printf("Trash: %s\n", trash.c_str());
    */

    instanciaG =  Instance(numClientes, numItens, numVeiculos);
    std::filesystem::path path(strFile);
    instanciaG.numDim = 3;
    instanciaG.maxPayload = maxPayload;
    instanciaG.vetDimVeiculo = veicDim;
    instanciaG.nome = path.filename();

    Matrix<double> matCoord(numClientes, 2);

    int cust = -1;

    for(int i=0; i < numClientes; ++i)
    {
        file>>cust>>matCoord.get(i, 0)>>matCoord.get(i, 1);
        std::getline(file, trash);
    }

    for(int i=0; i < numClientes; ++i)
    {
        instanciaG.matDist.get(i, i) = INF_Double;
        for(int j=i+1; j < numClientes; ++j)
        {
            double dist = std::sqrt(std::pow(matCoord(i, 0) - matCoord(j, 0), 2) +
                                    std::pow(matCoord(i, 1) - matCoord(j, 1), 2));

            instanciaG.matDist.get(i, j) = dist;
            instanciaG.matDist.get(j, i) = dist;
        }
    }

    //std::cout<<instanciaG.matDist<<"\n\n";

    std::getline(file, trash);  // Empty line
    std::getline(file, trash);  // Type	Length	Width	Height	Mass	Fragility	LoadBearingStrength
    std::getline(file, trash);

    Array<int, 4> vetDimMass;
    std::map<int, Array<int, 4>> mapItem_id_to_ItemDimMass;
    std::string bt;
    int btInt;

    for(int i=0; i < numItens; ++i)
    {
        file>>bt>>vetDimMass[0]>>vetDimMass[1]>>vetDimMass[2]>>vetDimMass[3];
        std::getline(file, trash);

        bt.erase(0, 2);
        //std::printf("%d\n", std::stoi(bt));
        btInt = std::stoi(bt);

        //std::cout<<btInt<<" "<<vetDimMass<<"\n";

        mapItem_id_to_ItemDimMass[btInt] = vetDimMass;
    }

    std::getline(file, trash);   // new line
    std::getline(file, trash);   // DEMANDS PER CUSTOMER
    std::getline(file, trash);   // i	Type Quantity

    std::string line;

    int nextItem = 0;
    int maxItemsPerCust = 0;

    while(std::getline(file, line))
    {
        std::stringstream ss(line);

        int i;
        ss >> i;

        //std::cout << "i = " << i << std::endl;

        std::string type;
        int quantity;
        int numItems = 0;
        instanciaG.matCliItensIniFim.get(i, 0) = nextItem;
        instanciaG.vetDemandaCliente[i] = 0;

        while(ss>>type>>quantity)
        {

            // remove "Bt"
            numItems += 1;
            int id = std::stoi(type.substr(2));

            //std::cout << "Type: " << id
            //          << " Quantity: " << quantity << std::endl;
            Array<int, 4> &arrayDimMass = mapItem_id_to_ItemDimMass[id];

            for(int k=0; k < quantity; ++k)
            {
                instanciaG.matCliItensIniFim.get(k, 1) = nextItem;

                Item item;// = instanciaG.vetItens[nextItem];
                //std::cout<<arrayDimMass<<"\n";
                item.oroloc3D_item_id = id;
                item.set((double)arrayDimMass[0], (double)arrayDimMass[1], (double)arrayDimMass[2]);
                item.weight = arrayDimMass[3];
                //std::cout<<item.weight<<"\n";
                item.weightForce = item.weight*Gravity;
                item.customer = i;

                instanciaG.vetItens.push_back(item);
                instanciaG.mapItem_IdItem[id] = nextItem;
                instanciaG.vetPesoItens[nextItem] = arrayDimMass[3];
                instanciaG.vetDemandaCliente[i] += arrayDimMass[3];

                nextItem += 1;
            }
        }

        maxItemsPerCust = std::max(maxItemsPerCust, numItems);

        //std::printf("\n");

    }

    instanciaG.maxNumItensPorClie = maxItemsPerCust;

    file.close();
    //EXIT_PRINT();

}
