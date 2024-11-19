//
// Created by igor on 15/11/24.
//
#include <iostream>
#include "Instancia.h"
#include <fstream>

InstanciaNS::InstVRP_TW::InstVRP_TW(int numClie_)
{
    numClientes = numClie_;
//    numVeic     = numVeic_;

    matDist     = Eigen::MatrixXd(numClientes, numClientes);
    vetClieTime = Eigen::VectorX<ClieTime>(numClientes);
    vetClieDem  = Eigen::VectorXi(numClientes);
}


void InstanciaNS::leInstanciaSalomon(const std::string &strFile, InstVRP_TW &instVrpTw)
{

    std::ifstream file(strFile);
    if(!file.is_open())
    {
        std::cout<< "Nao foi possivel abrir o arquivo: "<<strFile<<"\n\n";
        throw "ERRO";
    }

    std::string instName;
    file>>instName;

    std::cout<<"instName: "<<instName<<"\n";

    std::string strLixo;

    for(int i=0; i < 4; ++i)
    {
        getline(file, strLixo);
        std::cout<<i<<": "<<strLixo;
        strLixo = "";
    }

    int numClie, cap;
    file>>numClie>>cap;
    numClie += 1;

    std::cout<<"numeClie("<<numClie<<"); cap("<<cap<<")\n\n";

    instVrpTw = InstVRP_TW(numClie);
    instVrpTw.instName = instName;
    instVrpTw.capVeic = cap;

    Eigen::MatrixXi matData(numClie, 7);


    for(int i=0; i < 5; ++i)
        getline(file, strLixo);


    for(int i=0; i < numClie; ++i)
    {
        for(int j=0; j < 7; ++j)
            file>>matData(i,j);
    }

    std::cout<<matData<<"\n\n";

    int totalDem = 0;

    for(int i=0; i < numClie; ++i)
    {
        for(int j=i+1; j < numClie; ++j)
        {

            double distI_J = calculateDistance(matData(i, XCoord), matData(i,XCoord),
                                               matData(j, XCoord), matData(j, YCoord));

            instVrpTw.matDist(i,j) = distI_J;
            instVrpTw.matDist(j,i) = distI_J;

        }

        instVrpTw.vetClieTime[i].readyTime = matData(i, ReadyTime);
        instVrpTw.vetClieTime[i].dueTime   = matData(i, DueTime);
        instVrpTw.vetClieTime[i].servTime  = matData(i, ServTime);

        instVrpTw.vetClieDem[i] = matData(i, Dem);
        totalDem += instVrpTw.vetClieDem[i];

        std::cout<<i<<": ["<<instVrpTw.vetClieTime[i].readyTime<<";"<<instVrpTw.vetClieTime[i].dueTime<<
                      "]; servTime("<<instVrpTw.vetClieTime[i].servTime<<"); dem("<<instVrpTw.vetClieDem[i]<<")\n";
    }

    instVrpTw.numVeic = std::ceil(double(totalDem)/instVrpTw.capVeic);

    std::cout<<"numVeic: "<<instVrpTw.numVeic<<"\n\n";

    std::cout<<"\n"<<instVrpTw.matDist<<"\n";

}


void InstanciaNS::leInstanciaAugerat(const std::string &strFile, InstVRP_TW &instCvrp)
{


    std::ifstream file(strFile);
    if(!file.is_open())
    {
        std::cout<< "Nao foi possivel abrir o arquivo: "<<strFile<<"\n\n";
        throw "ERRO";
    }

    std::string line;
    int numClientes, capVeic;

    for(int i=0;i<3;++i)
        getline(file, line);

    file>>line>>line>>numClientes;
    getline(file, line);
    getline(file, line);

    file>>line>>line>>capVeic;
    getline(file, line);
    getline(file, line);

std::cout<<"numClientes: "<<numClientes<<"; cap: "<<capVeic<<"\n\n";

//    instancia->matrixDistancia.resize(numClientes, numClientes, false);
    instCvrp = InstVRP_TW(numClientes);
    instCvrp.capVeic = capVeic;

    int *vetorX = new int[numClientes];
    int *vetorY = new int[numClientes];

    int aux;

    for(int i=0; i < numClientes; ++i)
    {

        file>>aux>>vetorX[i]>>vetorY[i];

        //cout<<i<<" "<<vetorX[i]<<" "<<vetorY[i]<<"\n";
    }


    getline(file, line);
    getline(file, line);

    for(int i=0; i<numClientes; ++i)
    {
        instCvrp.matDist(i, i) = 0;

        for(int j=i+1; j < numClientes; ++j)
        {
            double temp = (sqrt(pow(vetorX[i] - vetorX[j], 2) + pow(vetorY[i] - vetorY[j], 2)));
            instCvrp.matDist(i, j) = static_cast<int>(std::round(temp));
            instCvrp.matDist(j, i) = instCvrp.matDist(i, j);
        }

    }

    //instancia->demanda = new int[numClientes];
    int demandaTotal = 0;

    for(int i=0; i < numClientes; ++i)
    {
        file >>aux>>instCvrp.vetClieDem[i];
        demandaTotal += instCvrp.vetClieDem[i];
    }

    instCvrp.numVeic = ceil(demandaTotal/double(instCvrp.capVeic));
    std::cout<<instCvrp.numVeic<<"\n";

    file.close();

    delete []vetorX;
    delete []vetorY;

}


double InstanciaNS::calculateDistance(double x1, double y1, double x2, double y2)
{
    // Calculate Euclidean distance
    double distance = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

    // Round down to the first lower decimal
    distance = std::floor(distance * 10) / 10.0;

    return distance;
}

double InstanciaNS::somaDist(const InstVRP_TW &instVrpTw)
{
    double dist = 0.0;

    for(int i=0; i < instVrpTw.numClientes; ++i)
    {
        for(int j=0; j < instVrpTw.numClientes; ++j)
        {
            if(i == j)
                continue;

            dist += instVrpTw.matDist(i, j);
        }
    }

    return dist;
}
