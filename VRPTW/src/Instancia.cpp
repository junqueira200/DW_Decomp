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

void InstanciaNS::leInstancia(const std::string &strFile, InstVRP_TW &instVrpTw)
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

        std::cout<<i<<": ["<<instVrpTw.vetClieTime[i].readyTime<<";"<<instVrpTw.vetClieTime[i].dueTime<<
                      "]; servTime("<<instVrpTw.vetClieTime[i].servTime<<"); dem("<<instVrpTw.vetClieDem[i]<<")\n";
    }

    std::cout<<"\n"<<instVrpTw.matDist<<"\n";

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
