//
// Created by igor on 15/11/24.
//

#ifndef DW_INSTANCIA_H
#define DW_INSTANCIA_H

#include <Eigen/Eigen>
#include <string>
#include "Aux.h"

namespace InstanceVRPTW_NS
{

    struct ClieTime
    {
        double readyTime = 0.0;
        double dueTime = 0.0;
        double servTime = 0.0;
    };

    class InstanceVRPTW
    {
    public:

        int numClientes = 0;
        int numVeic = 0;
        int capVeic = 0;

        EigenMatrixRowD matDist;
        Eigen::VectorX<ClieTime> vetClieTime;
        Eigen::VectorXi vetClieDem;
        std::string instName;
        bool subInstancia = false;

        explicit InstanceVRPTW(int numClie);
        InstanceVRPTW()=default;
        double sumDist();
        int sumDem();

    };

    void leInstanciaAugerat(const std::string &strFile, InstanceVRPTW &instVrpTw);
    void leInstanciaSalomon(const std::string &strFile, InstanceVRPTW &instCvrp);
    double calculateDistance(double x1, double y1, double x2, double y2);
    double somaDist(const InstanceVRPTW &instVrpTw);

    void getSubInstancia(int numClientes, InstanceVRPTW& instVrpTw);

    enum DataIndex
    {
        CustNo = 0,
        XCoord,
        YCoord,
        Dem,
        ReadyTime,
        DueTime,
        ServTime
    };

    inline InstanceVRPTW* ptr_instVrpG = nullptr;
}
#endif //DW_INSTANCIA_H
