//
// Created by igor on 15/11/24.
//

#ifndef DW_INSTANCIA_H
#define DW_INSTANCIA_H

#include <Eigen/Eigen>
#include <string>

namespace InstanciaNS
{

    struct ClieTime
    {
        double readyTime = 0.0;
        double dueTime = 0.0;
        double servTime = 0.0;
    };

    class InstVRP_TW
    {
    public:

        int numClientes = 0;
        int numVeic = 0;
        int capVeic = 0;

        Eigen::MatrixX<double> matDist;
        Eigen::VectorX<ClieTime> vetClieTime;
        Eigen::VectorXi vetClieDem;
        std::string instName;

        explicit InstVRP_TW(int numClie);
        InstVRP_TW()=default;

    };

    void leInstanciaSolomon(const std::string &strFile, InstVRP_TW &instVrpTw);
    void leInstanciaAugerat(const std::string &strFile, InstVRP_TW &instCvrp);
    double calculateDistance(double x1, double y1, double x2, double y2);
    double somaDist(const InstVRP_TW &instVrpTw);


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

}
#endif //DW_INSTANCIA_H
