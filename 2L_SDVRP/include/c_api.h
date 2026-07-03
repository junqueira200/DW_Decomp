#ifndef C_API_H
#define C_API_H

#include "Instancia.h"
#include "InputOutput.h"

using namespace ParseInputNS;


extern "C"
{
    void ini_3D_Packing(char* strInst_c, int oroloc3D);
    int testRoute(int* vet_c, int vetSize);

    int getNumberOfCustoms();
    int getNumberOfTrucks();
    int getDemandFromCustomr(int i);
    double getVolumeFromCustomr(int i);
    double getDistance(int i, int j);
    int getVehicleCapacity();
    double getVehicleVolume();

}

#endif // C_API_H
