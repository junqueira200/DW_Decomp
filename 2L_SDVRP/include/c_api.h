#ifndef C_API_H
#define C_API_H

#include "Instancia.h"
#include "InputOutput.h"

using namespace ParseInputNS;


extern "C"
{
    void ini_3D_Packing(char* strInst_c, int oroloc3D);
    int testRoute(int* vet_c, int vetSize, int onlyHeuristic=0, int doInverseRoute=0);
    void setClassical3DPackingProblem();
    void setOroloc3DProblem();

    int getNumberOfCustoms();
    int getNumberOfTrucks();
    int getDemandFromCustomr(int i);
    double getVolumeFromCustomr(int i);
    double getDistance(int i, int j);
    int getVehicleCapacity();
    double getVehicleVolume();

    // TODO: Criar uma estrutura para as rotas;
    // TODO: Criar uma funcao para recuperar as solucoes do binpacking
}

#endif // C_API_H
