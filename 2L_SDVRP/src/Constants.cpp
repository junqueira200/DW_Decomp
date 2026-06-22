#include "Constants.h"
#include "Instancia.h"

using namespace InstanceNS;

void startConstGlobalVaribles()
{
    DistCompactnessFront = PercentDistCompactnessFront * instanciaG.vetDimVeiculo[0];
    //std::printf("DistCompactnessFront: %.2f\n", DistCompactnessFront);
}
