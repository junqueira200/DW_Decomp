#ifndef C_API_H
#define C_API_H

#include "Instancia.h"
#include "InputOutput.h"

using namespace ParseInputNS;


extern "C"
{
    void ini_3D_Packing(char* strInst_c);
    bool testRoute(int* vet_c, int vetSize);


}

#endif // C_API_H
