/* ****************************************
 * ****************************************
 *  Data:    06/05/26
 *  Arquivo: AuxGeometry.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/


#ifndef AUXGEOMETRY_H
#define AUXGEOMETRY_H

#include "Instancia.h"
#include "AuxT.h"

namespace AuxGeometryNS
{

    INLINE
    double collision1D(double a0, double a1, double b0, double b1)
    {
        return std::max(0.0, (std::min(a1, b1) - std::max(a0, b0)));
    }

}

#endif // AUXGEOMETRY_H
