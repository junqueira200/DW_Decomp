/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: Constates.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_CONSTANTES_H
#define INC_2L_SDVRP_CONSTANTES_H

#include <limits>

constexpr int    NumItensPorBin    = 300;
constexpr int    NumEpPorBin       = 5*NumItensPorBin;
constexpr double INF_Double        = HUGE_VAL;
constexpr int    INF_Int           = std::numeric_limits<int>::infinity();
constexpr int    TamRota           = 60;
constexpr double DiffPermetidaDist = 1E-5;
constexpr double Gravity           = 10.0;//9.81;

constexpr bool PrintEP    = false;
constexpr bool PrintConst = false;

#endif //INC_2L_SDVRP_CONSTANTES_H
