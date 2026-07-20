/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: ParseInput.h
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_INPUTOUTPUT_H
#define INC_2L_SDVRP_INPUTOUTPUT_H

#include "Commit.h"
#include "Solucao.h"
#include "rand.h"
#include "string"
#include <chrono>
#include <iostream>

namespace ParseInputNS
{
class File
{
  public:
    std::string fileSol;      // X
    std::string fileNum;      // X
    std::string fileResulCSV; // X
    std::string fileSolPrint;
    std::string fileSeed; // X
};

class Input
{
  public:
    std::string strInstCompleto;
    std::string strInst;
    bool        splitInstancia 				= false;
    bool        splitVrp 					= false;
    double      aphaBin 					= 0.1; // 0.15
    double      aphaBinEscolhaEp 			= 0.45; // 0.6
    double      alphaVrp 					= 0.45;
    int         numItIG 					= 500;
    double      gapIgReset 					= 0.2; // 0.2
    bool        comprimentoAlturaIguais1 	= false;
    bool        cpSat 						= true;
    double      cpSatTime 					= 60.0*60*2;
    // 0 Forward, 1 Backard, 2 Bidirectional
    int         labelingType 				= 0;
    bool        lifo 						= true;
    bool        mlifo 						= true;
    bool        removeFromShortSide 		= false;
    bool        inst2d 						= false;
    bool        instOroloc3D 				= false;
    bool        instOroloc3D_2 				= false;
    std::string strSolOroloc3D_2;
    std::string strSolOroloc3D_output;
    bool        support                     = true;
    double      minSupportArea 				= 0.75;
    bool        compactness					= true;
    double      minLeftSupportArea 			= 0.2;
    bool        axleWights 					= true;
    int         supportLimit 				= 0;
    double      balancedLoadingD 			= 0.74;
    bool        balancedLoading             = true;
    double 		maxTimePackingHeuristic 	= 10.0;
    bool        fragility				    = true;
    // bool		inst3d                   = true;

    File file;
    // std::string commit = "54f84fa9027eeb1a17566368b19204f726b1e4d0";
};

class Output
{
  public:
    std::string  strMsg = "CP MODEL WITHOUT ROTATION!";
    std::string  data;
    std::string  fileSaida;
    unsigned int semente;
    double       tempoCpu = 0.0;
    double       tempoRelogio = 0.0;

    void setup()
    {
        semente = RandNs::estado_;
        std::time_t result = duration_cast<std::chrono::seconds>(
                                 std::chrono::system_clock::now().time_since_epoch())
                                 .count();
        data = std::string(std::asctime(std::localtime(&result)));
    }
};

inline Input  input;
inline Output output;

void parseInput(int argc, const char *argv[]);

void escreveFileNum(const std::string &fileNum);
void escreveFileResulCsv(const SolucaoNS::Solucao &sol,
                         const std::string        &fileResulCSV,
                         double                    tempoCpu);
void escreveFileSeed(const std::string &fileSeed);
void escreveFileSol(const SolucaoNS::Solucao &sol, const std::string &fileSol);
void escreveFileSolPrint(const SolucaoNS::Solucao &sol, const std::string &fileSol);
void escreveSaidas(const SolucaoNS::Solucao &best, double tempoCpu);

std::string getNomeInstancia(std::string str);

} // namespace ParseInputNS

#endif // INC_2L_SDVRP_INPUTOUTPUT_H
