/* ****************************************
 * ****************************************
 *  Data:    06/11/24
 *  Arquivo: ConstrutivoBin.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_CONSTRUTIVOBIN_H
#define INC_2L_SDVRP_CONSTRUTIVOBIN_H

#include "Instancia.h"
#include "Solucao.h"

namespace ConstrutivoBinNS
{

    struct EpRot
    {
        int epId = 0;
        InstanceNS::Rotation r = InstanceNS::Rot0;
        double atributo = 0.0;

        EpRot()=default;

        bool operator<(const EpRot &outro) const
        {
            return atributo < outro.atributo;
        }
    };

    bool canInsert(const SolucaoNS::Ponto &ep, const int itemId, const SolucaoNS::Bin &bin, InstanceNS::Rotation r);
    bool epColideItem(const SolucaoNS::Ponto &ep, const SolucaoNS::Ponto &ponto, const int itemId);

    int construtivoBinPacking(Vector<SolucaoNS::Bin> &vetBin,
                               const int vetBinTam,
                               const VectorI &vetItensC,
                               const int vetItensTam,
                               const double alpha);

    bool construtivoBinPacking(SolucaoNS::Bin &bin,
                               VectorI &vetItens,
                               const int vetItensTam,
                               const double alpha,
                               const int numRepeticoes);

    double computeXY_Overlap(InstanceNS::Item& item0, InstanceNS::Rotation r0, const SolucaoNS::Ponto& p0,
                             InstanceNS::Item& item1, InstanceNS::Rotation r1, const SolucaoNS::Ponto& p1);

}

#endif //INC_2L_SDVRP_CONSTRUTIVOBIN_H
