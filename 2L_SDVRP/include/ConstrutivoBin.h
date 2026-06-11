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
    int                  epId   = 0;
    InstanceNS::Rotation r      = InstanceNS::Rot0;
    SolucaoNS::Ponto     point;
    double               maxDif = 0.0;

    EpRot() = default;

    INLINE
    bool operator<(const EpRot &outro) const
    {
        if(doubleEqual(point.vetDim[2], outro.point.vetDim[2], 1E-5))
        {
            /*
            if(doubleEqual(point.vetDim[0], outro.point.vetDim[0], 1E-5))
                return (point.vetDim[1] < outro.point.vetDim[1]);
            else
                return point.vetDim[0] < outro.point.vetDim[0];
            */
            return maxDif < outro.maxDif;
        }
        else
            return point.vetDim[2] < outro.point.vetDim[2];

    }
};

struct ItemRandom
{
    int itemId    = -1;
    int randomKey =  0;

    //ItemRandom(int itemId_, int randomKey_):itemId(itemId_), randomKey(randomKey_){}

    INLINE
    bool operator<(const ItemRandom& other)const{return randomKey < other.randomKey;}
};

bool canInsert(const SolucaoNS::Ponto &ep,
               const int               itemId,
               const SolucaoNS::Bin   &bin,
               InstanceNS::Rotation    r,
               double                 &maxDif,
               double				   wightLimit);

bool epColideItem(const SolucaoNS::Ponto &ep,
                  const SolucaoNS::Ponto &ponto,
                  const int               itemId);

int construtivoBinPacking(Vector<SolucaoNS::Bin> &vetBin,
                          const int               vetBinTam,
                          const VectorI          &vetItensC,
                          const int               vetItensTam,
                          const double            alpha);

bool construtivoBinPacking(SolucaoNS::Bin  &bin,
                           VectorI         &vetItens,
                           const int        vetItensTam,
                           const double     alpha,
                           const int        numRepeticoes,
                           SolucaoNS::Rota *rota = nullptr);

double computeXY_Overlap(InstanceNS::Item       &item0,
                         InstanceNS::Rotation    r0,
                         const SolucaoNS::Ponto &p0,
                         InstanceNS::Item       &item1,
                         InstanceNS::Rotation    r1,
                         const SolucaoNS::Ponto &p1);

void sortVetItemsByCustomer(VectorI& vetItems, int size);

} // namespace ConstrutivoBinNS

#endif // INC_2L_SDVRP_CONSTRUTIVOBIN_H
