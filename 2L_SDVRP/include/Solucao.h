/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: Solucao.h
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_SOLUCAO_H
#define INC_2L_SDVRP_SOLUCAO_H

#include "AuxT.h"
#include "Instancia.h"
#include "safe_matrix.h"
#include "safe_vector.h"
#include <Eigen/Eigen>
// #include "InputOutput.h"

namespace SolucaoNS
{

struct Ponto
{
    Array<double, 3> vetDim;

    Ponto() { vetDim.setAll(0.0); }
    Ponto(double d0, double d1, double d2)
    {
        vetDim[0] = d0;
        vetDim[1] = d1;
        vetDim[2] = d2;
    }
    void set(double d0 = 0.0, double d1 = 0.0, double d2 = 0.0)
    {
        vetDim[0] = d0;
        vetDim[1] = d1;
        vetDim[2] = d2;
    }

    std::string print() const;
};

struct PontoRot
{
    Ponto                ponto;
    InstanceNS::Rotation rot;

    PontoRot()
    {
        ponto.set();
        rot = InstanceNS::Rot0;
    }
    PontoRot(double px, double py, double pz, InstanceNS::Rotation r)
    {
        ponto.vetDim[0] = px;
        ponto.vetDim[1] = py;
        ponto.vetDim[2] = pz;
        rot = r;
    }
    void set(double               px = 0.0,
             double               py = 0.0,
             double               pz = 0.0,
             InstanceNS::Rotation r = InstanceNS::Rot0)
    {
        ponto.vetDim[0] = px;
        ponto.vetDim[1] = py;
        ponto.vetDim[2] = pz;
        rot = r;
    }
};

enum Face
{
    Top = 0, // z
    Bottom,  // z
    Front,   // y
    Back,    // y
    Left,    // x
    Right    // x
};

inline std::vector<std::string> vetFaceStr =
    {"Top", "Bottom", "Front", "Back", "Left", "Right"};

class Rota;

struct Bin
{

    // VectorD vetX, vetY;          // Indica aa posicoes do canto inferior esquerdo onde
    // o item eh posicionado no bin VectorD vetEpX, vetEpY;      // Indica as posicoes dos
    // pontos extremos

    Vector<Ponto> vetPosItem; // Indica as posicoes do canto inferior esquerdo onde o item
                              // eh posicionado no bin
    Vector<Ponto>                vetEp;     // Indica as posicoes dos pontos extremos
    VectorI                      vetItemId; // Indica o id do iº item
    Vector<InstanceNS::Rotation> vetRotacao;
    Vector<int8_t>               vetItens; // Indica se o bin empacota o i° item
    Array<double, 3>             binDim;
    double                       volumeTotal = 0.0;
    double                       volumeOcupado = 0.0;
    double                       demandaTotal = 0.0;
    int                          numItens = 0;
    int                          numEps = 0; // Numero de pontos extremos
    double                       sumLeftBalancedLoading = 0.0;
    double                       sumRightBalancedLoading = 0.0;

    void        addItem(int idEp, int idItem, InstanceNS::Rotation r);
    void        addEp(const Ponto &ep);
    std::string printPlot();

    inline __attribute__((always_inline)) void
    setItem(int pos, int itemId, double x, double y, double z)
    {

#if VAR_SOLUTION_BORROW_CHECKER
        if(pos >= numItens)
            throw std::out_of_range("");
#endif

        vetPosItem[pos].vetDim[0] = x;
        vetPosItem[pos].vetDim[1] = y;
        vetPosItem[pos].vetDim[2] = z;
        vetItemId[pos] = itemId;
    };

    inline __attribute__((always_inline)) void
    setEp(int pos, double x, double y, double z)
    {

#if VAR_SOLUTION_BORROW_CHECKER
        if(pos >= numEps)
            throw std::out_of_range("");
#endif

        vetEp[pos].vetDim[0] = x;
        vetEp[pos].vetDim[1] = y;
        vetEp[pos].vetDim[2] = z;
    };

    inline __attribute__((always_inline)) bool vazio() const { return numItens == 0; }

    bool rmI_Item(int i);
    bool rmI_Ep(int i);

    void reset();
    int  getEpComMenorCoord(const VectorI &vetIdEp, int tam);
    bool checkFeasibility(Rota *rota = nullptr, bool fromCp = false, bool print=false);

    void rmItens(const VectorI &vetItensRm, const int tam);

    double getPorcentagemUtilizacao() const;
    void   computeLoadingBalancing();

    Bin();
    Bin(const Bin &bin) = delete;
    Bin &operator=(const Bin &bin) = delete;
};

int  getBinVazio(const Vector<Bin> &vetBin, int tam);
void copiaBin(const Bin &binFonte, Bin &bin);

class Rota
{
  public:
    Rota();
    Rota(const Rota &rota);
    void        reset();
    std::string printRota(bool printDist=true);
    void        computeDistance();

    VectorI        vetRota;
    VectorD        vetTempoSaida;
    VectorD        vetDemClie;
    Vector<int8_t> vetItens; // Indica se o bin empacota o i° item
    int            numPos = 2;
    // double  demTotal        = 0.0;
    double distTotal = 0.0;
    Bin   *binPtr = nullptr;
    size_t hashVal = 0;

    INLINE
    bool operator==(const Rota& r) const
    {
        if(hashVal != r.hashVal)
            return false;

        if(numPos != r.numPos)
            return false;

        for(int i=0; i < numPos; ++i)
        {
            if(vetRota[i] != r.vetRota[i])
                return false;
        }

        return true;
    }

};
#define __PRETTYFILE__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
class HashRoute
{
  public:
    INLINE
    std::size_t operator()(const Rota& r) const
    {
        //std::cout<<"FILE: " << __PRETTYFILE__            \
        //    << "  FUNC: " << __PRETTY_FUNCTION__ << "  LINHA: " << __LINE__ << "\n";

        //std::printf("numPos: %d\n", r.numPos);
        //https://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-16.html
        // ELF Hash algorithms
        size_t valHash = 0;
        for(int i=0; i < r.numPos; ++i)
        {   //         valHash * 16
            valHash = (valHash<<4) + r.vetRota[i];
            uint64_t g = valHash & 0xF0000000L;

            if(g != 0)
                valHash ^= g >> 24;
            valHash &= ~g;
        }

        void* ptr = (void*)&r.hashVal;

        size_t* size_tPtr = (size_t*)ptr;
        *size_tPtr = valHash;

        return valHash;
    }
};

class Solucao
{
  public:
    Vector<Bin>  vetBin;
    Vector<Rota> vetRota;
    double       distTotal = 0.0;

    Solucao();
    explicit Solucao(const InstanceNS::Instance &instancia);

    bool   verificaSol(std::string &error);
    int    getBinVazio();
    void   copiaSolucao(const Solucao &sol);
    void   reset();
    double getUtilizacaoMediaBins() const;
    double getUtilizacaoMedianaBins() const;
    double getTamMedianaRota() const;
    std::string printSol();
};

void          copiaRota(const Rota &rotaFonte, Rota &rota);
std::ostream &operator<<(std::ostream &os, const Solucao &sol);
std::ostream &operator<<(std::ostream &os, const Bin &bin);

std::string printBinEps(const Bin &bin);

bool verificaColisaoDoisItens(int                  item0,
                              int                  item1,
                              const Ponto         &p0,
                              const Ponto         &p1,
                              InstanceNS::Rotation r0,
                              InstanceNS::Rotation r1);

std::string printPonto(const Ponto &ponto, int dim);

int    calculaNumBinOcupados(const Solucao &solucao);
double calculaVolumeOcupado(const Solucao &solucao);
double calculaVolumeLivre(const Solucao &solucao);
double calculaMenorAreaLivre(const Solucao &solucao);

inline __attribute__((always_inline)) bool pontosIguais(const Ponto &p0, const Ponto &p1)
{
#pragma GCC unroll 3
    for(int i = 0; i < 3; ++i)
    {
        if(p0.vetDim[i] != p0.vetDim[i])
            return false;
    }

    return true;
}

bool checkUnloadingSequence(Bin 										&bin,
                            Rota 										&rota,
                            Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &matSupportItems);

inline __attribute__((always_inline)) int findPos(Rota &rota, int itemId)
{

    int customer = InstanceNS::instanciaG.vetItens[itemId].customer;
    return std::distance(
        rota.vetRota.begin(),
        std::find(rota.vetRota.begin(), rota.vetRota.begin() + rota.numPos, customer));
}

inline __attribute__((always_inline)) bool isBehind(InstanceNS::Item    &item0,
                                                    Ponto                p0,
                                                    InstanceNS::Rotation r0,
                                                    InstanceNS::Item    &item1,
                                                    Ponto                p1,
                                                    InstanceNS::Rotation r1)
{
    // double maxX0, maxY0, maxZ0, maxX1, maxY1, maxZ1;

    double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
    double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);
    double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

    // double maxX1 = p1.vetDim[0] + item1.vetDim[0];
    double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);
    double maxZ1 = p1.vetDim[2] + item1.getDimRotacionada(2, r1);

    return maxX0 <= p1.vetDim[0] && maxZ0 > p1.vetDim[2] && p0.vetDim[0] < maxZ1 &&
           maxY0 > p1.vetDim[1] && p0.vetDim[1] < maxY1;
}

inline __attribute__((always_inline)) bool isBelow(InstanceNS::Item    &item0,
                                                   Ponto                p0,
                                                   InstanceNS::Rotation r0,
                                                   InstanceNS::Item    &item1,
                                                   Ponto                p1,
                                                   InstanceNS::Rotation r1,
                                                   bool                 touch)
{
    double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
    double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);

    double maxX1 = p1.vetDim[0] + item1.getDimRotacionada(0, r1);
    double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);

    double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

    return ((touch && maxZ0 == p1.vetDim[2]) || !touch && maxZ0 <= p1.vetDim[2]) &&
           p0.vetDim[1] < maxY1 && p1.vetDim[1] < maxY0 && p0.vetDim[0] < maxX1 &&
           p1.vetDim[0] < maxX0;
}

inline __attribute__((always_inline))
bool lifo(InstanceNS::Item    							&item0,
          Ponto                							p0,
          InstanceNS::Rotation 							r0,
          InstanceNS::Item    							&item1,
          Ponto                							p1,
          InstanceNS::Rotation 							r1,
          bool                 						   	mlifo,
          bool						                   	removeFromShortSide,
          Eigen::Matrix<int, -1, -1, Eigen::RowMajor> 	&matSupportItems)
{
    // Item0 is delevery first

    double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
    double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);
    double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

    double maxX1 = p1.vetDim[0] + item1.getDimRotacionada(0, r1);
    double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);
    double maxZ1 = p1.vetDim[2] + item1.getDimRotacionada(2, r1);

    bool overlapX = !(maxX0 <= p1.vetDim[0] || maxX1 <= p0.vetDim[0]);
    bool overlapY = !(maxY0 <= p1.vetDim[1] || maxY1 <= p0.vetDim[1]);

           // matSupportItems[i][j] = 1 → j supports i
    bool item1_supports_item0 = matSupportItems(item1.itemId, item0.itemId);
    bool item0_supports_item1 = matSupportItems(item0.itemId, item1.itemId);

    if(mlifo && !removeFromShortSide)
    {
        if (overlapX)
        {
            bool in_front_block = (maxY1 < p0.vetDim[1]);
            bool vertical_block = item1_supports_item0;
            return !(in_front_block && vertical_block);
        }

        return true;
    }
    else if (!mlifo && removeFromShortSide)
    {
        if (overlapY)
        {
            // <=
            bool infront = (maxX1 <= p0.vetDim[0]);
            bool above  =  (maxZ1 <= p0.vetDim[2]);
            return infront || above;
        }

        return true;
    }
    else if(!mlifo && !removeFromShortSide)
    {
        if (overlapX)
        {
            bool left  =  (maxY1 <= p0.vetDim[1]);
            bool below = (maxZ1 <= p0.vetDim[2]);
            return left || below;
        }

        return true;
    }
    else if(mlifo && removeFromShortSide)
    {
        if (overlapY)
        {
            bool behind_block = (maxX1 < p0.vetDim[0]);
            bool vertical_block = item1_supports_item0;
            return !(behind_block && vertical_block);
        }
        return true;
    }
    else
    {
        std::printf("ERROR in LIFO or removeFromShortSide parameters!\n");
        PRINT_THROW();
    }


}

inline const Ponto PontoZero(0.0, 0.0, 0.0);

INLINE
double computeLeftBalancedLoading(double y, double width, int mass)
{
    double center = (InstanceNS::instanciaG.vetDimVeiculo[1] / 2.0);
    // center -
    return (mass / width) *
           (std::max(0.0, center - y) - std::max(0.0, center - (y + width)));
}

INLINE
double computeRightBalancedLoading(double y, double width, int mass)
{
    double center = (InstanceNS::instanciaG.vetDimVeiculo[1] / 2.0);
    return (mass / width) *
           (std::max(0.0, (y + width) - center) - std::max(0.0, y - center));
}

double getIntercetion(int                  item0,
                      Ponto                p0,
                      InstanceNS::Rotation r0,
                      Face                 f0,
                      int                  item1,
                      Ponto                p1,
                      InstanceNS::Rotation r1,
                      Face                 f1);

INLINE
bool tochRightSideOfTruck(int item, Ponto p, InstanceNS::Rotation r)
{
    double maxY =
        p.vetDim[0] + InstanceNS::instanciaG.vetItens[item].getDimRotacionada(0, r);

    return doubleEqual(maxY, InstanceNS::instanciaG.vetDimVeiculo[0], DifDistColision);
}

INLINE
bool tochLeftSideOfTruck(int item, Ponto p, InstanceNS::Rotation r)
{
    return doubleEqual(0.0, p.vetDim[1], DifDistColision);
}

INLINE
bool tochBackSideOfTruck(int item, Ponto p, InstanceNS::Rotation r)
{
    return doubleEqual(0.0, p.vetDim[0], DifDistColision);
}

INLINE
bool tochFrontSideOfTruck(int item, Ponto p, InstanceNS::Rotation r)
{
    double maxX =
        p.vetDim[1] + InstanceNS::instanciaG.vetItens[item].getDimRotacionada(1, r);

    return doubleEqual(maxX, InstanceNS::instanciaG.vetDimVeiculo[1], DifDistColision);
}

bool checkCompactness(Bin &bin, const VectorI &vetTop, std::string *strError = nullptr);
} // namespace SolucaoNS

#endif // INC_2L_SDVRP_SOLUCAO_H
