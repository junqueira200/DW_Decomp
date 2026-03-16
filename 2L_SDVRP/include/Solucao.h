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

#include "safe_vector.h"
#include "safe_matrix.h"
#include "Instancia.h"
#include "AuxT.h"
//#include "InputOutput.h"


namespace SolucaoNS
{

    struct Ponto
    {
        Array<double,3> vetDim;

        Ponto(){vetDim.setAll(0.0);}
        Ponto(double d0, double d1, double d2)
        {
            vetDim[0] = d0;
            vetDim[1] = d1;
            vetDim[2] = d2;
        }
        void set(double d0=0.0, double d1=0.0, double d2=0.0)
        {
            vetDim[0] = d0;
            vetDim[1] = d1;
            vetDim[2] = d2;
        }

        std::string print()const;
    };

    struct PontoRot
    {
        Ponto ponto;
        InstanceNS::Rotation rot;

        PontoRot(){ponto.set(); rot=InstanceNS::Rot0;}
        PontoRot(double px, double py, double pz, InstanceNS::Rotation r)
        {
            ponto.vetDim[0] = px;
            ponto.vetDim[1] = py;
            ponto.vetDim[2] = pz;
            rot = r;
        }
        void set(double px=0.0, double py=0.0, double pz=0.0, InstanceNS::Rotation r=InstanceNS::Rot0)
        {
            ponto.vetDim[0] = px;
            ponto.vetDim[1] = py;
            ponto.vetDim[2] = pz;
            rot = r;
        }
    };

    struct Bin
    {

        //VectorD vetX, vetY;          // Indica aa posicoes do canto inferior esquerdo onde o item eh posicionado no bin
        //VectorD vetEpX, vetEpY;      // Indica as posicoes dos pontos extremos

        Vector<Ponto>                   vetPosItem;             // Indica as posicoes do canto inferior esquerdo onde o item eh posicionado no bin
        Vector<Ponto>                   vetEp;                  // Indica as posicoes dos pontos extremos
        VectorI                         vetItemId;              // Indica o id do iº item
        Vector<InstanceNS::Rotation>    vetRotacao;
        Vector<int8_t>                  vetItens;               // Indica se o bin empacota o i° item
        Array<double,3>                 binDim;
        double                          volumeTotal     = 0.0;
        double                          volumeOcupado   = 0.0;
        double                          demandaTotal    = 0.0;
        int                             numItens        = 0;
        int                             numEps          = 0;    // Numero de pontos extremos

        void addItem(int idEp, int idItem, InstanceNS::Rotation r);
        void addEp(const Ponto &ep);

        inline __attribute__((always_inline))
        void setItem(int pos, int itemId, double x, double y, double z)
        {

            #if VAR_SOLUTION_BORROW_CHECKER
                if(pos >= numItens)
                    throw std::out_of_range("");
            #endif

            vetPosItem[pos].vetDim[0]   = x;
            vetPosItem[pos].vetDim[1]   = y;
            vetPosItem[pos].vetDim[2]   = z;
            vetItemId[pos]              = itemId;

        };

        inline __attribute__((always_inline))
        void setEp(int pos, double x, double y, double z)
        {

            #if VAR_SOLUTION_BORROW_CHECKER
                if(pos >= numEps)
                    throw std::out_of_range("");
            #endif

            vetEp[pos].vetDim[0] = x;
            vetEp[pos].vetDim[1] = y;
            vetEp[pos].vetDim[2] = z;
        };

        inline __attribute__((always_inline))
        bool vazio() const
        {
            return numItens == 0;
        }

        bool rmI_Item(int i);
        bool rmI_Ep(int i);

        void reset();
        int getEpComMenorCoord(const VectorI &vetIdEp, int tam);
        bool verificaViabilidade();

        void rmItens(const VectorI &vetItensRm, const int tam);

        double getPorcentagemUtilizacao()const;

        Bin();
        Bin(const Bin &bin)=delete;
        Bin& operator=(const Bin &bin)=delete;
    };

    int getBinVazio(const Vector<Bin> &vetBin, int tam);
    void copiaBin(const Bin &binFonte, Bin &bin);

    class Rota
    {
    public:

        Rota();
        Rota(const Rota &rota)=delete;
        void reset();
        std::string printRota();
        void computeDistance();

        VectorI        vetRota;
        VectorD        vetTempoSaida;
        VectorD        vetDemClie;
        Vector<int8_t> vetItens;        // Indica se o bin empacota o i° item
        int     numPos          = 2;
        //double  demTotal        = 0.0;
        double  distTotal       = 0.0;
        Bin*    binPtr          = nullptr;
    };

    class Solucao
    {
    public:

        Vector<Bin>     vetBin;
        Vector<Rota>    vetRota;
        double          distTotal = 0.0;

        Solucao();
        explicit Solucao(const InstanceNS::Instance &instancia);

        bool verificaSol(std::string &error);
        int getBinVazio();
        void copiaSolucao(const Solucao &sol);
        void reset();
        double getUtilizacaoMediaBins()const;
        double getUtilizacaoMedianaBins()const;
        double getTamMedianaRota()const;


    };

    void copiaRota(const Rota &rotaFonte, Rota &rota);
    std::ostream& operator<<(std::ostream &os, const Solucao &sol);
    std::ostream& operator<<(std::ostream &os, const Bin& bin);

    std::string printBinEps(const Bin &bin);

    bool verificaColisaoDoisItens(int item0,
                                  int item1,
                                  const Ponto &p0,
                                  const Ponto &p1,
                                  InstanceNS::Rotation r0,
                                  InstanceNS::Rotation r1);

    std::string printPonto(const Ponto &ponto, int dim);

    int calculaNumBinOcupados(const Solucao &solucao);
    double calculaVolumeOcupado(const Solucao &solucao);
    double calculaVolumeLivre(const Solucao &solucao);
    double calculaMenorAreaLivre(const Solucao &solucao);

    inline __attribute__((always_inline))
    bool pontosIguais(const Ponto &p0, const Ponto &p1)
    {
        #pragma GCC unroll 3
        for(int i=0; i < 3; ++i)
        {
            if(p0.vetDim[i] != p0.vetDim[i])
                return false;
        }

        return true;
    }

    bool checkUnloadingSequence(Bin& bin, Rota& rota);

    inline __attribute__((always_inline))
    int findPos(Rota& rota, int itemId)
    {

        int customer = InstanceNS::instanciaG.vetItens[itemId].customer;
        return std::distance(rota.vetRota.begin(), std::find(rota.vetRota.begin(), rota.vetRota.begin()+rota.numPos,
                                                             customer));

    }

    inline __attribute__((always_inline))
    bool isBehind(InstanceNS::Item& item0, Ponto p0, InstanceNS::Rotation r0,
                  InstanceNS::Item& item1, Ponto p1, InstanceNS::Rotation r1)
    {
        //double maxX0, maxY0, maxZ0, maxX1, maxY1, maxZ1;

        double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
        double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);
        double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

        //double maxX1 = p1.vetDim[0] + item1.vetDim[0];
        double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);
        double maxZ1 = p1.vetDim[2] + item1.getDimRotacionada(2, r1);

        return maxX0        <= p1.vetDim[0] &&
               maxZ0        > p1.vetDim[2]  &&
               p0.vetDim[0] < maxZ1         &&
               maxY0        > p1.vetDim[1]  &&
               p0.vetDim[1] < maxY1;

    }

    inline __attribute__((always_inline))
    bool isBelow(InstanceNS::Item& item0, Ponto p0, InstanceNS::Rotation r0,
                 InstanceNS::Item& item1, Ponto p1, InstanceNS::Rotation r1, bool touch)
    {
        double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
        double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);

        double maxX1 = p1.vetDim[0] + item1.getDimRotacionada(0, r1);
        double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);

        double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

        return ((touch && maxZ0 == p1.vetDim[2]) || !touch && maxZ0 <= p1.vetDim[2])
               && p0.vetDim[1] < maxY1
               && p1.vetDim[1] < maxY0
               && p0.vetDim[0] < maxX1
               && p1.vetDim[0] < maxX0;

    }

    inline __attribute__((always_inline))
    bool lifo(InstanceNS::Item& item0, Ponto p0, InstanceNS::Rotation r0,
              InstanceNS::Item& item1, Ponto p1, InstanceNS::Rotation r1, bool mlifo)
    {
            double maxX0 = p0.vetDim[0] + item0.getDimRotacionada(0, r0);
            double maxY0 = p0.vetDim[1] + item0.getDimRotacionada(1, r0);
            double maxZ0 = p0.vetDim[2] + item0.getDimRotacionada(2, r0);

            double maxX1 = p1.vetDim[0] + item1.getDimRotacionada(0, r1);
            double maxY1 = p1.vetDim[1] + item1.getDimRotacionada(1, r1);
            double maxZ1 = p1.vetDim[2] + item1.getDimRotacionada(2, r1);

            if((maxX1 <= p0.vetDim[0]) ||   // The end of j is less then the begining of i. It's correct for LIFO
                (maxY0 <= p1.vetDim[1]) ||   // Item i is complete at left of item j. It's correct for LIFO
                (maxY1 <= p0.vetDim[1]) ||   // Item i is complete at right of item j. It's correct for LIFO
                (maxZ1 <= p0.vetDim[2]) ||   // Item i is above item j. It's correct for LIFO
                (maxZ0 <= p1.vetDim[2]))     // Item i is below item j. It's correct for LIFO
                return true;

            return mlifo;
    }

    inline const Ponto PontoZero(0.0, 0.0, 0.0);


}

#endif //INC_2L_SDVRP_SOLUCAO_H
