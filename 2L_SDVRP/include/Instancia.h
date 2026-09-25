/* ****************************************
 * ****************************************
 *  Data:    05/11/24
 *  Arquivo: Instance.h
 * ****************************************
 * ****************************************/

#ifndef INC_2L_SDVRP_INSTANCIA_H
#define INC_2L_SDVRP_INSTANCIA_H

#include "Constants.h"
#include "safe_matrix.h"
#include "safe_vector.h"
#include "sefe_array.h"
#include "string"
#include "AuxT.h"
#include <iostream>
#include <map>
#include <string>

namespace InstanceNS
{ // 0->0
// 2 -> 1
// 2 -> 2
enum Rotation
{
    Rot0 = 0, // LWH
    Rot1,     // WLH
    Rot2,     // LHW
    Rot3,     // WHL
    Rot4,     // HLW
    Rot5      // HWL
};

class Item
{
  public:
    Array<double, 3> vetDim; // Length, width, height
    Array<double, 3> vetDimNor;
    double           volume           = 0.0;
    double           weight           = 0.0;
    double           weightForce      = 0.0;
    bool             fragility        = false;
    int              customer 		  = -1;
    int              oroloc3D_item_id = -1;
    int              itemId           = -1;
    double			 volNormDual      = 0;
    std::string      oroloc3D_item_id_str;

    Item() = default;
    Item(double x, double y, double z, double peso_, int itemId_);
    void        set(double x, double y, double z, int itemId_);
    std::string print(Rotation r=Rot0, bool printVol = false);
    double      getDimRotacionada(int d, Rotation r);
    INLINE
    void setNorm(double term0, double term1, double term2)
    {
        //for(int d=0; d < 3; ++d)
        vetDimNor[0] = vetDim[0]/term0;
        vetDimNor[1] = vetDim[1]/term1;
        vetDimNor[2] = vetDim[2]/term2;
        volNormDual = vetDimNor[0]*vetDimNor[1]*vetDimNor[2];
    }
};

struct TW
{
    double ini = 0.0, fim = INF_Double;
};

class Instance
{
  public:
    std::string nome;
    int         numClientes   		= 0;
    int         numItens      		= 0;
    Vector<int> vetNumItensPorCli;
    int         numVeiculos   		= 0;
    int         numDim        		= 2;
    int         numRotation   		= 6;
    double      maxPayload    		= 0;
    double      maxItemVolume 		= 0.0;

    Array<double, 3> vetDimVeiculo;
    // double veicAltura = 0.0;
    // double veicLargura = 0.0;

    bool split = false;
    bool packing = true;

    Matrix<double> matDist;
    Matrix<double> matTempo;
    Vector<TW>     vetTw;

    Vector<Item> vetItens;
    VectorD      vetPesoItens;
    VectorD      vetMinDimItens;
    VectorD 	 vetVolNormDualCust;

    int maxNumItensPorClie = 0;
    double volNormal       = 0.0;

    // Vector<double> vetItemAltura;
    // Vector<double> vetItemLargura;
    // Vector<double> vetItemArea;

    Vector<double> vetDemandaCliente;
    Vector<double> vetVolumeCliente;

    Vector<int> vetItemCliente;    // Indica o cliente dado um item; mat[itemId] = cliente
    Matrix<int> matCliItensIniFim; // Indica o id Inicial e id Final do primeiro e ultimo
                                   // item de um cliente; mat[clienteId,0] = iten0;
                                   // mat[clienteId,1] = itenFim
    Vector<int>                vetOrderId;     // Indicates for every item, its orderId
    std::map<int, Vector<int>> mapOrderIdItem; // Maps orderId to its items
    std::map<int, int>         mapOrderIdCust;
    std::map<int, int> mapItem_IdItem; // Maps original item id to itemId of the instance
    std::map<int, int> mapCustomer_idToCustomer; // Maps the original customer_id to the
                                                 // customer of the instance
    std::map<int, int> mapCustomerToCustomer_id; // Maps the  customer of the instance to
                                                 // the original customer_id

    Instance();
    Instance(int numClientes_, int numItens_, int numVeiculos_);
    void atualizaVetMinDimItens();
    void setMaxItemVolume();
};

void   read2dInstance(const std::string &strFile);
void   read3dInstance(const std::string &strFile);
void   readOroloc3D(const std::string &strFile);
void   readOroloc3D2(const std::string &strFile);
void   convertInstanceToCm(Instance& instance);
int    copiaItensCliente(int cliente, VectorI &vetItens);
int    copiaItensClientes(VectorI &vetClientes,
                          int      tam,
                          VectorI &vetItens,
                          bool     push = false);
double calculaDistancia(VectorI &vet, int tam);
int    generateRandomListOfItems(int numItens, VectorI &vetItems);


inline Instance                                    instanciaG;
inline static const Array<InstanceNS::Rotation, 1> vetRot = {Rot0};//, Rot1}; //, Rot2};
// std::string printItem(int itemId);
} // namespace InstanceNS

#endif // INC_2L_SDVRP_INSTANCIA_H
