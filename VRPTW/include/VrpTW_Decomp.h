//
// Created by igor on 15/11/24.
//

#ifndef DW_VRPTW_DECOMP_H
#define DW_VRPTW_DECOMP_H

#include "DW_Decomp.h"
#include "Instancia.h"

namespace VrpTW_DecompNS
{

    void criaMestre(const InstanciaNS::InstVRP_TW &instVrpTw, GRBModel &model);
    void criaVRP_TW_CompleteModel(const InstanciaNS::InstVRP_TW &instVrpTw, GRBModel &model);
    int getIndex(int i, int j, int numClie);
    double geraSolHeuristica(const InstanciaNS::InstVRP_TW &instVrpTw,
                             Eigen::VectorXd &vetRedCost,
                             Eigen::VectorXi &rota,
                             int &routeTam);

    class VrpSubProb : public DW_DecompNS::SubProb
    {
    public:
        void buildSubProbModel();

        std::unique_ptr<GRBModel> subProb;
        GRBVar* grbVarX       = nullptr;
        GRBVar* grbVarF       = nullptr;
        GRBVar* grbVarT_Cheg  = nullptr;
        GRBVar* grbVarT_Saida = nullptr;
        GRBVar* grbVarU       = nullptr;

        int numSubProb = 1;
        bool convConstIni = false;
        InstanciaNS::InstVRP_TW *instVrpTw;                  // NAO DELETAR

        int getNumConvConstr() override {return 0;}
        VrpSubProb(GRBEnv &e, InstanciaNS::InstVRP_TW &instVrpTw);
        //int64_t getNumberOfConvConstr() override;// {return numSubProb;}
        ~VrpSubProb() override;
        void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) override;
        int resolveSubProb(const Eigen::VectorXd &vetC,
                           const Eigen::RowVectorXd &vetRowPi,
                           GRBModel &mestre,
                           Eigen::VectorXd &vetX,
                           int itCG,
                           bool &custoRedNeg,
                           void* data,
                           const int iniConv,
                           int indSubProb,
                           Eigen::VectorXd &vetCooefRestConv,
                           const std::pair<int, int> &pairSubProb) override;

    }; // FIM MySubProb

}

#endif //DW_VRPTW_DECOMP_H
