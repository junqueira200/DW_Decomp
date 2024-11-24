//
// Created by igor on 24/11/24.
//

#ifndef DW_VRPTW_DECOMPLABELING_H
#define DW_VRPTW_DECOMPLABELING_H

#include "DW_Decomp.h"
#include "Instancia.h"
#include "LabelingAlgorithm.h"

namespace VrpTW_DecompLabelingNS
{


    class VrpLabelingSubProb : public DW_DecompNS::SubProb
    {
    public:


        int numSubProb = 1;
        bool convConstIni = false;
        InstanciaNS::InstVRP_TW *instVrpTw;                  // NAO DELETAR

        Eigen::Vector<LabelingAlgorithmNS::Step, 2> vetStepSize;
        LabelingAlgorithmNS::LabelingData labelingData;
        LabelingAlgorithmNS::VetMatResCost vetMatResCost;
        LabelingAlgorithmNS::VetVetResBound vetVetResBound;
        LabelingAlgorithmNS::NgSet ngSet;

        int getNumConvConstr() override {return 0;}
        VrpLabelingSubProb()=default;
        explicit VrpLabelingSubProb(InstanciaNS::InstVRP_TW &instVrpTw);
        int64_t getNumberOfConvConstr() override{return 0;};// {return numSubProb;}
        ~VrpLabelingSubProb() override =default;
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

#endif //DW_VRPTW_DECOMPLABELING_H
