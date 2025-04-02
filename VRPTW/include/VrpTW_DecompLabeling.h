/*  *****************************************************************
 *  *****************************************************************
 *  File:    VrpTW_DecompLabeling.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    20/11/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef DW_VRPTW_DECOMPLABELING_H
#define DW_VRPTW_DECOMPLABELING_H

#include "DW_Decomp.h"
#include "Instancia.h"
#include "LabelingAlgorithm.h"

namespace VrpTW_DecompLabelingNS
{
    // Linear index for a nxn matrix
    inline __attribute__((always_inline))
    int getIndex(int i, int j, int numClie){return i*numClie+j;}


    class VrpLabelingSubProb : public DW_DecompNS::SubProb
    {
    public:


        int numSubProb = 1;
        bool convConstIni = false;
        InstanciaNS::InstVRP_TW* instVrpTw = nullptr;                  // NAO DELETAR

        Eigen::Vector<LabelingAlgorithmNS::Step, 2> vetStepSize;
        LabelingAlgorithmNS::LabelingData           labelingData;
        LabelingAlgorithmNS::Vet3D_ResCost          vetMatResCost;
        LabelingAlgorithmNS::MatBoundRes            vetVetResBound;
        LabelingAlgorithmNS::NgSet                  ngSet;

        int getNumConvConstr() override {return 1;}
        VrpLabelingSubProb()=default;
        explicit VrpLabelingSubProb(InstanciaNS::InstVRP_TW &instVrpTw, double startDis);
        //int64_t getNumberOfConvConstr() override{return 0;};// {return numSubProb;}
        ~VrpLabelingSubProb() override =default;
        void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) override;

        int resolveSubProb(const Eigen::VectorXd &vetC,
                           Eigen::RowVectorXd &vetRowPi,
                           GRBModel &mestre,
                           int itCG,
                           bool &custoRedNeg,
                           void *data,
                           const int iniConv,
                           int indSubProb,
                           Eigen::VectorXd &vetCooefRestConv,
                           const std::pair<int, int> &pairSubProb,
                           Eigen::MatrixXd &matColX,
                           int &numSol,
                           double &redCost,
                           double constPiValue,
                           const VectorI &vetVar0,
                           const VectorI &vetVar1,
                           DW_DecompNS::PhaseStatus phaseStatus,
                           bool exact) override;

    }; // FIM MySubProb

    bool exactPricing(const LabelingAlgorithmNS::Vet3D_ResCost&          vetMatResCost,
                      const FloatType                                    startVal,
                      Eigen::Matrix<double, -1, 1, Eigen::ColMajor>&     vetColSolX,
                      InstanciaNS::InstVRP_TW&                           instVrpTw,
                      double&                                            cost);



}

#endif //DW_VRPTW_DECOMPLABELING_H
