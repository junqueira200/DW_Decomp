//
// Created by igor on 26/09/24.
//

#ifndef DW_MNFPDECOMP_H
#define DW_MNFPDECOMP_H
#include "MNFP_Inst.h"
#include "DW_Decomp.h"
#include "Grafo.h"


namespace MnfpDecompNS
{

    void criaSubProbFlow(const MNFP::MNFP_Inst &mnfp, GRBModel &model, int k);
    void criaMestre(const MNFP::MNFP_Inst &mnfp, GRBModel &model);

    class MySubProbFlow : public DW_DecompNS::SubProb
    {
    public:

        Eigen::VectorX<std::unique_ptr<GRBModel>> vetSubProb;
        const int numSubProb = 2;
        bool convConstIni = false;

        MySubProbFlow(GRBEnv &e, const MNFP::MNFP_Inst &mnfp)
        {
            vetSubProb = Eigen::VectorX<std::unique_ptr<GRBModel>>(numSubProb);
            for(int i=0; i < numSubProb; ++i)
            {
                vetSubProb(i) = std::make_unique<GRBModel>(e);
                MnfpDecompNS::criaSubProbFlow(mnfp, *(vetSubProb(i)), i);
            }
        }


        ~MySubProbFlow() override {};
        void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) override;
        int resolveSubProb(Eigen::VectorXd &subProbCooef,
                           GRBModel &mestre,
                           Eigen::VectorXd &vetX,
                           int itCG,
                           bool &custoRedNeg,
                           void *data,
                           const int iniConv,
                           int indSubProb,
                           Eigen::VectorXd &vetCooefRestConv,
                           const std::pair<int, int> &pairSubProb) override;

    }; // FIM MySubProb

    enum VerticeType{TypeZero, TypeSource, TypeSink};


    class MySubProbPath : public DW_DecompNS::SubProb
    {
    public:

        const int numSubProb;
        bool convConstIni = false;
        MNFP::MNFP_Inst mnfp;
        Eigen::VectorX<GraphNS::Graph<double>> vetGraphCost;
        GraphNS::Graph<double> vetGraphModCost;             // Armazena os custos modificados para cada aresta
        Eigen::MatrixX<int64_t> matCovConst;                // Armazenam para cada (k, vertice) o indice da const de conv
        Eigen::MatrixX<VerticeType> matVerticeType;         // Armazenam para cada (k, vertice) o seu tipo


        MySubProbPath(GRBEnv &e, const MNFP::MNFP_Inst &mnfp_);
        ~MySubProbPath() override {};
        void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) override;
        int resolveSubProb(Eigen::VectorXd &subProbCooef,
                           GRBModel &rmlp,
                           Eigen::VectorXd &vetX,
                           int itCG,
                           bool &custoRedNeg,
                           void *data,
                           const int iniConv,
                           int indSubProb,
                           Eigen::VectorXd &vetCooefRestConv,
                           const std::pair<int, int> &pairSubProb) override;

        void restoreGraphModCost();

    }; // FIM MySubProb
}


#endif //DW_MNFPDECOMP_H
