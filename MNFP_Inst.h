#include <Eigen/Eigen>
#include "gurobi_c++.h"

#ifndef DW_DECOMP_MNFP_INST_H
#define DW_DECOMP_MNFP_INST_H


namespace MNFP
{

    class MNFP_Inst
    {
    public:

        int K = 1;
        int N = 0;

        Eigen::VectorX<Eigen::MatrixXi> vetArcCost;
        Eigen::MatrixXi matVertexDem;
        Eigen::MatrixXi matCapacidade;

        MNFP_Inst() = default;

        MNFP_Inst(int k_, int n_) : K(k_), N(n_)
        {
            vetArcCost = Eigen::VectorX<Eigen::MatrixXi>(K);
            for(int k = 0; k < K; ++k)
            {
                vetArcCost[k] = Eigen::MatrixXi(N, N);
                vetArcCost[k].setZero();
            }

            matVertexDem = Eigen::MatrixXi(K, N);
            matVertexDem.setZero();

            matCapacidade = Eigen::MatrixXi(N, N);
            matCapacidade.setZero();
        }

    };

    MNFP_Inst criaToyInstance();

    void criaSubProb(const MNFP_Inst &mnfp, GRBModel &model, int k);
    void criaMestre(const MNFP::MNFP_Inst &mnfp, GRBModel &model);

    int getId(int i, int j, int n);
}

#endif //DW_DECOMP_MNFP_INST_H