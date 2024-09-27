#include "gurobi_c++.h"
#include <Eigen/Eigen>
#include <memory>
#include "Aux.h"

#ifndef DW_DECOMP_DW_DECOMP_H
#define DW_DECOMP_DW_DECOMP_H


namespace DW_DecompNS
{
    inline const double TolObjSubProb = 1E-5;

    Eigen::MatrixXd getMatA_Model(GRBModel &mestre);
    Eigen::VectorXd getVetC_Model(GRBModel &mestre);
    Eigen::VectorX<char> getConstSenseModel(GRBModel &model);
    Eigen::VectorXd getRhsModel(GRBModel &model);
    void recuperaX(GRBVar *var, Eigen::VectorXd &vetX, int numVar);


    enum StatusSubProb
    {
        StatusSubProb_Otimo    = 0,
        StatusSubProb_Inviavel,
        StatusSubProb_Unbounded,
        StatusSubProb_Outro
    };


    class SubProb
    {
        public:

        virtual ~SubProb(){};

        // Responsavel por adicionar as restricoes de convexidade
        virtual void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) =0;

        /** ************************************************************************
         *  ************************************************************************
         *
         *  Resolve o subproblema e add cooef das restricoes de conv
         *
         *  @param subProbCooef
         *  @param rmlp
         *  @param vetX
         *  @param itCG
         *  @param custoRedNeg
         *  @param data
         *  @param iniConv
         *  @param indSubProb
         *  @param vetCooefRestConv
         *  @param pairSubProb           pair representa o indice do inicio das variaveis do sub problema e o seu tamanho
         *
         *  ************************************************************************
         *  ************************************************************************
         */
        virtual int resolveSubProb(Eigen::VectorXd &subProbCooef,
                                   GRBModel &rmlp,
                                   Eigen::VectorXd &vetX,
                                   int itCG,
                                   bool &custoRedNeg,
                                   void *data,
                                   const int iniConv,
                                   int indSubProb,
                                   Eigen::VectorXd &vetCooefRestConv,
                                   const std::pair<int, int> &pairSubProb)=0;
    };


    void dwDecomp(GRBEnv &env,
                  GRBModel &mestre,
                  double custoVarA,
                  const std::vector<std::pair<int,int>> &&vetPairSubProb,
                  SubProb *subProb,
                  void *data,
                  const int numConstrsConv,
                  const int numSubProb);
}


#endif //DW_DECOMP_DW_DECOMP_H
