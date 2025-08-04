/*  *****************************************************************
 *  *****************************************************************
 *  File:    DW_Decomp.h
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    18/09/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "gurobi_c++.h"
#include <Eigen/Eigen>
#include <memory>
#include "Aux.h"
#include "Sparse.h"
#include <boost/unordered_set.hpp>
#include <set>


#ifndef DW_DECOMP_DW_DECOMP_H
#define DW_DECOMP_DW_DECOMP_H

using namespace SparseNS;

namespace DW_DecompNS
{

    constexpr double TolObjSubProb       = 1E-3; // 1E-3
    // TODO Retornar para 25!
    constexpr int    NumMaxSolSubProb    = 25;//25;
    constexpr double StabilizationAlpha  = 0.6; // .9 // .45  // 0.60
    inline bool      Stabilization       = true;
    constexpr double gapLimit            = 1E-3;
    constexpr int    NumCandidatesBranch = 3;
    constexpr bool   PrintDebug          = false;
    constexpr int    BigM_maxMult        = 10;
    constexpr double GapExactPricing     = 1E-1; // gap  <= ExactPricingTol%
    constexpr double GapTolStop          = 1E-1; // gap* <= GapTolStop%



    Eigen::MatrixXd getMatA_Model(GRBModel &mestre);
    void getSparseMatModel(GRBModel &model, Eigen::SparseMatrix<double, Eigen::RowMajor> &matA);

    Eigen::VectorXd getVetC_Model(GRBModel &mestre);
    void getVetC_Model(GRBModel &model, Eigen::RowVectorXd &vetC);

    Eigen::VectorX<char> getConstSenseModel(GRBModel &model);
    Eigen::VectorXd getRhsModel(GRBModel &model);
    void recuperaX(GRBVar *var, Eigen::VectorXd &vetX, int numVar);


    enum StatusProb
    {
        StatusSubProb_Otimo    = 0,
        StatusSubProb_Inviavel,
        StatusSubProb_Unbounded,
        StatusSubProb_Outro
    };

    enum class PhaseStatus
    {
        PhaseStatusColGen   = 0,
        PhaseStatusBigM     = 1,
        PhaseStatusTwoPhase = 2
    };


    class SubProb
    {
        public:

        Eigen::MatrixXd *matA = nullptr;

        SubProb()=default;
        virtual ~SubProb(){};

        virtual int getNumConvConstr()=0;

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
         *  @param pairSubProb           Pair representa o indice do inicio das
         *                             variaveis do sub problema e o seu tamanho
         *
         *  ************************************************************************
         *  ************************************************************************
         */

        virtual int
        resolveSubProb(const Eigen::VectorXd &vetC,
                       Eigen::RowVectorXd &vetRowPi,
                       GRBModel &rmlp,
                       int itCG,
                       bool &custoRedNeg,
                       void *data,
                       const int iniConv,
                       int indSubProb,
                       Eigen::VectorXd &vetCooefRestConv,
                       const std::pair<int, int> &pairSubProb,
                       Eigen::MatrixXd &matColX,
                       int &numSol,
                       Eigen::Array<double, 1, NumMaxSolSubProb>& vetRedCost,
                       double constPiValue,
                       const VectorI &vetVar0,
                       const VectorI &vetVar1,
                       PhaseStatus phaseStatus,
                       bool exact) =0;


        //virtual int64_t getNumberOfConvConstr() = 0;
    };


    void dwDecomp(GRBEnv &env,
                  GRBModel &mestre,
                  double custoVarA,
                  const std::vector<std::pair<int, int>> &&vetPairSubProb,
                  SubProb *subProb,
                  void *data,
                  const int numSubProb);

    template<typename T>
    class DelVetFunctor
    {
    public:
        void operator()(T *ptr){ delete []ptr;}
    };

    class DW_DecompNode;

    class StabilizationData
    {
    public:
        Eigen::Matrix<double, -1, -1, Eigen::RowMajor> matPi;
        int numRow = 0;


        explicit StabilizationData(int numVarRmlpPi);
        StabilizationData()=default;
        void start(int numVarRmlpPi);

        void addPi(const Eigen::RowVectorXd &vetRowRmlpPi);
        void getWeightSum(const double alpha, Eigen::RowVectorXd &vetRowRmlpPi);


    };

    struct AuxData
    {
        /// For each subproblem keeps the start and the size of the variables
        Vector<std::pair<int, int>> vetPairSubProb;

        Eigen::RowVectorXd          vetRowRmlpPi;
        Eigen::RowVectorXd          vetRowRmlpSmoothPi;
        Eigen::RowVectorXd          vetRowC;
        Eigen::VectorXd             vetColSubProbCooef;
        Eigen::MatrixXd             matColX_solSubProb;
        Eigen::VectorXd             vetColConvCooef;
        Eigen::VectorXd             vetColCooef;
        double*                     auxVetCooef = nullptr;

        StabilizationData           stabD;

        void updateSizes(DW_DecompNS::DW_DecompNode &e);
        AuxData()=default;
        ~AuxData();

        inline
        void updateAuxVetCooef()
        {
            for(int i=0; i < vetColCooef.size(); ++i)
                auxVetCooef[i] = vetColCooef[i];
        }
    };

    struct Info
    {
        int     numSubProb               = 0;        // Number of sub problems
        double  costA_Var                = 0.0;
        int64_t numConstrsConv           = 0;        // Number of constants of convexity
        int     numConstrsMaster         = 0;        // Number of the master's constrants
        int     numVarMaster             = 0;        // Number of the variables into the master mip
        int     numVarRmlpPi             = 0;        // Number of dual variables into the rmlp
        int     numConstrsOrignalProblem = 0;
    };


    class SolXHash
    {
    public:

        const Eigen::VectorXd &vetX;
        uint64_t hashVal;

        explicit SolXHash(const Eigen::VectorXd &vetX_);

    };

    inline
    bool operator == (const SolXHash &sol0, const SolXHash &sol1)
    {
        if(sol0.hashVal != sol1.hashVal)
            return false;

        if(sol0.vetX.size() != sol1.vetX.size())
            return false;

        for(int i=0; i < sol0.vetX.size(); ++i)
        {
            if(sol0.vetX[i] != sol1.vetX[i])
                return false;
        }

        return true;
    }

    inline
    std::size_t  hash_value(const SolXHash& solXHash){return solXHash.hashVal;}




    class DW_DecompNode
    {
    public:
        SubProb* ptrSubProb = nullptr;               // Ptr for class SubProb (DO NOT DELETE)
        std::unique_ptr<GRBModel>  uRmlp;            // Ptr for rmlp model
        int itCG = 0;
        Info info;
        double funcObj = 0.0;
        double rhsConv = 0.0;

        Eigen::SparseMatrix<double, Eigen::RowMajor> matA;
        Vector<std::unique_ptr<Eigen::VectorXd>> vetVarLambdaCol;         // Vectors generated by sub problems
        boost::unordered_set<SolXHash> setVarLamdaCol;
        Vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> vetSubMatA;
        Eigen::VectorXd vetSolX;

        GRBConstr* vetRmlpConstr = nullptr;
        GRBVar* vetVarArtifRmlp  = nullptr;

        VectorI vetVar0;   // set of variables with have x_i <= 0
        VectorI vetVar1;   // set of variables with have x_i >= 1

        VectorD vetObjCof;

        PhaseStatus phaseStatus = PhaseStatus::PhaseStatusBigM;

        DW_DecompNode(GRBEnv &env_,
                      GRBModel &master_,
                      double costA_Var_,
                      SubProb *ptrSubProb_,
                      const int numSubProb_,
                      AuxData &auxVect);

        StatusProb columnGeneration(AuxData &auxVect);
        void getSubProbCooef(int k, AuxData &auxVect);
        void addColumn(const double cost, int k, AuxData &auxVect);
        void addColumnX(const double cost, int k, AuxData& auxVect, Eigen::VectorXd& vetX);
        void updateRmlpPi(Eigen::RowVectorXd &vetRowRmlpPi);
        double getLagrangeDualBound(double objRmlp, double redCost);
        void getSolX();


        DW_DecompNode(const DW_DecompNode &decomp);

        ~DW_DecompNode()
        {
            //std::cout<<"~DW_DecompNode\n";
            delete []vetRmlpConstr;
            delete []vetVarArtifRmlp;
        }

    };


}


#endif //DW_DECOMP_DW_DECOMP_H
