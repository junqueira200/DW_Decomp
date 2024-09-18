#include "gurobi_c++.h"
#include <Eigen/Eigen>
#include <memory>

#ifndef DW_DECOMP_DW_DECOMP_H
#define DW_DECOMP_DW_DECOMP_H


namespace DW_Decomp
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

    void dwDecomp(GRBEnv env,
                  GRBModel &mestre,
                  Eigen::VectorX<std::unique_ptr<GRBModel>> &vetSubProb,
                  double custoVarA,
                  int numVarSubProblema,
                  const std::vector<std::pair<int,int>> &&vetPairSubProb,
                  std::vector<int (*)(const Eigen::VectorXd &, int k, void *data, Eigen::VectorXd &vetX, int itCG, bool &custoRedNeg)> &&vetFuncPtrSubProb);



class SubProb
{
    virtual ~SubProb(){};
    virtual void iniSubProb()=0;
    virtual int resolveSubProb(const Eigen::VectorXd &, int k, Eigen::VectorXd &vetX, int itCG, bool &custoRedNeg)=0;
    virtual void addRestConv(const Eigen::VectorXd &vetX, int k, GRBModel &master, int iniConv)=0;
};


/* ****************************************************************************************************
 * ****************************************************************************************************
 *
 * @param env                   Env do Gurobi
 * @param mestre                Modelo com as restricoes do mestre e a func obj do problema original
 * @param vetSubProb            Vetor com os modelos de subproblemas
 * @param custoVarA             Custo das variaveis artificiais, para inicializar a CG
 * @param numVarSubProblema     Numero de variaveis de cada subproblema
 * @param vetVarSubProb         Vetor de pair por subproblema, cada pair representa o indice do inicio das variaveis e o seu tamanho
 * @param vetFuncPtr            Vetor de ptr de funcoes
 *
 * ****************************************************************************************************
 * *****************************************************************************************************/
    void DW_Decomp::dwDecomp(GRBEnv env,
                             GRBModel &mestre,
                             Eigen::VectorX<std::unique_ptr<GRBModel>> &vetSubProb,
                             double custoVarA,
                             int numVarSubProblema,
                             const std::vector<std::pair<int,int>> &&vetPairSubProb,
                             std::vector<int (*)(const Eigen::VectorXd &, int k, void *data, Eigen::VectorXd &vetX, int itCG,
                                                 bool &custoRedNeg)> &&vetFuncPtrSubProb)
    {

        // TODO Alterar numVarSubProblema para vetor de pairs
        // TODO Add vetor de void*, para ser os dados de cada subproblema, ou so receber void*

        std::cout<<"********************************DW DECOMP********************************\n";

        Eigen::MatrixXd matA                = getMatA_Model(mestre);
        Eigen::VectorXd vetC                = getVetC_Model(mestre);
        Eigen::VectorX<char> vetMestreSense = getConstSenseModel(mestre);
        Eigen::VectorXd vetMestreRhs        = getRhsModel(mestre);
        const int numConstrsMestre          = mestre.get(GRB_IntAttr_NumConstrs);
        const int numConstrsRmlp            = numConstrsMestre + int(vetSubProb.size());
        const int numVarMestre              = mestre.get(GRB_IntAttr_NumVars);

        auto funcGetSubMatrix = [](const Eigen::MatrixXd &matA,
                                   Eigen::MatrixXd &matSaida,
                                   size_t numLin,
                                   size_t numCol,
                                   size_t iniCol,
                                   size_t iniLin=0)
        {
            for(size_t i=iniLin; (i-iniLin) < numLin; ++i)
            {
                for(size_t j = iniCol; (j-iniCol) < numCol; ++j)
                {
                    matSaida(i-iniLin, j-iniCol) = matA(i, j);
                }
            }
        };

        Eigen::VectorX<Eigen::MatrixXd> vetSubMatA(numVarSubProblema);
        for(int i=0; i < vetSubProb.size(); ++i)
        {
            vetSubMatA[i] = Eigen::MatrixXd(numConstrsMestre, numVarSubProblema);
            vetSubMatA[i].setZero();
            funcGetSubMatrix(matA, vetSubMatA[i], numConstrsMestre, vetPairSubProb.at(i).second,
                             vetPairSubProb.at(i).first, 0);
        }


        Eigen::MatrixXd matA0(numConstrsMestre, numVarSubProblema);
        Eigen::MatrixXd matA1(numConstrsMestre, numVarSubProblema);
        matA0.setZero();
        matA1.setZero();

        std::cout<<"numVarSubProblema: "<<numVarSubProblema<<"\n\n";
        std::cout<<"numConstrsMestre: "<<numConstrsMestre<<"\n";

        funcGetSubMatrix(matA, matA0, numVarSubProblema, numVarSubProblema, 0, 0);
        funcGetSubMatrix(matA, matA1, numVarSubProblema, numVarSubProblema, numVarSubProblema, 0);

        std::cout << "Mat A:\n" << matA << "\n\nMatA1:\n"<<matA1<<"\n\n\n";
        std::cout<<"vetMat[1]: \n"<<vetSubMatA[1]<<"\n\n";


        GRBModel rmlp(env);
        GRBLinExpr linExprObjRmlp;

        GRBVar* vetVarArtifRmlp        = rmlp.addVars(numConstrsRmlp);
        GRBLinExpr* vetLinExprRmlp     = new GRBLinExpr[numConstrsRmlp];
        char* vetRmlpConstrsSense      = new char[numConstrsRmlp];
        double* vetRmlpRhs             = new double[numConstrsRmlp];
        std::string* vetStrConstrs     = new std::string[numConstrsRmlp];
        int itCG                       = 0;
        bool subProbCustR_neg          = true;
        double* vetRmlpLambda          = nullptr;
        std::vector<std::unique_ptr<Eigen::VectorXd>> vetVarLambda;             // Vetores gerados pelos subproblemas
        Eigen::VectorXd vetRmlpDuals(numConstrsRmlp);
        Eigen::VectorXd subProbCooef(numVarSubProblema+1);
        Eigen::VectorXd vetX_solSubProb(numVarSubProblema);
        GRBConstr *rmlpConstrs;
        Eigen::VectorX<GRBVar*> vetVarSubProb(vetSubProb.size());

        vetRmlpDuals.setZero();
        subProbCooef.setZero();


        std::cout << "Mat A:\n" << matA << "\n\nvetC: " << vetC.transpose() << "\n\n" << vetMestreRhs.transpose() << "\n\n";

        rmlp.setObjective(linExprObjRmlp, mestre.get(GRB_IntAttr_ModelSense));

        auto funcAdd1VerVarLambda = [&]()
        {
            vetVarLambda.push_back(std::make_unique<Eigen::VectorXd>(numVarMestre));
            vetVarLambda[vetVarLambda.size()-1]->setZero();

        };

        // Adiciona as variaveis artificiais a vetVarLambda
        for(int i=0; i < numConstrsRmlp; ++i)
            funcAdd1VerVarLambda();


        // Adiciona variaveis artificiais ao RMLP
        for(int i=0; i < numConstrsRmlp; ++i)
        {
            vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, custoVarA);
            vetLinExprRmlp[i] = vetVarArtifRmlp[i];

            if(i < numConstrsMestre)
            {
                vetRmlpConstrsSense[i] = vetMestreSense(i);
                vetRmlpRhs[i]          = vetMestreRhs(i);
                vetStrConstrs[i]       = "mestre_" + std::to_string(i);
            }
            else
            {
                vetRmlpConstrsSense[i] = GRB_EQUAL;
                vetRmlpRhs[i]          = 1;
                vetStrConstrs[i]       = "conv_"+std::to_string(i-numConstrsMestre);
            }
        }

        // Adiciona as restricoes do RMLP
        rmlpConstrs = rmlp.addConstrs(vetLinExprRmlp, vetRmlpConstrsSense, vetRmlpRhs, vetStrConstrs, numConstrsRmlp);

/*    // Pega as variaveis dos subprobelemas
    for(int p=0; p < vetSubProb.size(); ++p)
    {
        vetSubProb[p]->addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, "const");
        vetSubProb[p]->update();

        vetVarSubProb[p] = vetSubProb[p]->getVars();
        (vetVarSubProb[p])[numVarSubProblema].set(GRB_DoubleAttr_UB, 1.0);
        (vetVarSubProb[p])[numVarSubProblema].set(GRB_DoubleAttr_LB, 1.0);
    }*/


        auto funcUpdateRmlpDuals = [&]()
        {
            for(int i=0; i < numConstrsRmlp; ++i)
            {
                vetRmlpDuals(i) = rmlpConstrs[i].get(GRB_DoubleAttr_Pi);

            }

        };

        auto funcGetSubProbCooef = [&](int k)
        {
            //subProbCooef.segment(0, numVarSubProblema) = vetC.segment(k*numVarSubProblema, numVarSubProblema) - vetRmlpDuals.segment(0, numConstrsMestre);
            subProbCooef.segment(0, numVarSubProblema) = (vetC.segment(k*numVarSubProblema, numVarSubProblema).transpose() -
                                                          vetRmlpDuals.segment(0, numVarSubProblema).transpose()*vetSubMatA[k]).transpose();


            subProbCooef(numVarSubProblema) = -vetRmlpDuals(numConstrsMestre+k);

            std::cout<<"subProbCooef: "<<subProbCooef.transpose()<<"\n\n";
            std::cout<<"subProbCooef: "<<subProbCooef.transpose()<<"\n\n";
        };

        auto funcSetSubProbObj = [&](int k)
        {
            GRBModel &model = *vetSubProb[k];
            model.set(GRB_DoubleAttr_Obj, vetVarSubProb[k], &subProbCooef(0), vetSubProb[k]->get(GRB_IntAttr_NumVars));
            model.update();
            model.write("colGen_subProb_" + std::to_string(k)+"_it_"+std::to_string(itCG)+".lp");

        };

        auto funcAddColumn = [&](const Eigen::VectorXd &vetCooef, const double cost, int k)
        {
            GRBColumn grbColumn;
            grbColumn.addTerms(&vetCooef[0], rmlpConstrs, int(vetCooef.size()));
            rmlp.addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS, grbColumn,
                        "l_"+std::to_string(itCG)+"_"+std::to_string(k));
        };


        auto funcRecuperaX = [](std::vector<std::unique_ptr<Eigen::VectorXd>> &vetVarLambda,
                                const Eigen::VectorXd &vetEigenLambda,
                                Eigen::VectorXd &vetX)
        {
            if(vetVarLambda.size() != vetEigenLambda.size())
            {
                std::cout<<"Erro, vetVarLambda.size() != vetEigenLambda\n";
                std::cout<<"\nvetVarLambda.size(): "<<vetVarLambda.size()<<"\nvetEigenLambda.size(): "<<vetEigenLambda.size()<<"\n\n";
                throw "Erro";
            }

            vetX.setZero();
            for(int i=0; i < vetVarLambda.size(); ++i)
            {
                std::cout<<i<<": "<<vetVarLambda[i]->transpose()<<"\n";
                vetX = vetX + vetEigenLambda[i]*(*vetVarLambda[i]);
            }
        };

        rmlp.set(GRB_IntParam_Method, GRB_METHOD_PRIMAL);
        rmlp.set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);

        while(subProbCustR_neg)
        {
            std::cout<<"CG It: "<<itCG<<"\n\n";

            for(int i=0; i < vetVarLambda.size(); ++i)
            {
                std::cout<<"\t"<< i << ": " << vetVarLambda[i]->transpose() << "\n";
            }

            std::cout<<"\n\n";


            rmlp.update();
            rmlp.write("rmlp_"+std::to_string(itCG)+".lp");
            rmlp.optimize();


            funcUpdateRmlpDuals();

            GRBVar *vetVar = rmlp.getVars();

            // Recupera variavel lambda do RMLP
            vetRmlpLambda = rmlp.get(GRB_DoubleAttr_X, vetVar, rmlp.get(GRB_IntAttr_NumVars));

            std::cout << "PI: " << vetRmlpDuals.transpose() << "\n";
            std::cout << "Lambda: ";

            for(int i=0; i < rmlp.get(GRB_IntAttr_NumVars); ++i)
                std::cout<<vetRmlpLambda[i]<<" ";

            std::cout<<"\n\n";

            Eigen::VectorXd cgCooef(numConstrsRmlp);
            subProbCustR_neg = false;

            // Atualiza e resolve os subproblemas
            for(int k = 0; k < vetSubProb.size(); ++k)
            {

                funcGetSubProbCooef(k);
                //funcSetSubProbObj(k);
                //vetSubProb[k]->optimize();

/*            if(vetSubProb[k]->get(GRB_IntAttr_Status) != GRB_OPTIMAL)
            {
                std::cout << "Erro ao resolver o subproblema " << k << " na iteracao " << itCG << "\n";
                return;
            }*/

                //if(vetSubProb[k]->get(GRB_DoubleAttr_ObjVal) > -TolObjSubProb)
                //    continue;

                bool subProbK_CustoR_neg = false;

                vetX_solSubProb.setZero();
                vetFuncPtrSubProb.at(k)(subProbCooef, k, (void*)(&(*vetSubProb[k])), vetX_solSubProb, itCG, subProbK_CustoR_neg);

                subProbCustR_neg = subProbCustR_neg || subProbK_CustoR_neg;
                if(!subProbK_CustoR_neg)
                    continue;

                // Recupera Solucao do subproblema k
                funcAdd1VerVarLambda();

                Eigen::VectorXd &vetSol = *vetVarLambda[vetVarLambda.size()-1];

                for(int i = 0; i < numVarSubProblema; ++i)
                {
                    vetSol[k*numVarSubProblema+i] = vetX_solSubProb[i];

                }

                double cgCooefObj = vetSol.dot(vetC);
                std::cout<<"cgCooefObj "<<k<<": "<<cgCooefObj<<"\n";

                cgCooef.setZero();
                cgCooef.segment(0, numConstrsMestre) = matA*vetSol;
                cgCooef(numConstrsMestre+k) = 1;

                std::cout<<"Coef rmlp: "<<cgCooef.transpose()<<"\n\n";

                funcAddColumn(cgCooef, cgCooefObj, k);




                std::cout << "\n**************************************************\n\n\n";

            }

            std::cout<<"#####################################################\n\n";

            itCG += 1;

            if(subProbCustR_neg)
            {
                delete[]vetRmlpLambda;
                vetRmlpLambda = nullptr;
            }
        }

        // Recupera Solucao
        Eigen::VectorXd vetEigenLambda(rmlp.get(GRB_IntAttr_NumVars));
        vetEigenLambda.setZero();
        for(int i=0; i < rmlp.get(GRB_IntAttr_NumVars); ++i)
        {
            vetEigenLambda[i] = vetRmlpLambda[i];
        }

        Eigen::VectorXd vetX(numVarMestre);
        funcRecuperaX(vetVarLambda, vetEigenLambda, vetX);

        std::cout<<"Lambda: "<<vetEigenLambda.transpose()<<"\n\n";
        std::cout<<"X: "<<vetX.transpose()<<"\n\n";
        std::cout<<"z: "<<vetX.dot(vetC)<<"\n\n";


        delete []vetVarArtifRmlp;
        delete []vetLinExprRmlp;
        delete []vetRmlpConstrsSense;
        delete []vetRmlpRhs;
        delete []vetStrConstrs;
        delete []rmlpConstrs;
        delete []vetRmlpLambda;

/*    for(int p=0; p < vetSubProb.size(); ++p)
        delete []vetVarSubProb[p];*/
    }
}


#endif //DW_DECOMP_DW_DECOMP_H
