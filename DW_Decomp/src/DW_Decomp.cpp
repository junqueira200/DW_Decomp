#include "DW_Decomp.h"
#include "SparseOp.h"
#include <memory>
#include <utility>

using namespace SparseOpNS;

Eigen::MatrixXd DW_DecompNS::getMatA_Model(GRBModel &mestre)
{


    const int numVar   = mestre.get(GRB_IntAttr_NumVars);
    const int numConst = mestre.get(GRB_IntAttr_NumConstrs);

    //std::cout<<numConst<<" x "<<numVar<<"\n";

    Eigen::MatrixXd matA = Eigen::MatrixXd(numConst, numVar);
    matA.setZero();

    mestre.update();

    GRBVar *vetVar = mestre.getVars();
    GRBConstr *vetConstr = mestre.getConstrs();

    for(int i=0; i < numVar; ++i)
    {

        for(int c=0; c < numConst; ++c)
        {
            matA(c, i) = mestre.getCoeff(vetConstr[c], vetVar[i]);
        }
    }


    delete []vetVar;
    delete []vetConstr;

    return std::move(matA);

}

void DW_DecompNS::getSparseMatColModel(GRBModel &model, SparseMatColD &matA)
{

    const int numVar   = model.get(GRB_IntAttr_NumVars);
    const int numConst = model.get(GRB_IntAttr_NumConstrs);

    matA = SparseMatColD(numConst, numVar);
    //std::cout<<"matA.vetVetInd[0].size("<<matA.vetVetInd[0]->size()<<")\n";
    model.update();

    GRBVar *vetVar = model.getVars();
    GRBConstr *vetConstr = model.getConstrs();

    for(int i=0; i < numVar; ++i)
    {

        for(int c=0; c < numConst; ++c)
        {
            double coeff = model.getCoeff(vetConstr[c], vetVar[i]);

            if(coeff != 0.0)
            {   std::cout<<"("<<c<<","<<i<<")\n";
                matA(c, i) = coeff;
            }
        }
    }


    delete []vetVar;
    delete []vetConstr;



}

Eigen::VectorXd DW_DecompNS::getVetC_Model(GRBModel &mestre)
{

    const int numVar = mestre.get(GRB_IntAttr_NumVars);
    Eigen::VectorXd vetC(numVar);
    GRBVar *vetVar = mestre.getVars();

    for(int i=0; i < numVar; ++i)
        vetC[i] = vetVar[i].get(GRB_DoubleAttr_Obj);

    //std::cout<<vetC.transpose()<<"\n";

    delete []vetVar;

    return std::move(vetC);

}

void DW_DecompNS::getVetC_Model(GRBModel &model, SparseVectorD &vetC)
{
    GRBVar *vetVar = model.getVars();
    const int numVar = model.get(GRB_IntAttr_NumVars);
    if(vetC.n != numVar)
        vetC = SparseVectorD(numVar);

    for(int i=0; i < numVar; ++i)
    {
        double val = vetVar[i].get(GRB_DoubleAttr_Obj);
        if(val != 0.0)
            vetC(i) = val;
    }

    delete []vetVar;
}

Eigen::VectorX<char> DW_DecompNS::getConstSenseModel(GRBModel &model)
{

    const int NumConstrs = model.get(GRB_IntAttr_NumConstrs);
    Eigen::VectorX<char> vetSense(NumConstrs);

    auto vetConstr = model.getConstrs();

    for(int c=0; c < NumConstrs; ++c)
    {
        vetSense[c] = vetConstr[c].get(GRB_CharAttr_Sense);
    }

    delete []vetConstr;

    return std::move(vetSense);

}

Eigen::VectorXd DW_DecompNS::getRhsModel(GRBModel &model)
{


    const int NumConstrs = model.get(GRB_IntAttr_NumConstrs);
    Eigen::VectorXd vetRhs(NumConstrs);

    auto vetConstr = model.getConstrs();

    for(int c=0; c < NumConstrs; ++c)
    {
        vetRhs[c] = vetConstr[c].get(GRB_DoubleAttr_RHS);
    }

    delete []vetConstr;

    return std::move(vetRhs);

}


void DW_DecompNS::recuperaX(GRBVar* var, Eigen::VectorXd &vetX, int numVar)
{

    for(int i = 0; i < numVar; ++i)
        vetX[i] = var[i].get(GRB_DoubleAttr_X);

}



static void funcGetSubMatrix(const Eigen::MatrixXd &matA,
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


static void funcGetSubMatrix(const SparseMatColD &matA,
                             SparseMatColD &matSaida,
                             int64_t numLin,
                             int64_t numCol,
                             int64_t iniCol,
                             int64_t iniLin=0)
{
    matSaida = SparseMatColD(numLin, numCol);

    for(int64_t i=iniLin; (i-iniLin) < numLin; ++i)
    {
        for(int64_t j = iniCol; (j-iniCol) < numCol; ++j)
        {
            //std::cout<<"\t\tA("<<i<<", "<<j<<")\n";
            double val = matA.getVal(i, j);
            //std::cout<<"\t\t\tget\n";
            if(val != 0.0)
            {   //std::cout<<"\t\t\t\tacess: ("<<i-iniLin<<", "<<j-iniCol<<")\n";
                matSaida(i - iniLin, j - iniCol) = val;
            }
        }
    }

    std::cout<<matSaida<<"\n";
};


/** ****************************************************************************************************
 * ****************************************************************************************************
 *
 * @param env                   Env do Gurobi
 * @param mestre                Modelo com as restricoes do mestre e a func obj do problema original
 * @param custoVarA             Custo das variaveis artificiais, para inicializar a CG
 * @param vetPairSubProb        Vetor de pair por subproblema, cada pair representa o indice do inicio das variaveis e o seu tamanho
 * @param subProb               Ptr para class SubProb
 * @param data                  Ptr para void* eh passado para os metodos de @param subProb
 * @param numConstrsConv        Numero de restricoes de convexidade do RMP
 * @param numSubProb            Numero de subproblemas
 *
 * ****************************************************************************************************
 * *****************************************************************************************************/
void DW_DecompNS::dwDecomp(GRBEnv &env,
                           GRBModel &mestre,
                           double custoVarA,
                           const std::vector<std::pair<int, int>> &&vetPairSubProb,
                           SubProb *subProb,
                           void *data,
                           const int numSubProb)
{

    const int64_t numConstrsConv = subProb->getNumberOfConvConstr();

    std::cout<<"********************************DW DECOMP********************************\n";

    Eigen::MatrixXd matA                = getMatA_Model(mestre);
    Eigen::VectorXd vetC                = getVetC_Model(mestre);
    Eigen::VectorX<char> vetMestreSense = getConstSenseModel(mestre);
    Eigen::VectorXd vetMestreRhs        = getRhsModel(mestre);
    const int numConstrsMestre          = mestre.get(GRB_IntAttr_NumConstrs);
    int numConstrsRmlp                  = numConstrsMestre;// + int(vetSubProb.size());
    const int numVarMestre              = mestre.get(GRB_IntAttr_NumVars);

    int temp = 0;
    for(const auto &it:vetPairSubProb)
    {
        if(it.second > temp)
            temp = it.second;
    }

    const int maxNumVarSubProb         = temp;


    Eigen::VectorX<Eigen::MatrixXd> vetSubMatA(numSubProb);
    for(int i=0; i < numSubProb; ++i)
    {
        vetSubMatA[i] = Eigen::MatrixXd(numConstrsMestre, vetPairSubProb.at(i).second);
        vetSubMatA[i].setZero();
        funcGetSubMatrix(matA, vetSubMatA[i], numConstrsMestre, vetPairSubProb.at(i).second,
                         vetPairSubProb.at(i).first, 0);
    }




    std::cout<<"numVarSubProblema: "<<maxNumVarSubProb<<"\n\n";
    std::cout<<"numConstrsMestre: "<<numConstrsMestre<<"\n";

    std::cout << "Mat A:\n" << matA << "\n\n";
    std::cout<<"vetMat[1]: \n"<<vetSubMatA[1]<<"\n\n";


    GRBModel rmlp(env);
    GRBLinExpr linExprObjRmlp;

    GRBVar* vetVarArtifRmlp        = rmlp.addVars(numConstrsMestre);
    GRBLinExpr* vetLinExprRmlp     = new GRBLinExpr[numConstrsMestre];
    char* vetRmlpConstrsSense      = new char[numConstrsMestre];
    double* vetRmlpRhs             = new double[numConstrsMestre];
    std::string* vetStrConstrs     = new std::string[numConstrsMestre];
    int itCG                       = 0;
    bool subProbCustR_neg          = true;
    double* vetRmlpLambda          = nullptr;
    std::vector<std::unique_ptr<Eigen::VectorXd>> vetVarLambda;             // Vetores gerados pelos subproblemas
    Eigen::VectorXd vetRmlpDuals(numConstrsMestre);
    Eigen::VectorXd subProbCooef(maxNumVarSubProb+1);
    Eigen::VectorXd vetX_solSubProb(maxNumVarSubProb);
    GRBConstr *rmlpConstrs;
    //Eigen::VectorX<GRBVar*> vetVarSubProb(vetSubProb.size());
    Eigen::VectorXd vetConvCooef(numConstrsConv);

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
    for(int i=0; i < numConstrsMestre+numConstrsConv; ++i)
        funcAdd1VerVarLambda();


    // Cria as variaveis artificiais ao RMLP
    for(int i=0; i < numConstrsMestre; ++i)
    {
        vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, custoVarA);
        vetLinExprRmlp[i] = vetVarArtifRmlp[i];

        //if(i < numConstrsMestre)

        vetRmlpConstrsSense[i] = vetMestreSense(i);
        vetRmlpRhs[i]          = vetMestreRhs(i);
        vetStrConstrs[i]       = "mestre_" + std::to_string(i);
    }

    // Adiciona as restricoes do RMLP
    rmlpConstrs = rmlp.addConstrs(vetLinExprRmlp, vetRmlpConstrsSense, vetRmlpRhs, vetStrConstrs, numConstrsMestre);
    delete []rmlpConstrs;
    rmlp.update();

    subProb->iniConvConstr(rmlp, data, custoVarA);
    //return;

    rmlp.update();
    rmlpConstrs = rmlp.getConstrs();

    auto funcUpdateRmlpDuals = [&]()
    {
        for(int i=0; i < numConstrsMestre; ++i)
        {
            vetRmlpDuals(i) = rmlpConstrs[i].get(GRB_DoubleAttr_Pi);

        }

    };

    auto funcGetSubProbCooef = [&](int k)
    {
        const int iniVarSubProbK = vetPairSubProb.at(k).first;
        const int numVarSubProbK = vetPairSubProb.at(k).second;

        subProbCooef.segment(0, numVarSubProbK) =
                                                       (vetC.segment(iniVarSubProbK, numVarSubProbK) -
                                                       (vetRmlpDuals.segment(0, numVarSubProbK).transpose()*vetSubMatA[k]).transpose()).transpose();

    };

    auto funcAddColumn = [&](Eigen::VectorXd &vetCooef, const Eigen::VectorXd &vetConvCooef, const double cost, int k)
    {
        GRBColumn grbColumn;
        std::cout<<"ini funcAddColumn\n";
        std::cout<<"vetConvCoorf.size(): "<<vetConvCooef.size()<<"\n";

        for(int i=0; i < vetConvCooef.size(); ++i)
            vetCooef[numConstrsMestre+i] = vetConvCooef[i];

        std::cout<<"Fim for\n";
        std::cout<<"vetConvCoorf: "<<vetConvCooef.transpose()<<"\n";

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


        if(itCG == 8)
            exit(-1);

        Eigen::VectorXd cgCooef(numConstrsRmlp+numConstrsConv);
        subProbCustR_neg = false;

        // Atualiza e resolve os subproblemas
        for(int k = 0; k < numSubProb; ++k)
        {

            funcGetSubProbCooef(k);
            bool subProbK_CustoR_neg = false;

            vetX_solSubProb.setZero();
            vetConvCooef.setZero();
            subProb->resolveSubProb(subProbCooef,
                                    rmlp,
                                    vetX_solSubProb,
                                    itCG,
                                    subProbK_CustoR_neg,
                                    data,
                                    numConstrsMestre,
                                    k,
                                    vetConvCooef,
                                    vetPairSubProb[k]);


            subProbCustR_neg = subProbCustR_neg || subProbK_CustoR_neg;
            if(!subProbK_CustoR_neg)
                continue;

            // Add a var lambda da solucao do subproblema k
            funcAdd1VerVarLambda();

            Eigen::VectorXd &vetSol = *vetVarLambda[vetVarLambda.size()-1];

            for(int i = 0; i < vetPairSubProb[k].second; ++i)
            {
                vetSol[vetPairSubProb[k].first+i] = vetX_solSubProb[i];

            }

            double cgCooefObj = vetSol.dot(vetC);
            std::cout<<"cgCooefObj "<<k<<": "<<cgCooefObj<<"\n";

            cgCooef.setZero();
            cgCooef.segment(0, numConstrsMestre) = matA*vetSol;
            funcAddColumn(cgCooef, vetConvCooef, cgCooefObj, k);


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

    // Recupera Solucao X
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

} // FIM DW_Decomp::dwDecomp


DW_DecompNS::DW_DecompNode::DW_DecompNode(GRBEnv &env_,
                                          GRBModel &master_,
                                          double costA_Var_,
                                          Vector<std::pair<int, int>> vetPairSubProb_,
                                          SubProb *ptrSubProb_,
                                          const int numSubProb_): ptrEnv(&env_),
                                                                  ptrMaster(&master_),
                                                                  costA_Var(costA_Var_),
                                                                  vetPairSubProb(std::move(vetPairSubProb_)),
                                                                  ptrSubProb(ptrSubProb_),
                                                                  numSubProb(numSubProb_)
{
    numConstrsConv = ptrSubProb->getNumberOfConvConstr();
    getSparseMatColModel(*ptrMaster, matA);

std::cout<<matA<<"\n";

    getVetC_Model(*ptrMaster, vetC);

std::cout<<"vetC: "    <<vetC<<"\n";


    const int numConstrsMestre          = ptrMaster->get(GRB_IntAttr_NumConstrs);
    int numConstrsRmlp                  = numConstrsMestre;// + int(vetSubProb.size());
    const int numVarMaster              = ptrMaster->get(GRB_IntAttr_NumVars);

    int temp = 0;
    for(const auto &it:vetPairSubProb)
    {
        if(it.second > temp)
            temp = it.second;
    }

    const int maxNumVarSubProb          = temp;


    vetSubMatA = Vector<SparseMatColD>(numSubProb);
    for(int i=0; i < numSubProb; ++i)
    {   std::cout<<"ini for it: "<<i<<"\n";
        //vetSubMatA[i] = SparseMatColD(numConstrsMestre, vetPairSubProb.at(i).second);
        //std::cout<<"antes funcGetSubMatrix\n";

        funcGetSubMatrix(matA, vetSubMatA[i], numConstrsMestre, vetPairSubProb.at(i).second,
                         vetPairSubProb.at(i).first, 0);
    }


    uRmlp = std::make_unique<GRBModel>(env_);
    GRBLinExpr objRmlp;
    vetVarArtifRmlp = uRmlp->addVars(numConstrsMestre);
    vetLinExprRmlp = Vector<GRBLinExpr>(numConstrsMestre);

    vetRmlpPi       = SparseVectorD(numConstrsMestre);
    vetSubProbCooef = SparseVectorD(maxNumVarSubProb+1);
    vetX_solSubProb = SparseVectorD(maxNumVarSubProb);
    vetConvCooef    = SparseVectorD(numConstrsConv);

    GRBLinExpr linExprObjRmlp;

    uRmlp->setObjective(linExprObjRmlp, ptrMaster->get(GRB_IntAttr_ModelSense));

    // Add the artificial variables to vetVarLambda
    for(int i=0; i < numConstrsMestre+numConstrsConv; ++i)
        vetVarLambda.emplace_back(std::make_unique<SparseVectorD>(numVarMaster));

    // Creates the artificial variables to the RMLP
    char* vetRmlpConstrsSense      = new char[numConstrsMestre];
    double* vetRmlpRhs             = new double[numConstrsMestre];
    std::string* vetStrConstrs     = new std::string[numConstrsMestre];

    auto vetConstr = ptrMaster->getConstrs();

    for(int i=0; i < numConstrsMestre; ++i)
    {
        vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, costA_Var);
        vetLinExprRmlp[i] = vetVarArtifRmlp[i];

        //if(i < numConstrsMestre)

        vetRmlpConstrsSense[i] = vetConstr[i].get(GRB_CharAttr_Sense);
        vetRmlpRhs[i]          = vetConstr[i].get(GRB_DoubleAttr_RHS);
        vetStrConstrs[i]       = "master_" + std::to_string(i);
    }

    delete []vetConstr;

    // Adiciona as restricoes do RMLP
    vetRmlpConstr = uRmlp->addConstrs(&vetLinExprRmlp[0], vetRmlpConstrsSense, vetRmlpRhs, vetStrConstrs, numConstrsMestre);
    delete []vetRmlpConstr;
    uRmlp->update();


    ptrSubProb->iniConvConstr(*uRmlp, nullptr, costA_Var);
    //return;

    uRmlp->update();
    vetRmlpConstr = uRmlp->getConstrs();

    uRmlp->set(GRB_IntParam_Method, GRB_METHOD_PRIMAL);
    uRmlp->set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);

    std::cout<<"ini del!\n";
    delete []vetRmlpConstrsSense;
    delete []vetRmlpRhs;
    delete []vetStrConstrs;
    std::cout<<"fim del!\n";
}

DW_DecompNS::StatusProb DW_DecompNS::DW_DecompNode::columnGeneration()
{

    return DW_DecompNS::StatusSubProb_Otimo;


}

void DW_DecompNS::DW_DecompNode::funcUpdateRmlpPi()
{

    vetRmlpPi.resetVector();

    for(int i=0; i < uRmlp->get(GRB_IntAttr_NumConstrs); ++i)
    {
        double val = vetRmlpConstr[i].get(GRB_DoubleAttr_Pi);
        if(val != 0.0)
            vetRmlpPi(i) = val;
    }


}

void DW_DecompNS::DW_DecompNode::funcGetSubProbCooef(int k)
{

    multLinTimesSparseMatCol(vetRmlpPi, vetSubMatA[k], vetSubProbCooef);


}
