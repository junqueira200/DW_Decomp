#include "DW_Decomp.h"
#include "SparseOp.h"
#include <memory>
#include <utility>
#include <boost/unordered_set.hpp>


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

void DW_DecompNS::getSparseMatModel(GRBModel &model, Eigen::SparseMatrix<double, Eigen::RowMajor> &matA)
{
std::cout<<"getSparseMatModel ini\n\n";

    const int numVar   = model.get(GRB_IntAttr_NumVars);
    const int numConst = model.get(GRB_IntAttr_NumConstrs);

    //matA = SparseMatColD(numConst, numVar);
    matA.resize(numConst, numVar);
    matA.setZero();

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
            {   //std::cout<<"("<<c<<","<<i<<")\n";
                matA.coeffRef(c, i) = coeff;
            }
        }
    }


    delete []vetVar;
    delete []vetConstr;


std::cout<<"getSparseMatModel fim\n\n";

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

void DW_DecompNS::getVetC_Model(GRBModel &model, Eigen::RowVectorXd &vetC)
{
    model.update();

    GRBVar *vetVar = model.getVars();
    const int numVar = model.get(GRB_IntAttr_NumVars);

    vetC = Eigen::RowVectorXd(numVar);
    vetC.setZero();

    for(int i=0; i < numVar; ++i)
    {
        double val = vetVar[i].get(GRB_DoubleAttr_Obj);
        if(val != 0.0)
        {
            //std::cout<<i<<"("<<val<<"); ";
            vetC.coeffRef(0, i) = val;
        }
    }
    //std::cout<<"\n"<<vetC<<"\n\n";
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


static void funcGetSubMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor> &matA,
                             Eigen::SparseMatrix<double, Eigen::RowMajor> &matSaida,
                             int64_t numLin,
                             int64_t numCol,
                             int64_t iniCol,
                             int64_t iniLin=0)
{
    std::cout<<"funcGetSubMatrix\n";

    //matSaida = SparseMatColD(numLin, numCol);
    matSaida.resize(numLin, numCol);
    for(int64_t i=iniLin; (i-iniLin) < numLin; ++i)
    {
        for(int64_t j = iniCol; (j-iniCol) < numCol; ++j)
        {
            //std::cout<<"\t\tA("<<i<<", "<<j<<")\n";
            double val = matA.coeff(i, j);
            //std::cout<<"\t\t\tget\n";
            if(val != 0.0)
            {   //std::cout<<"\t\t\t\tacess: ("<<i-iniLin<<", "<<j-iniCol<<")\n";
                matSaida.coeffRef(i-iniLin, j-iniCol) = val;
            }
        }
    }

    //std::cout<<TempSpMatPrint(matSaida)<<"\n";
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

    const int64_t numConstrsConv = subProb->getNumConvConstr();

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
/*            subProb->resolveSubProb(vetC,
                                    vetRmlpDuals,
                                    rmlp,
                                    vetX_solSubProb,
                                    itCG,
                                    subProbK_CustoR_neg,
                                    data,
                                    numConstrsMestre,
                                    k,
                                    vetConvCooef,
                                    vetPairSubProb[k]);*/


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
                                          GRBModel &master,
                                          double costA_Var_,
                                          SubProb *ptrSubProb_,
                                          const int numSubProb_,
                                          AuxData &auxVect,
                                          Info &info):
                                                                         ptrSubProb(ptrSubProb_)
{

std::cout<<"ini DW_DecompNode\n";

    info.costA_Var  = costA_Var_;
    info.numSubProb = numSubProb_;
std::cout<<"antes\n";
    info.numConstrsConv = ptrSubProb->getNumConvConstr();
std::cout<<"depois\n";

    info.numConstrsMaster = master.get(GRB_IntAttr_NumConstrs);
    info.numVarMaster     = master.get(GRB_IntAttr_NumVars);

std::cout<<"info set\n";

    getSparseMatModel(master, matA);

//std::cout<<TempSpMatPrint(matA)<<"\n";

    getVetC_Model(master, auxVect.vetRowC);

std::cout << "vetC: "<<auxVect.vetRowC<< "\n";


    int numConstrsRmlp    = info.numConstrsMaster + ptrSubProb->getNumConvConstr();// + int(vetSubProb.size());
    info.numVarRmlpPi     = numConstrsRmlp;

std::cout<<"numVarRmlpPi: "<<info.numVarRmlpPi<<"\n";
std::cout<<"numConstrsMaster: "<<info.numConstrsMaster<<"\n\n";
std::cout<<"ptrSubProb->getNumConvConstr: "<<ptrSubProb->getNumConvConstr()<<"\n\n";

    int temp = 0;
    for(const auto &it:auxVect.vetPairSubProb)
    {
        if(it.second > temp)
            temp = it.second;
    }

    const int maxNumVarSubProb = temp;


    vetSubMatA = Vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>(info.numSubProb);
    for(int i=0; i < info.numSubProb; ++i)
    {   //std::cout<<"ini for it: "<<i<<"\n";
        //vetSubMatA[i] = SparseMatColD(numConstrsMestre, vetPairSubProb.at(i).second);
        //std::cout<<"antes funcGetSubMatrix\n";

        funcGetSubMatrix(matA,
                         vetSubMatA[i],
                         info.numConstrsMaster,
                         auxVect.vetPairSubProb.at(i).second,
                         auxVect.vetPairSubProb.at(i).first,
                         0);
    }


    uRmlp = std::make_unique<GRBModel>(env_);
    GRBLinExpr objRmlp;
    vetVarArtifRmlp = uRmlp->addVars(info.numConstrsMaster);
    vetLinExprRmlp = Vector<GRBLinExpr>(info.numConstrsMaster);

    auxVect.vetRowRmlpPi.resize(1, info.numVarRmlpPi);
    auxVect.vetRowRmlpPi.setZero();
    auxVect.vetRowRmlpSmoothPi.resize(1, info.numVarRmlpPi);
    auxVect.vetRowRmlpSmoothPi.setZero();

    auxVect.stabD.start(info.numVarRmlpPi);

    auxVect.vetColSubProbCooef.resize(maxNumVarSubProb+1, 1);
    auxVect.vetColSubProbCooef.setZero();



    //vetColX_solSubProb = SparseVectorD(maxNumVarSubProb);
    auxVect.matColX_solSubProb.resize(maxNumVarSubProb, NumMaxSolSubProb);
    auxVect.matColX_solSubProb.setZero();

    //vetRowConvCooef    = SparseVectorD(numConstrsConv);
    auxVect.vetColConvCooef.resize(info.numConstrsConv, 1);
    auxVect.vetColConvCooef.setZero();

    auxVect.vetColCooef.resize(numConstrsRmlp, 1);
    auxVect.vetColCooef.setZero();

    GRBLinExpr linExprObjRmlp;

    uRmlp->setObjective(linExprObjRmlp, master.get(GRB_IntAttr_ModelSense));

    // Add the artificial variables to vetVarLambda
    for(int i=0; i < info.numConstrsMaster + info.numConstrsConv; ++i)
        vetVarLambdaCol.emplace_back(std::make_unique<Eigen::VectorXd>(info.numVarMaster));

    for(int i=0; i < info.numConstrsMaster + info.numConstrsConv; ++i)
        vetVarLambdaCol[i]->setZero();

    // Creates the artificial variables to the RMLP
    char* vetRmlpConstrsSense      = new char[info.numConstrsMaster];
    double* vetRmlpRhs             = new double[info.numConstrsMaster];
    std::string* vetStrConstrs     = new std::string[info.numConstrsMaster];

    auto vetConstr = master.getConstrs();

    for(int i=0; i < info.numConstrsMaster; ++i)
    {
        char sense = vetConstr[i].get(GRB_CharAttr_Sense);

        if(sense != '<')
        {
            vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, info.costA_Var);
            vetLinExprRmlp[i] = vetVarArtifRmlp[i];
        }
        //if(i < numConstrsMestre)

        vetRmlpConstrsSense[i] = sense;
        vetRmlpRhs[i]          = vetConstr[i].get(GRB_DoubleAttr_RHS);
        vetStrConstrs[i]       = "master_" + std::to_string(i);
    }

    delete []vetConstr;

    ptrSubProb->iniConvConstr(*uRmlp, nullptr, info.costA_Var);

    // Adiciona as restricoes do RMLP
    vetRmlpConstr = uRmlp->addConstrs(&vetLinExprRmlp[0],
                                      vetRmlpConstrsSense,
                                      vetRmlpRhs,
                                      vetStrConstrs,
                                      info.numConstrsMaster);
    delete []vetRmlpConstr;

    //uRmlp->set(GRB_DoubleParam_OptimalityTol, 1E-9);
    uRmlp->update();


    //ptrSubProb->iniConvConstr(*uRmlp, nullptr, info.costA_Var);
    //return;

    //uRmlp->update();
    uRmlp->write("rmlp_"+std::to_string(-1)+".lp");
    vetRmlpConstr = uRmlp->getConstrs();

    uRmlp->set(GRB_IntParam_Method, GRB_METHOD_PRIMAL);
    uRmlp->set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);

    std::cout<<"ini del!\n";
    delete []vetRmlpConstrsSense;
    delete []vetRmlpRhs;
    delete []vetStrConstrs;
    std::cout<<"fim del!\n";
}

DW_DecompNS::StatusProb DW_DecompNS::DW_DecompNode::columnGeneration(AuxData &auxVect, const Info &info)
{

    bool subProbCustR_neg = true;
    itCG = 0;

    const int iniConv = 0;
    bool setAllVarAZero = false;
    int numSol = 0;
    double redCost = 0.0;
    double lagrangeDualBound = std::numeric_limits<double>::max();
    double gap = std::numeric_limits<double>::max();
    double privObjRmlp = gap;
    int numLimit = 0;

    static Eigen::Matrix<double, 2, -1, Eigen::RowMajor> matDual;
    if(matDual.cols() < info.numVarRmlpPi)
        matDual.resize(2, info.numVarRmlpPi);

    matDual.setConstant(std::numeric_limits<double>::infinity());


    while(subProbCustR_neg && numLimit < 5)
    {
        subProbCustR_neg = false;

        std::cout<<"CG It: "<<itCG<<"\n\n";

        uRmlp->update();
        // TODO del
        //uRmlp->write("rmlp_"+std::to_string(itCG)+".lp");
        uRmlp->optimize();

        std::cout<<"Val fun OBJ: "<<uRmlp->get(GRB_DoubleAttr_ObjVal)<<"\n";
        double objRmlp = uRmlp->get(GRB_DoubleAttr_ObjVal);

        if(itCG > 50)
        {
            gap = (std::abs(objRmlp-privObjRmlp)/objRmlp)*100.0;
            std::cout<<"GAP("<<gap<<"%)\n";
            if(gap <= gapLimit)
                numLimit += 1;
            else
                numLimit = 0;

            std::cout<<"numLimit: "<<numLimit<<"\n";

        }


        // TODO Remove
        GRBVar *vetVar        = uRmlp->getVars();
        double *vetRmlpLambda = uRmlp->get(GRB_DoubleAttr_X, vetVar, uRmlp->get(GRB_IntAttr_NumVars));

        std::cout<<"vetRmlpLambda: ";
        for(int i=0; i < uRmlp->get(GRB_IntAttr_NumVars); ++i)
            std::cout<<vetRmlpLambda[i]<<" ";

        std::cout<<"\n\n";


        updateRmlpPi(auxVect.vetRowRmlpPi, info);

        if(Stabilization)
        {
            //auxVect.stabD.getWeightSum(StabilizationAlpha, auxVect.vetRowRmlpSmoothPi);
            auxVect.vetRowRmlpSmoothPi *= StabilizationAlpha;
            auxVect.vetRowRmlpSmoothPi += (1.0-StabilizationAlpha) * auxVect.vetRowRmlpPi;


            std::cout << "SPI: " << auxVect.vetRowRmlpSmoothPi << "\n\n";
        }
        else
            auxVect.vetRowRmlpSmoothPi = auxVect.vetRowRmlpPi;

        // Check if the dual solution has changed

        /*
        for(int64_t i=0; i < info.numVarRmlpPi; ++i)
            matDual((itCG%2), i) = auxVect.vetRowRmlpPi[i];

        bool equal = true;

        for(int64_t i=0; i < info.numVarRmlpPi; ++i)
        {
            if(!doubleEqual(matDual((itCG%2), i), matDual(((itCG+1)%2), i)))
            {
                equal = false;
                break;
            }
        }

        if(equal)
        {
            std::cout<<"Dual solution is the same of the previous iteration\n";
            break;
        }
        */

        // Update and solve the subproblems
        for(int k=0; k < info.numSubProb; ++k)
        {
            bool subProbK_CustoR_neg = false;
            auxVect.matColX_solSubProb.setZero();
            auxVect.vetColConvCooef.setZero();

            //getSubProbCooef(k, auxVect);
            numSol = 0;
            Eigen::VectorXd vetX;

            ptrSubProb->resolveSubProb(auxVect.vetRowC,
                                       auxVect.vetRowRmlpSmoothPi,
                                       *uRmlp,
                                       itCG,
                                       subProbK_CustoR_neg,
                                       nullptr,
                                       iniConv,
                                       k,
                                       auxVect.vetColConvCooef,
                                       auxVect.vetPairSubProb[k],
                                       auxVect.matColX_solSubProb,
                                       numSol,
                                       redCost);

            if(!subProbK_CustoR_neg)
                continue;

            subProbCustR_neg = true;


            for(int l=0; l < numSol; ++l)
            {
                // Create a new lambda variable
                vetVarLambdaCol.push_back(std::make_unique<Eigen::VectorXd>(info.numVarMaster));
                Eigen::VectorXd &vetSol = *vetVarLambdaCol[vetVarLambdaCol.size()-1];
                vetSol.setZero();

                // Shift the x solution
                for(int i = 0; i < auxVect.vetPairSubProb[k].second; ++i)
                    vetSol[auxVect.vetPairSubProb[k].first + i] = auxVect.matColX_solSubProb(i, l);

                if(setVarLamdaCol.count(SolXHash(vetSol)) == 1)
                {
                    std::cout<<"miss pricing!\n";
                    PRINT_DEBUG("", "");
                    throw "MISS_PRICING";
                }
                else
                    setVarLamdaCol.emplace(vetSol);


                double cgCooefObj = vetSol.dot(auxVect.vetRowC);
                //std::cout<<"cgCooefObj "<<k<<": "<<cgCooefObj<<"\n";
                //std::cout<<vetSol.transpose()<<"\n";

                auxVect.vetColCooef.setZero();

                auxVect.vetColCooef.segment(0, info.numConstrsConv) = auxVect.vetColConvCooef;
                auxVect.vetColCooef.segment(info.numConstrsConv, info.numConstrsMaster) = matA * vetSol;

                addColumn(cgCooefObj, l, auxVect, info);
            }

            break;
        }

        objRmlp = uRmlp->get(GRB_DoubleAttr_ObjVal);
        privObjRmlp = objRmlp;
        //lagrangeDualBound = getLagrangeDualBound(objRmlp, redCost);
        //gap = (std::abs(redCost)/objRmlp)*100.0;
        //std::cout<<"GAP("<<gap<<"%)\n";



        //if(Stabilization)
        //    auxVect.stabD.addPi(auxVect.vetRowRmlpPi);

        delete []vetVar;
        delete []vetRmlpLambda;
        itCG += 1;

/*        if(itCG == 1150)
        {
            std::cout<<"CG atingiu numero max de iteracoes!\n";
            break;
        }*/

    }


    uRmlp->update();
    uRmlp->optimize();

    std::cout<<"FIM CG!\n";
    std::cout<<"Val fun OBJ: "<<uRmlp->get(GRB_DoubleAttr_ObjVal)<<"\n";
    uRmlp->write("rmlp_"+std::to_string(itCG)+".lp");

    GRBVar *vetVar        = uRmlp->getVars();
    double *vetRmlpLambda = uRmlp->get(GRB_DoubleAttr_X, vetVar, uRmlp->get(GRB_IntAttr_NumVars));

    std::cout<<"vetRmlpLambda: ";
    for(int i=0; i < uRmlp->get(GRB_IntAttr_NumVars); ++i)
        std::cout<<vetRmlpLambda[i]<<" ";

    std::cout<<"\n\n";
    delete []vetRmlpLambda;
    delete []vetVar;


    return DW_DecompNS::StatusSubProb_Otimo;


}

void DW_DecompNS::DW_DecompNode::updateRmlpPi(Eigen::RowVectorXd &vetRowRmlpPi, const Info &info)
{

std::cout<<"updateRmlpPi\n";

    //vetRowRmlpPi.resetVector();

    vetRowRmlpPi.setZero();

    for(int i=0; i < info.numVarRmlpPi; ++i)
    {
        double val = vetRmlpConstr[i].get(GRB_DoubleAttr_Pi);
        //if(val != 0.0)
            vetRowRmlpPi.coeffRef(0, i) = val;
    }

std::cout<<"PI: "<<vetRowRmlpPi<<"\n\n";

}

void DW_DecompNS::DW_DecompNode::getSubProbCooef(int k, AuxData &auxVect)
{

    const int iniVarSubProbK = auxVect.vetPairSubProb[k].first;
    const int numVarSubProbK = auxVect.vetPairSubProb[k].second;

std::cout<<"iniVarSubProbK: "<<iniVarSubProbK<<"\nnumVarSubProbK: "<<numVarSubProbK<<"\n\n";
std::cout<<"size(auxVect.vetColSubProbCooef)("<<auxVect.vetColSubProbCooef.size()<<"\n\n";
std::cout<<"size(auxVect.vetRowC)("<<auxVect.vetRowC.size()<<")\n";
std::cout<<"size(auxVect.vetRowRmlpPi)("<<auxVect.vetRowRmlpPi.size()<<"\n";
std::cout<<"size(vetSubMatA[k].numLin)("<<vetSubMatA[k].rows()<<")\n";
std::cout<<"size(vetSubMatA[k].numCol)("<<vetSubMatA[k].cols()<<")\n";
std::cout<<"vetLin: "<<auxVect.vetRowRmlpPi.segment(0, numVarSubProbK)<<"\n";

auto resul = (auxVect.vetRowRmlpPi.segment(0, numVarSubProbK)*vetSubMatA[k]);
std::cout<<"resul: \n"<<resul.transpose()<<"\n";

    auxVect.vetColSubProbCooef.segment(0, numVarSubProbK).noalias() =
            (auxVect.vetRowC.segment(iniVarSubProbK, numVarSubProbK) -
             (auxVect.vetRowRmlpPi.segment(0, numVarSubProbK)*vetSubMatA[k]));

std::cout<<"vetColSubProbCooef: "<<auxVect.vetColSubProbCooef.segment(0, numVarSubProbK).transpose()<<"\n\n";


}

void DW_DecompNS::DW_DecompNode::addColumn(const double cost, int k, AuxData &auxVect, const Info &info)
{

    GRBColumn grbColumn;
    grbColumn.addTerms(&auxVect.vetColCooef[0], vetRmlpConstr, int(auxVect.vetColCooef.size()));
    uRmlp->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS, grbColumn,
                "l_"+std::to_string(itCG)+"_"+std::to_string(k));


}

double DW_DecompNS::DW_DecompNode::getLagrangeDualBound(double objRmlp, double redCost)
{
   return objRmlp + redCost;
}

// TODO
void DW_DecompNS::AuxData::updateSizes(DW_DecompNS::DW_DecompNode &e)
{

    assertm(true, "Nao implementado!");

}

DW_DecompNS::StabilizationData::StabilizationData(int numVarRmlpPi)
{
    //matPi.resize(1000, numVarRmlpPi);
    start(numVarRmlpPi);
}

void DW_DecompNS::StabilizationData::addPi(const Eigen::RowVectorXd &vetRowRmlpPi)
{
    if((numRow+1) > matPi.rows())
        matPi.conservativeResize(2*numRow, matPi.cols());

    matPi.row(numRow) = vetRowRmlpPi;
    numRow += 1;
}

void DW_DecompNS::StabilizationData::getWeightSum(const double alpha, Eigen::RowVectorXd &vetRowRmlpPi)
{
    vetRowRmlpPi.setZero();

    if(numRow == 0)
        return;

    const int numIt = numRow+1;

    const double alpha_1 = (1.0-alpha)*alpha;

    for(int i=0; i < numRow; ++i)
    {
        const double coef = alpha_1;// * //std::pow(alpha, numIt-i);
        vetRowRmlpPi += coef*matPi.row(i);
    }

}

void DW_DecompNS::StabilizationData::start(int numVarRmlpPi)
{

    matPi.resize(1000, numVarRmlpPi);
    matPi.setZero();
    numRow = 0;

}

DW_DecompNS::SolXHash::SolXHash(const Eigen::VectorXd &vetX_):vetX(vetX_)
{

    // https://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-16.html
    // ELF Hash algorithm
    hashVal = 0;
    for(int i=0; i < vetX.size(); ++i)
    {   //         valHash * 16
        hashVal = (hashVal<<4) + vetX[i];
        uint64_t g = hashVal & 0xF0000000L;

        if(g != 0)
            hashVal ^= g >> 24;
        hashVal &= ~g;
    }

}
