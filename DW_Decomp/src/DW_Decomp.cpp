/*  *****************************************************************
 *  *****************************************************************
 *  File:    DW_Decomp.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    18/09/24
 *
 *  *****************************************************************
 *  *****************************************************************/

#include "DW_Decomp.h"
#include "SparseOp.h"
#include <memory>
#include <utility>
#include <boost/unordered_set.hpp>
#include <format>

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

    return matA;

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

    return vetC;

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

    return vetSense;

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

    return vetRhs;

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


    //std::cout << "Mat A:\n" << matA << "\n\nvetC: " << vetC.transpose() << "\n\n" << vetMestreRhs.transpose() << "\n\n";

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
        if((int)vetVarLambda.size() != (int)vetEigenLambda.size())
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

        for(int i=0; i < (int)vetVarLambda.size(); ++i)
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


DW_DecompNS::DW_DecompNode::DW_DecompNode(GRBEnv&   env_,
                                          GRBModel& master,
                                          double    costA_Var_,
                                          SubProb*  ptrSubProb_,
                                          const int numSubProb_,
                                          AuxData&  auxVect):
                                                                         ptrSubProb(ptrSubProb_)
{
    vetVar0 = VectorI();

std::cout<<"ini DW_DecompNode\n";

    info.costA_Var                = costA_Var_;
    info.numSubProb               = numSubProb_;
    info.numConstrsConv           = ptrSubProb->getNumConvConstr();
    info.numConstrsMaster         = master.get(GRB_IntAttr_NumConstrs);
    info.numConstrsOrignalProblem = info.numConstrsMaster;
    info.numVarMaster             = master.get(GRB_IntAttr_NumVars);

    getSparseMatModel(master, matA);
    getVetC_Model(master, auxVect.vetRowC);

//std::cout << "vetC: "<<auxVect.vetRowC<< "\n";


    int numConstrsRmlp    = info.numConstrsMaster + ptrSubProb->getNumConvConstr();// + int(vetSubProb.size());
    info.numVarRmlpPi     = numConstrsRmlp;

//std::cout<<"numVarRmlpPi: "<<info.numVarRmlpPi<<"\n";
//std::cout<<"numConstrsMaster: "<<info.numConstrsMaster<<"\n\n";
//std::cout<<"ptrSubProb->getNumConvConstr: "<<ptrSubProb->getNumConvConstr()<<"\n\n";

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
    Vector<GRBLinExpr> vetLinExprRmlp(info.numConstrsMaster);

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
    auxVect.auxVetCooef = new double[numConstrsRmlp];

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
            if(this->phaseStatus == PhaseStatus::PhaseStatusBigM)
                vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, info.costA_Var);
            else
                vetVarArtifRmlp[i].set(GRB_DoubleAttr_Obj, 1);

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

    vetSolX.resize(info.numVarMaster);
    vetSolX.setZero();

    //uRmlp->update();
    uRmlp->write("rmlp_"+std::to_string(-1)+".lp");
    vetRmlpConstr = uRmlp->getConstrs();

    uRmlp->set(GRB_IntParam_Method, GRB_METHOD_PRIMAL);
    //uRmlp->set(GRB_IntParam_Presolve, GRB_PRESOLVE_OFF);
    uRmlp->set(GRB_IntParam_OutputFlag, 0);
    uRmlp->set(GRB_IntParam_Quad, 1);


    delete []vetRmlpConstrsSense;
    delete []vetRmlpRhs;
    delete []vetStrConstrs;
    delete []vetRmlpConstr;
    delete []vetVarArtifRmlp;

    vetRmlpConstr = nullptr;
    vetVarArtifRmlp = nullptr;


    uRmlp->update();
    vetRmlpConstr = uRmlp->getConstrs();


}

DW_DecompNS::StatusProb DW_DecompNS::DW_DecompNode::columnGeneration(AuxData &auxVect)
{

std::cout<<"*******************Column Generation*******************\n\n";


    std::cout<<"rhsConv: "<<rhsConv<<"\n";
    bool subProbCustR_neg = true;
    itCG = 0;

    const int iniConv = 0;
    //bool setAllVarAZero = false;
    int numSol = 0;
    //double redCost = 0.0;
    Eigen::Array<double, 1, NumMaxSolSubProb> vetRedCost;
    double lagrangeDualBound = std::numeric_limits<double>::infinity();
    double gap = std::numeric_limits<double>::infinity();
    double privObjRmlp = gap;
    int numLimit = 0;
    bool missPricing = false;
    bool exactPi = false;
    bool exactPricing = false;

    uRmlp->update();

    if(vetRmlpConstr != nullptr)
        delete []vetRmlpConstr;


    vetRmlpConstr = uRmlp->getConstrs();

    //PhaseStatus phaseStatus = PhaseStatus::PhaseStatusBigM;
    VectorD vetMult(info.numConstrsOrignalProblem+info.numConstrsConv, 1.0);
    VectorI vetLastIt(info.numConstrsOrignalProblem+info.numConstrsConv, -1);    // The last iteration the artificial variable is different from 0


    if(uRmlp->get(GRB_IntAttr_NumVars) != (int)vetVarLambdaCol.size())
    {
        std::cout<<"ERROR\n";
        std::cout<<"LP: "<<uRmlp->get(GRB_IntAttr_NumVars)<<"; vet: "<<vetVarLambdaCol.size()<<"\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }


    while(subProbCustR_neg)
    {
        //std::cout<<"ItCG("<<itCG<<")\n\n";
        if(!missPricing)
            uRmlp->update();
        subProbCustR_neg = false;

        if(uRmlp->get(GRB_IntAttr_NumVars) != (int)vetVarLambdaCol.size())
        {
            std::cout<<"ERROR\n";
            std::cout<<"LP: "<<uRmlp->get(GRB_IntAttr_NumVars)<<"; vet: "<<vetVarLambdaCol.size()<<"\n";
            PRINT_DEBUG("", "");
            throw "ERROR";
        }

        //std::cout<<"LP: "<<uRmlp->get(GRB_IntAttr_NumVars)<<"; vet: "<<vetVarLambdaCol.size()<<"\n";

        if(!missPricing)
            uRmlp->optimize();

        //uRmlp->write("rmlp_"+std::to_string(itCG)+".lp");

        if(uRmlp->get(GRB_IntAttr_Status) != GRB_OPTIMAL)
            return StatusSubProb_Inviavel;

        //std::cout<<"Val fun OBJ: "<<uRmlp->get(GRB_DoubleAttr_ObjVal)<<"\n";
        double objRmlp = uRmlp->get(GRB_DoubleAttr_ObjVal);

            //std::cout<<"GAP("<<gap<<"%)\n";
        /*
            if(gap <= gapLimit)
                numLimit += 1;
            else
                numLimit = 0;
        */
            //std::cout<<"numLimit: "<<numLimit<<"\n";


        GRBVar *vetVar = uRmlp->getVars();

        if(phaseStatus == PhaseStatus::PhaseStatusBigM && itCG > 0)
        {
            bool allZero  = true;
            bool increase = true;

            int numTotal   = 0;
            int numReached = 0;

            for(int i = 0; i < info.numConstrsOrignalProblem + info.numConstrsConv; ++i)
            {
                if(!doubleEqual(vetVar[i].get(GRB_DoubleAttr_X), 0.0, 1E-5))
                {
                        numTotal += 1;

                        if(vetMult[i] == BigM_maxMult)
                        {
                            numReached += 1;
                            continue;
                        }

                        vetMult[i] += 1.0;

                        std::cout<<"set "<<i<<"_"<<vetVar[i].get(GRB_StringAttr_VarName)<<"("<<vetVar[i].get(GRB_DoubleAttr_X)<<") variable to: "<<vetMult[i]*info.costA_Var<<"\n";
                        vetVar[i].set(GRB_DoubleAttr_Obj, vetMult[i] * info.costA_Var);
                        break;

                }
                else
                {
                    vetVar[i].set(GRB_DoubleAttr_LB, 0.0);
                    vetVar[i].set(GRB_DoubleAttr_UB, 0.0);
                }

            }

            //std::printf("numTotal(%d); numReached(%d)\n\n", numTotal, numReached);
            if(numTotal > 0 && numTotal == numReached)
            {
                std::cout<<"numTotal("<<numTotal<<"), numReached("<<numReached<<")\n";
                phaseStatus = PhaseStatus::PhaseStatusTwoPhase;
                std::cout << "Start PhaseStatusTwoPhase\n";
            }
        }

        else if(phaseStatus == PhaseStatus::PhaseStatusTwoPhase)
        {
            GRBLinExpr linExpr;

            int numZero = 0;
            for(int i = 0; i < info.numConstrsOrignalProblem + info.numConstrsConv; ++i)
            {
                if(!doubleEqual(vetVar[i].get(GRB_DoubleAttr_X), 0.0, 1E-5))
                {
                    linExpr += vetVar[i];
                    std::cout<<i<<"("<<vetVar[i].get(GRB_DoubleAttr_X)<<"); ";
                }
                else
                {
                    vetVar[i].set(GRB_DoubleAttr_LB, 0.0);
                    vetVar[i].set(GRB_DoubleAttr_UB, 0.0);
                    numZero += 1;
                }
            }

            if(numZero == info.numConstrsOrignalProblem + info.numConstrsConv)
            {
                std::cout<<"All artificial variables are equal to zero!\nChanging to Column Genration phase!\n\n";
                delete []vetVar;

                uRmlp->update();
                vetVar = uRmlp->getVars();

                // Recreate the obj function with the distance of the columns
                GRBLinExpr newObj;
                int ini = (int) info.numConstrsOrignalProblem + info.numConstrsConv;

                if(uRmlp->get(GRB_IntAttr_NumVars) != (int) vetVarLambdaCol.size())
                {
                    std::cout << "ERROR\n";
                    std::cout << "LP: " << uRmlp->get(GRB_IntAttr_NumVars) << "; vet: " << vetVarLambdaCol.size()
                              << "\n";
                    PRINT_DEBUG("", "");
                    throw "ERROR";
                }

                for(int i = ini; i < (int) vetVarLambdaCol.size(); ++i)
                {
                    Eigen::VectorXd &vetSol = *vetVarLambdaCol[i];
                    double cgCooefObj = vetSol.dot(auxVect.vetRowC);
                    newObj += cgCooefObj * vetVar[i];

                }

                uRmlp->setObjective(newObj, GRB_MINIMIZE);
                phaseStatus = PhaseStatus::PhaseStatusColGen;
                subProbCustR_neg = true;
                //auxVect.vetRowRmlpSmoothPi.setZero();
                exactPi = true;
                continue;

            }

            //std::cout<<"\n";
            uRmlp->setObjective(linExpr, GRB_MINIMIZE);
        }

        delete []vetVar;


        //GRBVar *vetVar        = uRmlp->getVars();
        //double *vetRmlpLambda = uRmlp->get(GRB_DoubleAttr_X, vetVar, uRmlp->get(GRB_IntAttr_NumVars));
        if(PrintDebug)
            std::cout<<"UpdatePi\n";

        updateRmlpPi(auxVect.vetRowRmlpPi);

        if(PrintDebug)
            std::cout<<"~UpdatePi\n";

        if(Stabilization && !exactPi)// && phaseStatus == PhaseStatus::PhaseStatusColGen)
        {
            auxVect.vetRowRmlpSmoothPi = (1.0-StabilizationAlpha)*auxVect.vetRowRmlpSmoothPi +
                                          StabilizationAlpha * (auxVect.vetRowRmlpPi);
            //std::cout << "SPI: " << auxVect.vetRowRmlpSmoothPi << "\n\n";
        }
        else
        {
            auxVect.vetRowRmlpSmoothPi = auxVect.vetRowRmlpPi;
            numLimit = 0;
            //std::cout << "PI: " << auxVect.vetRowRmlpSmoothPi << "\n\n";

            //if(exactPi)
            //    std::cout<<"exactPi\n";
        }

        // **********************************************************************************************

        /*
        phaseStatus = PhaseStatus::PhaseStatusTwoPhase;
        PRINT_DEBUG("", "");
        std::cout<<"Seting Pi\n";
        auxVect.vetRowRmlpSmoothPi.setZero();
        auxVect.vetRowRmlpSmoothPi[6] = 1.0;
        std::cout << "PI: " << auxVect.vetRowRmlpSmoothPi << "\n\n";
        //exit(-1);
        */
        // **********************************************************************************************
        double constVal = 0;

        if(PrintDebug)
            std::cout<<"ConstVal conv\n";

        for(int i=0; i < (int)info.numConstrsConv; ++i)
        {
            constVal += -auxVect.vetRowRmlpSmoothPi[i];
            if(PrintDebug)
                std::cout<<auxVect.vetRowRmlpSmoothPi[i]<<" ";
        }

        if(PrintDebug)
            std::cout<<"\n\n~ConstVal conv\n";

        if(!doubleEqual(constVal, 0.0) && PrintDebug)
            std::cout<<"ConstVal: "<<constVal<<"\n";

        if(PrintDebug)
        {
            std::cout << "\n\n~ConstVal Branching\n\n";
            if(!doubleEqual(constVal, 0.0))
                std::cout << "*ConstVal: " << constVal << "\n";
        }

        // Update and solve the subproblems
        for(int k=0; k < info.numSubProb; ++k)
        {
            bool subProbK_CustoR_neg = false;
            auxVect.matColX_solSubProb.setZero();
            auxVect.vetColConvCooef.setZero();

            //getSubProbCooef(k, auxVect);
            numSol = 0;
            Eigen::VectorXd vetX;

            ptrSubProb->resolveSubProb(auxVect.vetRowC, auxVect.vetRowRmlpSmoothPi, *uRmlp, itCG, subProbK_CustoR_neg,
                                       nullptr, iniConv, k, auxVect.vetColConvCooef, auxVect.vetPairSubProb[k],
                                       auxVect.matColX_solSubProb, numSol, vetRedCost, constVal, vetVar0, vetVar1,
                                       phaseStatus, exactPricing);

            if(!subProbK_CustoR_neg)
                continue;

            double minRedCost = std::numeric_limits<double>::infinity();

            subProbCustR_neg = true;
            int numSolRep = 0;

            for(int l=0; l < numSol; ++l)
            {
                if(PrintDebug)
                    std::cout<<"Processando "<<l<<" Coluna\n";
                // Create a new lambda variable
                vetVarLambdaCol.push_back(std::make_unique<Eigen::VectorXd>(info.numVarMaster));
                Eigen::VectorXd& vetSol = *vetVarLambdaCol[vetVarLambdaCol.size()-1];
                vetSol.setZero();

                if(PrintDebug)
                    std::cout<<"Shift X\n";

                // Shift the x solution
                for(int i = 0; i < auxVect.vetPairSubProb[k].second; ++i)
                    vetSol[auxVect.vetPairSubProb[k].first+i] = auxVect.matColX_solSubProb(i, l);

                if(PrintDebug)
                    std::cout<<"~Shift X\n";

                if(setVarLamdaCol.count(SolXHash(vetSol)) == 1)
                {
                    numSolRep += 1;
                    vetVarLambdaCol.pop_back();

                    if(!Stabilization)
                    {
                        uRmlp->write("missPricing.lp");
                        PRINT_DEBUG("", "");
                        throw "MISS_PRICING";
                    }

                    continue;
                }
                else
                    setVarLamdaCol.emplace(vetSol);

                minRedCost = std::min(minRedCost, vetRedCost[l]);

                if(PrintDebug)
                    std::cout<<"cgCooefObj\n";
                double cgCooefObj = 0.0;//vetSol.dot(auxVect.vetRowC);

                if(phaseStatus == PhaseStatus::PhaseStatusTwoPhase)
                    cgCooefObj = 0.0;
                else
                    cgCooefObj = vetSol.dot(auxVect.vetRowC);

                if(PrintDebug)
                    std::cout<<"~cgCooefObj\n";

                auxVect.vetColCooef.setZero();
                //std::cout<<"numConstrsConv: "<<info.numConstrsConv<<"\n";
                if(info.numConstrsConv > 0)
                    auxVect.vetColCooef.segment(0, info.numConstrsConv) = auxVect.vetColConvCooef;

                if(PrintDebug)
                {
                    std::cout << "A*X\n";
                    std::cout << "\tA(" << matA.rows() << "," << matA.cols() << ")\n";
                    std::cout << "\tvetSol: " << vetSol.size() << "\n";
                    std::cout << "\tinfo.numConstrsMaster: " << (info.numConstrsMaster) << "\n";
                    std::cout << "\tnumConstrsConv: " << info.numConstrsConv << "\n";
                    std::cout << "\tauxVect.vetColCooef.size: " << auxVect.vetColCooef.size() << "\n";
                }

                auxVect.vetColCooef.segment(info.numConstrsConv, (info.numConstrsMaster)) = matA * vetSol;

                /*
                if(phaseStatus != PhaseStatus::PhaseStatusColGen)
                {
                    //uRmlp->update();

                    vetVar = uRmlp->getVars();

                    for(int i=1; i < info.numVarRmlpPi; ++i)
                    {
                        if(!doubleEqual(auxVect.vetColCooef[i], 0.0))
                        {
                            vetVar[i-1].set(GRB_DoubleAttr_UB, 0.0);
                            vetVar[i-1].set(GRB_DoubleAttr_LB, 0.0);

                            std::cout<<"seting var("<<vetVar[i-1].get(GRB_StringAttr_VarName)<<") to 0\n";
                        }
                    }
                }
                */

                //std::cout<<"\n";

                if(PrintDebug)
                    std::cout<<"~A*X\n\n";
                //std::cout<<"cooef: \n"<<auxVect.vetColCooef<<"\n\n";


                addColumn(cgCooefObj, l, auxVect);
            }

            if(numSolRep == numSol)
            {
                //std::cout<<"MISS PRICING\n";
                missPricing = true;
                numLimit += 1;
                lagrangeDualBound = std::numeric_limits<double>::infinity();
                gap =  std::numeric_limits<double>::infinity();
            }
            else
            {
                missPricing = false;

                lagrangeDualBound = objRmlp + rhsConv*minRedCost;
                gap = (std::abs(rhsConv*minRedCost)/objRmlp)*100.0;

                /*
                if(exactPricing && gap <= GapTolStop)
                {
                    subProbCustR_neg = false;
                    std::cout<<"Stoping CG by gap tolerance\n";
                }
                */

                if(exactPricing)
                    std::cout<<"*";

                if(gap <= GapExactPricing && phaseStatus == PhaseStatus::PhaseStatusColGen)
                {
                    //exactPricing = true;
                    //std::cout<<"exactPricing\n";
                }
                else
                    exactPricing = false;
            }

            break;
        }

        if(!missPricing && !subProbCustR_neg)
            std::cout<<"END!\n";

        if(!subProbCustR_neg)
            std::cout<<"Sub Problem have a positive value!\n";

        //uRmlp->update();

        bool allZero = true;


        //if(!missPricing)
        {
            objRmlp = uRmlp->get(GRB_DoubleAttr_ObjVal);
            privObjRmlp = objRmlp;
        }


/*
        if(!subProbCustR_neg && Stabilization && !exactPi)
        {
            exactPi = true;
            subProbCustR_neg = true;
        }
        else
            exactPi = false;

*/
        /*
        if(phaseStatus != PhaseStatus::PhaseStatusColGen)
        {

            vetVar = uRmlp->getVars();

            for(int i = 0; i < info.numConstrsOrignalProblem + info.numConstrsConv; ++i)
            {
                if(!doubleEqual(vetVar[i].get(GRB_DoubleAttr_X), 0.0, 1E-5))
                {
                    allZero = false;
                    break;
                }
            }


            if(phaseStatus == PhaseStatus::PhaseStatusTwoPhase && allZero)
            {
                delete []vetVar;
                uRmlp->update();
                vetVar = uRmlp->getVars();

                // Recreate the obj function with the distance of the columns
                GRBLinExpr newObj;
                int ini = (int) info.numConstrsOrignalProblem + info.numConstrsConv;

                if(uRmlp->get(GRB_IntAttr_NumVars) != (int) vetVarLambdaCol.size())
                {
                    std::cout << "ERROR\n";
                    std::cout << "LP: " << uRmlp->get(GRB_IntAttr_NumVars) << "; vet: " << vetVarLambdaCol.size()
                              << "\n";
                    PRINT_DEBUG("", "");
                    throw "ERROR";
                }

                for(int i = ini; i < (int) vetVarLambdaCol.size(); ++i)
                {
                    Eigen::VectorXd &vetSol = *vetVarLambdaCol[i];
                    double cgCooefObj = vetSol.dot(auxVect.vetRowC);
                    newObj += cgCooefObj * vetVar[i];

                }

                // "delete" all artificial variables
                for(int i = 0; i < ini; ++i)
                {
                    vetVar[i].set(GRB_DoubleAttr_LB, 0.0);
                    vetVar[i].set(GRB_DoubleAttr_UB, 0.0);
                }

                uRmlp->setObjective(newObj, GRB_MINIMIZE);
                phaseStatus = PhaseStatus::PhaseStatusColGen;
                std::cout<<"Changing phase to colGen!\n";
                subProbCustR_neg = true;
                //auxVect.vetRowRmlpSmoothPi.setZero();
                exactPi = true;
                continue;

                std::cout<<"Change to second phase!\n";
            }

            if(phaseStatus == PhaseStatus::PhaseStatusBigM && allZero)
            {
                int ini = (int) info.numConstrsOrignalProblem + info.numConstrsConv;

                // "delete" all artificial variables
                for(int i = 0; i < ini; ++i)
                {
                    vetVar[i].set(GRB_DoubleAttr_LB, 0.0);
                    vetVar[i].set(GRB_DoubleAttr_UB, 0.0);
                }

                phaseStatus = PhaseStatus::PhaseStatusColGen;
                std::cout<<"Changing phase to colGen!\n";

                auxVect.vetRowRmlpSmoothPi.setZero();
                subProbCustR_neg = true;
            }

            delete[]vetVar;

            if(!subProbCustR_neg && phaseStatus == PhaseStatus::PhaseStatusTwoPhase && !allZero)
                return StatusSubProb_Inviavel;

        }
        else
        */
        {

            if(!exactPi && Stabilization && !missPricing && !subProbCustR_neg)
            {
                //exactPi = true;
                //subProbCustR_neg = true;
                std::cout<<"Seting exactPi\n";
            }
            else if(exactPi && Stabilization && !missPricing && subProbCustR_neg)
                exactPi = false;
            else if(!exactPi && Stabilization && missPricing && subProbCustR_neg)
            {
                //exactPi = true;
                //std::cout<<"Set exactPi\n";
            }
            /*
            if(missPricing && exactPi)
            {
                std::cout<<"ERROR!m miss pricing with exact duals!\n";
                PRINT_DEBUG("", "");
                throw "ERROR";
            }
             */
        }


        //if(!missPricing)
        /*
        {
            objRmlp = uRmlp->get(GRB_DoubleAttr_ObjVal);
            privObjRmlp = objRmlp;
        }
        */


        if(missPricing && !Stabilization)
        {
            uRmlp->write("missPricing.lp");
            std::cout<<"miss pricing\n";
            std::cout<<auxVect.vetRowRmlpPi<<"\n";
            PRINT_DEBUG("", "");

            std::cout<<"vet0: "<<vetVar0<<"\n";
            std::cout<<"vet1: "<<vetVar1<<"\n";

            std::cout<<"phaseStatus: "<<(int)phaseStatus<<"\n";

            throw "MISS_PRICING";
        }


        //lagrangeDualBound = getLagrangeDualBound(objRmlp, redCost);

        //std::cout<<"GAP("<<gap<<"%)\n";

        if((itCG%5) == 0)
        {
            //std::cout<<"\t"<<itCG<<"\t"<<uRmlp->get(GRB_DoubleAttr_ObjVal)<<"\t\""<<gap<<"%\"\n";
            std::cout<<std::format("\t{0}\t{1:.1f}\t{2:.1f}\t{3:.1f}%\n", itCG, objRmlp, lagrangeDualBound, gap);
        }

        itCG += 1;

        //std::cout<<"PI: "<<auxVect.vetRowRmlpPi<<"\n\n";

        //if(itCG == 40)
        //    return DW_DecompNS::StatusProb::StatusSubProb_Inviavel;
        //delete []vetRmlpConstr;


    }

    //std::cout<<"PI: "<<auxVect.vetRowRmlpPi<<"\n";

    uRmlp->update();
    uRmlp->optimize();

    funcObj = uRmlp->get(GRB_DoubleAttr_ObjVal);
    //std::cout<<"\t"<<itCG-1<<"\t"<<uRmlp->get(GRB_DoubleAttr_ObjVal)<<"\t\""<<gap<<"\"\n\n\n";
    std::cout<<std::format("\t{0}\t{1:.2f}\t{2:.2f}%\n", itCG-1, uRmlp->get(GRB_DoubleAttr_ObjVal), gap);
    //uRmlp->write("rmlp_"+std::to_string(itCG)+".lp");



    getSolX();
    //std::cout<<"vetX: \n"<<<<"\n";

    StatusProb status = DW_DecompNS::StatusProb::StatusSubProb_Otimo;
//std::cout<<"numConstrsOrignalProblem: "<<info.numConstrsOrignalProblem+info.numConstrsConv<<"\n\n";
    GRBVar *vetVar        = uRmlp->getVars();
    for(int i=0; i < info.numConstrsOrignalProblem+info.numConstrsConv; ++i)
    {
        if(!doubleEqual(vetVar[i].get(GRB_DoubleAttr_X), 0.0))
        {
            //std::cout<<"!0\n";
            status = DW_DecompNS::StatusProb::StatusSubProb_Inviavel;
            break;
        }
    }

    delete []vetVar;
    delete []vetRmlpConstr;
    vetRmlpConstr = nullptr;

    //uRmlp->write("lpFinal.lp");


    std::cout<<"PI: "<<auxVect.vetRowRmlpPi<<"\n\n";
    uRmlp->write("lastRmlp.lp");

    std::cout<<"*******************~Column Generation*******************\n\n";

    return status;


}

void DW_DecompNS::DW_DecompNode::updateRmlpPi(Eigen::RowVectorXd &vetRowRmlpPi)
{

//std::cout<<"updateRmlpPi\n";

    //vetRowRmlpPi.resetVector();

    vetRowRmlpPi.setZero();

    for(int i=0; i < info.numVarRmlpPi; ++i)
    {
        double val = vetRmlpConstr[i].get(GRB_DoubleAttr_Pi);
        //if(val != 0.0)
            vetRowRmlpPi.coeffRef(0, i) = val;
    }

    //if(PrintDebug)
    //    std::cout<<"PI: "<<vetRowRmlpPi<<"\n\n";

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


void DW_DecompNS::DW_DecompNode::addColumn(const double cost, int k, AuxData &auxVect)
{

    auxVect.updateAuxVetCooef();
    GRBColumn grbColumn;
    if(vetRmlpConstr == nullptr)
    {
        std::cout<<"vetRmlpConstr is null";
        PRINT_EXIT();
    }
    grbColumn.addTerms(auxVect.auxVetCooef, vetRmlpConstr, info.numConstrsMaster+info.numConstrsConv);

    uRmlp->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS, grbColumn,
                "l_"+std::to_string(itCG)+"_"+std::to_string(k));
}

void  DW_DecompNS::DW_DecompNode::addColumnX(const double cost, int k, AuxData& auxVect, Eigen::VectorXd& vetX)
{
    vetVarLambdaCol.push_back(std::make_unique<Eigen::VectorXd>(vetX));
    Eigen::VectorXd& vetSol = *vetVarLambdaCol[vetVarLambdaCol.size()-1];

    if(setVarLamdaCol.count(SolXHash(vetSol)) == 1)
    {
        vetVarLambdaCol.pop_back();
        return;
    }

    //std::cout<<"Add: "<<vetSol.transpose()<<"\n";

    setVarLamdaCol.emplace(vetSol);
    auxVect.vetColCooef.segment(0, info.numConstrsConv) = auxVect.vetColConvCooef;

    auxVect.vetColCooef.segment(info.numConstrsConv, info.numConstrsMaster) = matA*vetSol;

    addColumn(cost, k, auxVect);
}

double DW_DecompNS::DW_DecompNode::getLagrangeDualBound(double objRmlp, double redCost)
{
   return objRmlp + redCost;
}

DW_DecompNS::DW_DecompNode::DW_DecompNode(const DW_DecompNS::DW_DecompNode &decomp)
{

    ptrSubProb = decomp.ptrSubProb;
    //decomp.uRmlp->update();
    uRmlp      = std::make_unique<GRBModel>(*decomp.uRmlp);
    uRmlp->update();
    info       = decomp.info;
    itCG       = 0;
    matA       = decomp.matA;
    vetVarArtifRmlp = nullptr;
    vetVar0  = decomp.vetVar0;
    rhsConv  = decomp.rhsConv;


    vetVarLambdaCol.reserve(decomp.vetVarLambdaCol.size());
    for(int i=0; i < int(decomp.vetVarLambdaCol.size()); ++i)
        vetVarLambdaCol.push_back(std::make_unique<Eigen::VectorXd>(*decomp.vetVarLambdaCol[i]));

    if(vetVarLambdaCol.size() != decomp.vetVarLambdaCol.size())
    {
        std::cout<<"ERROR in copy vetVarLambaCol\n";
        PRINT_DEBUG("", "");
        throw "ERROR";

    }

    if(decomp.uRmlp->get(GRB_IntAttr_NumVars) != uRmlp->get(GRB_IntAttr_NumVars))
    {
        std::cout<<"ERROR in copy model\n";
        PRINT_DEBUG("", "");
        throw "ERROR";
    }

    for(int i=0; i < int(decomp.vetVarLambdaCol.size()); ++i)
        setVarLamdaCol.insert(SolXHash(*vetVarLambdaCol[i]));


    vetSolX = decomp.vetSolX;

    vetSubMatA = Vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>(info.numSubProb);
    for(int i=0; i < info.numSubProb; ++i)
        vetSubMatA[i] = decomp.vetSubMatA[i];

    vetRmlpConstr = nullptr;


}

void DW_DecompNS::DW_DecompNode::getSolX()
{

    GRBVar *vetVar        = uRmlp->getVars();
    double *vetRmlpLambda = uRmlp->get(GRB_DoubleAttr_X, vetVar, uRmlp->get(GRB_IntAttr_NumVars));
    double sum = 0.0;

    //std::cout<<"Lambda: ";
    for(int i=0; i < uRmlp->get(GRB_IntAttr_NumVars); ++i)
    {
        //std::cout << vetRmlpLambda[i] << " ";
        sum += vetRmlpLambda[i];
    }
    //std::cout<<"\nsum: "<<sum<<"\n\n";
    //std::cout<<"vetVarLambdaCol.size(): "<<vetVarLambdaCol.size()<<"\n";

    vetSolX.setZero();
    for(int i=0; i < uRmlp->get(GRB_IntAttr_NumVars); ++i)
        vetSolX += vetRmlpLambda[i]*(*vetVarLambdaCol[i]);

    delete []vetRmlpLambda;
    delete []vetVar;
}

void DW_DecompNS::AuxData::updateSizes(DW_DecompNS::DW_DecompNode &e)
{
    Info &info = e.info;
    int numConstrsMaster = info.numConstrsMaster + (int)info.numConstrsConv;


    if(vetRowRmlpPi.size() >= numConstrsMaster)
        return;

    vetRowRmlpPi.resize(numConstrsMaster);
    vetRowRmlpSmoothPi.resize(numConstrsMaster);
    vetColCooef.resize(numConstrsMaster);

    vetRowRmlpPi.setZero();
    vetRowRmlpSmoothPi.setZero();
    vetColCooef.setZero();

    delete []auxVetCooef;
    auxVetCooef = new double[numConstrsMaster];

}

DW_DecompNS::AuxData::~AuxData()
{
    delete []auxVetCooef;
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
