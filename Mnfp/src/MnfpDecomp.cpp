//
// Created by igor on 26/09/24.
//
#include "MnfpDecomp.h"
using namespace MNFP;


void MnfpDecompNS::criaSubProbFlow(const MNFP::MNFP_Inst &mnfp, GRBModel &model, int k)
{
    const int N = mnfp.N;
    std::cout<<"N*N: "<<N*N<<"\n";

    model.reset(1);
    GRBVar* vetX;

    // Cria variavel x
    //for(int k=0; k < K; ++k)
    vetX = model.addVars(N*N+1);

    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            vetX[getId(i, j, N)].set(GRB_StringAttr_VarName, "x_" + std::to_string(k)+"_"+
                                                             std::to_string(i) + "_" +
                                                             std::to_string(j));

            std::cout << "(" <<k<<","<< i << "," << j << "): " << getId(i, j, N) << "\n";
        }

        std::cout << "\n";
    }

    std::cout << "\n\n";

    for(int i = 0; i < N; ++i)
    {
        GRBLinExpr linExpr;
        for(int j=0; j < N; ++j)
        {
            if(j == i)
                continue;

            if(mnfp.matCapacidade(i, j) != 0)
            {
                linExpr += -vetX[getId(i, j, N)];
                vetX[getId(i, j, N)].set(GRB_DoubleAttr_UB, mnfp.matCapacidade(i, j));
            }

            if(mnfp.matCapacidade(j, i) != 0)
            {
                linExpr += vetX[getId(j, i, N)];
                vetX[getId(j, i, N)].set(GRB_DoubleAttr_UB, mnfp.matCapacidade(j, i));
            }
        }

        model.addConstr(linExpr == mnfp.matVertexDem(k, i));

    }

    vetX[N*N].set(GRB_DoubleAttr_UB, 1.0);
    vetX[N*N].set(GRB_DoubleAttr_LB, 1.0);
    vetX[N*N].set(GRB_StringAttr_VarName, "const");

    delete []vetX;

}

void MnfpDecompNS::criaMestre(const MNFP::MNFP_Inst &mnfp, GRBModel &model)
{


    const int N = mnfp.N;
    const int K = mnfp.K;

    model.reset(1);
    Eigen::VectorX<GRBVar*> vetX(K);


    for(int k=0; k < K; ++k)
    {
        vetX(k) = model.addVars(N*N);

        for(int i = 0; i < N; ++i)
        {
            for(int j = 0; j < N; ++j)
            {
                (vetX(k))[getId(i, j, N)].set(GRB_StringAttr_VarName, "x_" + std::to_string(k)+"_"+
                                                                      std::to_string(i) + "_" +
                                                                      std::to_string(j));

                std::cout << "(" <<k<<","<< i << "," << j << "): " << getId(i, j, N) << "\n";
            }

            std::cout << "\n";
        }

        std::cout << "\n\n";
    }




    std::cout<<"\n\n";

    GRBLinExpr obj;

    for(int i=0; i < N; ++i)
    {
        for(int j=0; j < N; ++j)
        {

            GRBLinExpr linExpr;
            for(int k = 0; k < K; ++k)
            {
                linExpr += (vetX(k))[getId(i, j, N)];
                obj += (mnfp.vetArcCost(k))(i, j) * (vetX(k))[getId(i, j, N)];
            }

            model.addConstr(linExpr <= mnfp.matCapacidade(i, j));
        }
    }

    model.setObjective(obj, GRB_MINIMIZE);

    for(int k=0; k < K; ++k)
    {
        delete []vetX(k);
    }

}


void MnfpDecompNS::MySubProbFlow::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{
    if(convConstIni)
        return;

    rmlp.update();
    int numVar = rmlp.get(GRB_IntAttr_NumVars);

    GRBVar a0 = rmlp.addVar(0.0, GRB_INFINITY, custoVarA, GRB_CONTINUOUS, "mestre_" + std::to_string(numVar));
    numVar += 1;
    GRBVar a1 = rmlp.addVar(0.0, GRB_INFINITY, custoVarA, GRB_CONTINUOUS, "mestre_" + std::to_string(numVar));

    rmlp.addConstr(a0, GRB_EQUAL, 1.0, "conv_0");
    rmlp.addConstr(a1, GRB_EQUAL, 1.0, "conv_1");
    rmlp.update();

    convConstIni = true;

}



/** ************************************************************************
 *  ************************************************************************
 *
 *  Resolve o subproblema e add cooef das restricoes de conv
 *
 *  @param subProbCooef
 *  @param mestre
 *  @param vetX
 *  @param itCG
 *  @param custoRedNeg
 *  @param data
 *  @param iniConv
 *  @param indSubProb
 *  @param vetCooefRestConv
 *
 *  ************************************************************************
 *  ************************************************************************
 */
int MnfpDecompNS::MySubProbFlow::resolveSubProb(Eigen::VectorXd &subProbCooef,
                                                GRBModel &mestre,
                                                Eigen::VectorXd &vetX,
                                                int itCG,
                                                bool &custoRedNeg,
                                                void *data,
                                                const int iniConv,
                                                int indSubProb,
                                                Eigen::VectorXd &vetCooefRestConv,
                                                const std::pair<int, int> &pairSubProb)
{
    if(!convConstIni)
    {
        std::cout<<"Restricoes de convexidade nao foram inicializadas!\n";
        PRINT_DEBUG("", "");
        exit(-1);
    }

    const int k = indSubProb;

    MNFP::MNFP_Inst &mnfp = *((MNFP::MNFP_Inst*)data);
    Eigen::VectorXi vetStatus(numSubProb);
    vetX.setZero();

//for(int k=0; k < numSubProb; ++k)



    GRBConstr constr = mestre.getConstr(iniConv+k);
    subProbCooef[pairSubProb.second] = -constr.get(GRB_DoubleAttr_Pi);

    std::cout << "Funcao resolveSubProb\n\n\n";
    std::cout<<"subProbCooef: "<<subProbCooef.segment(0, (pairSubProb.second+1)).transpose()<<"\n\n";
    std::cout<<"vetX.size(): "<<vetX.size()<<"\n";

    DW_DecompNS::StatusSubProb status = DW_DecompNS::StatusSubProb_Otimo;
    custoRedNeg = false;

    GRBModel &model = (*vetSubProb[k]);
    model.update();
    GRBVar *varX = model.getVars();

    std::cout<<"num var subprob "<<k<<": "<<model.get(GRB_IntAttr_NumVars)<<"\n";
    std::cout<<"num var subprob: "<<(pairSubProb.second+1)<<"\n";

    try
    {

        model.set(GRB_DoubleAttr_Obj, varX, &subProbCooef(0), model.get(GRB_IntAttr_NumVars));
        model.update();
        model.write("colGen_subProb_" + std::to_string(k) + "_it_" + std::to_string(itCG) + ".lp");
        model.optimize();

        int s = model.get(GRB_IntAttr_Status);
        if(s == GRB_OPTIMAL)
        {
            status = DW_DecompNS::StatusSubProb_Otimo;

            if(model.get(GRB_DoubleAttr_ObjVal) < -DW_DecompNS::TolObjSubProb)
            {
                std::cout<<"ini for\n";

                custoRedNeg = true;
                for(int i = 0; i < vetX.size(); ++i)
                {
                    vetX[i] = varX[i].get(GRB_DoubleAttr_X);
                }

                std::cout<<"fim for\n";

                vetCooefRestConv.setZero();
                vetCooefRestConv[k] = 1;

                std::cout<<"subProb("<<k<<") X: "<<vetX.segment(0, pairSubProb.second).transpose()<<"\n\n";

                std::cout<<"FIM resolve sub prob.\n";

                delete []varX;
                return s;
            }
        }
        else
        {

            if(s == GRB_UNBOUNDED)
            {
                status = DW_DecompNS::StatusSubProb_Unbounded;
            } else if(s == GRB_INFEASIBLE)
            {
                status = DW_DecompNS::StatusSubProb_Inviavel;
            } else
                status = DW_DecompNS::StatusSubProb_Outro;
        }

        std::cout<<"FIM resolveSubProb\n\n\n";

        delete []varX;
        vetStatus[k] = status;

    }
    catch (GRBException &e)
    {
        std::cout<<e.getMessage()<<"\n";
        exit(-1);
    }


    return status;

} // FIM resolveSubProb


void MnfpDecompNS::MySubProbPath::iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA)
{
    std::cout<<"iniConvConstr\nmatDem: \n";
    std::cout<<mnfp.matVertexDem<<"\n";

    rmlp.update();

    int64_t nextConstr = rmlp.get(GRB_IntAttr_NumConstrs);

    for(int64_t k=0; k < mnfp.K; ++k)
    {
        for(int64_t i=0; i < mnfp.N; ++i)
        {

            if(mnfp.matVertexDem(k, i) != 0)
            {
                try
                {
                    int neg = 1;
                    if(mnfp.matVertexDem(k, i) < 0)
                    {
                        neg = -1;
                    }


                    std::string name = std::to_string(k) + "_" + std::to_string(i);

                    GRBVar var = rmlp.addVar(0, GRB_INFINITY, custoVarA, GRB_CONTINUOUS, "a_"+name);
                    var.set(GRB_StringAttr_VarName, "a_" + name);
                    GRBLinExpr exp;
                    exp += var;
                    rmlp.addConstr(exp == neg * mnfp.matVertexDem(k, i), "conv_" + name);

                    matCovConst(k, i) = nextConstr;
                    nextConstr += 1;

                }
                catch(GRBException &e)
                {
                    std::cout<<e.getMessage()<<"\n";
                    exit(-1);
                }
            }
        }
    }



    rmlp.update();
    rmlp.write("rmlp.lp");

}

void MnfpDecompNS::MySubProbPath::restoreGraphModCost(int64_t k)
{

    for(int64_t i=0; i < vetGraphCost[k].numVertices; ++i)
    {
        for(auto &it:vetGraphCost[k].getArcsRange(i))
            vetGraphModCost[k].addArc(i, it.first, it.second);
    }
}



/**
 *
 *
 * @param pairSubProb pair representa o indice do inicio das variaveis do sub problema e o seu tamanho
 */
int MnfpDecompNS::MySubProbPath::resolveSubProb(Eigen::VectorXd &subProbCooef,
                                                GRBModel &rmlp,
                                                Eigen::VectorXd &vetX,
                                                int itCG,
                                                bool &custoRedNeg,
                                                void *data,
                                                const int iniConv,
                                                int indSubProb,
                                                Eigen::VectorXd &vetCooefRestConv,
                                                const std::pair<int, int> &pairSubProb)

{

    const int k = indSubProb;

    std::cout<<"resolveSubProb: "<<k<<"\n";

    restoreGraphModCost(k);
    GraphNS::Graph<double> &graphM = vetGraphModCost[k];
    GRBConstr *constr = rmlp.getConstrs();


    for(int64_t i=0; i < mnfp.N; ++i)
    {
        for(auto &it:graphM.getArcsRange(i))
        {
            if(it.first == idT)
                continue;

            //graphM.addArc(i, it.first, it.second- constr[getId(i, it.first, (int64_t)mnfp.N)].get(GRB_DoubleAttr_Pi));
            graphM.addArc(i, it.first, subProbCooef[getId(i, it.first, (int64_t)mnfp.N)]);
        }

        if(matVerticeType(k, i) == TypeSource)
        {
            int64_t idK_I = matCovConst(k, i);
            graphM.addArc(idS, i, -constr[matCovConst(k, i)].get(GRB_DoubleAttr_Pi));
        }
    }

    std::cout<<GraphNS::printGraph(graphM);

    Eigen::VectorXd vetDist(graphM.numVertices);
    Eigen::VectorX<int64_t> vetPredecessor(graphM.numVertices);
    int64_t verticeNegCycle;

    GraphNS::bellmanFord(graphM, idS, vetDist, vetPredecessor, verticeNegCycle);

    exit(-1);

    delete constr;
    return 0;

}

MnfpDecompNS::MySubProbPath::MySubProbPath(GRBEnv &e, const MNFP::MNFP_Inst &mnfp_):mnfp(mnfp_), numSubProb(mnfp_.K)
{
    std::cout << "MySubProbPath:\n\n";
    std::cout << mnfp.matVertexDem << "\n";

    matCovConst = Eigen::MatrixX<int64_t>(mnfp.K, mnfp.N);
    matCovConst.setConstant(-1);
    std::cout << "vetCovConst: \n" << matCovConst << "\n\n";

    matVerticeType = Eigen::MatrixX<VerticeType>(mnfp.K, mnfp.N);
    matVerticeType.setConstant(TypeZero);

    vetGraphCost = Eigen::VectorX<GraphNS::Graph<double>>(mnfp.K);
    vetGraphModCost = Eigen::VectorX<GraphNS::Graph<double>>(mnfp.K);

    // Start vetGraphCost
    for(int64_t k = 0; k < mnfp.K; ++k)
    {
        vetGraphCost[k].reset(mnfp.N);

        for(int i = 0; i < mnfp.N; ++i)
        {

            if(mnfp.matVertexDem(k, i) != 0)
            {


                matVerticeType(k, i) = TypeSink;

                if(mnfp.matVertexDem(k, i) < 0)
                    matVerticeType(k, i) = TypeSource;

            }

            for(int j = 0; j < mnfp.N; ++j)
            {
                if(i == j)
                    continue;
                double val = (mnfp.vetArcCost[k])(i, j);

                if(val != 0.0)
                {
                    vetGraphCost[k].addArc(i, j, val);
                }

            }
        }
    }


    idS = mnfp.N;
    idT = mnfp.N+1;

    std::cout<<"N: "<<mnfp.N<<"\n\n";
    std::cout<<"idS: "<<idS<<"\nidT: "<<idT<<"\n\n";


    for(int64_t k = 0; k < mnfp.K; ++k)
    {
        GraphNS::copyGraph(vetGraphCost[k], vetGraphModCost[k]);

        std::cout<<"antes("<<k<<"): "<<vetGraphModCost[k].numVertices<<"\n";

        vetGraphModCost[k].addVertice();
        vetGraphModCost[k].addVertice();


        std::cout<<"antes("<<k<<"): "<<vetGraphModCost[k].numVertices<<"\n";

        for(int64_t i=0; i < mnfp.N; ++i)
        {
            if(matVerticeType(k, i) == TypeSource)
                vetGraphModCost[k].addArc(idS, i, 0.0);
            else if(matVerticeType(k, i) == TypeSink)
                vetGraphModCost[k].addArc(i, idT, 0.0);
        }

        restoreGraphModCost(k);
        vetGraphModCost[k].loadVetArcs();
        std::cout<<vetGraphModCost[k].printVetArcs()<<"\n";

    }
}