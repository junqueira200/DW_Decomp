#include <iostream>
#include <memory>
#include "MNFP_Inst.h"
#include "DW_Decomp.h"
#include "Grafo.h"

class MySubProb : public DW_DecompNS::SubProb
{
public:

    Eigen::VectorX<std::unique_ptr<GRBModel>> vetSubProb;
    const int numSubProb = 2;
    bool convConstIni = false;

    MySubProb(GRBEnv &e, const MNFP::MNFP_Inst &mnfp)
    {
        vetSubProb = Eigen::VectorX<std::unique_ptr<GRBModel>>(numSubProb);
        for(int i=0; i < numSubProb; ++i)
        {
            vetSubProb(i) = std::make_unique<GRBModel>(e);
            criaSubProb(mnfp, *(vetSubProb(i)), i);
        }
    }


    ~MySubProb() override {};
    void iniConvConstr(GRBModel &rmlp, void *data, const double custoVarA) override
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
    int resolveSubProb(Eigen::VectorXd &subProbCooef,
                       GRBModel &mestre,
                       Eigen::VectorXd &vetX,
                       int itCG,
                       bool &custoRedNeg,
                       void *data,
                       const int iniConv,
                       int indSubProb,
                       Eigen::VectorXd &vetCooefRestConv,
                       const std::pair<int, int> &pairSubProb) override
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

}; // FIM MySubProb

int resolveSubProb(const Eigen::VectorXd &subProbCooef, int k, void *data, Eigen::VectorXd &vetX, int itCG, bool &custoRedNeg)
{
    try
    {

        std::cout << "Funcao resolveSubProb\n\n\n";
        std::cout<<"subProbCooef: "<<subProbCooef.transpose()<<"\n\n";

        DW_DecompNS::StatusSubProb status = DW_DecompNS::StatusSubProb_Otimo;
        custoRedNeg = false;

        GRBModel &model = (*(GRBModel *) data);
        GRBVar *varX = model.getVars();


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
                custoRedNeg = true;
                for(int i = 0; i < model.get(GRB_IntAttr_NumVars)-1; ++i)
                {
                    vetX[i] = varX[i].get(GRB_DoubleAttr_X);
                }
            }
        } else
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
        return status;

    }
    catch (GRBException &e)
    {
        std::cout<<e.getMessage()<<"\n";
        exit(-1);
    }

} // FIM resolveSubProb

int main()
{

    GraphNS::Graph<int> graph(5);

    graph.addArc(0, 1, -1);
    graph.addArc(0, 2, 4);

    graph.addArc(1, 2, 3);
    graph.addArc(1, 3, 2);
    graph.addArc(1, 4, 2);

    graph.addArc(3, 2, 5);
    graph.addArc(3, 1, 1);

    graph.addArc(4, 3, -3);


    for(int i=0; i < 5; ++i)
    {
        auto pair = graph.getArcs(i);

        for(auto it=pair.first; it != pair.second; ++it)
        {
            std::cout<<"("<<i<<", "<<it->first<<"): "<<it->second<<"\n";
        }
    }

    if(graph.arcExist(0, 1))
        std::cout<<"(0, 1) exist!\n";

    float val = graph.getArc(0, 1);
    std::cout<<"val: "<<val<<"\n\n";

    std::cout << "arcExist(0,0): " << graph.arcExist(0, 1) << "\n";


    for(int i=0; i < 5; ++i)
    {

        for(auto &it:graph.getArcsRange(i))
        {
            std::cout<<"("<<i<<", "<<it.first<<"): "<<it.second<<"\n";
        }
    }


    GraphNS::bellmanFord(graph, 0);


    return 0;





    MNFP::MNFP_Inst mnfp = MNFP::criaToyInstance();
    const int K = mnfp.K;
    const int N = mnfp.N;

    GRBEnv env;
    GRBModel mestre(env);
    MySubProb subProb(env, mnfp);

/*    Eigen::VectorX<std::unique_ptr<GRBModel>> vetSubProb(2);

    for(int k=0; k < K; ++k)
        vetSubProb(k) = std::make_unique<GRBModel>(env);

    criaSubProb(mnfp, *(vetSubProb(0)), 0);
    criaSubProb(mnfp, *(vetSubProb(1)), 1);*/
    criaMestre(mnfp, mestre);

    mestre.update();
/*    for(int k=0; k < K; ++k)
    {
        vetSubProb(k)->update();
        vetSubProb(k)->write("subProb_"+std::to_string(k)+".lp");
    }*/

    mestre.write("mestre.lp");
    auto vetPairSubProb = std::vector<std::pair<int,int>>{std::make_pair(0, N*N), std::make_pair(N*N, N*N)};

    DW_DecompNS::dwDecomp(env,
                          mestre,
                          99999.0,
                          std::forward<std::vector<std::pair<int,int>>>(vetPairSubProb),
                          (DW_DecompNS::SubProb*)&subProb,
                          (void*)&mnfp,
                          2,
                          2);


    return 0;
} // FIM main
