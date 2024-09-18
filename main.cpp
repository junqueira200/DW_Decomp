#include <iostream>
#include <memory>
#include "MNFP_Inst.h"
#include "DW_Decomp.h"

int resolveSubProb(const Eigen::VectorXd &subProbCooef, int k, void *data, Eigen::VectorXd &vetX, int itCG, bool &custoRedNeg)
{
    try
    {

        std::cout << "Funcao resolveSubProb\n\n\n";
        std::cout<<"subProbCooef: "<<subProbCooef.transpose()<<"\n\n";

        DW_Decomp::StatusSubProb status = DW_Decomp::StatusSubProb_Otimo;
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
            status = DW_Decomp::StatusSubProb_Otimo;

            if(model.get(GRB_DoubleAttr_ObjVal) < -DW_Decomp::TolObjSubProb)
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
                status = DW_Decomp::StatusSubProb_Unbounded;
            } else if(s == GRB_INFEASIBLE)
            {
                status = DW_Decomp::StatusSubProb_Inviavel;
            } else
                status = DW_Decomp::StatusSubProb_Outro;
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

}

int main()
{

    MNFP::MNFP_Inst mnfp = MNFP::criaToyInstance();
    const int K = mnfp.K;
    const int N = mnfp.N;

    GRBEnv env;
    GRBModel mestre(env);
    Eigen::VectorX<std::unique_ptr<GRBModel>> vetSubProb(2);

    for(int k=0; k < K; ++k)
        vetSubProb(k) = std::make_unique<GRBModel>(env);

    criaSubProb(mnfp, *(vetSubProb(0)), 0);
    criaSubProb(mnfp, *(vetSubProb(1)), 1);
    criaMestre(mnfp, mestre);

    mestre.update();
    for(int k=0; k < K; ++k)
    {
        vetSubProb(k)->update();
        vetSubProb(k)->write("subProb_"+std::to_string(k)+".lp");
    }

    mestre.write("mestre.lp");


    DW_Decomp::dwDecomp(env, mestre, vetSubProb,99999.0,N*N,
                        {std::make_pair(0, N*N), std::make_pair(N*N, N*N)}, {&resolveSubProb, &resolveSubProb});


    return 0;
}