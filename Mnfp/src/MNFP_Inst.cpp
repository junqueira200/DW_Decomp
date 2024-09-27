//
// Created by igor on 07/09/24.
//

#include <iostream>
#include "MNFP_Inst.h"


MNFP::MNFP_Inst MNFP::criaToyInstance()
{
    MNFP_Inst mnfp(2, 4);

    (mnfp.vetArcCost(0))(0, 1) = (mnfp.vetArcCost(0))(1, 0) = 3;
    (mnfp.vetArcCost(0))(0, 2) = (mnfp.vetArcCost(0))(2, 0) = 2;
    (mnfp.vetArcCost(0))(1, 2) = (mnfp.vetArcCost(0))(2, 1) = 2;
    (mnfp.vetArcCost(0))(1, 3) = (mnfp.vetArcCost(0))(3, 1) = 4;
    (mnfp.vetArcCost(0))(2, 3) = (mnfp.vetArcCost(0))(3, 2) = 1;

    (mnfp.vetArcCost(1))(0, 1) = (mnfp.vetArcCost(1))(1, 0) = 4;
    (mnfp.vetArcCost(1))(0, 2) = (mnfp.vetArcCost(1))(2, 0) = 2;
    (mnfp.vetArcCost(1))(1, 2) = (mnfp.vetArcCost(1))(2, 1) = 1;
    (mnfp.vetArcCost(1))(1, 3) = (mnfp.vetArcCost(1))(3, 1) = 2;
    (mnfp.vetArcCost(1))(2, 3) = (mnfp.vetArcCost(1))(3, 2) = 5;

    //std::cout<<mnfp.vetArcCost(0)<<"\n\n"<<mnfp.vetArcCost(1)<<"\n";


    mnfp.matVertexDem(0, 0) = -6;
    mnfp.matVertexDem(0, 1) = -2;
    mnfp.matVertexDem(0, 2) =  0;
    mnfp.matVertexDem(0, 3) =  8;


    mnfp.matVertexDem(1, 0) = -7;
    mnfp.matVertexDem(1, 1) =  0;
    mnfp.matVertexDem(1, 2) =  3;
    mnfp.matVertexDem(1, 3) =  4;

    //std::cout<<mnfp.matVertexDem<<"\n";


    mnfp.matCapacidade(0, 1) = 10;
    mnfp.matCapacidade(0, 2) = 5;
    mnfp.matCapacidade(1, 2) = 15;
    mnfp.matCapacidade(1, 3) = 6;
    mnfp.matCapacidade(2, 3) = 7;

    //std::cout<<"\n\n"<<mnfp.matCapacidade<<"\n\n";

    return std::move(mnfp);
}

/*
void MNFP::criaSubProbFlow(const MNFP_Inst &mnfp, GRBModel &model, int k)
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

void MNFP::criaMestreFlow(const MNFP::MNFP_Inst &mnfp, GRBModel &model)
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


    Eigen::VectorX<std::pair<int, int>> arcos(N*N);
    int numArcos = 0;

    // \sum\limits_{k \in K} x_{a}^{k} <= u_{a}   \forall a \in A

    for(int i=0; i < N; ++i)
    {
        for(int j=0; j < N; ++j)
        {
            if(mnfp.matCapacidade(i, j) == 0)
                continue;

            arcos(numArcos) = std::make_pair(i, j);
            std::cout<<"("<<arcos(numArcos).first<<","<<arcos(numArcos).second<<"); ";
            numArcos += 1;
        }
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
*/

int MNFP::getId(int i, int j, int n){return (i*n+j);}
int64_t MNFP::getId(int64_t i, int64_t j, int64_t n){return (i*n+j);}
