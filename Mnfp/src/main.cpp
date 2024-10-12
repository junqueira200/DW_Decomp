#include <iostream>
#include <memory>
#include "MNFP_Inst.h"
#include "DW_Decomp.h"
#include "Grafo.h"
#include "MnfpDecomp.h"
#include "Aux.h"


int main()
{


    SparseMatNS::SparseMatLin<double> sparseMat(3, 3);
    sparseMat(0, 1) = 2.0;
    sparseMat(0, 0) = 1.0;

    std::cout<<sparseMat<<"\n";

    sparseMat.addCollums(1);

    std::cout<<sparseMat<<"\n";

    std::cout<<"sizeof: "<<sizeof(DW_DecompNS::DW_DecompNode)<<"\n";
    std::cout<<"alignof: "<< alignof(DW_DecompNS::DW_DecompNode)<<"\n";

    return 0;


    MNFP::MNFP_Inst mnfp = MNFP::criaToyInstance();
    const int K = mnfp.K;
    const int N = mnfp.N;

    GRBEnv env;
    GRBModel mestre(env);
    MnfpDecompNS::MySubProbFlow subProb(env, mnfp);
    std::cout<<"subProbFlow\n";
    MnfpDecompNS::MySubProbPath subProbPath(env, mnfp);

    //return 0;

/*    Eigen::VectorX<std::unique_ptr<GRBModel>> vetSubProb(2);

    for(int k=0; k < K; ++k)
        vetSubProb(k) = std::make_unique<GRBModel>(env);

    criaSubProb(mnfp, *(vetSubProb(0)), 0);
    criaSubProb(mnfp, *(vetSubProb(1)), 1);*/
    MnfpDecompNS::criaMestre(mnfp, mestre);

    mestre.update();
/*    for(int k=0; k < K; ++k)
    {
        vetSubProb(k)->update();
        vetSubProb(k)->write("subProb_"+std::to_string(k)+".lp");
    }*/

    mestre.write("mestre.lp");
    auto vetPairSubProb = std::vector<std::pair<int,int>>{std::make_pair(0, N*N), std::make_pair(N*N, N*N)};

    DW_DecompNS::dwDecomp(env, mestre, 9999.0, std::forward<std::vector<std::pair<int, int>>>(vetPairSubProb),
                          (DW_DecompNS::SubProb *) &subProbPath, (void *) &mnfp, 2);


    return 0;
} // FIM main
