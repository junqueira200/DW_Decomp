#include <iostream>
#include <memory>
#include "MNFP_Inst.h"
#include "DW_Decomp.h"
#include "Grafo.h"
#include "MnfpDecomp.h"



int main()
{

/*    GraphNS::Graph<int> graph(5);

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


    return 0;*/





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

    DW_DecompNS::dwDecomp(env,
                          mestre,
                          99999.0,
                          std::forward<std::vector<std::pair<int,int>>>(vetPairSubProb),
                          (DW_DecompNS::SubProb*)&subProbPath,
                          (void*)&mnfp,
                          2,
                          2);


    return 0;
} // FIM main
