#include "SCIP.h"


SCIP_NS::MatrixVar::MatrixVar(scippp::Model &model, size_t n_, size_t m_, std::string&& prex, scippp::VarType type)
{
    std::printf("Start MatrixVar\n");
    n = n_;
    m = m_;
    vetVar = model.addVars(prex, n*m, scippp::COEFF_ZERO, type);
    SCIP *scip = model.scip();
    for(size_t i=0; i < n; ++i)
    {
        for(size_t j=0; j < m; ++j)
        {
            size_t index = m*i+j;
            SCIP_Var* varPtr = vetVar[index].getVar();
            std::string name = prex + std::format("_[{}][{}]", i, j);
            SCIPchgVarName(scip, varPtr, name.c_str());
        }
    }

    std::printf("END MatrixVar\n");
}

SCIP_NS::Scip3dPacking::Scip3dPacking(const VectorI &vetItems_, const int numItems_, const SolucaoNS::Rota &rota_):
                        vetItems(vetItems_),numItems(numItems_), rota(rota_), model("Scip3dPacking")
{
    std::printf("Scip3dPacking Start!\n\n");
    ptrScip = model.scip();

    MatrixVar matVar(model, 5, 5, "Test", scippp::VarType::BINARY);
    Matrix3DVar mat3dVar(model, 3, 3, 3, "t", scippp::VarType::BINARY);

    scippp::LinExpr lin = matVar(0, 0);
    model.addConstr(lin == 1, "const0");
    model.solve();

    std::printf("Scip3dPacking END!\n\n");

}

SCIP_NS::Matrix3DVar::Matrix3DVar(scippp::Model &model, size_t n_, size_t m_, size_t p_, std::string &&prex, scippp::VarType type)
{

    std::printf("Start Matrix3DVar\n\n");
    n = n_;
    m = m_;
    p = p_;

    vetVar = model.addVars(prex, n*m*p, scippp::COEFF_ZERO, type);
    SCIP *scip = model.scip();
    for(size_t i=0; i < n; ++i)
    {
        for(size_t j=0; j < m; ++j)
        {
            for(size_t k=0; k < p; ++k)
            {
                size_t index = (i*m*p)+(j*p)+k;
                std::printf("%ld,%ld,%ld: %ld\n", i, j, k, index);
                SCIP_Var* varPtr = vetVar[index].getVar();
                std::string name = prex + std::format("_[{}][{}]", i, j);
                SCIPchgVarName(scip, varPtr, name.c_str());
            }
        }
    }

    std::printf("\n\nEND Matrix3DVar\n");

}
