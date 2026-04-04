#include "SCIP.h"


using namespace scippp;
using namespace InstanceNS;

SCIP_NS::MatrixVar::MatrixVar(scippp::Model &model, size_t n_, size_t m_, std::string&& prex, scippp::VarType type)
{
    std::printf("Start MatrixVar\n");
    n = n_;
    m = m_;
    vetVar = model.addVars(prex, n*m, scippp::COEFF_ZERO, type);

    if(prex.empty())
        return;

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

SCIP_NS::Scip3dPacking::Scip3dPacking(const VectorI &vetItems_,
                                      const int numItems_,
                                      const SolucaoNS::Rota &rota_,
                                      SolucaoNS::Bin& bin):
                                                            vetItems(vetItems_),
                                                            numItems(numItems_),
                                                            rota(rota_),
                                                            model("Scip3dPacking")
{
    std::printf("Scip3dPacking Start!\n\n");
    ptrScip = model.scip();

    MY_SCIP_CALL(SCIPsetEmphasis(ptrScip, SCIP_PARAMEMPHASIS_CPSOLVER, TRUE));

    /*
    MatrixVar matVar(model, 5, 5, "Test", scippp::VarType::BINARY);
    Matrix3DVar mat3dVar(model, 3, 3, 3, "t", scippp::VarType::BINARY);

    scippp::LinExpr lin = matVar(0, 0);
    //scippp::LinIneq expr = lin == 0;

    model.createIndicatorConstraint(matVar(0,1) == 1, matVar(0, 0), "teste");

    model.addConstr(lin == 1, "const0");
    std::printf("Before solve()\n");
    model.solve();
    std::printf("After solve()\n");
    std::filesystem::directory_entry file("model.lp");
    model.writeOrigProblem(file);

    */

    //auto sol = model.getBestSol();
    //std::cout<<"(0,0): "<<matVar(0,0).getSolVal(sol);
    //std::cout<<"\n(0,1): "<<matVar(0,1).getSolVal(sol);


    createVariables();
    createConstraints();

    model.solve();

    if(model.getStatus() != SCIP_STATUS_OPTIMAL)
        return;
    //auto iis = model.generateIIS();


    //std::filesystem::directory_entry file("model.lp");
    //model.writeOrigProblem(file);

    auto sol = model.getBestSol();

    for(int i=0; i < numItems; ++i)
    {
        double px, py, pz;
        px = mStartPositionsX[i].getSolVal(sol);
        py = mStartPositionsY[i].getSolVal(sol);
        pz = mStartPositionsZ[i].getSolVal(sol);

        double length = mLengths[i].getSolVal(sol);
        double width = mWidths[i].getSolVal(sol);
        double height = mHeights[i].getSolVal(sol);

        int indRot = -1;

        for(int k=0; k < vetRot.size(); ++k)
        {
            auto r = vetRot[k];
            if(mOrientation(i, (int)r).getSolValAsInt(sol) >= 1)
            {
                indRot = (int)r;
                break;
            }
        }

        if(indRot == -1)
        {
            std::printf("Error in rotation!\n");
            EXIT_PRINT();
        }

        bin.vetPosItem[i].set(px, py, pz);
        bin.vetRotacao[i] = (InstanceNS::Rotation)indRot;

        /*
        for(int j=0; j < numItems; ++j)
        {
            for(int k=0; k < 6; ++k)
            {
                int mRel = mRelativeDirections(i, j, k).getSolVal(sol);
                std::printf("mRel[%d][%d][%d]: %d\n", i, j, k, mRel);
            }

            std::printf("\n");
        }

        bin.vetPosItem[i].set(px, py, pz);
        std::printf("%f %f %f; rot: %d; Length: %f; width: %f; Height: %f\n", px, py, pz, indRot, length, width, height);

        Item& item = instanciaG.vetItens[vetItems[i]];
        length = item.vetDim[0];
        width = item.vetDim[1];
        height = item.vetDim[2];

        std::printf("Length: %f; width: %f; Height: %f\n\n\n", length, width, height);
        */

    }

    if(bin.verificaViabilidade())
        std::printf("MODEL IS FEASIBLE!\n");

    else
    {
        std::printf("MODEL IS INFEASIBLE!\n");
        EXIT_PRINT();
    }

    std::printf("\n\nScip3dPacking END!\n\nIIS:\n");
    //for(std::string& name:iis.consIds)
    //    std::printf("%s\n", name.c_str());

}

void SCIP_NS::Scip3dPacking::createVariables()
{

    auto typePosition = VarType::CONTINUOUS;

    VectorI vetZero(numItems);
    vetZero.setAll(0);

    mStartPositionsX = model.addVars("starX", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[0]);
    mEndPositionsX   = model.addVars("endX", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[0]);

    mStartPositionsY = model.addVars("starY", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[1]);
    mEndPositionsY   = model.addVars("endY", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[1]);

    mStartPositionsZ = model.addVars("starZ", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[2]);
    mEndPositionsZ   = model.addVars("endZ", numItems, vetZero, typePosition, 0, instanciaG.vetDimVeiculo[2]);

    mLengths         = model.addVars("length", numItems, vetZero, typePosition, 0, 0);
    mWidths          = model.addVars("width", numItems, vetZero, typePosition, 0, 0);
    mHeights         = model.addVars("height", numItems, vetZero, typePosition, 0, 0);

    mPlacedOnFloor   = model.addVars("placedOnFloor", numItems, vetZero, VarType::BINARY, 0, 1);

    for(size_t i=0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];

        int min = std::min(item.vetDim[0], item.vetDim[1]);
        int max = std::max(item.vetDim[0], item.vetDim[1]);

        MY_SCIP_CALL(SCIPchgVarLb(ptrScip, mLengths[i].getVar(), min));
        MY_SCIP_CALL(SCIPchgVarUb(ptrScip, mLengths[i].getVar(), max));

        MY_SCIP_CALL(SCIPchgVarLb(ptrScip, mWidths[i].getVar(), min));
        MY_SCIP_CALL(SCIPchgVarUb(ptrScip, mWidths[i].getVar(), max));

        MY_SCIP_CALL(SCIPchgVarLb(ptrScip, mHeights[i].getVar(), item.vetDim[2]));
        MY_SCIP_CALL(SCIPchgVarUb(ptrScip, mHeights[i].getVar(), item.vetDim[2]));


    }

    mRelativeDirections = Matrix3DVar(model, numItems, numItems, 6, "", VarType::BINARY);
    mOrientation = MatrixVar(model, numItems, 2, "", VarType::BINARY);

    mItemsOverlapsXY = MatrixVar(model, numItems, numItems, "over", VarType::BINARY);
    mSupportXY = MatrixVar(model, numItems, numItems, "supp", VarType::BINARY);
    mOverlapAreasXY = MatrixVar(model, numItems, numItems, "overA", VarType::CONTINUOUS);

    for(int i=0; i < numItems; ++i)
    {
        for(int j=0; j < numItems; ++j)
        {
            for(auto dir:arrayRelaDirections)
            {
                std::string dirStr = mapRelaDirections.at(dir);
                SCIP_Var* var = mRelativeDirections(i, j, dir).getVar();
                std::string name = "relDir"+std::format("[{}][{}][{}]", i, j, dirStr);
                MY_SCIP_CALL(SCIPchgVarName(ptrScip, var, name.c_str()));
            }
        }

        std::string nameNoRot = "orien"+std::format("[{}][noRot]", i);
        std::string nameRotZ  = "orien"+std::format("[{}][rotZ]", i);

        MY_SCIP_CALL(SCIPchgVarName(ptrScip, mOrientation(i, 0).getVar(), nameNoRot.c_str()));
        MY_SCIP_CALL(SCIPchgVarName(ptrScip, mOrientation(i, 1).getVar(), nameRotZ.c_str()));

    }

}

void SCIP_NS::Scip3dPacking::createConstraints()
{

    CreateItemOrientations();
    CreateEnd();
    CreateNoOverlap();
    CreateOnFloorConstraints();

}

void SCIP_NS::Scip3dPacking::CreateNoOverlap()
{

    const double eps = 1E-5;

    for(int i=0; i < numItems; ++i)
    {
        for(int j=0; j < numItems; ++j)
        {
            if(i == j)
                continue;

            scippp::LinExpr sum;
            std::string str6 = std::format("R6_({},{})", i, j);

            for (size_t d = 0; d < mDimensions.size(); ++d)
            {
                const Dimension& dimension = mDimensions[d];
                auto startPositionI = getIntVars(dimension.Type, true, i);
                auto startPositionJ = getIntVars(dimension.Type, true, j);

                auto endPositionI   = getIntVars(dimension.Type, false, i);
                auto endPositionJ   = getIntVars(dimension.Type, false, j);

                sum += mRelativeDirections(i, j, dimension.FirstDirection) +
                       mRelativeDirections(i, j, dimension.SecondDirection);

                std::string str1 = std::format("R1_({},{},{})", i, j, d);
                std::string str3 = std::format("R3_({},{},{})", i, j, d);
                std::string str5 = std::format("R5_({},{},{})", i, j, d);
                std::string str5_ = std::format("R5^_({},{},{})", i, j, d);
                std::string str9 = std::format("R9_({},{},{})", i, j, d);

                // endPositionJ <= startPositionI
                model.createIndicatorConstraintLessThan(endPositionJ - startPositionI, 0.0,
                                                        mRelativeDirections(i, j, dimension.FirstDirection), str1);
                //model.createIndicatorConstraint(endPositionJ >= startPositionI-eps,
                //                                mRelativeDirections(i, j, dimension.FirstDirection), "2", true);

                // endPositionI <= startPositionJ
                model.createIndicatorConstraintLessThan(endPositionI - startPositionJ, 0.0,
                                                        mRelativeDirections(i, j, dimension.SecondDirection), str3);
                //model.createIndicatorConstraint(endPositionI >= startPositionJ-eps,
                //                                mRelativeDirections(i, j, dimension.SecondDirection), "4", true);



                model.addConstr(mRelativeDirections(i, j, dimension.FirstDirection) ==
                          mRelativeDirections(j, i, dimension.SecondDirection), str5);

                model.addConstr(mRelativeDirections(i,j,dimension.FirstDirection) +
                                mRelativeDirections(i, j, dimension.SecondDirection) <= 1, str9);

                model.addConstr(mRelativeDirections(i, j, dimension.SecondDirection) ==
                                    mRelativeDirections(j, i, dimension.FirstDirection), str5_);

                // Simetric break constraint
                //model.add(mRelativeDirections[i][j][dimension.SecondDirection] ==
                //          mRelativeDirections[j][i][dimension.FirstDirection]);

                //model.add(mRelativeDirections[i][j][dimension.FirstDirection] +
                //          mRelativeDirections[i][j][dimension.SecondDirection] <= 1);

            }

            model.addConstr(sum >= 1, str6);

        }
    }

}

void SCIP_NS::Scip3dPacking::CreateItemOrientations()
{

    for (int i = 0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        int dx, dy, dz;

        scippp::LinExpr sum;

        std::string strName15 = std::format("R15_({})", i);

        for(InstanceNS::Rotation r:vetRot)
        {
            dx = item.getDimRotacionada(0, r);
            dy = item.getDimRotacionada(1, r);
            dz = item.getDimRotacionada(2, r);

            std::string strName7 = std::format("R7_({},{})", i, (int)r);
            std::string strName7_ = std::format("R7!_({},{})", i, (int)r);
            std::string strName8 = std::format("R8_({},{})", i, (int)r);
            std::string strName8_ = std::format("R8!_({},{})", i, (int)r);
            std::string strName9 = std::format("R9_({},{})", i, (int)r);
            std::string strName9_ = std::format("R9!_({},{})", i, (int)r);
//            std::string strName2 = std::format("R7_[{}][{}]_2", i, (int)r);

            // mLengths[i] <= dx
            model.createIndicatorConstraintLessThan(mLengths[i], dx, mOrientation(i, r), strName7);
            // mLengths[i] >= dx => -mLengths[i] <= -dx
            model.createIndicatorConstraintLessThan(-1*mLengths[i], -dx, mOrientation(i, r), "^"+strName7);


            model.createIndicatorConstraintLessThan(mWidths[i], dy, mOrientation(i, r), strName8);
            model.createIndicatorConstraintLessThan(-1*mWidths[i], -dy, mOrientation(i, r), "^"+strName8);

            model.createIndicatorConstraintLessThan(mHeights[i], dz, mOrientation(i, r), strName9);
            model.createIndicatorConstraintLessThan(-1*mHeights[i], -dz, mOrientation(i, r), "^"+strName9);

            sum +=mOrientation(i, r);
        }

        model.addConstr(sum == 1, strName15);
    }
}

void SCIP_NS::Scip3dPacking::CreateEnd()
{
    std::printf("Start CreateEnd\n");
    std::string empty = "";

    for(int i=0; i < numItems; ++i)
    {
        std::string strI = std::format("_({})", i);
        std::string str10 = "R10"+strI;
        std::string str11 = "R11"+strI;
        std::string str12 = "R12"+strI;

        model.addConstr(mEndPositionsX[i] == mStartPositionsX[i] + mLengths[i], str10);
        //std::printf("1\n");
        model.addConstr(mEndPositionsY[i] == mStartPositionsY[i] + mWidths[i], str11);
        //std::printf("2\n");
        model.addConstr(mEndPositionsZ[i] == mStartPositionsZ[i] + mHeights[i], str12);
        //std::printf("3\n");
    }

    std::printf("END CreateEnd\n");
}

void SCIP_NS::Scip3dPacking::CreateOnFloorConstraints()
{
    double minZ = INF_Double;
    for (size_t i = 0; i < numItems; ++i)
    {
        double z = instanciaG.vetItens[vetItems[i]].vetDim[2];

        if(z < minZ)
            minZ = z;
    }

    for (size_t i = 0; i < numItems; ++i)
    {
        std::string strI = std::format("_({})", i);
        std::string str13 = "R13"+strI;
        std::string str13_ = "R^13"+strI;
        std::string str14 = "R14"+strI;

       // model.addConstr(mStartPositionsZ[i] == 0, str13);
        // If mPlacedOnFloor[i] -> mStartPositionsZ[i] == 0

        // mStartPositionsZ[i] >= 0 => -mStartPositionsZ[i] <= 0
        model.createIndicatorConstraintLessThan(-1*mStartPositionsZ[i], 0.0, mPlacedOnFloor[i], str13);
        // mStartPositionsZ[i] <= 0
        model.createIndicatorConstraintLessThan(mStartPositionsZ[i], 0.0, mPlacedOnFloor[i], str13_);

        // mStartPositionsZ[i] >= minZ
        // -mStartPositionsZ[i] <= -minZ
        model.createIndicatorConstraintLessThan(-1*mStartPositionsZ[i], -minZ, mPlacedOnFloor[i], str14);

        //model.createIndicatorConstraint(mStartPositionsZ[i] >= minZ, mPlacedOnFloor[i], str14, true);
    }
}

Var &SCIP_NS::Scip3dPacking::getIntVars(DimensionType dimension, bool first, int i)
{

    switch (dimension)
    {
    case AxisX:
        if(first)
            return mStartPositionsX[i];
        else
            return mEndPositionsX[i];

    case AxisY:
        if(first)
            return mStartPositionsY[i];
        else
            return mEndPositionsY[i];

    case AxisZ:
        if(first)
            return mStartPositionsZ[i];
        else
            return mEndPositionsZ[i];
    default:
        throw std::runtime_error("DimensionType not implemented.");

    }

}

SCIP_NS::Matrix3DVar::Matrix3DVar(scippp::Model &model, size_t n_, size_t m_, size_t p_, std::string &&prex,
                                  scippp::VarType type)
{

    n = n_;
    m = m_;
    p = p_;

    vetVar = model.addVars(prex, n*m*p, scippp::COEFF_ZERO, type);

    if(prex.empty())
        return;

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


}
