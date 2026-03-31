#include "IBM_CpOptimizer.h"
#include "OPP_CP_3D.h"
#include <ilcplex/ilocplex.h>
#include "AuxT.h"
#include <scippp/model.hpp>

using namespace IBM_CpOptimizerNS;
using namespace ParseInputNS;
using namespace InstanceNS;

using namespace scippp;

IloEnv* CpOptimizer::envPtr = nullptr;

/** ************************************************************************
 *  ************************************************************************
 * @param vetItems_   [in]  Vector containing the indices of the items
 *                   to be packed.
 * @param numItems_   [in]  Number of items to consider from vetItems.
 * @param route_      [in]  Route associated with the items (define
 *                   unloading order).
 *
 * ******************************************************************************
 * ******************************************************************************
 */


CpOptimizer::CpOptimizer(const VectorI &vetItems_, int numItems_, const SolucaoNS::Rota &rota_) :vetItems(vetItems_),
                                                                                                numItems(numItems_),
                                                                                                rota(rota_)
{
}

bool CpOptimizer::solve(SolucaoNS::Bin &bin)
{
    std::printf("\n\n**************************************\n");
    std::printf("***********IBM_CP_OPTIMIZER***********\n\n");

    std::printf("numItems: %d\n", numItems);


    int numItems_ = 4;
    numItems = numItems_;
    bin.numItens = numItems_;


    createVariables();
    CreateNoOverlap();
    CreateItemOrientations();



    CreateXYIntersectionBool();
    CreateSupportItem();
    CreateXYIntersectionArea();
    CreateSupportArea();


    CreateLifo();
    CreateOnFloorConstraints();
    CreateAxleWeights();


    //for(int i=0; i < 3; ++i)
    //{
    //    model.add(mStartPositionsZ[i] == 0);
    //}

    IloCP cp(model);
    cp.setParameter(IloCP::Workers, 8);
    //cp.setParameter(IloCP::SearchType, IloCP::MultiPoint);
    cp.setParameter(IloCP::TimeLimit, 120);
    cp.setParameter(IloCP::LogVerbosity, IloCP::Verbose);
    cp.setParameter(IloCP::LogPeriod, 1);

    if(cp.solve())
    {
        std::printf("IloCP Sucess!\n");
        for(int i=0; i < numItems; ++i)
        {
            //int px, py, pz, rot;
            Array<int, 3> start, end;
            start[0] = cp.getValue(mStartPositionsX[i]);
            start[1] = cp.getValue(mStartPositionsY[i]);
            start[2] = cp.getValue(mStartPositionsZ[i]);

            end[0] = cp.getValue(mEndPositionsX[i]);
            end[1] = cp.getValue(mEndPositionsY[i]);
            end[2] = cp.getValue(mEndPositionsZ[i]);


            int rot = cp.getValue(mOrientation[i][Rot1]);

            Array<int, 3> array;

            array[0] = cp.getValue(mLengths[i]);
            array[1] = cp.getValue(mWidths[i]);
            array[2] = cp.getValue(mHeights[i]);

            //std::printf("Widths: %d; Length: %d\n", array[0], array[1]);

            for(int d=0; d < 3; ++d)
            {
                int dimD = instanciaG.vetItens[vetItems[i]].getDimRotacionada(d, (InstanceNS::Rotation)rot);
                if(dimD != array[d])
                {
                    std::printf("Error: dimensios are diferents\n%d != %d\n", dimD, array[d]);
                    EXIT_PRINT();
                }

                if(end[d] != start[d] + array[d])
                {
                    std::printf("end(%d) != start(%d) + array(%d) (%d)\n", end[d], start[d], array[d], start[d] + array[d]);
                    EXIT_PRINT();
                }
            }

            //std::printf("Widths: %d; Length: %d\n", array[0], array[1]);

            bin.vetPosItem[i].set(start[0], start[1], start[2]);
            std::printf("%d\n", start[2]);
            bin.vetRotacao[i] = (InstanceNS::Rotation)((int)Rot1*rot);
        }

        if(bin.verificaViabilidade())
        {
            std::printf("BIN is feaseble!\n");
        }
        else
        {
            std::printf("BIN is infeaseble\n");
            EXIT_PRINT();
        }
    }
    else
    {
        std::printf("IloCP Faild!\n");
        cp.refineConflict();
        cp.writeConflict(std::cout);
    }



    std::printf("***********IBM_CP_OPTIMIZER***********\n");
    std::printf("\n**************************************\n\n");
    EXIT_PRINT();

    return false;
}

void IBM_CpOptimizerNS::CpOptimizer::createVariables()
{
    if(!envPtr)
        envPtr = new IloEnv();
    IloEnv& env = *envPtr;

    model = IloModel(env);

    mStartPositionsX = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[0]);
    NameVars(mStartPositionsX, "startX");

    mEndPositionsX   = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[0]);
    NameVars(mEndPositionsX, "endX");

    mStartPositionsY = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[1]);
    NameVars(mStartPositionsY, "startY");

    mEndPositionsY   = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[1]);
    NameVars(mEndPositionsY, "endY");


    mStartPositionsZ = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[2]);
    NameVars(mStartPositionsZ, "startZ");

    mEndPositionsZ   = IloIntVarArray(env, numItems, 0, (int)instanciaG.vetDimVeiculo[2]);
    NameVars(mEndPositionsZ, "endZ");

    mHeights         = IloIntVarArray(env, numItems, 0, 0);
    NameVars(mHeights, "height");


    mLengths		 = IloIntVarArray(env, numItems, 0, 0);
    NameVars(mLengths, "length");

    mWidths			 = IloIntVarArray(env, numItems, 0, 0);
    NameVars(mWidths, "width");

    mRelativeDirections = IloArray<IloArray<IloBoolVarArray>>(env, numItems);

    mItemsOverlapsXY    = IloArray<IloBoolVarArray>(env, numItems);
    mSupportXY          = IloArray<IloBoolVarArray>(env, numItems);
    mOverlapAreasXY     = IloArray<IloIntVarArray>(env, numItems);
    mOrientation        = IloArray<IloBoolVarArray>(env, numItems);

    mPlacedOnFloor = IloBoolVarArray(env, numItems);
    NameVars(mPlacedOnFloor, "placedOnFloor");

    for(int i=0; i < numItems; ++i)
    {
        mRelativeDirections[i] = IloArray<IloBoolVarArray>(env, numItems);
        char name[100];


        mItemsOverlapsXY[i] = IloBoolVarArray(env, numItems);
        snprintf(name, 100, "over[%i]", i);
        NameVars(mItemsOverlapsXY[i], name);


        mSupportXY[i] = IloBoolVarArray(env, numItems);
        snprintf(name, 100, "sup[%i]", i);
        NameVars(mSupportXY[i], name);


        mOverlapAreasXY[i] = IloIntVarArray(env, numItems, 0, 0);
        snprintf(name, 100, "overA[%i]", i);
        NameVars(mOverlapAreasXY[i], name);

        mOrientation[i] = IloBoolVarArray(env, vetRot.size());

        for(int r:vetRot)
        {
            char name[100];
            snprintf(name, 100, "orient[%i][%i]", i, r);
            mOrientation[i][r].setName(name);
        }

        //return;

        for(int j=0; j < numItems; ++j)
        {
            if(i == j)
                continue;
            char name[100];


            mRelativeDirections[i][j] = IloBoolVarArray(env, 6);

            snprintf(name, 100, "relaD[%i][%i][RightY]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::RightY].setName(name);

            snprintf(name, 100, "relaD[%i][%i][LeftY]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::LeftY].setName(name);

            snprintf(name, 100, "relaD[%i][%i][InFrontX]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::InFrontX].setName(name);

            snprintf(name, 100, "relaD[%i][%i][BehindX]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::BehindX].setName(name);

            snprintf(name, 100, "relaD[%i][%i][BelowZ]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::BelowZ].setName(name);

            snprintf(name, 100, "relaD[%i][%i][AboveZ]", i, j);
            mRelativeDirections[i][j][ContainerLoading::Model::AboveZ].setName(name);
        }
    }

    for(int i=0; i < numItems; ++i)
    {
        const Item& item = instanciaG.vetItens[vetItems[i]];

        IloIntArray arrayXy(env, 2);
        arrayXy[0] = item.vetDim[0];
        arrayXy[1] = item.vetDim[1];

        IloIntArray arrayZ(env, 1);
        arrayZ[0] = item.vetDim[2];

        int dimMax = std::max(arrayXy[0], arrayXy[1]);
        int dimMin = std::min(arrayXy[0], arrayXy[1]);

        mLengths[i].setBounds(dimMin, dimMax);
        mWidths[i].setBounds(dimMin, dimMax);
        mHeights[i].setBounds(arrayZ[0], arrayZ[0]);

        mLengths[i].setPossibleValues(arrayXy);
        mWidths[i].setPossibleValues(arrayXy);

        arrayZ[0] = item.vetDim[2];
        mHeights[i].setPossibleValues(arrayZ);

    }

    //IloConstraint a = mWidths[0] + floatVet[0] <= instanciaG.vetItens[vetItems[0]].vetDim[0];
    //a.setName("Name");
    //model.add(mWidths[0] + floatVet[0] <= instanciaG.vetItens[vetItems[0]].vetDim[0]);
    //model.add(IloIfThen(env, mWidths[0] <= 5000, mWidths[1] == 800));

    //IloCP cp(model);

    /*
    if(cp.solve())
        std::printf("Sucess!\n");
    else
        std::printf("!Sucess\n");
    */

    std::printf("END createVariables\n");
}

// Checked 5
void CpOptimizer::CreateNoOverlap()
{

    std::printf("Beging CreateNoOverlap");

    for (size_t i = 0; i < numItems; ++i)
    {
        for (size_t j = 0; j < numItems; ++j)
        {
            if (i == j)
            {
                continue;
            }

            IloExpr sum(*envPtr);

            for (size_t d = 0; d < mDimensions.size(); ++d)
            {
                const Dimension& dimension = mDimensions[d];
                //if(dimension.Type != AxisZ)
                {
                    sum += mRelativeDirections[i][j][dimension.FirstDirection] +
                           mRelativeDirections[i][j][dimension.SecondDirection];
                }

                auto startPositionI = getIntVars(dimension.Type, true, i);
                auto startPositionJ = getIntVars(dimension.Type, true, j);

                auto endPositionI   = getIntVars(dimension.Type, false, i);
                auto endPositionJ   = getIntVars(dimension.Type, false, j);

                model.add(((mRelativeDirections[i][j][dimension.FirstDirection]==1) && (endPositionJ <= startPositionI)) ||
                          (mRelativeDirections[i][j][dimension.FirstDirection]==0) && (endPositionJ > startPositionI));
                //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][dimension.FirstDirection]==1,
                //                    endPositionJ <= startPositionI));


                //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][dimension.FirstDirection] ==0,
                //                    endPositionJ > startPositionI));


                model.add(((mRelativeDirections[i][j][dimension.SecondDirection]==1) && (endPositionI <= startPositionJ)) ||
                          ((mRelativeDirections[i][j][dimension.SecondDirection]==0) && (endPositionI > startPositionJ)));
                //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][dimension.SecondDirection]==1,
                //                    endPositionI <= startPositionJ));


                //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][dimension.SecondDirection]==0,
                //                    endPositionI > startPositionJ));

                model.add(mRelativeDirections[i][j][dimension.FirstDirection] ==
                          mRelativeDirections[j][i][dimension.SecondDirection]);

                /*
                model.add(mRelativeDirections[i][j][dimension.FirstDirection] ==
                          mRelativeDirections[j][i][dimension.SecondDirection]);

                model.add(mRelativeDirections[i][j][dimension.SecondDirection] ==
                          mRelativeDirections[j][i][dimension.FirstDirection]);
                */
                //model.add(mRelativeDirections[i][j][dimension.FirstDirection] +
                //          mRelativeDirections[i][j][dimension.SecondDirection] <= 1);

            }

            // No overlap constraints
            model.add(sum >= 1);
        }
    }

    std::printf("END CreateNoOverlap");
}


/*
void CpOptimizer::CreateNoOverlap()
{
    std::printf("Begin CreateNoOverlap\n");

    int Mx = instanciaG.vetDimVeiculo[0];
    int My = instanciaG.vetDimVeiculo[1];
    int Mz = instanciaG.vetDimVeiculo[2];

    for (size_t i = 0; i < numItems; ++i)
    {
        for (size_t j = i + 1; j < numItems; ++j)
        {
            {
                auto left  = mRelativeDirections[i][j][BehindX];   // i before j
                auto right = mRelativeDirections[i][j][InFrontX];  // i after j

                model.add(mEndPositionsX[i] <= mStartPositionsX[j] + Mx * (1 - left));
                model.add(mEndPositionsX[j] <= mStartPositionsX[i] + Mx * (1 - right));

                model.add(left + right <= 1);

                // symmetry
                model.add(left  == mRelativeDirections[j][i][InFrontX]);
                model.add(right == mRelativeDirections[j][i][BehindX]);
            }

            {
                auto left  = mRelativeDirections[i][j][LeftY];
                auto right = mRelativeDirections[i][j][RightY];

                model.add(mEndPositionsY[i] <= mStartPositionsY[j] + My * (1 - left));
                model.add(mEndPositionsY[j] <= mStartPositionsY[i] + My * (1 - right));

                model.add(left + right <= 1);

                // symmetry
                model.add(left  == mRelativeDirections[j][i][RightY]);
                model.add(right == mRelativeDirections[j][i][LeftY]);
            }

            {
                auto below = mRelativeDirections[i][j][BelowZ];
                auto above = mRelativeDirections[i][j][AboveZ];

                model.add(mEndPositionsZ[i] <= mStartPositionsZ[j] + Mz * (1 - below));
                model.add(mEndPositionsZ[j] <= mStartPositionsZ[i] + Mz * (1 - above));

                model.add(below + above <= 1);

                // symmetry
                model.add(below == mRelativeDirections[j][i][AboveZ]);
                model.add(above == mRelativeDirections[j][i][BelowZ]);
            }


            IloExpr sepXY(*envPtr);
            sepXY += mRelativeDirections[i][j][LeftY];
            sepXY += mRelativeDirections[i][j][RightY];
            sepXY += mRelativeDirections[i][j][BehindX];
            sepXY += mRelativeDirections[i][j][InFrontX];

            IloExpr sepZ(*envPtr);
            sepZ += mRelativeDirections[i][j][AboveZ];
            sepZ += mRelativeDirections[i][j][BelowZ];

            model.add(sepXY + sepZ >= 1);

            sepXY.end();
            sepZ.end();
        }
    }

    std::printf("End CreateNoOverlap\n");
}
*/


void CpOptimizer::CreateItemOrientations()
{

    for (int i = 0; i < numItems; ++i)
    {
        Item& item = instanciaG.vetItens[vetItems[i]];
        int dx, dy, dz;

        IloExpr sum(*envPtr);

        for(InstanceNS::Rotation r:vetRot)
        {
            dx = item.getDimRotacionada(0, r);
            dy = item.getDimRotacionada(1, r);
            dz = item.getDimRotacionada(2, r);

            model.add(mOrientation[i][r]==0 || mLengths[i]==dx);
            model.add(mOrientation[i][r]==0 || mWidths[i]==dy);
            model.add(mOrientation[i][r]==0 || mHeights[i]==dz);

            //model.add(IloIfThen(*envPtr, mOrientation[i][r]==1, mLengths[i]==dx));
            //model.add(IloIfThen(*envPtr, mOrientation[i][r]==1, mWidths[i]==dy));
            //model.add(IloIfThen(*envPtr, mOrientation[i][r]==1, mHeights[i]==dz));


            sum +=mOrientation[i][r];
        }

        model.add(sum == 1);

        model.add(mStartPositionsX[i]+mLengths[i] == mEndPositionsX[i]);
        model.add(mStartPositionsY[i]+mWidths[i] == mEndPositionsY[i]);
        model.add(mStartPositionsZ[i]+mHeights[i] == mEndPositionsZ[i]);

        /*
            mModelCP.AddEquality(mLengths[i], itemLength).OnlyEnforceIf(mOrientation[i][o]);
            mModelCP.AddEquality(mWidths[i], itemWidth).OnlyEnforceIf(mOrientation[i][o]);
            mModelCP.AddEquality(mHeights[i], itemHeight).OnlyEnforceIf(mOrientation[i][o]);

        mModelCP.AddExactlyOne(mOrientation[i]);
        */
    }

}

// Checked 3
void CpOptimizer::CreateSupportItem()
{
    for (size_t i = 0; i < numItems; ++i)
    {
        model.add(mSupportXY[i][i] == 0);
        for (size_t j = 0; j < numItems; ++j)
        {
            if (i == j)
            {
                continue;
            }

            IloBoolVar isVerticallyAdjacent(*envPtr);
            //model.add(IloIfThen(*envPtr, isVerticallyAdjacent==1, mEndPositionsZ[j]==mStartPositionsZ[i]));
            model.add(((isVerticallyAdjacent==1) && (mEndPositionsZ[j]==mStartPositionsZ[i])) ||
                       (isVerticallyAdjacent==0) && (mEndPositionsZ[j]!=mStartPositionsZ[i]));

            //model.add(IloIfThen(*envPtr, isVerticallyAdjacent==0, mEndPositionsZ[j]!=mStartPositionsZ[i]));

            model.add(((isVerticallyAdjacent==0) && (mSupportXY[i][j]==0)) || (isVerticallyAdjacent==1));

            //model.add(IloIfThen(*envPtr, isVerticallyAdjacent==0, mSupportXY[i][j]==0));

            if (i < j)
            {
                auto position = j - i - 1;

                model.add(mSupportXY[i][j] + (1-isVerticallyAdjacent) + (1-mItemsOverlapsXY[i][position]) >= 1);

                model.add(((mItemsOverlapsXY[i][position]==0) && mSupportXY[i][j]==0) || (mItemsOverlapsXY[i][position]==1));
//                model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][position]==0, mSupportXY[i][j]==0));

            }
            else
            {
                auto position = i - j - 1;

                model.add(mSupportXY[i][j] + (1-isVerticallyAdjacent) + (1-mItemsOverlapsXY[j][position]) >= 1);
                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[j][position]==0, mSupportXY[i][j]==0));
                 model.add(((mItemsOverlapsXY[j][position]==0) && mSupportXY[i][j]==0) || (mItemsOverlapsXY[j][position]==1));
            }
        }
    }
}


/*
void CpOptimizer::CreateSupportItem()
{
    int M = instanciaG.vetDimVeiculo[2]; // Big-M for Z dimension

    for (size_t i = 0; i < numItems; ++i)
    {
        model.add(mSupportXY[i][i] == 0);

        for (size_t j = 0; j < numItems; ++j)
        {
            if (i == j)
                continue;

            IloBoolVar isVerticallyAdjacent(*envPtr);

            model.add(mEndPositionsZ[j] - mStartPositionsZ[i] <= M * (1 - isVerticallyAdjacent));
            model.add(mStartPositionsZ[i] - mEndPositionsZ[j] <= M * (1 - isVerticallyAdjacent));

            IloBoolVar overlapVar;

            if (i < j)
            {
                auto position = j - i - 1;
                overlapVar = mItemsOverlapsXY[i][position];
            }
            else
            {
                auto position = i - j - 1;
                overlapVar = mItemsOverlapsXY[j][position];
            }

            model.add(mSupportXY[i][j] <= isVerticallyAdjacent);
            model.add(mSupportXY[i][j] <= overlapVar);
            model.add(mSupportXY[i][j] >= isVerticallyAdjacent + overlapVar - 1);

        }
    }
}
*/
// Checked 3
void CpOptimizer::CreateSupportArea()
{

    for (size_t i = 0; i < numItems; ++i)
    {
        IloExpr supportedAreaExpr(*envPtr);
        int areaI = instanciaG.vetItens[vetItems[i]].vetDim[0] * instanciaG.vetItens[vetItems[i]].vetDim[1];
        for (size_t j = 0; j < numItems; ++j)
        {
            if (i == j)
            {
                continue;
            }


            int areaJ = instanciaG.vetItens[vetItems[j]].vetDim[0] * instanciaG.vetItens[vetItems[j]].vetDim[1];
            int minArea = std::min(areaI, areaJ);

            IloIntVar usableArea(*envPtr, 0, minArea);
            if (i < j)
            {
                auto position = j - i - 1;
                model.add(usableArea == mOverlapAreasXY[i][position] * mSupportXY[i][j]);
            }
            else
            {
                auto position = i - j - 1;
                model.add(usableArea == mOverlapAreasXY[j][position] * mSupportXY[i][j]);
            }

//            model.add(IloIfThen(*envPtr, mSupportXY[i][j]==0, usableArea==0));
            model.add(((mSupportXY[i][j]==0) && (usableArea==0)) || (mSupportXY[i][j] == 1));
            supportedAreaExpr += usableArea;

        }

        //operations_research::sat::IntVar supportedArea = mModelCP.NewIntVar({0, areaI});

        //IloIntVar supportedArea(*envPtr, 0, areaI);
        //model.add(IloIfThen(*envPtr, mPlacedOnFloor[i] == 0, supportedArea == supportedAreaExpr));

        int minArea = (int)std::ceil(input.minSupportArea*instanciaG.vetItens[vetItems[i]].vetDim[0]*
                                     instanciaG.vetItens[vetItems[i]].vetDim[1]);

        model.add(((mPlacedOnFloor[i] == 0) && (supportedAreaExpr >= minArea)) || (mPlacedOnFloor[i] == 1));

//        model.add(IloIfThen(*envPtr, mPlacedOnFloor[i]==0, supportedAreaExpr >= minArea));
        //model.add(supportedAreaExpr <= areaI * (1 - mPlacedOnFloor[i]));
        //model.add(IloIfThen(*envPtr, supportedAreaExpr<minArea, mPlacedOnFloor[i]==1));
    }
}

// Checked 3
void CpOptimizer::CreateXYIntersectionBool()
{
    for (size_t i = 0; i < numItems - 1; ++i)
    {
        for (size_t j = i + 1; j < numItems; ++j)
        {
            auto positionJ = j - i - 1;

            // TODO Posible wrong!
            model.add(mItemsOverlapsXY[i][positionJ] + mRelativeDirections[i][j][LeftY] +
                      mRelativeDirections[i][j][RightY] + mRelativeDirections[i][j][BehindX] +
                      mRelativeDirections[i][j][InFrontX] >= 1);


            //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][LeftY]==1, mItemsOverlapsXY[i][positionJ]==0));
            model.add(((mRelativeDirections[i][j][LeftY]==1) && (mItemsOverlapsXY[i][positionJ]==0)) ||
                        (mRelativeDirections[i][j][LeftY]==0));
            model.add(((mRelativeDirections[i][j][RightY]==1) && (mItemsOverlapsXY[i][positionJ]==0)) ||
                      (mRelativeDirections[i][j][RightY]==0));
            model.add(((mRelativeDirections[i][j][BehindX]==1) && (mItemsOverlapsXY[i][positionJ]==0)) ||
                      (mRelativeDirections[i][j][BehindX]==0));
            model.add(((mRelativeDirections[i][j][InFrontX]==1) && (mItemsOverlapsXY[i][positionJ]==0)) ||
                      (mRelativeDirections[i][j][InFrontX]==0));
            //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][RightY]==1, mItemsOverlapsXY[i][positionJ]==0));
            //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][BehindX]==1, mItemsOverlapsXY[i][positionJ]==0));
            //model.add(IloIfThen(*envPtr, mRelativeDirections[i][j][InFrontX]==1, mItemsOverlapsXY[i][positionJ]==0));
        }
    }
}

// Checked 2
void CpOptimizer::CreateXYIntersectionArea()
{

    for (size_t i = 0; i < numItems - 1; ++i)
    {
        for (size_t j = i + 1; j < numItems; ++j)
        {
            if (instanciaG.vetItens[vetItems[i]].vetDim[2] + instanciaG.vetItens[vetItems[j]].vetDim[2] <=
                instanciaG.vetDimVeiculo[2])
            {
                // Overlap in x
                auto positionJ = j - i - 1;
                IloIntVar diffXij(*envPtr, 0, instanciaG.vetDimVeiculo[0]);
                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ]==1, diffXij ==
                //                                                              (mEndPositionsX[i]-mStartPositionsX[j])));
                model.add(((mItemsOverlapsXY[i][positionJ]==1) && (diffXij==mEndPositionsX[i]-mStartPositionsX[j])) ||
                           (mItemsOverlapsXY[i][positionJ]==0));

                IloIntVar diffXji(*envPtr, 0, instanciaG.vetDimVeiculo[0]);

                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ]==1, diffXji ==
                //                                                                (mEndPositionsX[j]-mStartPositionsX[i])));
                model.add(((mItemsOverlapsXY[i][positionJ]==1) && (diffXji==mEndPositionsX[j]-mStartPositionsX[i])) ||
                          (mItemsOverlapsXY[i][positionJ]==0));

                IloExpr min = IloMin(IloMin(diffXij, diffXji),  IloMin(mLengths[i], mLengths[j]));

                IloIntVar xOverlap(*envPtr, 0, instanciaG.vetDimVeiculo[0]);
                model.add(((mItemsOverlapsXY[i][positionJ] == 1) && (xOverlap == min)) ||
                           (mItemsOverlapsXY[i][positionJ] == 0));

                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ] == 1,
                //                    xOverlap == IloMin(IloMin(diffXij, diffXji),  IloMin(mLengths[i], mLengths[j]))));

                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ] == 0, xOverlap == 0));
                model.add(((mItemsOverlapsXY[i][positionJ] == 0) && (xOverlap == 0)) ||
                           (mItemsOverlapsXY[i][positionJ] == 1));

                // Overlap in y

                IloIntVar diffYij(*envPtr, 0, instanciaG.vetDimVeiculo[1]);
                model.add(((mItemsOverlapsXY[i][positionJ]==1) && (diffYij == (mEndPositionsY[i]-mStartPositionsY[j]))) ||
                           (mItemsOverlapsXY[i][positionJ]==0));
  //                model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ]==1, diffYij ==
  //                                                                              (mEndPositionsY[i]-mStartPositionsY[j])));

                IloIntVar diffYji(*envPtr, 0, instanciaG.vetDimVeiculo[1]);
                model.add(((mItemsOverlapsXY[i][positionJ]==1) && (diffYji == (mEndPositionsY[j]-mStartPositionsY[i]))) ||
                          (mItemsOverlapsXY[i][positionJ]==0));

                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ]==1, diffYji ==
                //                                                                (mEndPositionsY[j]-mStartPositionsY[i])));

                IloIntVar yOverlap(*envPtr, 0, instanciaG.vetDimVeiculo[1]);
                IloExpr min2 = IloMin(IloMin(diffYij, diffYji),  IloMin(mWidths[i], mWidths[j]));
                model.add(((mItemsOverlapsXY[i][positionJ] == 1) && (yOverlap == min2)) || (mItemsOverlapsXY[i][positionJ] == 0));

//                model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ] == 1,
//                                    yOverlap == IloMin(IloMin(diffYij, diffYji),  IloMin(mWidths[i], mWidths[j]))));
                //model.add(yOverlap == IloMin(IloMin(diffYij, diffYji),  IloMin(mWidths[i], mWidths[j])));
                //model.add(IloIfThen(*envPtr, mItemsOverlapsXY[i][positionJ] == 0, yOverlap == 0));
                model.add(((mItemsOverlapsXY[i][positionJ] == 0) && (yOverlap == 0)) || mItemsOverlapsXY[i][positionJ] == 1);

                model.add(mOverlapAreasXY[i][positionJ] == xOverlap*yOverlap);


            }
        }
    }
}

void CpOptimizer::CreateLifo()
{

}


void CpOptimizer::CreateOnFloorConstraints()
{

    for (size_t i = 0; i < numItems; ++i)
    {
        model.add(IloIfThen(*envPtr, mPlacedOnFloor[i] == 1, mStartPositionsZ[i] == 0));        
        model.add(IloIfThen(*envPtr, mPlacedOnFloor[i] == 0, mStartPositionsZ[i] > 0));

        model.add(((mPlacedOnFloor[i] == 1) && (mStartPositionsZ[i] == 0)) ||
                   (mPlacedOnFloor[i] == 0) && (mStartPositionsZ[i] != 0));
    }

}

void CpOptimizer::CreateAxleWeights()
{

}

IloIntVar& CpOptimizer::getIntVars(DimensionType dimension, bool first, int i)
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

void IBM_CpOptimizerNS::testSCIP()
{

    scippp::Model model("Simple");
    auto x1 = model.addVar("x_1", 1);
    auto x2 = model.addVar("x_2", 1);
    model.addConstr(3 * x1 + 2 * x2 <= 1, "capacity");
    model.setObjsense(Sense::MAXIMIZE);
    model.solve();

}
