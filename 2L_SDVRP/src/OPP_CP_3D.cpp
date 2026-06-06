#include "OPP_CP_3D.h"

#include "AxleWeights.h"
#include "InputOutput.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <ostream>
#include <ranges>
#include <stdexcept>
#include <string>

namespace ContainerLoading
{
using namespace Model;
using namespace ParseInputNS;
using namespace AxleWeightsNS;

namespace Algorithms
{
void ContainerLoadingCP::BuildModel()
{
    GravityMM = 1;

    CreateVariables();

    AddConstraints();

    ////AddObjective();
}

LoadingStatus ContainerLoadingCP::Solve()
{
    BuildModel();

    operations_research::sat::SatParameters parameters;
    SetParameters(parameters);

    operations_research::sat::Model model = operations_research::sat::Model();
    model.Add(operations_research::sat::NewSatParameters(parameters));

    operations_research::sat::CpModelProto protoModel = mModelCP.Build();
    ////auto validationResponse = operations_research::sat::ValidateCpModel(protoModel);
    ////LOG(INFO) << validationResponse;

    mResponse = operations_research::sat::SolveCpModel(protoModel, &model);

    ////LOG(INFO) << operations_research::sat::CpSolverResponseStats(mResponse);

    GravityMM = GravityMM_const;

    ////PrintSolution(Items, mResponse);
    switch(mResponse.status())
    {
    case operations_research::sat::OPTIMAL:
        [[fallthrough]];
    case operations_research::sat::FEASIBLE:
        return LoadingStatus::FeasOpt;
    case operations_research::sat::INFEASIBLE:
        return LoadingStatus::Infeasible;
    case operations_research::sat::UNKNOWN:
        return LoadingStatus::Unknown;
    default:
        return LoadingStatus::Invalid;
    }
}

void ContainerLoadingCP::WriteProtoModel(
    const operations_research::sat::CpModelProto &protoModel) const
{
    std::string   protoModelString = protoModel.DebugString();
    std::ofstream file("protoModel_basicDomain.txt");
    file << protoModelString;
    file.close();
}

void ContainerLoadingCP::SetParameters(
    operations_research::sat::SatParameters &parameters) const
{
    parameters.set_num_search_workers(mParams.Threads);
    parameters.set_log_search_progress(mParams.LogFlag);
    ////parameters.set_search_branching(parameters.PORTFOLIO_SEARCH);
    parameters.set_max_time_in_seconds(mMaxRuntime);
    parameters.set_stop_after_first_solution(true);
    // Setting seed value is without effect for parallel mode
    // https://github.com/google/or-tools/issues/2793
    ////parameters.set_random_seed(mParams.Seed);

    ////parameters.set_cp_model_presolve(false);

    // Example how to set decision strategies, but not helpful here.
    /*
    cp_model.AddDecisionStrategy(mStartPositionsY,
        operations_research::sat::DecisionStrategyProto::VariableSelectionStrategy::DecisionStrategyProto_VariableSelectionStrategy_CHOOSE_LOWEST_MIN,
        operations_research::sat::DecisionStrategyProto::DomainReductionStrategy::DecisionStrategyProto_DomainReductionStrategy_SELECT_MIN_VALUE);
    mModelCP.AddDecisionStrategy(mStartPositionsX,
        operations_research::sat::DecisionStrategyProto::VariableSelectionStrategy::DecisionStrategyProto_VariableSelectionStrategy_CHOOSE_LOWEST_MIN,
        operations_research::sat::DecisionStrategyProto::DomainReductionStrategy::DecisionStrategyProto_DomainReductionStrategy_SELECT_MIN_VALUE);
    cp_model.AddDecisionStrategy(mStartPositionsZ,
        operations_research::sat::DecisionStrategyProto::VariableSelectionStrategy::DecisionStrategyProto_VariableSelectionStrategy_CHOOSE_LOWEST_MIN,
        operations_research::sat::DecisionStrategyProto::DomainReductionStrategy::DecisionStrategyProto_DomainReductionStrategy_SELECT_MIN_VALUE);
        */

    // Example how to set solution observer
    /*
    auto observer = operations_research::sat::NewFeasibleSolutionObserver([&](const
    operations_research::sat::CpSolverResponse mResponse) { LOG(INFO) << "Length " <<
    operations_research::sat::SolutionIntegerValue(mResponse, mMaxLength); LOG(INFO) <<
    "LB " << mResponse.best_objective_bound();
        });
    model.Add(observer);
    */
}

/*
void ContainerLoadingCP::PrintSolution()
{
    for (size_t i = 0; i < mItems.size(); ++i)
    {
        LOG(INFO) << "Item " << std::to_string(i) << ": ("
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mStartPositionsX[i]) << ","
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mStartPositionsY[i]) << ","
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mStartPositionsZ[i]) << ") | ("
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mLengths[i]) << ","
                  << operations_research::sat::SolutionIntegerValue(mResponse, mWidths[i])
<< ","
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mHeights[i]) << ") | ("
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mEndPositionsX[i]) << ","
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mEndPositionsY[i]) << ","
                  << operations_research::sat::SolutionIntegerValue(mResponse,
mEndPositionsZ[i]) << ")";
    }

    for (size_t i = 0; i < mItems.size(); ++i)
    {
        mRelativeDirections.emplace_back();
        for (size_t j = 0; j < mItems.size(); ++j)
        {
            mRelativeDirections[i].emplace_back();
            std::stringstream vars;
            vars << std::to_string(i) << "_" << std::to_string(j) << ": ";
            for (const auto& dimension: mDimensions)
            {
                vars << operations_research::sat::SolutionBooleanValue(
                    mResponse, mRelativeDirections[i][j][dimension.FirstDirection])
                     << " ";
                vars << operations_research::sat::SolutionBooleanValue(
                    mResponse, mRelativeDirections[i][j][dimension.SecondDirection])
                     << " ";
            }

            vars << " | ";
            vars << operations_research::sat::SolutionIntegerValue(mResponse,
mOverlapAreasXY[i][j]); LOG(INFO) << vars.str();
        }
    }
}
*/

void ContainerLoadingCP::PrintSolution(std::vector<Array<int, 4>> &vetPos)
{
    vetPos = std::vector<Array<int, 4>>();
    vetPos.reserve(mItems.size());

    Array<int, 4> array;

    int sumLeft =
        operations_research::sat::SolutionIntegerValue(mResponse, sumLeftBalancedLoading);
    std::printf("sumLeft: %f\n", sumLeft / (double)scaleBalancedLoading);

    int sumRight = operations_research::sat::SolutionIntegerValue(
        mResponse, sumRightBalancedLoading);
    std::printf("sumRight: %f\n", sumRight / (double)scaleBalancedLoading);

    static std::map<Orientation, int> mapOritentationRotation = {
        {NoRotation, 0}, {RotationZ, 1}, {RotationX, 2}};

    if(input.axleWights)
    {
        int fK = operations_research::sat::SolutionIntegerValue(mResponse, forceK);
        int fFA = operations_research::sat::SolutionIntegerValue(mResponse, forceFA);
        int fRA = operations_research::sat::SolutionIntegerValue(mResponse, forceRA);
        int fTA = operations_research::sat::SolutionIntegerValue(mResponse, forceTA);
        std::printf("fK: %d; fFA: %d; fRA: %d; fTA: %d\n", fK, fFA, fRA, fTA);
    }

    for(size_t i = 0; i < mItems.size(); ++i)
    {
        bool top = operations_research::sat::SolutionBooleanValue(mResponse, mTopBool[i]);
        if(top)
            std::printf("Item(%d) is topItem\n", (int)mItems[i].ExternId);

        // if (input.axleWights) {
        //     std::cout<<"item: "<<i<<"\n";
        //     std::cout<<"mStartPositionsX:
        //     "<<operations_research::sat::SolutionIntegerValue(mResponse,
        //     mStartPositionsX[i])<<std::endl; std::cout<<"mR:
        //     "<<operations_research::sat::SolutionIntegerValue(mResponse,
        //     mR[i])<<std::endl; std::cout<<"MR
        //     CHECK"<<2*operations_research::sat::SolutionIntegerValue(mResponse, mR[i])
        //     << " " << 2*semiTrailer.distanceCargoSpaceTrailerAxle
        //     -2*operations_research::sat::SolutionIntegerValue(mResponse,
        //     mStartPositionsX[i]) -
        //     operations_research::sat::SolutionIntegerValue(mResponse,
        //     mLengths[i])<<std::endl;
        // }
        array[0] = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsX[i]);
        array[1] = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsY[i]);
        array[2] = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsZ[i]);

        Array<int, 3> arrayEnd;
        arrayEnd[0] =
            operations_research::sat::SolutionIntegerValue(mResponse, mEndPositionsX[i]);
        arrayEnd[1] =
            operations_research::sat::SolutionIntegerValue(mResponse, mEndPositionsY[i]);
        arrayEnd[2] =
            operations_research::sat::SolutionIntegerValue(mResponse, mEndPositionsZ[i]);

        int dx, dy, dz;
        dx = operations_research::sat::SolutionIntegerValue(mResponse, mLengths[i]);
        dy = operations_research::sat::SolutionIntegerValue(mResponse, mWidths[i]);
        dz = operations_research::sat::SolutionIntegerValue(mResponse, mHeights[i]);

        if((array[0] + dx) != arrayEnd[0])
        {
            std::cout << "Erro: start0(" << array[0] << ") + dx(" << dx
                      << "): " << array[0] + dx << " != " << arrayEnd[0] << "\n";
        }

        if((array[1] + dy) != arrayEnd[1])
        {
            std::cout << "Erro: start0(" << array[1] << ") + dy(" << dy
                      << "): " << array[1] + dy << " != " << arrayEnd[1] << "\n";
        }

        if((array[2] + dz) != arrayEnd[2])
        {
            std::cout << "Erro: start0(" << array[2] << ") + dz(" << dz
                      << "): " << array[2] + dz << " != " << arrayEnd[2] << "\n";
        }

        for(size_t o = 0; o < mItemOrientations.size(); ++o)
        {
            bool util = operations_research::sat::SolutionBooleanValue(
                mResponse, mOrientation[i][o]);
            if(util)
            {
                array[3] = (int)mapOritentationRotation[mItemOrientations[o]];
                break;
            }
        }

        vetPos.push_back(array);

        /*
        for(int j=0; j < mItems.size(); ++j)
        {
            if(i == j)
                std::cout<<"0 ";

            bool sup = operations_research::sat::SolutionBooleanValue(mResponse,
        mSupportXY[i][j]); std::cout<<sup<<" ";
        }

        std::cout<<"\n";
        */
    }

    std::printf("mMinX: %d",
                (int)operations_research::sat::SolutionIntegerValue(mResponse, mMinX));

    /*
    size_t numberOfItems = mItems.size();
    for (size_t i = 0; i < numberOfItems; ++i)
    {
        for (size_t j = 0; j < numberOfItems; ++j)
        {
            if (i == j || mItems[i].pos >= mItems[j].pos)
            {
                continue;
            }

            for (size_t d = 0; d < mDimensions.size(); ++d)
            {
                const Dimension& dimension = mDimensions[d];

                bool fisrtDir =  operations_research::sat::SolutionBooleanValue(mResponse,
    mRelativeDirections[i][j][dimension.FirstDirection]); bool secondDir =
    operations_research::sat::SolutionBooleanValue(mResponse,
    mRelativeDirections[i][j][dimension.SecondDirection]);


                switch (d)
                {
                case AxisX:
                    // InFrontX, BehindX
                    std::printf("\tmRelativeDirections[%d][%d][InFrontX]: %d", (int)i,
    (int)j, (int)fisrtDir);
                    //std::printf("\tmRelativeDirections[%d][%d][BehindX]: %d", (int)i,
    (int)j, (int)secondDir); break;

                case AxisY:
                    // AxisY, RightY, LeftY
                    std::printf("\tmRelativeDirections[%d][%d][RightY]: %d", (int)i,
    (int)j, (int)fisrtDir); std::printf("\tmRelativeDirections[%d][%d][LeftY]: %d",
    (int)i, (int)j, (int)secondDir); break;

                 case AxisZ:
                    // AxisZ, AboveZ, BelowZ
                    std::printf("\tmRelativeDirections[%d][%d][AboveZ]: %d", (int)i,
    (int)j, (int)fisrtDir);
                    //std::printf("\tmRelativeDirections[%d][%d][LeftY]: %d", (int)i,
    (int)j, (int)secondDir); break; default: break;
                }

            }

            std::cout<<"\n";

        }
        std::cout<<"\n\n";
    }


    for(size_t i=0; i < numberOfItems; ++i)
    {
        int z = operations_research::sat::SolutionIntegerValue(mResponse,
    mStartPositionsZ[i]); std::printf("%d: %d\n", (int)i, z);
    }

    std::cout<<"\n";
  */
}

void ContainerLoadingCP::ExtractPacking(std::vector<Cuboid> &items) const
{
    for(size_t i = 0; i < items.size(); ++i)
    {
        auto &item = items[i];

        item.Rotated = (Rotation)operations_research::sat::SolutionBooleanValue(
            mResponse, mOrientation[i][RotationZ]);

        item.X = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsX[i]);
        item.Y = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsY[i]);
        item.Z = operations_research::sat::SolutionIntegerValue(
            mResponse, mStartPositionsZ[i]);
    }
}

std::vector<int> ContainerLoadingCP::ExtractSequence() const
{
    std::vector<std::tuple<int, int>> assignments;
    assignments.reserve(mNumberCustomers);
    for(size_t i = 0; i < mNumberCustomers; ++i)
    {
        assignments.emplace_back(operations_research::sat::SolutionIntegerValue(
                                     mResponse, mCustomerPosition[i]),
                                 i);
    }

    std::ranges::sort(assignments);

    std::vector<int> sequence = std::vector<int>();
    sequence.reserve(assignments.size());
    for(const auto &[position, id] : assignments)
    {
        sequence.push_back(id);
    }

    return sequence;
}

void ContainerLoadingCP::CreateVariables()
{
    size_t numberOfItems = mItems.size();

    if(input.axleWights)
    {
        // cleaned this up, shouldnt make a difference
        int maxK = semiTrailer.maxMassRearAxle + semiTrailer.maxMassFrontAxle;
        forceK =
            mModelCP.NewIntVar({-scale * GravityMM * maxK, scale * GravityMM * maxK});
        forceRA = mModelCP.NewIntVar({-scale * GravityMM * semiTrailer.maxMassRearAxle,
                                      scale * GravityMM * semiTrailer.maxMassRearAxle});
        forceFA = mModelCP.NewIntVar({-scale * GravityMM * semiTrailer.maxMassFrontAxle,
                                      scale * GravityMM * semiTrailer.maxMassFrontAxle});
        forceTA =
            mModelCP.NewIntVar({-scale * GravityMM * semiTrailer.maxMassTrailerAxle,
                                scale * GravityMM * semiTrailer.maxMassTrailerAxle});

        // std::printf("max FK: %d\nmax FRA: %d\nmax FA: %d\nmax FTA: %d", (int)maxFK,
        // GravityMM*semiTrailer.maxMassRearAxle,
        //             GravityMM*semiTrailer.maxMassFrontAxle,
        //             GravityMM*semiTrailer.maxMassTrailerAxle);
    }

    int limitWight =
        InstanceNS::instanciaG.maxPayload * input.balancedLoadingD * scaleBalancedLoading;
    sumRightBalancedLoading = mModelCP.NewIntVar({0, limitWight});
    sumLeftBalancedLoading = mModelCP.NewIntVar({0, limitWight});

    std::vector<Cuboid> itemCopy;
    itemCopy.reserve(mItems.size());
    for(const auto &item : mItems)
    {
        itemCopy.push_back(item);
    }

    auto placementPointsPerType =
        PlacementPointGenerator::GeneratePlacementPatterns(mContainer,
                                                           itemCopy,
                                                           mPlacementPatternTypeX,
                                                           mPlacementPatternTypeY,
                                                           mPlacementPatternTypeZ);

    mStartPositionsX.reserve(numberOfItems);
    mEndPositionsX.reserve(numberOfItems);
    mStartPositionsY.reserve(numberOfItems);
    mEndPositionsY.reserve(numberOfItems);
    mStartPositionsZ.reserve(numberOfItems);
    mEndPositionsZ.reserve(numberOfItems);

    for(size_t i = 0; i < numberOfItems; i++)
    {
        const Cuboid &item = mItems[i];
        int           minLength =
            item.EnableHorizontalRotation ? std::min(item.Dx, item.Dy) : item.Dx;
        int minWidth =
            item.EnableHorizontalRotation ? std::min(item.Dx, item.Dy) : item.Dy;

        // operations_research::Domain::FromValues(placementPointsPerType[item].X
        mStartPositionsX.emplace_back(mModelCP.NewIntVar({0, mContainer.Dx}));
        mEndPositionsX.emplace_back(mModelCP.NewIntVar({minLength, mContainer.Dx}));

        // operations_research::Domain::FromValues(placementPointsPerType[item].Y
        mStartPositionsY.emplace_back(mModelCP.NewIntVar({0, mContainer.Dy}));
        mEndPositionsY.emplace_back(mModelCP.NewIntVar({minWidth, mContainer.Dy}));

        // operations_research::Domain::FromValues(placementPointsPerType[item].Z
        mStartPositionsZ.emplace_back(mModelCP.NewIntVar({0, mContainer.Dz}));
        mEndPositionsZ.emplace_back(mModelCP.NewIntVar({item.Dz, mContainer.Dz}));

        if(input.axleWights)
            mR.emplace_back(mModelCP.NewIntVar({-mContainer.Dx, mContainer.Dx}));
    }

    mLengths.reserve(numberOfItems);
    mWidths.reserve(numberOfItems);
    mHeights.reserve(numberOfItems);

    for(size_t i = 0; i < numberOfItems; i++)
    {
        const Cuboid &item = mItems[i];

        int minLength =
            item.EnableHorizontalRotation ? std::min(item.Dx, item.Dy) : item.Dx;
        int minWidth =
            item.EnableHorizontalRotation ? std::min(item.Dx, item.Dy) : item.Dy;

        int maxLength =
            item.EnableHorizontalRotation ? std::max(item.Dx, item.Dy) : item.Dx;
        int maxWidth =
            item.EnableHorizontalRotation ? std::max(item.Dx, item.Dy) : item.Dy;

        mLengths.emplace_back(
            mModelCP.NewIntVar(operations_research::Domain::FromIntervals(
                {{minLength, minLength}, {maxLength, maxLength}})));

        mWidths.emplace_back(
            mModelCP.NewIntVar(operations_research::Domain::FromIntervals(
                {{minWidth, minWidth}, {maxWidth, maxWidth}})));

        mHeights.emplace_back(mModelCP.NewIntVar({item.Dz, item.Dz}));
    }

    mPlacedOnFloor.reserve(numberOfItems);
    mPlacedOnLeft.reserve(numberOfItems);

    for(size_t i = 0; i < numberOfItems; i++)
    {
        mPlacedOnFloor.emplace_back(mModelCP.NewBoolVar());
        mPlacedOnLeft.emplace_back(mModelCP.NewBoolVar());
    }

    mIntervalsX.reserve(numberOfItems);
    mIntervalsY.reserve(numberOfItems);
    mIntervalsZ.reserve(numberOfItems);

    for(size_t i = 0; i < numberOfItems; i++)
    {
        mIntervalsX.emplace_back(
            mModelCP.NewIntervalVar(mStartPositionsX[i], mLengths[i], mEndPositionsX[i]));
        mIntervalsY.emplace_back(
            mModelCP.NewIntervalVar(mStartPositionsY[i], mWidths[i], mEndPositionsY[i]));
        mIntervalsZ.emplace_back(
            mModelCP.NewIntervalVar(mStartPositionsZ[i], mHeights[i], mEndPositionsZ[i]));
    }

    mOrientation.reserve(numberOfItems);
    for(size_t i = 0; i < numberOfItems; i++)
    {
        const Cuboid &itemI = mItems[i];
        mOrientation.emplace_back();
        mOrientation[i].reserve(mItemOrientations.size());

        for(size_t o = 0; o < mItemOrientations.size(); ++o)
        {
            mOrientation[i].emplace_back(mModelCP.NewBoolVar());
        }

        if(!itemI.EnableHorizontalRotation)
        {
            // TODO: consider not creating the variable instead of fixing it to zero.
            mModelCP.FixVariable(mOrientation[i][RotationZ], false);
        }
    }

    mRelativeDirections.reserve(numberOfItems);
    mSupportXY.reserve(numberOfItems);
    mLeftYZ.reserve(numberOfItems);

    for(size_t i = 0; i < numberOfItems; i++)
    {
        mRelativeDirections.emplace_back();
        mRelativeDirections[i].reserve(numberOfItems);

        if(mEnableFragility || mEnableSupport)
        {
            mSupportXY.emplace_back();
            mSupportXY[i].reserve(numberOfItems);
        }

        mLeftYZ.emplace_back();
        mLeftYZ[i].reserve(numberOfItems);

        for(size_t j = 0; j < numberOfItems; ++j)
        {
            mRelativeDirections[i].emplace_back();
            mRelativeDirections.reserve(mDimensions.size());
            for(size_t d = 0; d < mDimensions.size(); ++d)
            {
                mRelativeDirections[i][j].emplace_back(mModelCP.NewBoolVar());
                mRelativeDirections[i][j].emplace_back(mModelCP.NewBoolVar());
            }

            if(mEnableFragility || mEnableSupport)
            {
                mSupportXY[i].emplace_back(mModelCP.NewBoolVar());
            }

            mLeftYZ[i].emplace_back(mModelCP.NewBoolVar());
        }
    }

    mOverlapAreasXY.reserve(numberOfItems);
    mItemsOverlapsXY.reserve(numberOfItems);
    mItemsOverlapsYZ.reserve(numberOfItems);
    mOverlapAreasYZ.reserve(numberOfItems);
    mTopSum.reserve(numberOfItems);
    mTopBool.reserve(numberOfItems);

    for(int i=0; i < numberOfItems; ++i)
    {
        mTopSum.emplace_back(mModelCP.NewIntVar({0, (int)numberOfItems}));
        mTopBool.emplace_back(mModelCP.NewBoolVar());
    }


    for(size_t i = 0; i < numberOfItems - 1; i++)
    {
        const Cuboid &itemI = mItems[i];

        mItemsOverlapsXY.emplace_back();
        mItemsOverlapsXY.reserve(numberOfItems - i);

        mItemsOverlapsYZ.emplace_back();
        mItemsOverlapsYZ.reserve(numberOfItems - 1);



        if(mEnableSupport)
        {
            mOverlapAreasXY.emplace_back();
            mOverlapAreasXY.reserve(numberOfItems - i);
        }

        mOverlapAreasYZ.emplace_back();
        mOverlapAreasYZ.reserve(numberOfItems-i);

        for(size_t j = i + 1; j < numberOfItems; j++)
        {
            const Cuboid &itemJ = mItems[j];
            int maxIntersection = std::max(itemI.Dx * itemI.Dy, itemJ.Dx * itemJ.Dy);
            int maxIntersection2 = std::max(std::max(itemI.Dz*itemI.Dx, itemI.Dz*itemI.Dy),
                                            std::max(itemJ.Dz*itemJ.Dx, itemJ.Dz*itemJ.Dy));

            mItemsOverlapsXY[i].emplace_back(mModelCP.NewBoolVar());
            mItemsOverlapsYZ[i].emplace_back(mModelCP.NewBoolVar());

            if(mEnableSupport)
            {
                mOverlapAreasXY[i].emplace_back(mModelCP.NewIntVar({0, maxIntersection}));
            }

            mOverlapAreasYZ[i].emplace_back(mModelCP.NewIntVar({0, maxIntersection2}));
        }
    }

    if(mEnableLifo && !mFixedSequence)
    {
        mCustomerPosition.reserve(mNumberCustomers);
        for(size_t i = 0; i < mNumberCustomers; ++i)
        {
            mCustomerPosition.emplace_back(
                mModelCP.NewIntVar({1, static_cast<int>(mNumberCustomers)}));
        }

        mSuccessionMatrix.reserve(mNumberCustomers);
        for(size_t i = 0; i < mNumberCustomers - 1; i++)
        {
            mSuccessionMatrix.emplace_back();
            mSuccessionMatrix.reserve(mNumberCustomers - i);
            for(size_t j = i + 1; j < mNumberCustomers; j++)
            {
                mSuccessionMatrix[i].emplace_back(mModelCP.NewBoolVar());
            }
        }
    }

    mMaxLength = mModelCP.NewIntVar({0, mContainer.Dx});
    mMinX      = mModelCP.NewIntVar({0, mContainer.Dx});
}

void ContainerLoadingCP::CreateTopItem()
{

    for(int i=0; i < mItems.size(); ++i)
    {
        operations_research::sat::LinearExpr linExpr;

        for(int t=0; t < mItems.size(); ++t)
        {
            if(i == t)
                continue;

            linExpr += mSupportXY[t][i];
        }

        mModelCP.AddEquality(mTopSum[i], linExpr);
        mModelCP.AddGreaterOrEqual(mTopSum[i], 1).OnlyEnforceIf(mTopBool[i].Not());
        mModelCP.AddLessOrEqual(mTopSum[i], 0).OnlyEnforceIf(mTopBool[i]);
    }

}

std::tuple<ORIntVars1D, ORIntVars1D>
ContainerLoadingCP::GetIntVars(DimensionType dimension) const
{
    switch(dimension)
    {
    case AxisX:
        return std::make_tuple(mStartPositionsX, mEndPositionsX);
    case AxisY:
        return std::make_tuple(mStartPositionsY, mEndPositionsY);
    case AxisZ:
        return std::make_tuple(mStartPositionsZ, mEndPositionsZ);
    default:
        throw std::runtime_error("DimensionType not implemented.");
    }
}

void ContainerLoadingCP::AddConstraints()
{
    CreateNoOverlap();

    CreateItemOrientations();

    if(mEnableFragility || mEnableSupport)
    {
        CreateXYIntersectionBool();
        CreateSupportItem();
    }

    if(mEnableFragility)
    {
        CreateFragility();
    }

    CreateOnFloorConstraints();

    if(mEnableSupport)
    {
        CreateXYIntersectionArea();
        CreateSupportArea();
    }

    if(mEnableLifo)
    {
        if(mFixedSequence)
        {
            CreateLifoSequence();
        }
        else
        {
            CreatePositioningConstraints();
            CreateLifoNoSequence();
        }
    }

    if(ParseInputNS::input.axleWights)
        CreateAxleWeights();

    CreateBalancedLoading();


    CreateTopItem();


    CreateCompactnessItem();
    CreateOnLeftConstraints();
    CreateCompactnessArea();
    CreateYZIntersectionBool();
    CreateYZIntersectionArea();

}

void ContainerLoadingCP::CreateAxleWeights()
{
    int tolerance = 1;
    //GravityMM = 1;

    operations_research::sat::LinearExpr sumMoments;
    int                                  sumForces = 0;
    for(int i = 0; i < mItems.size(); ++i)
    {
        // x = y
        // x - y = 0
        // x - y >= -1
        // x - y <= 1
        // mModelCP.AddEquality(2*mR[i], 2*semiTrailer.distanceCargoSpaceTrailerAxle
        // -2*mStartPositionsX[i] - mLengths[i]);
        mModelCP
            .AddGreaterOrEqual(2 * mR[i],
                               2 * semiTrailer.distanceCargoSpaceTrailerAxle -
                                   2 * mStartPositionsX[i] - mLengths[i])
            .WithName("R0");

        mModelCP
            .AddLessThan(2 * mR[i],
                         2 * semiTrailer.distanceCargoSpaceTrailerAxle -
                             2 * mStartPositionsX[i] - mLengths[i] + 2)
            .WithName("R1");

        int itemF = mItems[i].Weight * GravityMM;
        sumMoments += itemF * mR[i];
        sumForces += itemF;

        /*
        std::cout << "i=" << i
                  << " R in [" << mR[i].Domain() << "]"
                  << " X in [" << mStartPositionsX[i].Domain() << "]"
                  << " W in [" << mWidths[i].Domain()<<"]"
                  << " D=" << semiTrailer.distanceCargoSpaceTrailerAxle
                  << std::endl;
        */
    }

    // EQ: 10
    mModelCP
        .AddGreaterOrEqual(
            semiTrailer.distanceKingpinTrailerAxle * forceK,
            scale * sumMoments + scale * semiTrailer.massTrailer * GravityMM *
                                     semiTrailer.distanceMassTrailerTrailerAxle)
        .WithName("EQ10_0");
    mModelCP
        .AddLessThan(semiTrailer.distanceKingpinTrailerAxle * forceK,
                     scale * sumMoments +
                         scale * semiTrailer.massTrailer * GravityMM *
                             semiTrailer.distanceMassTrailerTrailerAxle +
                         semiTrailer.distanceKingpinTrailerAxle)
        .WithName("EQ10_1");

    // EQ: 11
    // mModelCP.AddGreaterOrEqual(forceFA + forceRA - forceK -
    // scale*semiTrailer.massTractor*GravityMM, -tolerance).WithName("EQ11_0");
    // mModelCP.AddLessOrEqual(forceFA + forceRA - forceK -
    // scale*semiTrailer.massTractor*GravityMM, tolerance).WithName("EQ11_1");

    // EQ: 12
    mModelCP
        .AddGreaterOrEqual(semiTrailer.wheelBase * forceFA,
                           semiTrailer.distanceKingpinRearAxle * forceK +
                               scale * semiTrailer.massTractor * GravityMM *
                                   semiTrailer.distanceMassTractorRearAxle)
        .WithName("EQ12_0");
    mModelCP
        .AddLessThan(semiTrailer.wheelBase * forceFA,
                     semiTrailer.distanceKingpinRearAxle * forceK +
                         scale * semiTrailer.massTractor * GravityMM *
                             semiTrailer.distanceMassTractorRearAxle +
                         semiTrailer.wheelBase)
        .WithName("EQ12_1");

    // EQ: 9
    // mModelCP.AddGreaterOrEqual(forceTA, scale*sumForces +
    // scale*semiTrailer.massTrailer*GravityMM - forceK - tolerance)
    //                           .WithName("EQ9_0");
    // mModelCP.AddLessOrEqual(forceTA, scale*sumForces +
    // scale*semiTrailer.massTrailer*GravityMM - forceK + tolerance)
    //                         .WithName("EQ9_1");

    // EQ: 10
    // mModelCP.AddEquality(semiTrailer.distanceKingpinTrailerAxle*forceK,
    // scale*sumMoments
    // +
    // scale*semiTrailer.massTrailer*GravityMM*semiTrailer.distanceMassTrailerTrailerAxle);
    // EQ: 11 - can use this as there are no divisions
    mModelCP.AddEquality(
        forceFA + forceRA - forceK - scale * semiTrailer.massTractor * GravityMM, 0);
    // EQ: 12
    // mModelCP.AddEquality(semiTrailer.wheelBase*forceFA,
    // semiTrailer.distanceKingpinRearAxle*forceK +
    // scale*semiTrailer.massTractor*GravityMM*semiTrailer.distanceMassTractorRearAxle);
    // EQ: 9 - can use this as there are no divisions
    mModelCP.AddEquality(
        forceTA,
        scale * sumForces + scale * semiTrailer.massTrailer * GravityMM - forceK);
}

void ContainerLoadingCP::CreateBalancedLoading()
{

    operations_research::sat::LinearExpr exp0;
    operations_research::sat::LinearExpr exp1;

    const int w = InstanceNS::instanciaG.vetDimVeiculo[1];
    const int wDiv2 = InstanceNS::instanciaG.vetDimVeiculo[1] / 2;

    for(int i = 0; i < mItems.size(); ++i)
    {
        {
            auto tempLeftf0 = mModelCP.NewIntVar({0, 10 * w});
            mModelCP.AddMaxEquality(tempLeftf0, {wDiv2 - mStartPositionsY[i], 0});

            auto tempLeftf1 = mModelCP.NewIntVar({0, 10 * w});
            mModelCP.AddMaxEquality(
                tempLeftf1, {wDiv2 - (mStartPositionsY[i] + mWidths[i]), 0});

            const int m = mItems[i].Weight;
            int       ub = m * scaleBalancedLoading * 10 * w;

            auto resLeft = mModelCP.NewIntVar({-ub, ub});

            mModelCP.AddDivisionEquality(
                resLeft,
                scaleBalancedLoading * m * (tempLeftf0 - tempLeftf1),
                mWidths[i]);
            exp0 += resLeft;
        }

        {
            auto tempRight0 = mModelCP.NewIntVar({0, 10 * w});
            mModelCP.AddMaxEquality(
                tempRight0, {((mStartPositionsY[i] + mWidths[1]) - wDiv2), 0});

            auto tempRight1 = mModelCP.NewIntVar({0, 10 * w});
            mModelCP.AddMaxEquality(tempRight1, {(mStartPositionsY[i] - wDiv2), 0});

            const int m = mItems[i].Weight;
            int       ub = m * scaleBalancedLoading * 10 * w;

            auto resRight = mModelCP.NewIntVar({-ub, ub});

            mModelCP.AddDivisionEquality(
                resRight,
                scaleBalancedLoading * m * (tempRight0 - tempRight1),
                mWidths[i]);
            exp1 += resRight;
        }
    }

    mModelCP.AddLessOrEqual(
        exp0,
        (int)(scaleBalancedLoading *
              (input.balancedLoadingD * InstanceNS::instanciaG.maxPayload - 10)));
    mModelCP.AddLessOrEqual(
        exp1,
        (int)(scaleBalancedLoading *
              (input.balancedLoadingD * InstanceNS::instanciaG.maxPayload - 10)));

    mModelCP.AddEquality(sumRightBalancedLoading, exp1);
    mModelCP.AddEquality(sumLeftBalancedLoading, exp0);
}

/// Relative directions of items. Necessary for non overlapping items.
void ContainerLoadingCP::CreateNoOverlap()
{
    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems; ++i)
    {
        for(size_t j = 0; j < numberOfItems; ++j)
        {
            if(i == j)
            {
                continue;
            }

            for(size_t d = 0; d < mDimensions.size(); ++d)
            {
                const Dimension &dimension = mDimensions[d];
                const auto [startPosition, endPosition] = GetIntVars(dimension.Type);

                mModelCP.AddLessOrEqual(endPosition[j], startPosition[i])
                    .OnlyEnforceIf(mRelativeDirections[i][j][dimension.FirstDirection]);
                mModelCP.AddLessThan(startPosition[i], endPosition[j])
                    .OnlyEnforceIf(
                        mRelativeDirections[i][j][dimension.FirstDirection].Not());

                mModelCP.AddLessOrEqual(endPosition[i], startPosition[j])
                    .OnlyEnforceIf(mRelativeDirections[i][j][dimension.SecondDirection]);
                mModelCP.AddLessThan(startPosition[j], endPosition[i])
                    .OnlyEnforceIf(
                        mRelativeDirections[i][j][dimension.SecondDirection].Not());

                mModelCP.AddEquality(
                    mRelativeDirections[i][j][dimension.FirstDirection],
                    mRelativeDirections[j][i][dimension.SecondDirection]);
            }

            // No overlap constraints
            mModelCP.AddAtLeastOne(mRelativeDirections[i][j]);
        }
    }

    /*
    for (size_t i = 0; i < numberOfItems; ++i)
    {
        for (size_t j = i+1; j < numberOfItems; ++j)
        {
            std::vector<operations_research::sat::BoolVar> dirs;

            for (size_t d = 0; d < mDimensions.size(); ++d)
            {
                const Dimension& dimension = mDimensions[d];
                const auto [startPosition, endPosition] = GetIntVars(dimension.Type);

                mModelCP.AddLessOrEqual(endPosition[j], startPosition[i])
                    .OnlyEnforceIf(mRelativeDirections[i][j][dimension.FirstDirection]);

                mModelCP.AddLessOrEqual(endPosition[i], startPosition[j])
                    .OnlyEnforceIf(mRelativeDirections[i][j][dimension.SecondDirection]);

                dirs.push_back(mRelativeDirections[i][j][dimension.FirstDirection]);
                dirs.push_back(mRelativeDirections[i][j][dimension.SecondDirection]);

                mModelCP.AddEquality(
                    mRelativeDirections[i][j][dimension.FirstDirection],
                    mRelativeDirections[j][i][dimension.SecondDirection]);
            }

            mModelCP.AddBoolOr(dirs);
        }
    }
    */
}

/// Set dimensions of items based on orientation
void ContainerLoadingCP::CreateItemOrientations()
{
    size_t numberOfItems = mItems.size();
    for(size_t i = 0; i < numberOfItems; ++i)
    {
        const Cuboid &item = mItems[i];
        for(size_t o = 0; o < mItemOrientations.size(); ++o)
        {
            // Dimensions of item depending on orientation
            auto [itemLength, itemWidth, itemHeight] =
                item.DetermineDimensions(mItemOrientations[o]);

            mModelCP.AddEquality(mLengths[i], itemLength)
                .OnlyEnforceIf(mOrientation[i][o]);
            mModelCP.AddEquality(mWidths[i], itemWidth).OnlyEnforceIf(mOrientation[i][o]);
            mModelCP.AddEquality(mHeights[i], itemHeight)
                .OnlyEnforceIf(mOrientation[i][o]);
        }

        mModelCP.AddExactlyOne(mOrientation[i]);
    }
}

/// Fragility of items
void ContainerLoadingCP::CreateFragility()
{
    size_t numberOfItems = mItems.size();

    // Variant 2 - fragile items can be stacked on fragile items.
    // Fragile items can be stacked onto other fragile or non-fragile items, whereas
    // non-fragile items must not touch fragile items from above.
    for(size_t i = 0; i < numberOfItems; ++i)
    {
        for(size_t j = 0; j < numberOfItems; ++j)
        {
            if(mItems[j].Fragility == Fragility::Fragile &&
               mItems[i].Fragility == Fragility::None)
            {
                if(i != j)
                {
                    // Item j cannot support item i, if j is fragile and i non fragile.
                    mModelCP.FixVariable(mSupportXY[i][j], false);
                }
            }
        }
    }
}

/// Determines supported area of items
void ContainerLoadingCP::CreateSupportArea()
{
    for(size_t i = 0; i < mItems.size(); ++i)
    {
        ORLinExpr supportedAreaExpr;
        int       areaI = mItems[i].Dy * mItems[i].Dx;
        for(size_t j = 0; j < mItems.size(); ++j)
        {
            if(i == j)
            {
                continue;
            }

            if(!mEnableFragility || mItems[j].Fragility == Fragility::None ||
               (mItems[j].Fragility == Fragility::Fragile &&
                mItems[i].Fragility == Fragility::Fragile))
            {
                int                              areaJ = mItems[j].Dy * mItems[j].Dx;
                int                              minArea = std::min(areaI, areaJ);
                operations_research::sat::IntVar usableArea =
                    mModelCP.NewIntVar({0, minArea});
                if(i < j)
                {
                    auto position = j - i - 1;
                    mModelCP.AddMultiplicationEquality(
                        usableArea, {mOverlapAreasXY[i][position], mSupportXY[i][j]});
                }
                else
                {
                    auto position = i - j - 1;
                    mModelCP.AddMultiplicationEquality(
                        usableArea, {mOverlapAreasXY[j][position], mSupportXY[i][j]});
                }

                mModelCP.AddEquality(usableArea, 0).OnlyEnforceIf(mSupportXY[i][j].Not());
                supportedAreaExpr += usableArea;
            }
        }

        operations_research::sat::IntVar supportedArea = mModelCP.NewIntVar({0, areaI});
        mModelCP.AddEquality(supportedArea, supportedAreaExpr)
            .OnlyEnforceIf(mPlacedOnFloor[i].Not());

        mModelCP
            .AddGreaterOrEqual(
                supportedArea,
                static_cast<int>(std::ceil(mSupportArea * mItems[i].Dx * mItems[i].Dy)))
            .OnlyEnforceIf(mPlacedOnFloor[i].Not());
    }
}

/// Set bool variable, if items i and j intersect in XY
void ContainerLoadingCP::CreateXYIntersectionBool()
{
    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems - 1; ++i)
    {
        for(size_t j = i + 1; j < numberOfItems; ++j)
        {
            auto positionJ = j - i - 1;

            mModelCP.AddAtLeastOne({mItemsOverlapsXY[i][positionJ],
                                    mRelativeDirections[i][j][LeftY],
                                    mRelativeDirections[i][j][RightY],
                                    mRelativeDirections[i][j][BehindX],
                                    mRelativeDirections[i][j][InFrontX]});

            mModelCP.AddImplication(
                mRelativeDirections[i][j][LeftY], mItemsOverlapsXY[i][positionJ].Not());
            mModelCP.AddImplication(
                mRelativeDirections[i][j][RightY], mItemsOverlapsXY[i][positionJ].Not());
            mModelCP.AddImplication(
                mRelativeDirections[i][j][BehindX], mItemsOverlapsXY[i][positionJ].Not());
            mModelCP.AddImplication(mRelativeDirections[i][j][InFrontX],
                                    mItemsOverlapsXY[i][positionJ].Not());
        }
    }
}

/// Determine intersection area of all items
void ContainerLoadingCP::CreateXYIntersectionArea()
{
    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems - 1; ++i)
    {
        for(size_t j = i + 1; j < numberOfItems; ++j)
        {
            if(mItems[i].Dz + mItems[j].Dz <= mContainer.Dz)
            {
                // Variant 2
                auto positionJ = j - i - 1;

                // Overlap in x
                operations_research::sat::IntVar diffXij =
                    mModelCP.NewIntVar({0, mContainer.Dx});
                mModelCP
                    .AddEquality(diffXij,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsX[i], mStartPositionsX[j]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsXY[i][positionJ]});

                operations_research::sat::IntVar diffXji =
                    mModelCP.NewIntVar({0, mContainer.Dx});
                mModelCP
                    .AddEquality(diffXji,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsX[j], mStartPositionsX[i]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsXY[i][positionJ]});

                operations_research::sat::IntVar xOverlap =
                    mModelCP.NewIntVar({0, mContainer.Dx});
                mModelCP.AddMinEquality(
                    xOverlap, {diffXij, diffXji, mLengths[i], mLengths[j]});
                mModelCP.AddEquality(xOverlap, 0)
                    .OnlyEnforceIf(mItemsOverlapsXY[i][positionJ].Not());

                // Overlap in y
                operations_research::sat::IntVar diffYij =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP
                    .AddEquality(diffYij,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsY[i], mStartPositionsY[j]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsXY[i][positionJ]});

                operations_research::sat::IntVar diffYji =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP
                    .AddEquality(diffYji,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsY[j], mStartPositionsY[i]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsXY[i][positionJ]});

                operations_research::sat::IntVar yOverlap =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP.AddMinEquality(
                    yOverlap, {diffYij, diffYji, mWidths[i], mWidths[j]});
                mModelCP.AddEquality(yOverlap, 0)
                    .OnlyEnforceIf(mItemsOverlapsXY[i][positionJ].Not());

                // Area
                mModelCP.AddMultiplicationEquality(
                    mOverlapAreasXY[i][positionJ], {xOverlap, yOverlap});
            }
        }
    }
}

/// Set bool variable, if item j supports item i -> adjacent in Z AND items intersect in
/// XY
void ContainerLoadingCP::CreateSupportItem()
{
    for(size_t i = 0; i < mItems.size(); ++i)
    {
        mModelCP.FixVariable(mSupportXY[i][i], false);
        for(size_t j = 0; j < mItems.size(); ++j)
        {
            if(i == j)
            {
                continue;
            }

            operations_research::sat::BoolVar isVerticallyAdjacent =
                mModelCP.NewBoolVar();
            mModelCP.AddEquality(mEndPositionsZ[j], mStartPositionsZ[i])
                .OnlyEnforceIf(isVerticallyAdjacent);
            mModelCP.AddNotEqual(mEndPositionsZ[j], mStartPositionsZ[i])
                .OnlyEnforceIf(isVerticallyAdjacent.Not());

            mModelCP.AddImplication(isVerticallyAdjacent.Not(), mSupportXY[i][j].Not());

            if(i < j)
            {
                auto position = j - i - 1;

                mModelCP.AddAtLeastOne({mSupportXY[i][j],
                                        isVerticallyAdjacent.Not(),
                                        mItemsOverlapsXY[i][position].Not()});
                mModelCP.AddImplication(
                    mItemsOverlapsXY[i][position].Not(), mSupportXY[i][j].Not());
            }
            else
            {
                auto position = i - j - 1;

                mModelCP.AddAtLeastOne({mSupportXY[i][j],
                                        isVerticallyAdjacent.Not(),
                                        mItemsOverlapsXY[j][position].Not()});
                mModelCP.AddImplication(
                    mItemsOverlapsXY[j][position].Not(), mSupportXY[i][j].Not());
            }
        }
    }
}

/// LIFO unloading with given customer sequence
void ContainerLoadingCP::CreateLifoSequence()
{
    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems; ++i)
    {
        for(size_t j = 0; j < numberOfItems; ++j)
        {
            if(i != j && mItems[i].pos < mItems[j].pos)
            {

                // i is delevered first

                // std::printf("%d %d\n", (int)mItems[i].ExternId,
                // (int)mItems[j].ExternId);
                //  Item i must be placed behind or below item j if
                //  - item i is unloaded after item j (smaller group id) and
                //  - item i is not placed left or right of item j -> in way to rear end
                //  of container

                // Original Constraint
                // mModelCP.AddAtLeastOne({mRelativeDirections[i][j][BehindX],
                // mRelativeDirections[i][j][BelowZ]})
                //    .OnlyEnforceIf({mRelativeDirections[i][j][LeftY].Not(),
                //    mRelativeDirections[i][j][RightY].Not()});
                // Funciona

                operations_research::sat::LinearExpr linExp;

                /*
                if(input.mlifo)
                {
                    mModelCP.AddAtLeastOne({mRelativeDirections[i][j][InFrontX],
                mRelativeDirections[i][j][BehindX], mRelativeDirections[i][j][RightY],
                mRelativeDirections[i][j][AboveZ], mRelativeDirections[i][j][LeftY]});
                }
                else*/
                if(input.mlifo) // input.removeFromShortSide)
                {
                    operations_research::sat::BoolVar varStarEndEq =
                        mModelCP.NewBoolVar();
                    operations_research::sat::BoolVar var = mModelCP.NewBoolVar();

                    mModelCP.AddEquality(mStartPositionsZ[j], mEndPositionsZ[i])
                        .OnlyEnforceIf(varStarEndEq.Not());
                    mModelCP.AddNotEqual(mStartPositionsZ[j], mEndPositionsZ[i])
                        .OnlyEnforceIf(varStarEndEq);

                    mModelCP.AddEquality(var, 1).OnlyEnforceIf(
                        {varStarEndEq.Not(), mRelativeDirections[i][j][BelowZ]});

                    // var = mRelativeDirections[i][j][BelowZ] ^ (mStartPositionZ[j] !=
                    // mEndPositionZ[i])
                    // mModelCP.

                    mModelCP.AddAtLeastOne({mRelativeDirections[i][j][RightY],
                                            mRelativeDirections[i][j][LeftY],
                                            mRelativeDirections[i][j][InFrontX],
                                            mRelativeDirections[i][j][AboveZ],
                                            var});
                }
                else
                {
                    mModelCP.AddAtLeastOne({mRelativeDirections[i][j][InFrontX],
                                            mRelativeDirections[i][j][BehindX],
                                            mRelativeDirections[i][j][RightY],
                                            mRelativeDirections[i][j][AboveZ]});
                }
            }
        }
    }

    // std::printf("\n\n");
}

/// LIFO unloading without given customer sequence
void ContainerLoadingCP::CreateLifoNoSequence()
{
    for(size_t i = 0; i < mItems.size(); ++i)
    {
        auto customerI = mItems[i].GroupId;
        for(size_t j = 0; j < mItems.size(); ++j)
        {
            auto customerJ = mItems[j].GroupId;

            if(i != j && customerI != customerJ)
            {
                // Item i must be placed behind or below item j if
                // - item i is unloaded after item j (customer i succeeds customer j) and
                // - item i is not placed left or right of item j -> in way to rear end of
                // container.

                if(customerI < customerJ)
                {
                    auto position = customerJ - customerI - 1;
                    mModelCP
                        .AddAtLeastOne({mRelativeDirections[i][j][BehindX],
                                        mRelativeDirections[i][j][BelowZ]})
                        .OnlyEnforceIf({mRelativeDirections[i][j][LeftY].Not(),
                                        mRelativeDirections[i][j][RightY].Not(),
                                        mSuccessionMatrix[customerI][position]});
                }
                else
                {
                    auto position = customerI - customerJ - 1;
                    mModelCP
                        .AddAtLeastOne({mRelativeDirections[i][j][BehindX],
                                        mRelativeDirections[i][j][BelowZ]})
                        .OnlyEnforceIf({mRelativeDirections[i][j][LeftY].Not(),
                                        mRelativeDirections[i][j][RightY].Not(),
                                        mSuccessionMatrix[customerJ][position].Not()});
                }
            }
        }
    }
}

/// Positioning of customers in route if customer sequence is not given; needed for LIFO
void ContainerLoadingCP::CreatePositioningConstraints()
{
    mModelCP.AddAllDifferent(mCustomerPosition);

    for(size_t i = 0; i < mNumberCustomers - 1; i++)
    {
        for(size_t j = i + 1; j < mNumberCustomers; j++)
        {
            auto positionInVector = j - i - 1;

            // If customer i is visited after customer j => position of i is greater than
            // position of j Sequence of customers in ascending order
            mModelCP.AddGreaterThan(mCustomerPosition[i], mCustomerPosition[j])
                .OnlyEnforceIf(mSuccessionMatrix[i][positionInVector]);
            mModelCP.AddLessThan(mCustomerPosition[i], mCustomerPosition[j])
                .OnlyEnforceIf(mSuccessionMatrix[i][positionInVector].Not());
        }
    }
}

/// Determine which items are placed on the floor
void ContainerLoadingCP::CreateOnFloorConstraints()
{
    for(size_t i = 0; i < mItems.size(); ++i)
    {
        mModelCP.AddEquality(mStartPositionsZ[i], 0).OnlyEnforceIf(mPlacedOnFloor[i]);
        mModelCP.AddGreaterThan(mStartPositionsZ[i], 0)
            .OnlyEnforceIf(mPlacedOnFloor[i].Not());
    }
}

void ContainerLoadingCP::AddObjective()
{
    mModelCP.AddMaxEquality(mMaxLength, mEndPositionsX);

    mModelCP.Minimize(mMaxLength);
}

void ContainerLoadingCP::CreateCompactnessItem()
{
    mModelCP.AddMinEquality(mMinX, mStartPositionsX);

    for(size_t i = 0; i < mItems.size(); ++i)
    {
        mModelCP.FixVariable(mLeftYZ[i][i], false);
        for(size_t j = 0; j < mItems.size(); ++j)
        {
            if(i == j)
            {
                continue;
            }

            operations_research::sat::BoolVar isAdjacent = mModelCP.NewBoolVar();
            mModelCP.AddEquality(mEndPositionsX[j], mStartPositionsX[i])
                .OnlyEnforceIf(isAdjacent);
            mModelCP.AddNotEqual(mEndPositionsX[j], mStartPositionsX[i])
                .OnlyEnforceIf(isAdjacent.Not());

            mModelCP.AddImplication(isAdjacent.Not(), mLeftYZ[i][j].Not());

            if(i < j)
            {
                auto position = j - i - 1;

                mModelCP.AddAtLeastOne({mLeftYZ[i][j],
                                        isAdjacent.Not(),
                                        mItemsOverlapsYZ[i][position].Not()});
                mModelCP.AddImplication(
                    mItemsOverlapsYZ[i][position].Not(), mLeftYZ[i][j].Not());
            }
            else
            {
                auto position = i - j - 1;

                mModelCP.AddAtLeastOne({mLeftYZ[i][j],
                                        isAdjacent.Not(),
                                        mItemsOverlapsYZ[j][position].Not()});
                mModelCP.AddImplication(
                    mItemsOverlapsYZ[j][position].Not(), mLeftYZ[i][j].Not());
            }
        }
    }

}

void ContainerLoadingCP::CreateCompactnessArea()
{
    for(size_t i = 0; i < mItems.size(); ++i)
    {
        ORLinExpr supportedAreaExpr;
        //int       areaI = mItems[i].Dy * mItems[i].Dx;
        //int minAreaI = mItems[i].Dz * std::min(mItems[i].Dx, mItems[i].Dy);
        int maxAreaI = mItems[i].Dz * std::max(mItems[i].Dx, mItems[i].Dy);

        //operations_research::sat::IntVar areaI = mModelCP.NewIntVar({minAreaI, maxAreaI});
        //mModelCP.AddMultiplicationEquality(areaI, {mHeights[i], mLengths[i]});

        for(size_t j = 0; j < mItems.size(); ++j)
        {
            if(i == j)
            {
                continue;
            }

                   //if(!mEnableFragility || mItems[j].Fragility == Fragility::None ||
                   //   (mItems[j].Fragility == Fragility::Fragile &&
                   //    mItems[i].Fragility == Fragility::Fragile))
            {
                //int minAreaJ = mItems[j].Dz * std::min(mItems[j].Dx, mItems[j].Dy);
                int maxAreaJ = mItems[j].Dz * std::max(mItems[j].Dx, mItems[j].Dy);


               // mModelCP.AddMultiplicationEquality(areaJ, {mHeights[j], mLengths[j]});

                operations_research::sat::IntVar usableArea =
                    mModelCP.NewIntVar({0, std::max(maxAreaI, maxAreaJ)});
                if(i < j)
                {
                    auto position = j - i - 1;
                    mModelCP.AddMultiplicationEquality(
                        usableArea, {mOverlapAreasYZ[i][position], mLeftYZ[i][j]});
                }
                else
                {
                    auto position = i - j - 1;
                    mModelCP.AddMultiplicationEquality(
                        usableArea, {mOverlapAreasYZ[j][position], mLeftYZ[i][j]});
                }

                mModelCP.AddEquality(usableArea, 0).OnlyEnforceIf(mLeftYZ[i][j].Not());
                supportedAreaExpr += usableArea;
            }
        }

        operations_research::sat::IntVar supportedArea = mModelCP.NewIntVar({0, maxAreaI});
        mModelCP.AddEquality(supportedArea, supportedAreaExpr)
            .OnlyEnforceIf(mPlacedOnLeft[i].Not());


        mModelCP
            .AddGreaterOrEqual(
                supportedArea,
                static_cast<int>(std::ceil(msupportAreaLeft * mItems[i].Dy * mItems[i].Dz)))
            .OnlyEnforceIf({mPlacedOnLeft[i].Not(), mOrientation[i][NoRotation],
                            mTopBool[i].Not()});

        mModelCP
            .AddGreaterOrEqual(
                supportedArea,
                static_cast<int>(std::ceil(msupportAreaLeft * mItems[i].Dx * mItems[i].Dz)))
            .OnlyEnforceIf({mPlacedOnLeft[i].Not(), mOrientation[i][RotationZ],
                            mTopBool[i].Not()});

    }
}

void ContainerLoadingCP::CreateYZIntersectionBool()
{

    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems - 1; ++i)
    {
        for(size_t j = i + 1; j < numberOfItems; ++j)
        {
            auto positionJ = j - i - 1;

            mModelCP.AddAtLeastOne({mItemsOverlapsYZ[i][positionJ],
                                    mRelativeDirections[i][j][LeftY],
                                    mRelativeDirections[i][j][RightY],
                                    mRelativeDirections[i][j][BelowZ],
                                    mRelativeDirections[i][j][AboveZ]});

            mModelCP.AddImplication(
                mRelativeDirections[i][j][LeftY], mItemsOverlapsYZ[i][positionJ].Not());
            mModelCP.AddImplication(
                mRelativeDirections[i][j][RightY], mItemsOverlapsYZ[i][positionJ].Not());

            mModelCP.AddImplication(
                mRelativeDirections[i][j][BelowZ], mItemsOverlapsYZ[i][positionJ].Not());
            mModelCP.AddImplication(
                mRelativeDirections[i][j][AboveZ], mItemsOverlapsYZ[i][positionJ].Not());
        }
    }

}

void ContainerLoadingCP::CreateYZIntersectionArea()
{

    size_t numberOfItems = mItems.size();

    for(size_t i = 0; i < numberOfItems - 1; ++i)
    {
        for(size_t j = i + 1; j < numberOfItems; ++j)
        {
            //if(mItems[i].Dz + mItems[j].Dz <= mContainer.Dz)
            {
                // Variant 2
                auto positionJ = j - i - 1;

                       // Overlap in x
                operations_research::sat::IntVar diffZij =
                    mModelCP.NewIntVar({0, mContainer.Dz});
                mModelCP
                    .AddEquality(diffZij,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsZ[i], mStartPositionsZ[j]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsYZ[i][positionJ]});

                operations_research::sat::IntVar diffZji =
                    mModelCP.NewIntVar({0, mContainer.Dz});
                mModelCP
                    .AddEquality(diffZji,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsZ[j], mStartPositionsZ[i]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsYZ[i][positionJ]});

                operations_research::sat::IntVar zOverlap =
                    mModelCP.NewIntVar({0, mContainer.Dz});
                mModelCP.AddMinEquality(
                    zOverlap, {diffZij, diffZji, mHeights[i], mHeights[j]});
                mModelCP.AddEquality(zOverlap, 0)
                    .OnlyEnforceIf(mItemsOverlapsYZ[i][positionJ].Not());

                       // Overlap in y
                operations_research::sat::IntVar diffYij =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP
                    .AddEquality(diffYij,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsY[i], mStartPositionsY[j]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsYZ[i][positionJ]});

                operations_research::sat::IntVar diffYji =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP
                    .AddEquality(diffYji,
                                 ORLinExpr::LinearExpr::WeightedSum(
                                     {mEndPositionsY[j], mStartPositionsY[i]}, {1, -1}))
                    .OnlyEnforceIf({mItemsOverlapsYZ[i][positionJ]});

                operations_research::sat::IntVar yOverlap =
                    mModelCP.NewIntVar({0, mContainer.Dy});
                mModelCP.AddMinEquality(
                    yOverlap, {diffYij, diffYji, mWidths[i], mWidths[j]});
                mModelCP.AddEquality(yOverlap, 0)
                    .OnlyEnforceIf(mItemsOverlapsYZ[i][positionJ].Not());

                       // Area
                mModelCP.AddMultiplicationEquality(
                    mOverlapAreasYZ[i][positionJ], {zOverlap, yOverlap});
            }
        }
    }

}

void ContainerLoadingCP::CreateOnLeftConstraints()
{

    for(size_t i = 0; i < mItems.size(); ++i)
    {
        mModelCP.AddEquality(mStartPositionsX[i], mMinX).OnlyEnforceIf(mPlacedOnLeft[i]);
        mModelCP.AddGreaterThan(mStartPositionsX[i], mMinX)
            .OnlyEnforceIf(mPlacedOnLeft[i].Not());
    }

}

} // namespace Algorithms
} // namespace ContainerLoading
