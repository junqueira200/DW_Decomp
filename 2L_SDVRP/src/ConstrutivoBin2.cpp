/* ****************************************
 * ****************************************
 *  Data:    14/09/26
 *  Arquivo: ConstrutivoBin.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L/3L-SDVRP
 * ****************************************
 * ****************************************/

#include "ConstrutivoBin2.h"
#include "AxleWeights.h"
#include "ConstrutivoBin.h"
#include "AuxGeometry.h"
#include "AuxT.h"
#include "InputOutput.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include "c_api.h"
#include "Alarm.h"
#include <set>

using namespace InstanceNS;
using namespace SolucaoNS;
using namespace AuxGeometryNS;
using namespace ConstrutivoBinNS;

namespace ConstrutivoBin2NS
{

void getRouteItems(const SolucaoNS::Rota &rota, VectorI &outItemIds)
{
    outItemIds.clear();

           // Iterate through all customer stops (excluding depot at 0 and numPos - 1)
    for(int pos = 1; pos < rota.numPos - 1; ++pos)
    {
        int customer = rota.vetRota[pos];
        if(customer <= 0 || customer >= instanciaG.numClientes)
            continue;

               // Collect items from instance customer range
        int firstItem = instanciaG.matCliItensIniFim(customer, 0);
        int lastItem  = instanciaG.matCliItensIniFim(customer, 1);

        for(int itemId = firstItem; itemId <= lastItem; ++itemId)
        {
            if(itemId >= 0 && itemId < instanciaG.numItens)
                outItemIds.push_back(itemId);
        }
    }
}


/**
 * @brief Evaluates whether placing itemId at point p with rotation r is locally feasible.
 */
static bool isPlacementFeasible(const Bin                 &bin,
                                Rota                      &rota,
                                int                        itemId,
                                const Ponto               &p,
                                InstanceNS::Rotation       r,
                                Eigen::Matrix<int, -1, -1, Eigen::RowMajor> &matSupportItems,
                                double 					   sumM,
                                double					   sumF)
{
    Item &item = instanciaG.vetItens[itemId];

           // 1. Boundary check
    for(int d = 0; d < instanciaG.numDim; ++d)
    {
        double dimR = item.getDimRotacionada(d, r);
        if(p.vetDim[d] + dimR > bin.binDim[d] + 1e-5)
            return false;
    }

    /*
    double f = instanciaG.vetItens[itemId].weight * GravityCm;
    sumF += f;

    double r_ = (double)AxleWeightsNS::semiTrailer.distanceCargoSpaceTrailerAxle -
               p.vetDim[0] -
               instanciaG.vetItens[itemId].getDimRotacionada(0, r) / 2.0;

    sumM += f * r_;

    double fK = (1.0 / (double)AxleWeightsNS::semiTrailer.distanceKingpinTrailerAxle) *
         (sumM + (double)AxleWeightsNS::semiTrailer.massTrailer * GravityCm *
                     AxleWeightsNS::semiTrailer.distanceMassTrailerTrailerAxle);

    double fFA = (1.0 / (double)AxleWeightsNS::semiTrailer.wheelBase) *
          (fK * (double)AxleWeightsNS::semiTrailer.distanceKingpinRearAxle +
           (double)AxleWeightsNS::semiTrailer.massTractor *
                      GravityCm * AxleWeightsNS::semiTrailer.distanceMassTractorRearAxle);

    double fRA = fK + (double)AxleWeightsNS::semiTrailer.massTractor * GravityCm - fFA;
    double fTA = sumF + (double)AxleWeightsNS::semiTrailer.massTrailer * GravityCm - fK;

    if(fFA > (double)AxleWeightsNS::semiTrailer.maxMassFrontAxle* GravityCm*1.1)
        return false;


    if(fRA > (double)AxleWeightsNS::semiTrailer.maxMassRearAxle * GravityCm*1.1)
        return false;


    if(fTA > (double)AxleWeightsNS::semiTrailer.maxMassTrailerAxle * GravityCm*1.1) // ||
        return false;
    */

    // 2. Collision check with all already-packed items
    for(int k = 0; k < bin.numItens; ++k)
    {
        if(verificaColisaoDoisItens(itemId, bin.vetItemId[k], p, bin.vetPosItem[k], r,
                                    bin.vetRotacao[k]))
            return false;
    }

    // 3. Bottom Support & Fragility check (if 3D and item is not resting on the floor)
    if(instanciaG.numDim == 3 && ParseInputNS::input.support && p.vetDim[2] > 1e-4)
    {
        double supportedArea = 0.0;
        double itemBaseArea  = item.getDimRotacionada(0, r) * item.getDimRotacionada(1, r);

        for(int k = 0; k < bin.numItens; ++k)
        {
            int placedId = bin.vetItemId[k];
            Item &placedItem = instanciaG.vetItens[placedId];
            double topZ = bin.vetPosItem[k].vetDim[2] + placedItem.getDimRotacionada(2, bin.vetRotacao[k]);

            if(doubleEqual(topZ, p.vetDim[2], 1e-4))
            {
                double overlap = computeXY_Overlap(item, r, p, placedItem, bin.vetRotacao[k], bin.vetPosItem[k]);
                if(overlap > 0.0)
                {
                    // Non-fragile item cannot sit on a fragile item
                    if(ParseInputNS::input.fragility && placedItem.fragility && !item.fragility)
                        return false;

                    supportedArea += overlap;
                }
            }
        }

        if((supportedArea / itemBaseArea) < ParseInputNS::input.minSupportArea)
            return false;
    }

           // 4. LIFO / Unloading sequence check with already packed items
    if(ParseInputNS::input.lifo)
    {
        int posNew = findPos(rota, itemId);

        for(int k = 0; k < bin.numItens; ++k)
        {
            int placedId = bin.vetItemId[k];
            Item &placedItem = instanciaG.vetItens[placedId];

            if(item.customer == placedItem.customer)
                continue;

            int posPlaced = findPos(rota, placedId);
            if(posPlaced == posNew)
                continue;

            if(posNew < posPlaced)
            {
                Item tempNew = item;
                if(!lifo(tempNew, p, r, placedItem, bin.vetPosItem[k], bin.vetRotacao[k],
                         ParseInputNS::input.mlifo, ParseInputNS::input.removeFromShortSide, matSupportItems))
                    return false;
            }
            else
            {
                Item tempNew = item;
                if(!lifo(placedItem, bin.vetPosItem[k], bin.vetRotacao[k], tempNew, p, r,
                         ParseInputNS::input.mlifo, ParseInputNS::input.removeFromShortSide, matSupportItems))
                    return false;
            }
        }
    }

    // 5. Lateral load balance check
    if(ParseInputNS::input.balancedLoading && !ParseInputNS::input.comprimentoAlturaIguais1)
    {
        double width = item.getDimRotacionada(1, r);
        double left  = computeLeftBalancedLoading(p.vetDim[1], width, item.weight);
        double right = item.weight - left;

        double limit = std::ceil(ParseInputNS::input.balancedLoadingD * (bin.demandaTotal + item.weight));
        if((bin.sumLeftBalancedLoading + left) > limit || (bin.sumRightBalancedLoading + right) > limit)
            return false;
    }

    // 6. Left Compactness Check (Prevent floating gaps on the left side)
    if(ParseInputNS::input.compactness && p.vetDim[0] > 1e-4) // Not touching the left truck wall
    {
        double dy = item.getDimRotacionada(1, r);
        double dz = item.getDimRotacionada(2, r);
        double leftAreaTotal = dy * dz;
        double leftSupport = 0.0;

        for(int k = 0; k < bin.numItens; ++k)
        {
            int placedId = bin.vetItemId[k];
            Item &placedItem = instanciaG.vetItens[placedId];

            double placedDx = placedItem.getDimRotacionada(0, bin.vetRotacao[k]);
            double placedMaxX = bin.vetPosItem[k].vetDim[0] + placedDx;

                   // If placed item touches the left face of the new candidate item
            if(doubleEqual(placedMaxX, p.vetDim[0], 1e-4))
            {
                double placedDy = placedItem.getDimRotacionada(1, bin.vetRotacao[k]);
                double placedDz = placedItem.getDimRotacionada(2, bin.vetRotacao[k]);

                double pY = p.vetDim[1];
                double pkY = bin.vetPosItem[k].vetDim[1];
                double pZ = p.vetDim[2];
                double pkZ = bin.vetPosItem[k].vetDim[2];

                       // Calculate overlap area in the Y-Z plane
                double overlapY = std::max(0.0, std::min(pY + dy, pkY + placedDy) - std::max(pY, pkY));
                double overlapZ = std::max(0.0, std::min(pZ + dz, pkZ + placedDz) - std::max(pZ, pkZ));

                leftSupport += overlapY * overlapZ;
            }
        }

        if((leftSupport / leftAreaTotal) < ParseInputNS::input.minLeftSupportArea)
            return false;
    }

    return true;
}

/**
 * @brief Computes placement score using Deepest-Bottom-Left (DBL) and contact surface area.
 * Lower score is preferred.
 */
static double computePlacementScore(const Bin &bin, int itemId, const Ponto &p, InstanceNS::Rotation r)
{
    Item &item = instanciaG.vetItens[itemId];

    double dx = item.getDimRotacionada(0, r);
    double dy = item.getDimRotacionada(1, r);
    double dz = item.getDimRotacionada(2, r);

           // Calculate contact area with container walls
    double contactArea = 0.0;
    if(p.vetDim[0] <= 1e-4) contactArea += dy * dz;                               // Left wall
    if(doubleEqual(p.vetDim[0] + dx, bin.binDim[0], 1e-4)) contactArea += dy * dz; // Right wall
    if(p.vetDim[1] <= 1e-4) contactArea += dx * dz;                               // Back wall
    if(doubleEqual(p.vetDim[1] + dy, bin.binDim[1], 1e-4)) contactArea += dx * dz; // Front wall
    if(p.vetDim[2] <= 1e-4) contactArea += dx * dy;                               // Floor

           // DBL weights: prioritize minimum Z (floor), then minimum Y (depth), then minimum X (lateral)
    double score = 1000.0 * p.vetDim[2] + 10.0 * p.vetDim[1] + 1.0 * p.vetDim[0] - 0.05 * contactArea;
    return score;
}

struct PlacementCandidate
{
    double score;
    int epIdx;
    InstanceNS::Rotation rot;
};

bool packItemsIntoBin(SolucaoNS::Bin          &bin,
                      SolucaoNS::Rota         &rota,
                      const VectorI           &itemIds,
                      int					   numItems,
                      bool					   sortItemIds,
                      bool                     print,
                      SolucaoNS::Penalty*	   penalty)
{

    static VectorI sortedItems(MaxNumItemsBin);
    assertm(numItems > MaxNumItemsBin, "Error, numItems is greater them MaxNumItemsBin\n");

    //std::printf("size itemIds: %d\nItems: ", numItems);

    //for(int i=0; i < numItems; ++i)
    //	std::printf("%d ", itemIds[i]);
    //std::printf("\n\n");

    uint64_t it=0;
    uint64_t itMax = 10000;


    if(input.maxTimePackingHeuristic > 0)
    {
        itMax = std::numeric_limits<uint64_t>::max();
        setOffAlarm();
        setAlarm(input.maxTimePackingHeuristic);
    }


    static const SortStrategy strategies[] =
    {
        SortStrategy::ReverseRoute_AreaDesc,
        SortStrategy::ReverseRoute_VolumeDesc,
        SortStrategy::ReverseRoute_WeightDesc
    };

    if(!sortItemIds)
        itMax = 1;

    for(it=0; it < itMax; ++it)
    {
        //std::printf("it: %lld\n", it);
        if(input.maxTimePackingHeuristic > 0 && doStop() && sortItemIds)
        {
            //std::printf("doStop: %d\n", doStop());
            break;
        }

        for(SortStrategy strategy : strategies)
        {

            bin.reset();
            // Initialize with the origin Extreme Point (0, 0, 0)
            bin.addEp(Ponto(0.0, 0.0, 0.0));

            // Working copy of items to sort
            copyVet(itemIds, sortedItems, numItems);
            //std::vector<int> sortedItems(itemIds.begin(), itemIds.begin()+numItems);

            // 1. Generate persistent random noise for this packing pass (+/- 10%)
            // This allows randomized sorting while maintaining strict weak ordering rules for std::sort.
            static std::vector<double> noise(MaxNumItemsBin);
            for(int i=0; i < numItems; ++i)
                noise[sortedItems[i]] = RandNs::getRandDouble();


            // Sort items: primary key is reverse delivery position along the route (deepest first)
            if(sortItemIds)
            {
                std::sort(sortedItems.begin(), sortedItems.begin()+numItems, [&](int a, int b)
                {
                    int posA = findPos(rota, a);
                    int posB = findPos(rota, b);
                    if(posA != posB)
                        return posA > posB; // Higher delivery pos -> pack earlier (deeper in truck)

                    Item &itemA = instanciaG.vetItens[a];
                    Item &itemB = instanciaG.vetItens[b];

                       // Non-fragile before fragile (fragile on top)
                    if(itemA.fragility != itemB.fragility)
                        return !itemA.fragility;

                       // Randomized tie-breaking for items in the same delivery tier
                    switch(strategy)
                    {
                    case SortStrategy::ReverseRoute_AreaDesc: {
                        double areaA = (itemA.vetDim[0] * itemA.vetDim[1]) * noise[a];
                        double areaB = (itemB.vetDim[0] * itemB.vetDim[1]) * noise[b];
                        if(!doubleEqual(areaA, areaB, 1e-4)) return areaA > areaB;
                        return (itemA.volume * noise[a]) > (itemB.volume * noise[b]);
                    }
                    case SortStrategy::ReverseRoute_VolumeDesc: {
                        double volA = itemA.volume * noise[a];
                        double volB = itemB.volume * noise[b];
                        if(!doubleEqual(volA, volB, 1e-4)) return volA > volB;

                        return (itemA.weight * noise[a]) > (itemB.weight * noise[b]);
                    }
                    case SortStrategy::ReverseRoute_WeightDesc: {
                        double wA = itemA.weight * noise[a];
                        double wB = itemB.weight * noise[b];
                        if(!doubleEqual(wA, wB, 1e-4)) return wA > wB;

                        return (itemA.volume * noise[a]) > (itemB.volume * noise[b]);
                    }
                    }
                    return a < b; // stable fallback
                });

            }

                   // Scratch support matrix for LIFO check
            static Eigen::Matrix<int, -1, -1, Eigen::RowMajor>
                matSupportItems(instanciaG.numItens, instanciaG.numItens);

            matSupportItems.setConstant(0);

                   // Number of rotations to try (2 if vertical stability only, or up to 6 if 3D rotation allowed)
            int numRotationsToTry = (instanciaG.numRotation >= 6) ? 6 : 2;
            bool isInfeasible = false;

            double fk(0.0), fFA(0.0), fRA(0.0), fTA(0.0), sumM(0.0), sumF(0.0);

            int lastItem = -1;

            for(int i=0; i < numItems; ++i)
            {
                if(i > 0)
                {
                    AxleWeightsNS::semiTrailer.checkAxleWeights(bin, false, &fk, &fFA,
                                                                &fRA, &fTA, &sumM, &sumF);
                }

                int itemId = sortedItems[i];
                static Vector<PlacementCandidate> candidates(NumEpPorBin);
                int numOfCandidates = 0;
                assertm(bin.numEps > NumEpPorBin,
                        "Error, bin.numEps is greter them NumEpPorBin\n")

                // Evaluate all active Extreme Points
                //std::printf("bin.numEps %d\n", bin.numEps);
                for(int epIdx = 0; epIdx < bin.numEps; ++epIdx)
                {
                    //std::printf("\tepIdx: %d\n", epIdx);
                    const Ponto &ep = bin.vetEp[epIdx];

                    for(int rIdx = 0; rIdx < 2; ++rIdx)
                    {
                        Rotation r = vetRot[rIdx];

                        if(isPlacementFeasible(bin, rota, itemId, ep, r, matSupportItems,
                                               sumM, sumF))
                        {
                            double score = computePlacementScore(bin, itemId, ep, r);
                            candidates[numOfCandidates] = {score, epIdx, r};
                            numOfCandidates += 1;
                        }
                    }
                }

                // If no feasible EP found, packing fails under this strategy
                if(numOfCandidates == 0)
                {
                    if(print)
                        std::printf("Failed to pack item %d (cust: %d)\n", itemId, instanciaG.vetItens[itemId].customer);

                    isInfeasible = true;
                    lastItem = i;
                    break;
                    //return false;
                }

                // 2. Restricted Candidate List (RCL)
                std::sort(candidates.begin(), candidates.begin()+numOfCandidates,
                          [](const PlacementCandidate &a, const PlacementCandidate &b)
                          {return a.score < b.score;});

                int rclSize = std::max(1, (int)(numOfCandidates*input.aphaBinEscolhaEp));
                //std::printf("\t\trclSize: %d\n", rclSize);
                int chosen  = RandNs::getRandInt(0, rclSize-1);

                if(!sortItemIds)
                    chosen = 0;

                //std::printf("\t\tchosen: %d\n", chosen);

                // Place the item (this automatically updates EPs, dimensions, and loads)
                bin.addItem(candidates[chosen].epIdx, itemId, candidates[chosen].rot);
            }

            if(isInfeasible && sortItemIds)
                continue;

            if(isInfeasible && !sortItemIds && penalty)
            {
                if(lastItem < 0 || lastItem >= sortedItems.size() )
                {
                    std::printf("Error, isInfeasible: %d; lastItem: %d\n\n", isInfeasible, lastItem);
                    PRINT_THROW();
                }


                double vol = 0.0;
                for(int i=lastItem; i < numItems; ++i)
                {
                    vol += instanciaG.vetItens[sortedItems[i]].volume;
                }

                penalizeSolution(rota, bin, *penalty);
                penalty->set(SolucaoNS::Penalty::TypePenalty::VolItemsNotPacked, vol);


                return false;
            }

            // Comprehensive final feasibility check (axle weights, compactness, support, LIFO)

            // TODO remover
            //Penalty penalty;

            if(bin.checkFeasibility(&rota, false, print))
            {
                if(numItems != bin.numItens)
                {
                    std::printf("Erro, \nitemIds: %s\n", itemIds.printN(numItems).c_str());
                    std::printf("bin.vetItens: %s", bin.vetItens.printN(bin.numItens).c_str());

                    PRINT_THROW();
                }

                std::set<int> setItems;
                for(int i=0; i < numItems; ++i)
                {
                    if(setItems.contains(itemIds[i]))
                    {
                        std::printf("Repeted item\n");
                        PRINT_THROW();
                    }
                    else
                        setItems.insert(itemIds[i]);
                }

                for(int i=0; i < numItems; ++i)
                {
                    if(!setItems.contains(bin.vetItemId[i]))
                    {
                        std::printf("Error, item(%d) was not packing\n", bin.vetItens[i]);

                        std::printf("Erro, \nitemIds: %s\n", itemIds.printN(numItems).c_str());
                        std::printf("bin.vetItens: %s", bin.vetItens.printN(bin.numItens).c_str());

                        PRINT_THROW();
                    }
                }


                return true;
            }
            else if(!sortItemIds && penalty)
            {

                penalizeSolution(rota, bin, *penalty);
                penalty->set(SolucaoNS::Penalty::TypePenalty::VolItemsNotPacked, 0.0);


                return false;
            }


            if(!sortItemIds)
                return false;
        }
    }

    return false;
}

bool packRoute(SolucaoNS::Rota &rota, VectorI _itemIds, bool print)
{
    print = true;
    bool capVolFeas = testRouteCapVol(&rota.vetRota[0], rota.numPos);
    if(!capVolFeas)
        return false;


    if(!rota.binPtr)
    {
        if(print) std::printf("Error: rota.binPtr is null\n");
        return false;
    }

    Bin &bin = *rota.binPtr;
    bin.reset();
           // Route has no customers (only depot visits 0 -> 0)
    if(rota.numPos <= 2)
    {
        bin.reset();
        return true;
    }

    // 1. Collect all items demanded by customers in this route
    VectorI itemIds;
    int numItems;

    if(!_itemIds.empty())
    {
        itemIds = _itemIds;
        numItems = itemIds.size();
    }
    else
        numItems = copiaItensClientes(rota.vetRota, rota.numPos, itemIds, true);

    if(itemIds.empty())
    {
        bin.reset();
        return true;
    }

    return packItemsIntoBin(bin, rota, itemIds, numItems, true, false);



} // namespace ConstrutivoBinNS

}
