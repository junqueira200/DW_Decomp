/* ****************************************
 * ****************************************
 *  Data:    14/09/26
 *  Arquivo: ConstrutivoBin.cpp
 *  Autor:   Igor de Andrade Junqueira
 *  Projeto: 2L/3L-SDVRP
 * ****************************************
 * ****************************************/

#include "ConstrutivoBin2.h"
#include "Alarm.h"
#include "AuxGeometry.h"
#include "AuxT.h"
#include "AxleWeights.h"
#include "ConstrutivoBin.h"
#include "InputOutput.h"
#include "c_api.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>

using namespace InstanceNS;
using namespace SolucaoNS;
using namespace AuxGeometryNS;
using namespace ConstrutivoBinNS;
using namespace AxleWeightsNS;
using namespace StopDataNS;

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
        int lastItem = instanciaG.matCliItensIniFim(customer, 1);

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
// Check where the last call to isPlacementFeasibleModel0 failed
// (PackingErro index; PackingErroSize = Boundary).
static int placementFailErr = PackingErroSize;

// Order in which the checks run: a larger value means the candidate got further.
static int failDepth(int err)
{
    switch(err)
    {
        case PackingErroSize:          return 0;   // Boundary
        case PackingErroGeometric:     return 1;
        case PackingErroSupport:
        case PackingErroFragility:     return 2;
        case PackingErroLifo:          return 3;
        case PackingErroLoadBalancing: return 4;
        case PackingErroAxleWights:    return 5;
        case PackingErroCompactness:   return 6;
    }
    return 0;
}

// Soft-placement fallback: minimum support (fraction of base area) a soft-placed item must
// still have, so it never floats. Below this, support stays a hard constraint.
constexpr double MinSoftSupport = 0.5;

static inline double softW(SolucaoNS::Penalty::TypePenalty t)
{
    return ParseInputNS::input.arrayPenaltyWeight[(std::size_t)t];
}

static bool isPlacementFeasibleModel0(const Bin           &bin,
                                      Rota                &rota,
                                      int                  itemId,
                                      const Ponto         &p,
                                      InstanceNS::Rotation r,
                                      double sumM,
                                      double sumF,
                                      StopData *stopData,
                                      double   *ptrViolation = nullptr)  // soft mode if != nullptr
{
    using TP = SolucaoNS::Penalty::TypePenalty;
    Item &item = instanciaG.vetItens[itemId];

           // 1. Boundary check
    for(int d = 0; d < instanciaG.numDim; ++d)
    {
        double dimR = item.getDimRotacionada(d, r);
        if(p.vetDim[d] + dimR > bin.binDim[d] + 1e-5)
        {
            placementFailErr = PackingErroSize;   // Boundary
            return false;
        }
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
        if(verificaColisaoDoisItens(
               itemId, bin.vetItemId[k], p, bin.vetPosItem[k], r, bin.vetRotacao[k]))
        {
            vetPackinErros[PackingErroGeometric] += 1;
            placementFailErr = PackingErroGeometric;
            return false;
        }
    }

    // 3. Bottom Support & Fragility check (if 3D and item is not resting on the floor)
    if(instanciaG.numDim == 3 && ParseInputNS::input.support && p.vetDim[2] > 1e-4)
    {
        double supportedArea = 0.0;
        double itemBaseArea = item.getDimRotacionada(0, r) * item.getDimRotacionada(1, r);

        for(int k = 0; k < bin.numItens; ++k)
        {

            int    placedId = bin.vetItemId[k];
            Item  &placedItem = instanciaG.vetItens[placedId];
            double topZ =
                bin.vetPosItem[k].vetDim[2] + placedItem.getDimRotacionada(2, bin.vetRotacao[k]);

            if(doubleEqual(topZ, p.vetDim[2], 1e-4))
            {
                double overlap =
                    computeXY_Overlap(item, r, p, placedItem, bin.vetRotacao[k], bin.vetPosItem[k]);

                if(overlap > 0.0)
                {
                    // Non-fragile item cannot sit on a fragile item
                    if(ParseInputNS::input.fragility && placedItem.fragility && !item.fragility)
                    {
                        vetPackinErros[PackingErroFragility] += 1;
                        placementFailErr = PackingErroFragility;
                        return false;
                    }

                    supportedArea += overlap;
                }
            }
        }

        const double reqSupport = std::ceil(ParseInputNS::input.minSupportArea * itemBaseArea);
        if(supportedArea < reqSupport)
        {
            if(ptrViolation && supportedArea >= MinSoftSupport * itemBaseArea)
                *ptrViolation += softW(TP::Support) * (reqSupport - supportedArea) / reqSupport;
            else
            {
                vetPackinErros[PackingErroSupport] += 1;
                placementFailErr = PackingErroSupport;
                return false;
            }
        }
    }

    static VectorI vetIndex(instanciaG.numItens);
    vetIndex.setAll(-1);
    for(int k = 0; k < bin.numItens; ++k)
        vetIndex[bin.vetItemId[k]] = k;

           // 4. LIFO / Unloading sequence check with already packed items
    if(ParseInputNS::input.lifo)
    {
        int    posNew = findPos(rota, itemId);
        double topNew = p.vetDim[2] + item.getDimRotacionada(2, r);

        for(int k = 0; k < bin.numItens; ++k)
        {
            int      placedId = bin.vetItemId[k];
            Item    &placedItem = instanciaG.vetItens[placedId];
            Ponto    pPlaced = bin.vetPosItem[k];
            Rotation rPlaced = bin.vetRotacao[k];

            if(item.customer == placedItem.customer)
                continue;

            int posPlaced = findPos(rota, placedId);
            if(posPlaced == posNew)
                continue;

                   // --- NEW: DYNAMICALLY CALCULATE VERTICAL SUPPORT ---
            bool   new_supports_placed = false;
            bool   placed_supports_new = false;
            double topPlaced = pPlaced.vetDim[2] + placedItem.getDimRotacionada(2, rPlaced);

                   // Does candidate support placedItem? (Candidate is directly under placedItem)
            if(doubleEqual(topNew, pPlaced.vetDim[2], 1e-4))
            {
                if(computeXY_Overlap(item, r, p, placedItem, rPlaced, pPlaced) > 0.0)
                    new_supports_placed = true;
            }

            // Does placedItem support candidate? (placedItem is directly under Candidate)
            if(doubleEqual(topPlaced, p.vetDim[2], 1e-4))
            {
                if(computeXY_Overlap(item, r, p, placedItem, rPlaced, pPlaced) > 0.0)
                    placed_supports_new = true;
            }
            // ---------------------------------------------------

            if(posNew < posPlaced)
            {
                Item tempNew = item;
                if(!lifo(tempNew,
                         p,
                         r,
                         placedItem,
                         pPlaced,
                         rPlaced,
                         ParseInputNS::input.mlifo,
                         ParseInputNS::input.removeFromShortSide,
                         placed_supports_new,  // item1_supports_item0
                         new_supports_placed)) // item0_supports_item1
                {
                    if(ptrViolation)
                        *ptrViolation += softW(TP::Lifo);
                    else
                    {
                        vetPackinErros[PackingErroLifo] += 1;
                        placementFailErr = PackingErroLifo;
                        return false;
                    }
                }
            }
            else
            {
                Item tempNew = item;
                if(!lifo(placedItem,
                         pPlaced,
                         rPlaced,
                         tempNew,
                         p,
                         r,
                         ParseInputNS::input.mlifo,
                         ParseInputNS::input.removeFromShortSide,
                         new_supports_placed,  // item1_supports_item0
                         placed_supports_new)) // item0_supports_item1
                {
                    if(ptrViolation)
                        *ptrViolation += softW(TP::Lifo);
                    else
                    {
                        vetPackinErros[PackingErroLifo] += 1;
                        placementFailErr = PackingErroLifo;
                        return false;
                    }
                }
            }
        }
    }

    //return true;

    // 5. Lateral load balance check
    if(ParseInputNS::input.balancedLoading && !ParseInputNS::input.comprimentoAlturaIguais1)
    {

        static VectorD vetLeft(10);
        static VectorD vetRight(10);

        vetLeft.setAll(0.0);
        vetRight.setAll(0.0);

        // Model 0 only checks the initial load (stop 1)
        const int lastStop = (input.modelType <= 0) ? 1 : (rota.numPos - 2);

        // Candidate item
        {
            double width = item.getDimRotacionada(1, r);
            double left  = computeLeftBalancedLoading(p.vetDim[1], width, item.weight);
            double right = item.weight - left;

            for(int s = 1; s <= lastStop; ++s)
            {
                if(!stopData->vetStopsSet[s].contains(itemId))
                    break;
                vetLeft[s]  += left;
                vetRight[s] += right;
            }
        }

        // Items already in the bin (own rotation and own weight)
        for(int k = 0; k < bin.numItens; ++k)
        {
            int    placedId   = bin.vetItemId[k];
            Item  &placedItem = instanciaG.vetItens[placedId];
            double width = placedItem.getDimRotacionada(1, bin.vetRotacao[k]);
            double left  = computeLeftBalancedLoading(bin.vetPosItem[k].vetDim[1], width,
                                                      placedItem.weight);
            double right = placedItem.weight - left;

            for(int s = 1; s <= lastStop; ++s)
            {
                if(!stopData->vetStopsSet[s].contains(placedId))
                    break;
                vetLeft[s]  += left;
                vetRight[s] += right;
            }
        }

        // Side loads only grow as items are added: once a side exceeds the limit
        // of the stop's final load, no later placement can repair it.
        for(int s = 1; s <= lastStop; ++s)
        {
            double limit = std::ceil(input.balancedLoadingD * stopData->vetStopWeight[s]);
            if(vetLeft[s] > limit || vetRight[s] > limit)
            {
                // Hard even in soft mode: the excess stays at stop s, so every later item on
                // board at s would fail too and cascade into the fallback.
                vetPackinErros[PackingErroLoadBalancing] += 1;
                placementFailErr = PackingErroLoadBalancing;
                return false;
            }
        }
    }


    // 5b. Axle weights check on the partial load, for every stop the candidate is on board.
    // Stops without the candidate are unchanged and were checked when their items were placed.
    // Adding an item increases the axle loads unless its centre is behind the trailer axle,
    // so an overload in a partial load almost always remains in the final load.
    if(ParseInputNS::input.axleWights)
    {
        using AxleWeightsNS::semiTrailer;
        static VectorD vetSumF(10);
        static VectorD vetSumM(10);

        vetSumF.setAll(0.0);
        vetSumM.setAll(0.0);

        const int lastStop = (input.modelType <= 0) ? 1 : (rota.numPos - 2);
        const double g = GravityCmConst;

        auto addLoad = [&](int id, double x, Rotation rot)
        {
            Item &it = instanciaG.vetItens[id];
            double f  = it.weight * g;
            double rr = (double)semiTrailer.distanceCargoSpaceTrailerAxle - x -
                        it.getDimRotacionada(0, rot) / 2.0;

            for(int s = 1; s <= lastStop; ++s)
            {
                if(!stopData->vetStopsSet[s].contains(id))
                    break;
                vetSumF[s] += f;
                vetSumM[s] += f * rr;
            }
        };

        addLoad(itemId, p.vetDim[0], r);
        for(int k = 0; k < bin.numItens; ++k)
            addLoad(bin.vetItemId[k], bin.vetPosItem[k].vetDim[0], bin.vetRotacao[k]);

        for(int s = 1; s <= lastStop; ++s)
        {
            if(!stopData->vetStopsSet[s].contains(itemId))
                break;

            double fK = (vetSumM[s] + (double)semiTrailer.massTrailer * g *
                                          semiTrailer.distanceMassTrailerTrailerAxle) /
                        (double)semiTrailer.distanceKingpinTrailerAxle;

            double fFA = (fK * (double)semiTrailer.distanceKingpinRearAxle +
                          (double)semiTrailer.massTractor * g *
                              semiTrailer.distanceMassTractorRearAxle) /
                         (double)semiTrailer.wheelBase;

            double fRA = fK + (double)semiTrailer.massTractor * g - fFA;
            double fTA = vetSumF[s] + (double)semiTrailer.massTrailer * g - fK;

            const double limFA = (double)semiTrailer.maxMassFrontAxle * g;
            const double limRA = (double)semiTrailer.maxMassRearAxle * g;
            const double limTA = (double)semiTrailer.maxMassTrailerAxle * g;

            if(fFA > limFA + 1e-6 || fRA > limRA + 1e-6 || fTA > limTA + 1e-6)
            {
                // Hard even in soft mode (same cascade as balance).
                vetPackinErros[PackingErroAxleWights] += 1;
                placementFailErr = PackingErroAxleWights;
                return false;
            }
        }
    }

    // 6. Left Compactness Check (Prevent floating gaps on the left side)
    if(ParseInputNS::input.compactness && p.vetDim[0] > 1e-4) // Not touching the left truck wall
    {
        double dy = item.getDimRotacionada(1, r);
        double dz = item.getDimRotacionada(2, r);
        double leftAreaTotal = dy * dz;
        double leftSupport = 0.0;

        // Model 0 only checks the initial load (stop 1)
        const int lastStop = (input.modelType <= 0) ? 1 : (rota.numPos - 2);

        for(int s = 1; s <= lastStop; ++s)
        {
            if(!stopData->vetStopsSet[s].contains(itemId))
                break;

            leftSupport = 0.0;

            //for(int k = 0; k < bin.numItens; ++k)
            for(int itemJ:stopData->vetStopsSet[s])
            {
                int k = vetIndex[itemJ];
                if(k < 0)
                    continue;

                int   placedId = bin.vetItemId[k];
                if(placedId != itemJ)
                {
                    std::printf("Error, placedId != itemJ\n");
                    PRINT_THROW();
                }

                Item &placedItem = instanciaG.vetItens[placedId];

                double placedDx = placedItem.getDimRotacionada(0, bin.vetRotacao[k]);
                double placedMaxX = bin.vetPosItem[k].vetDim[0] + placedDx;

                double dif = p.vetDim[0] - placedMaxX;
                if(dif >= 0.0 && dif <= DifDistColision)
                {
                    double placedDy = placedItem.getDimRotacionada(1, bin.vetRotacao[k]);
                    double placedDz = placedItem.getDimRotacionada(2, bin.vetRotacao[k]);

                    double pY = p.vetDim[1];
                    double pkY = bin.vetPosItem[k].vetDim[1];
                    double pZ = p.vetDim[2];
                    double pkZ = bin.vetPosItem[k].vetDim[2];

                           // Calculate overlap area in the Y-Z plane
                    double overlapY =
                        std::max(0.0, std::min(pY + dy, pkY + placedDy) - std::max(pY, pkY));
                    double overlapZ =
                        std::max(0.0, std::min(pZ + dz, pkZ + placedDz) - std::max(pZ, pkZ));

                    leftSupport += overlapY * overlapZ;
                }
            }

            const double reqLeft = std::ceil(ParseInputNS::input.minLeftSupportArea * leftAreaTotal);
            if(leftSupport < reqLeft)
            {
                if(ptrViolation)
                    *ptrViolation += softW(TP::Compactness) * (reqLeft - leftSupport) / reqLeft;
                else
                {
                    vetPackinErros[PackingErroCompactness] += 1;
                    placementFailErr = PackingErroCompactness;
                    return false;
                }
            }
        }

    }

    return true;
}

/**
 * @brief Computes placement score using Deepest-Bottom-Left (DBL) and contact surface area.
 * Lower score is preferred.
 */
static double
computePlacementScore(const Bin &bin, int itemId, const Ponto &p, InstanceNS::Rotation r)
{
    Item &item = instanciaG.vetItens[itemId];

    double dx = item.getDimRotacionada(0, r);
    double dy = item.getDimRotacionada(1, r);
    double dz = item.getDimRotacionada(2, r);

    /*
    // Calculate contact area with container walls
    double contactArea = 0.0;
    if(p.vetDim[0] <= 1e-4)
        contactArea += dy * dz; // Left wall
    if(doubleEqual(p.vetDim[0] + dx, bin.binDim[0], 1e-4))
        contactArea += dy * dz; // Right wall
    if(p.vetDim[1] <= 1e-4)
        contactArea += dx * dz; // Back wall
    if(doubleEqual(p.vetDim[1] + dy, bin.binDim[1], 1e-4))
        contactArea += dx * dz; // Front wall
    if(p.vetDim[2] <= 1e-4)
        contactArea += dx * dy; // Floor

     // DBL weights: prioritize minimum Z (floor), then minimum Y (depth), then minimum X (lateral)
     double score =
         1000.0 * p.vetDim[2] + 10.0 * p.vetDim[1] + 1.0 * p.vetDim[0] - 0.05 * contactArea;

      */

           // Ideal longitudinal position (X-axis) for the trailer axle sweet spot
    double targetX = AxleWeightsNS::semiTrailer.distanceCargoSpaceTrailerAxle * 0.5;
    double itemCenterX = p.vetDim[0] + dx / 2.0;
    double axleDistPenalty = std::abs(itemCenterX - targetX);
    // Lateral balance: penalize deviation from truck centerline (Y-axis)
    double truckCenterY = bin.binDim[1] / 2.0;
    double itemCenterY = p.vetDim[1] + dy / 2.0;
    double lateralPenalty = std::abs(itemCenterY - truckCenterY);
    // Combine DBL with stability & balance
    double score = 1000.0 * p.vetDim[2]          // Keep low (stability)
                   + 10.0 * p.vetDim[0]            // Depth
                   + 2.0 * lateralPenalty          // Lateral balance
                   + 1.0 * axleDistPenalty;        // Longitudinal distribution

    return score;
}

struct PlacementCandidate
{
    double               score;
    int                  epIdx;
    InstanceNS::Rotation rot;
    Ponto                pos;      // used when epIdx == -1 (fallback position, not an EP)
};

bool packItemsIntoBin(SolucaoNS::Bin     	&bin,
                      SolucaoNS::Rota    	&rota,
                      const VectorI      	&itemIds,
                      int                 	numItems,
                      bool                	sortItemIds,
                      bool                	print,
                      SolucaoNS::Penalty 	*penalty,
                      bool				  	resetBin,
                      StopDataNS::StopData	*stopData,
                      double				*randomKey)
{
    double totalVol = 0.0;
    double leftArea = 0.0;

    for(int i=0; i < numItems; ++i)
    {
        totalVol += instanciaG.vetItens[itemIds[i]].volume;
        double a0 = instanciaG.vetItens[itemIds[i]].vetDim[2]*
                    instanciaG.vetItens[itemIds[i]].vetDim[0];

        double a1 = instanciaG.vetItens[itemIds[i]].vetDim[2]*
                    instanciaG.vetItens[itemIds[i]].vetDim[1];

        if(a0 > a1)
            leftArea += a0;
        else
            leftArea += a1;
    }

    static VectorI sortedItems(MaxNumItemsBin);
    static VectorI vetNotPackedItems(instanciaG.numItens);
    int sizeVetNotPackedItems = 0;

    assertm(numItems > MaxNumItemsBin, "Error, numItems is greater them MaxNumItemsBin\n");

           //std::printf("size itemIds: %d\nItems: ", numItems);

           //for(int i=0; i < numItems; ++i)
           //    std::printf("%d ", itemIds[i]);
           //std::printf("\n\n");

    uint64_t it = 0;
    uint64_t itMax = 10000;

    if(input.maxTimePackingHeuristic > 0 && sortItemIds)
    {
        itMax = std::numeric_limits<uint64_t>::max();
        setOffAlarm();
        setAlarm(input.maxTimePackingHeuristic);
    }

    static const SortStrategy strategies[] = {SortStrategy::ReverseRoute_AreaDesc,
                                              SortStrategy::ReverseRoute_VolumeDesc,
                                              SortStrategy::ReverseRoute_WeightDesc};

    if(!sortItemIds)
        itMax = 1;

    for(it = 0; it < itMax; ++it)
    {
        //std::printf("it: %lld\n", it);
        if(input.maxTimePackingHeuristic > 0 && doStop() && sortItemIds)
        {
            // std::printf("doStop: %d\n", doStop());
            break;
        }

        for(SortStrategy strategy : strategies)
        {
            vetLastNotPacked.clear();
            vetLastSoftPlaced.clear();
            // Initialize with the origin Extreme Point (0, 0, 0)
            if(resetBin)
            {
                bin.reset();
                bin.addEp(Ponto(0.0, 0.0, 0.0));
            }

                   // Working copy of items to sort
            copyVet(itemIds, sortedItems, numItems);
            // std::vector<int> sortedItems(itemIds.begin(), itemIds.begin()+numItems);

                   // 1. Generate persistent random noise for this packing pass (+/- 10%)
                   // This allows randomized sorting while maintaining strict weak ordering rules for
                   // std::sort.
            static Vector<double> noise(instanciaG.numItens);
            for(int i = 0; i < numItems; ++i)
                noise[sortedItems[i]] = RandNs::getRandDouble();

                   // Sort items: primary key is reverse delivery position along the route (deepest first)
            if(sortItemIds)
            {
                std::sort(
                    sortedItems.begin(),
                    sortedItems.begin() + numItems,
                    [&](int a, int b)
                    {
                        int posA = findPos(rota, a);
                        int posB = findPos(rota, b);
                        if(posA != posB)
                            return posA >
                                   posB; // Higher delivery pos -> pack earlier (deeper in truck)

                        Item &itemA = instanciaG.vetItens[a];
                        Item &itemB = instanciaG.vetItens[b];

                               // Non-fragile before fragile (fragile on top)
                        if(itemA.fragility != itemB.fragility)
                            return !itemA.fragility;

                               // Randomized tie-breaking for items in the same delivery tier
                        switch(strategy)
                        {
                        case SortStrategy::ReverseRoute_AreaDesc:
                        {
                            double areaA = (itemA.vetDim[0] * itemA.vetDim[1]) * noise[a];
                            double areaB = (itemB.vetDim[0] * itemB.vetDim[1]) * noise[b];
                            if(!doubleEqual(areaA, areaB, 1e-4))
                                return areaA > areaB;
                            return (itemA.volume * noise[a]) > (itemB.volume * noise[b]);
                        }
                        case SortStrategy::ReverseRoute_VolumeDesc:
                        {
                            double volA = itemA.volume * noise[a];
                            double volB = itemB.volume * noise[b];
                            if(!doubleEqual(volA, volB, 1e-4))
                                return volA > volB;

                            return (itemA.weight * noise[a]) > (itemB.weight * noise[b]);
                        }
                        case SortStrategy::ReverseRoute_WeightDesc:
                        {
                            double wA = itemA.weight * noise[a];
                            double wB = itemB.weight * noise[b];
                            if(!doubleEqual(wA, wB, 1e-4))
                                return wA > wB;

                            return (itemA.volume * noise[a]) > (itemB.volume * noise[b]);
                        }
                        }
                        return a < b; // stable fallback
                    });
            }

                   // Scratch support matrix for LIFO check
                   // static Eigen::Matrix<int, -1, -1, Eigen::RowMajor>
                   //     matSupportItems(instanciaG.numItens, instanciaG.numItens);

                   // matSupportItems.setConstant(0);

                   // Number of rotations to try (2 if vertical stability only, or up to 6 if 3D rotation
                   // allowed)
            int  numRotationsToTry = (instanciaG.numRotation >= 6) ? 6 : 2;
            bool isInfeasible = false;

            double fk(0.0), fFA(0.0), fRA(0.0), fTA(0.0), sumM(0.0), sumF(0.0);

            int lastItem = -1;

            for(int i = 0; i < numItems; ++i)
            {
                /*
                if(i > 0)
                {
                    AxleWeightsNS::semiTrailer.checkAxleWeights(bin, false, &fk, &fFA,
                                                                &fRA, &fTA, &sumM, &sumF);
                }
                */

                int                               itemId = sortedItems[i];
                // EPs x up to 6 rotations x (EP + 2 fallback positions)
                static Vector<PlacementCandidate> candidates(NumEpPorBin * 18);
                int                               numOfCandidates = 0;
                int deepestErr   = PackingErroSize;
                int deepestDepth = -1;
                assertm(bin.numEps > NumEpPorBin, "Error, bin.numEps is greter them NumEpPorBin\n")

                // Evaluate all active Extreme Points
                // std::printf("bin.numEps %d\n", bin.numEps);
                for(int epIdx = 0; epIdx < bin.numEps; ++epIdx)
                {
                    //int epIdx = randomKey[numItems+instanciaG.numItens+itemId]*(bin.numEps-1);
                    //epIdx = epIdx%bin.numEps;
                    // std::printf("\tepIdx: %d\n", epIdx);
                    const Ponto &ep = bin.vetEp[epIdx];

                    for(int rIdx = 0; rIdx < vetRot.size(); ++rIdx)
                    {
                        //int rIdx = randomKey[numItems+itemId]*(vetRot.size()-1);
                        //rIdx = rIdx%vetRot.size();
                        Rotation r = vetRot[rIdx];

                        if(isPlacementFeasibleModel0(bin, rota, itemId, ep, r, sumM, sumF,
                                                     stopData))
                        {
                            double score = computePlacementScore(bin, itemId, ep, r);
                            candidates[numOfCandidates] = {score, epIdx, r, ep};
                            numOfCandidates += 1;
                        }
                        else if(failDepth(placementFailErr) > deepestDepth)
                        {
                            deepestDepth = failDepth(placementFailErr);
                            deepestErr   = placementFailErr;
                        }
                    }
                }

                // Fallback when no EP works: try the item centred across the width and against
                // the right wall, at each EP's x and z. EPs only give y = 0 or flush with a
                // neighbour, so a lone heavy item can never be balanced at a late stop.
                if(numOfCandidates == 0)
                {
                    for(int epIdx = 0; epIdx < bin.numEps; ++epIdx)
                    {
                        for(int rIdx = 0; rIdx < vetRot.size(); ++rIdx)
                        {
                            Rotation     r  = vetRot[rIdx];
                            double       dy = instanciaG.vetItens[itemId].getDimRotacionada(1, r);
                            const double vetY[] = {(bin.binDim[1] - dy) / 2.0,  // centred
                                                   bin.binDim[1] - dy};         // right wall

                            for(double y : vetY)
                            {
                                if(y < 0.0 || doubleEqual(y, bin.vetEp[epIdx].vetDim[1], 1e-3))
                                    continue;

                                Ponto pc = bin.vetEp[epIdx];
                                pc.vetDim[1] = y;

                                if(isPlacementFeasibleModel0(bin, rota, itemId, pc, r, sumM, sumF,
                                                             stopData))
                                {
                                    double score = computePlacementScore(bin, itemId, pc, r);
                                    candidates[numOfCandidates] = {score, -1, r, pc};
                                    numOfCandidates += 1;
                                }
                                else if(failDepth(placementFailErr) > deepestDepth)
                                {
                                    deepestDepth = failDepth(placementFailErr);
                                    deepestErr   = placementFailErr;
                                }
                            }
                        }
                    }
                }

                // Last resort (GA decoding only): no position satisfies every constraint.
                // Place the item where the weighted soft violation is smallest, so the GA can
                // compare how far each layout is from feasible instead of all being "item missing".
                // Geometry, fragility and a minimum support stay hard. The final checkFeasibility
                // still rejects the packing, and penalizeSolution measures the violation.
                if(numOfCandidates == 0 && penalty)
                {
                    double   bestViol = std::numeric_limits<double>::max();
                    int      bestEp   = -1;
                    Rotation bestRot  = vetRot[0];
                    Ponto    bestPos;

                    for(int epIdx = 0; epIdx < bin.numEps; ++epIdx)
                    {
                        for(int rIdx = 0; rIdx < vetRot.size(); ++rIdx)
                        {
                            Rotation     r  = vetRot[rIdx];
                            double       dy = instanciaG.vetItens[itemId].getDimRotacionada(1, r);
                            const double vetY[] = {bin.vetEp[epIdx].vetDim[1],       // the EP
                                                   (bin.binDim[1] - dy) / 2.0,       // centred
                                                   bin.binDim[1] - dy};              // right wall

                            for(int k = 0; k < 3; ++k)
                            {
                                if(vetY[k] < 0.0)
                                    continue;

                                Ponto pc = bin.vetEp[epIdx];
                                pc.vetDim[1] = vetY[k];

                                double viol = 0.0;
                                if(isPlacementFeasibleModel0(bin, rota, itemId, pc, r, sumM, sumF,
                                                             stopData, &viol) &&
                                   viol < bestViol)
                                {
                                    bestViol = viol;
                                    bestEp   = (k == 0) ? epIdx : -1;
                                    bestRot  = r;
                                    bestPos  = pc;
                                }
                            }
                        }
                    }

                    if(bestViol < std::numeric_limits<double>::max())
                    {
                        candidates[0]   = {bestViol, bestEp, bestRot, bestPos};
                        numOfCandidates = 1;
                        vetLastSoftPlaced.push_back({itemId, deepestErr});
                    }
                }

                       // If no feasible EP found, packing fails under this strategy
                if(numOfCandidates == 0)
                {
                    // The furthest check any candidate reached is what blocked this item.
                    vetDeadEndItem[deepestErr] += 1;
                    vetLastNotPacked.push_back({itemId, deepestErr});

                    if(print)
                        std::printf("Failed to pack item %d (cust: %d)\n",
                                    itemId,
                                    instanciaG.vetItens[itemId].customer);

                    vetNotPackedItems[sizeVetNotPackedItems] = itemId;
                    sizeVetNotPackedItems += 1;
                    isInfeasible = true;

                    if(penalty)
                        continue;  	// for(int i = 0; i < numItems; ++i)
                    else
                        break;		// for(int i = 0; i < numItems; ++i)
                    //lastItem = i;
                    //break;
                    // return false;
                }

                       // 2. Restricted Candidate List (RCL)
                std::sort(candidates.begin(),
                          candidates.begin() + numOfCandidates,
                          [](const PlacementCandidate &a, const PlacementCandidate &b)
                          { return a.score < b.score; });

                // numOfCandidates
                // Let the gene pick among all feasible candidates; the pow bias keeps the best likely.
                int rclSize = randomKey ? numOfCandidates
                                        : std::max(1, (int)(numOfCandidates * input.aphaBinEscolhaEp));
                // std::printf("\t\trclSize: %d\n", rclSize);
                int chosen = 0;//randomKey[numItems+itemId]*rclSize;

                //if(!sortItemIds)
                //    chosen = 0;

                if(randomKey)
                {
                    double gen = randomKey[numItems+itemId];
                    // Candidates are sorted by score (best first). Raising gen to a power > 1
                    // pushes values toward 0, so the best positions are picked more often.
                    constexpr double CandidateBias = 3.0;   // 1 = uniform
                    chosen = (int)(std::pow(gen, CandidateBias) * rclSize);
                    chosen = std::min(chosen, rclSize - 1);
                }
                else
                    chosen = RandNs::getRandInt(0, rclSize - 1);

                       // std::printf("\t\tchosen: %d\n", chosen);

                       // Place the item (this automatically updates EPs, dimensions, and loads)
                if(candidates[chosen].epIdx >= 0)
                    bin.addItem(candidates[chosen].epIdx, itemId, candidates[chosen].rot);
                else
                    bin.addItemAt(candidates[chosen].pos, itemId, candidates[chosen].rot);
            } // END for(int i = 0; i < numItems; ++i)

            if(isInfeasible && sortItemIds)
                continue;

            if(isInfeasible && !sortItemIds && penalty)
            {

                double vol = 0.0;
                //for(int i = lastItem; i < numItems; ++i)
                for(int i=0; i < sizeVetNotPackedItems; ++i)
                {
                    int item = vetNotPackedItems[i];
                    vol += instanciaG.vetItens[item].volume;
                }

                penalizeSolution(rota, bin, *penalty, *stopData);
                penalty->set(SolucaoNS::Penalty::TypePenalty::VolItemsNotPacked, vol);

                return false;
            }

            if(isInfeasible)
                return false;

                   // Comprehensive final feasibility check (axle weights, compactness, support, LIFO)

                   // TODO remover
                   // Penalty penalty;

            const auto errBeforeFinal = vetPackinErros;
            const bool finalOk = bin.checkFeasibility(&rota, false, print, nullptr, stopData);
            if(!finalOk)
            {
                for(int e = 0; e < PackingErroSize; ++e)
                    if(vetPackinErros[e] > errBeforeFinal[e])
                        vetFinalCheckFail[e] += 1;
            }

            if(finalOk)
            {
                /*
                if(numItems != bin.numItens)
                {
                    std::printf("Erro, \nitemIds: %s\n", itemIds.printN(numItems).c_str());
                    std::printf("bin.vetItens: %s", bin.vetItemId.printN(bin.numItens).c_str());

                     PRINT_THROW();
                 }

                  std::set<int> setItems;
                  for(int i = 0; i < numItems; ++i)
                  {
                      if(setItems.contains(itemIds[i]))
                      {
                          std::printf("Repeted item\n");
                          PRINT_THROW();
                      }
                      else
                          setItems.insert(itemIds[i]);
                  }

                   for(int i = 0; i < numItems; ++i)
                   {
                       if(!setItems.contains(bin.vetItemId[i]))
                       {
                           std::printf("Error, item(%d) was not packing\n", bin.vetItens[i]);

                            std::printf("Erro, \nitemIds: %s\n", itemIds.printN(numItems).c_str());
                            std::printf("bin.vetItens: %s", bin.vetItens.printN(bin.numItens).c_str());

                             PRINT_THROW();
                         }
                     }
                     */

                return true;
            }
            else if(!sortItemIds && penalty)
            {

                penalizeSolution(rota, bin, *penalty, *stopData);
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
        if(print)
            std::printf("Error: rota.binPtr is null\n");
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
    int     numItems;

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

}

double getTotalVolume(const SolucaoNS::Bin &bin, const VectorI &itemIds, int numItems)
{
    double totalVol = 0.0;

    for(int i=0; i < numItems; ++i)
        totalVol += instanciaG.vetItens[itemIds[i]].volume;

    for(int i=0; i < bin.numItens; ++i)
        totalVol += instanciaG.vetItens[bin.vetItemId[i]].volume;

    return totalVol;
}

double getTotalLeftArea(const SolucaoNS::Bin &bin, const VectorI &itemIds, int numItems)
{
    double leftArea = 0.0;

    for(int i=0; i < numItems; ++i)
    {
        double a0 = instanciaG.vetItens[itemIds[i]].vetDim[2]*
                    instanciaG.vetItens[itemIds[i]].vetDim[0];

        double a1 = instanciaG.vetItens[itemIds[i]].vetDim[2]*
                    instanciaG.vetItens[itemIds[i]].vetDim[1];

        leftArea += std::max(a0, a1);
    }

    for(int i=0; i < bin.numItens; ++i)
    {
        int item = bin.vetItemId[i];
        double a0 = instanciaG.vetItens[item].vetDim[2]*
                    instanciaG.vetItens[item].vetDim[0];

        double a1 = instanciaG.vetItens[item].vetDim[2]*
                    instanciaG.vetItens[item].vetDim[1];

        leftArea += std::max(a0, a1);
    }

    return leftArea;
}

double getTotalAreaSupp(const SolucaoNS::Bin &bin, const VectorI &itemIds, int numItems)
{
    double area = 0.0;

    for(int i=0; i < numItems; ++i)
    {
        const Array<double, 3> dim = instanciaG.vetItens[itemIds[i]].vetDim;
        area += dim[0]*dim[1];
    }

    for(int i=0; i < bin.numItens; ++i)
    {
        int item = bin.vetItemId[i];
        const Array<double, 3> dim = instanciaG.vetItens[item].vetDim;
        area += dim[0]*dim[1];
    }

    return area;
}

double getTotalWeight(const SolucaoNS::Bin &bin, const VectorI &itemIds, int numItems)
{

    double weight = 0.0;

    for(int i=0; i < numItems; ++i)
        weight += instanciaG.vetItens[itemIds[i]].weight;


    for(int i=0; i < bin.numItens; ++i)
        weight += instanciaG.vetItens[bin.vetItemId[i]].weight;

    return weight;
}

void penalizeSolution(Rota &route, Bin &bin, Penalty &penalty, StopDataNS::StopData &stopData)
{
    if(bin.vazio())
    {
        //std::printf("BIN IS EMPTY!\n");
        return;
    }

    double colisionP = 0.0;
    double supportP = 0.0;
    double axleWightsP = 0.0;
    double loadBalancingP = 0.0;
    double compactnessP = 0.0;
    double lifoP = 0.0;
    double fragilityP = 0.0;

    const Eigen::MatrixXi &matStopItem = stopData.matStopItem;
    const VectorI &sizeMatStopItem = stopData.sizeMatStopItem;
    const VectorI &vetCustStop = stopData.vetCustStop;

    static VectorI vetItemIdex(instanciaG.numItens);
    static Bin localBin;

    vetItemIdex.setAll(-1);

    double fk, fFA, fRA, fTA;

    if(input.modelType >= 1)
    {
        for(int i=0; i < bin.numItens; ++i)
        {
            int item = bin.vetItemId[i];
            vetItemIdex[item] = i;
        }

        /*
        for(int stop=0; stop < (route.numPos-1); ++stop)
        {
            for(int i=0; i < sizeMatStopItem[stop]; ++i)
            {
                int item = matStopItem(stop, i);
                std::printf("%d ", item);
            }

              std::printf("\n");
          }
          */

               //std::string plot = bin.printPlot();

               //std::printf("****************************************\n");
               //std::printf("****************************************\n\n");

               //std::printf("\n\nPlot: \n\n%s", plot.c_str());

               //if(!bin.vazio())
               //    PRINT_THROW();
    }

    bin.sumLeftBalancedLoading = 0.0;
    bin.sumRightBalancedLoading = 0.0;

    static Eigen::Matrix<int, -1, -1, Eigen::RowMajor> matSupportItems(
        instanciaG.numItens, instanciaG.numItens);
    matSupportItems.setConstant(0);

    for(int i = 0; i < bin.numItens; ++i)
    {
        int      itemI = bin.vetItemId[i];
        Rotation rI = bin.vetRotacao[i];
        Ponto   &pI = bin.vetPosItem[i];

        bool tochLeft = tochLeftSideOfTruck(itemI, pI, rI);

        double areaSuport = 0.0;
        double sumAreasLeft = 0.0;

        for(int j=0; j < bin.numItens; ++j)
        {
            if(i == j)
                continue;

            int      itemJ = bin.vetItemId[j];
            Rotation rJ = bin.vetRotacao[j];
            Ponto   &pJ = bin.vetPosItem[j];

            if(input.support && pI.vetDim[2] > 1E-5)
            {
                double zMaxJ = pJ.vetDim[2] + instanciaG.vetItens[itemJ].getDimRotacionada(2, rJ);

                if(doubleEqual(zMaxJ, pI.vetDim[2]))
                {
                    double sup = computeXY_Overlap(
                        instanciaG.vetItens[itemI], rI, pI, instanciaG.vetItens[itemJ], rJ, pJ);

                    if(doubleGreater(sup, 0.0))
                    {
                        areaSuport += sup;
                        matSupportItems(itemI, itemJ) = 1;

                        if(input.fragility && instanciaG.vetItens[itemJ].fragility &&
                           !instanciaG.vetItens[itemI].fragility)
                        {
                            fragilityP += instanciaG.vetItens[itemJ].volume +
                                          instanciaG.vetItens[itemI].volume;
                        }
                    }
                }
            }


            if(!tochLeft && input.compactness)
                sumAreasLeft += getIntercetion(itemI, pI, rI, Left, itemJ, pJ, rJ, Right);
        }

        double area = instanciaG.vetItens[itemI].getDimRotacionada(0, rI) *
                      instanciaG.vetItens[itemI].getDimRotacionada(1, rI);
        double support = areaSuport / area;

        if(input.support && support < input.minSupportArea && pI.vetDim[2] > 1E-5)
            supportP += instanciaG.vetItens[itemI].volume;

        if(input.balancedLoading)
        {
            double w = instanciaG.vetItens[itemI].getDimRotacionada(1, rI);
            double left =
                computeLeftBalancedLoading(pI.vetDim[1], w, instanciaG.vetItens[itemI].weight);
            bin.sumLeftBalancedLoading += left;
            bin.sumRightBalancedLoading += instanciaG.vetItens[itemI].weight - left;
        }

        if(input.compactness && !tochLeft)
        {
            double dz = instanciaG.vetItens[itemI].getDimRotacionada(2, rI);
            double dy = instanciaG.vetItens[itemI].getDimRotacionada(1, rI);

            double areaTotal = dz * dy;
            double ratio = sumAreasLeft / areaTotal;
            if(ratio < input.minLeftSupportArea)
                compactnessP += areaTotal * (input.minLeftSupportArea - ratio);
        }
    }

    if(input.lifo)
    {
        for(int i = 0; i < bin.numItens; ++i)
        {
            int itemI = bin.vetItemId[i];
            int posItemI = findPos(route, bin.vetItemId[i]);

            for(int j = 0; j < bin.numItens; ++j)
            {
                if(i == j)
                    continue;

                int itemJ = bin.vetItemId[j];

                if(instanciaG.vetItens[itemI].customer != instanciaG.vetItens[itemJ].customer)
                {
                    int posItemJ = findPos(route, bin.vetItemId[j]);
                    if(posItemI < posItemJ)
                    {
                        if(!lifo(instanciaG.vetItens[itemI],
                                 bin.vetPosItem[i],
                                 bin.vetRotacao[i],
                                 instanciaG.vetItens[itemJ],
                                 bin.vetPosItem[j],
                                 bin.vetRotacao[j],
                                 input.mlifo,
                                 input.removeFromShortSide,
                                 matSupportItems(itemI, itemJ),
                                 matSupportItems(itemJ, itemI)))
                        {
                            lifoP += instanciaG.vetItens[itemI].volume +
                                     instanciaG.vetItens[itemJ].volume;
                        }
                    }
                }
            }
        }
    }

    if(input.balancedLoading)
    {
        double limit = std::ceil(input.balancedLoadingD * bin.demandaTotal);
        if(bin.sumLeftBalancedLoading > limit)
            loadBalancingP += bin.sumLeftBalancedLoading - limit;

        if(bin.sumRightBalancedLoading > limit)
            loadBalancingP += bin.sumRightBalancedLoading - limit;
    }

    // Stop 1 (all items on board) is already penalized by the global block above.
    // Model 0 only checks the initial load, so it has no per-stop part.
    const int lastStop = (input.modelType >= 1) ? (route.numPos - 2) : 0;
    for(int stop = 2; stop <= lastStop; ++stop)
    {

        double sumLeft 		= 0.0;
        double sumWeight 	= 0.0;
        localBin.reset();
        int next = 0;

        for(int i=0; i < sizeMatStopItem[stop]; ++i)
        {
            int itemI = matStopItem(stop, i);
            int indexI = vetItemIdex[itemI];

            if(indexI < 0)
                continue;

            Rotation rI = bin.vetRotacao[indexI];
            Ponto   &pI = bin.vetPosItem[indexI];

            localBin.vetItemId[next] = itemI;
            localBin.vetPosItem[next] = pI;
            localBin.vetRotacao[next] = rI;
            next += 1;
            localBin.numItens += 1;

            double w = instanciaG.vetItens[itemI].getDimRotacionada(1, rI);

            if(input.balancedLoading)
            {
                sumWeight += instanciaG.vetItens[itemI].weight;
                sumLeft   += computeLeftBalancedLoading(pI.vetDim[1], w,
                                                      instanciaG.vetItens[itemI].weight);
            }

            if(input.compactness)
            {
                bool tochLeft = tochLeftSideOfTruck(itemI, pI, rI);
                if(!tochLeft)
                {
                    double sumAreasLeft = 0.0;

                    for(int j=0; j < sizeMatStopItem[stop]; ++j)
                    {
                        if(i == j)
                            continue;

                        int itemJ = matStopItem(stop, j);
                        int indexJ = vetItemIdex[itemJ];
                        if(indexJ < 0)
                            continue;
                        Rotation rJ = bin.vetRotacao[indexJ];
                        Ponto   &pJ = bin.vetPosItem[indexJ];

                        sumAreasLeft += getIntercetion(itemI, pI, rI, Left, itemJ, pJ, rJ, Right);

                    }

                    double dz = instanciaG.vetItens[itemI].getDimRotacionada(2, rI);
                    double dy = instanciaG.vetItens[itemI].getDimRotacionada(1, rI);

                    double areaTotal = dz * dy;
                    double ratio = sumAreasLeft / areaTotal;
                    if(ratio < input.minLeftSupportArea)
                        compactnessP += areaTotal * (input.minLeftSupportArea - ratio);
                }
            }

        }

        if(input.balancedLoading)
        {
            double sumRight = sumWeight - sumLeft;
            double limit = std::ceil(input.balancedLoadingD * sumWeight);

            if(sumRight > limit)
                loadBalancingP += sumRight - limit;
            if(sumLeft > limit)
                loadBalancingP += sumLeft - limit;
        }

        if(input.axleWights)
        {
            if(!semiTrailer.checkAxleWeights(localBin, false, &fk, &fFA, &fRA, &fTA))
            {
                double overFA = std::max(0.0, fFA - semiTrailer.maxMassFrontAxle * GravityCm);
                double overRA = std::max(0.0, fRA - semiTrailer.maxMassRearAxle * GravityCm);
                double overTA = std::max(0.0, fTA - semiTrailer.maxMassTrailerAxle * GravityCm);

                axleWightsP += overFA+overRA+overTA;
            }
        }

    }

    if(input.axleWights)
    {
        if(!semiTrailer.checkAxleWeights(bin, false, &fk, &fFA, &fRA, &fTA))
        {
            double overFA = std::max(0.0, fFA - semiTrailer.maxMassFrontAxle   * GravityCm);
            double overRA = std::max(0.0, fRA - semiTrailer.maxMassRearAxle    * GravityCm);
            double overTA = std::max(0.0, fTA - semiTrailer.maxMassTrailerAxle * GravityCm);

            axleWightsP += overFA + overRA + overTA;
        }
    }

    penalty.set(Penalty::TypePenalty::Collision, colisionP);
    penalty.set(Penalty::TypePenalty::Support, supportP);
    penalty.set(Penalty::TypePenalty::AxleWeights, axleWightsP);
    penalty.set(Penalty::TypePenalty::LoadBalancing, loadBalancingP);
    penalty.set(Penalty::TypePenalty::Compactness, compactnessP);
    penalty.set(Penalty::TypePenalty::Lifo, lifoP);
    penalty.set(Penalty::TypePenalty::Fragility, fragilityP);
    penalty.set(Penalty::TypePenalty::VolItemsNotPacked, 0.0);
}

       // namespace ConstrutivoBinNS

} // namespace ConstrutivoBin2NS
