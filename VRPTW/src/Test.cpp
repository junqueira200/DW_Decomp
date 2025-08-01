/*  ********************************************************************************************************************
 *  ********************************************************************************************************************
 *  File:    Test.cpp
 *  Project: DW_Decomp
 *  Author:  Igor de Andrade Junqueira
 *  Date:    06/06/25
 *
 *  ********************************************************************************************************************
 *  ********************************************************************************************************************/

#include "Test.h"
#include <unordered_set>

using namespace std;
using namespace TestNS;
using namespace InstanceVRPTW_NS;

// Helper function to find all combinations
// of size r in an array of size n
void TestNS::combinationUtil(int ind, int r, VectorI &data, Vector<VectorI> &result, VectorI &arr)
{
    int n = arr.size();

    // If size of current combination is r
    if ((int)data.size() == r) {
        result.push_back(data);
        return;
    }

    // If no more elements are left to put in data
    if (ind >= n)
        return;

    // include the current element
    data.push_back(arr[ind]);

    // Recur for next elements
    combinationUtil(ind + 1, r, data, result, arr);

    // Backtrack to find other combinations
    data.pop_back();

    // exclude the current element and
    // move to the next unique element
    while(ind + 1 < n && arr[ind] == arr[ind + 1])
        ind++;

    combinationUtil(ind + 1, r, data, result, arr);
}

// Function to find all combinations of size r
// in an array of size n
Vector<VectorI> TestNS::findCombination(VectorI &arr, int r)
{
    // to store the result
    Vector<VectorI> result;

    // sort the array
    sort(arr.begin(), arr.end());

    // Temporary array to store current combination
    VectorI data;

    combinationUtil(0, r, data, result, arr);
    return result;
}

void TestNS::enumerateRoutes(InstanceVRPTW_NS::InstanceVRPTW& instVrp, int numMax, RouteHash& routeHash)
{
    VectorI vet;
    Vector<VectorI> vetResult;
    routeHash.reserve(1000);

    vet.reserve(instVrp.numClientes-2);

    for(int i=1; i < (instVrp.numClientes-1); ++i)
        vet.push_back(i);

    numMax = std::min(numMax, (int)vet.size());

    for(int i=1; i <= numMax; ++i)
    {
        Vector<VectorI> vetAux = findCombination(vet, i);
        std::move(vetAux.begin(), vetAux.end(), back_inserter(vetResult));
    }

    int numInsert = 0;
    for(VectorI& vetIntRoute:vetResult)
    {
        if(vetIntRoute.size() <= 0)
            continue;

        int demand = computeDemand(vetIntRoute);
        if(demand > instVrp.capVeic)
            continue;

        Route route(vetIntRoute.size()+2);

        for(int i=0; i < (int)vetIntRoute.size(); ++i)
            route.vetRoute[i+1] = vetIntRoute[i];

        computeDistance(route);
        computeHash(route);
        route.demand = computeDemand(route.vetRoute);
        routeHash.insert(std::move(route));
        numInsert += 1;
    }

    if(numInsert != (int)routeHash.size())
    {
        std::cout<<"ERROR\nnumInser("<<numInsert<<"); size("<<routeHash.size()<<")\n";
    }
    else
        std::cout<<"All routes("<<numInsert<<") were inserted in the hash\n";
}

void TestNS::computeDistance(Route& route)
{
    double dist = 0.0;
    for(int i=0; i < ((int)route.vetRoute.size()-1); ++i)
        dist += ptr_instVrpG->matDist(route.vetRoute[i], route.vetRoute[i+1]);

    route.dist = dist;
}

void TestNS::computeHash(Route& route)
{
    // https://cseweb.ucsd.edu/~kube/cls/100/Lectures/lec16/lec16-16.html
    // ELF Hash algorithms
    route.valHash = 0;
    for(int i=0; i < (int)route.vetRoute.size(); ++i)
    {   //               valHash * 16
        route.valHash = (route.valHash<<4) + route.vetRoute[i];
        uint64_t g = route.valHash & 0xF0000000L;

        if(g != 0)
            route.valHash ^= g >> 24;
        route.valHash &= ~g;
    }
}

bool TestNS::Route::operator == (const Route& route) const
{

    if(valHash != route.valHash)
        return false;

    if(vetRoute.size() != route.vetRoute.size())
        return false;

    for(int i=0; i < (int)route.vetRoute.size(); ++i)
    {
        if(vetRoute[i] != route.vetRoute[i])
            return false;
    }

    return true;
}

FloatType TestNS::computeReducedCost(const Route& route, const Eigen::RowVectorXd& vetRowPi)
{
    double redCost = -vetRowPi[0] + route.dist;

    for(int i=1; i < (int)route.vetRoute.size(); ++i)
        redCost += -vetRowPi[route.vetRoute[i]+1];

    return redCost;
}

int TestNS::computeDemand(VectorI& route)
{
    int demand = 0;
    for(int i=0; i < (int)route.size(); ++i)
        demand += ptr_instVrpG->vetClieDem[route[i]];

    return demand;
}
