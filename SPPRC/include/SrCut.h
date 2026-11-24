/*  *****************************************************************
 *  *****************************************************************
 *  File:    SrCut.h
 *  Project: SPPCR
 *  Author:  Igor de Andrade Junqueira
 *  Date:    23/11/26
 *
 *  *****************************************************************
 *  *****************************************************************/

#ifndef SRCUT_H
#define SRCUT_H

#include "LabelingConstants.h"   // NumMaxCust
#include "sefe_array.h"
#include <bitset>

namespace LabelingAlgorithmNS
{

    struct Rank1Cut
    {
        std::bitset<NumMaxCust> S;        // customers of the cut (|S| = 3, k = 2)
        Array<int, 5> vetS;
        int sizeVetS;
        std::bitset<NumMaxCust> memory;   // limited memory (later); for now = all customers
        double penalty;                   // = -sigma >= 0, already with the sign fixed
        int    colunaId;                  // which Coluna cut it is (to match the dual each call)
    };


}

#endif // SRCUT_H
