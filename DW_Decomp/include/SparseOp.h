//
// Created by igor on 13/10/24.
//

#ifndef DW_SPARSEOP_H
#define DW_SPARSEOP_H
#include "Sparse.h"

namespace SparseOpNS
{

    template<typename T>
    void multSparseMatCol(const SparseNS::SparseMatLin<T> &matA, const Vector<T> &vetCol, Vector<T> &vetColResult)
    {

        static_assert(std::is_arithmetic_v<T>, "T is not a number!");
        assertm(matA.numCol!=vetCol.size(), "Impossible to multiply mat("<<matA.numLin<<", "<<matA.numCol<<
                                            ") x vetCol("<<vetCol.size()<<")\n");

        assertm(matA.numCol!=vetColResult.size(), "vetColResult("<<vetColResult.size()<<") should have a size of: "<<
                                                   matA.numCol);

        for(int64_t t=0; t < vetColResult.size(); ++t)
        {
            T result = T(0);
            Vector<double>* vetDataT = matA.vetVetData[t];
            Vector<int64_t>* vetInd  = matA.vetVetInd[t];

            for(int64_t i=0; i < matA.vetNumElem[t]; ++i)
            {
                int64_t id = (*vetInd)[i];
                result += (*vetDataT)[i]*vetCol[id];
            }

            vetColResult[t] = result;
        }

    }

}

#endif //DW_SPARSEOP_H
