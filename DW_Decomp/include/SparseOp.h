//
// Created by igor on 13/10/24.
//

#ifndef DW_SPARSEOP_H
#define DW_SPARSEOP_H
#include "Sparse.h"

namespace SparseOpNS
{

    template<typename T>
    void multSparseMatLinTimesCol(const SparseNS::SparseMatLin<T> &matA, const Vector<T> &vetCol, Vector<T> &vetColResult)
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


    template<typename T>
    void multLinTimesSparseMatLinFast(const Vector<T> &vetLin, const SparseNS::SparseMatLin<T> &matA, Vector<T> &vetColResult)
    {

        static_assert(std::is_arithmetic_v<T>, "T is not a number!");
        assertm(matA.numCol!=vetLin.size(), "Impossible to multiply mat("<<matA.numLin<<", "<<matA.numCol<<
                                                                         ") x vetCol("<<matA.vetCol.size()<<")\n");

        assertm(matA.numCol!=vetColResult.size(), "vetColResult("<<vetColResult.size()<<") should have a size of: "<<
                                                                 matA.numCol);

        vetColResult.setAll(T(0));

        for(int64_t i=0; i < matA.numCol; ++i)
        {

            Vector<double>* vetDataT = matA.vetVetData[i];
            Vector<int64_t>* vetInd  = matA.vetVetInd[i];

            for(int64_t t =0; t < matA.vetNumElem[i]; ++t)
            {
                int64_t id = (*vetInd)[t];
                vetColResult[id] += (*vetDataT)[t]*vetLin[id];
            }
        }
    }

    template<typename T>
    void multLinTimesSparseMatCol(const SparseNS::SparseVector<T> &vetLin,
                                  const SparseNS::SparseMatCol<T> &matA,
                                  SparseNS::SparseVector<T> &vetColResult)
    {
        vetColResult.resetVector();

        static_assert(std::is_arithmetic_v<T>, "T is not a number!");


        assertm(matA.numCol!=vetLin.n, "Impossible to multiply mat("<<matA.numLin<<", "<<matA.numCol<<
                                                                         ") x vetCol("<<vetLin.n<<")\n");

        assertm(matA.numCol!=vetColResult.n, "vetColResult("<<vetColResult.n<<") should have a size of: "<<
                                                                 matA.numCol);

        for(int64_t j=0; j < matA.numCol; ++j)
        {
            T temp = T(0);
            const Vector<int64_t> &vetIndMatA = *matA.vetVetInd[j];
            const Vector<T> &vetDataColJ_MatA = *matA.vetVetData[j];

            const Vector<int64_t> &vetIndVetLin = vetLin.vetInd;
            const Vector<T> &vetDataVetLin      = vetLin.vetData;

            int64_t iMat = 0;
            int64_t iVet = 0;

            while(iMat < matA.vetNumElem[j] && iVet < vetLin.numElem)
            {
                if(vetIndMatA[iMat] == vetIndVetLin[iVet])
                {
                    temp += vetDataColJ_MatA[iMat] * vetDataVetLin[iVet];
                    iMat += 1;
                    iVet += 1;
                }
                else if(vetIndMatA[iMat] < vetIndVetLin[iVet])
                    iMat += 1;
                else
                    iVet += 1;
            }

            if(temp != 0.0)
                vetColResult(j) = temp;
        }
    }


    template<typename T>
    void multSparseMatColTimesVetCol(const SparseNS::SparseMatCol<T> &matA,
                                  const SparseNS::SparseVector<T> &vetCol,
                                  SparseNS::SparseVector<T> &vetLinResult)
    {
        static_assert(std::is_arithmetic_v<T>, "T is not a number!");

        vetLinResult.resetVector();
        assertm(matA.numCol != vetCol.n, "mat.numCol("<<matA.numCol<<") != vetCol.n("<<vetCol.n<<")");
        assertm(matA.numLin != vetLinResult.n, "mat.numLin("<<matA.numLin<<") != vetLinResult.n("<<vetLinResult.n<<")");

        for(int64_t j=0; j < vetLinResult.n; ++j)
        {
            T colJ = vetCol.getVal(j);
            if(colJ == T(0))
                continue;

            const Vector<T> &vetData = *matA.vetVetData[j];
            const Vector<int64_t> &vetInd = *matA.vetVetInd[j];

            for(int64_t t=0; t < matA.vetNumElem[j]; ++t)
            {
                int64_t id = vetInd[t];
                vetLinResult(id) += vetData[t]*colJ;
            }
        }
    }
}

#endif //DW_SPARSEOP_H
