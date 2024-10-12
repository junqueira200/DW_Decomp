//
// Created by igor on 12/10/24.
//

#ifndef DW_SPARSEMAT_H
#define DW_SPARSEMAT_H

#include "safe_vector.h"
#include "Aux.h"
#include <ranges>
#include <iostream>

namespace SparseMatNS
{

    template<typename T>
    int64_t binarySearch(Vector<T> &vet, T val)
    {
        if(vet.empty())
            return -1;

        int64_t ini = 0;
        int64_t mid = vet.size()/2;
        int64_t end = vet.size()-1;

        while(vet[mid] != val)
        {
            if(val > vet[mid])
                ini = mid;
            else
                end = mid;
            int64_t temp = (end-ini+1)/2;
            mid = ini + temp;

            if(ini == mid || mid == end)
            {
                if(vet[mid] == val)
                    return mid;
                else if(vet[ini] == val)
                    return ini;
                else if(vet[end] == val)
                    return end;
                else
                    return -1;

            }

            assertm(temp == 0, "temp nao deveria ser 0!");
        }

        if(vet[mid] == val)
            return mid;
        else
            return -1;

    }

    template<typename T>
    class SparseMatLin
    {
    public:

        SparseMatLin() = default;
        SparseMatLin(int64_t numLin_, int64_t numCol_):numLin(numLin_), numCol(numCol_)
        {
            assertm(numLin==0, "Number of lines have to be > 0");
            assertm(numCol==0, "number of columns have to be > 0");

            vetVetData = Vector<Vector<T>*>(numLin);
            vetVetInd  = Vector<Vector<int64_t>*>(numLin);
            vetNumElem = Vector<int64_t>(numLin, 0);

            for(int64_t i=0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(0);
                vetVetInd[i]  = new Vector<int64_t>(0);
            }
        };

        // Warning: adds a new element if elem (i,j) there isn't exist
        T& operator()(int64_t i, int64_t j)
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i>=numLin, "i("<<i<<") is out of range in numLine("<<numLin<<")");
                assertm(j>=numCol, "j("<<j<<") is out of range in numCol("<<numCol<<")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j);
            if(id >= 0)
            {
                Vector<T>* vetData = vetVetData[i];
                return (*vetData)[(*vetInd)[id]];
            }

            vetInd->push_back(j);
            Vector<T>* vetData = vetVetData[i];

            vetNumElem[i] += 1;

            if constexpr(std::is_arithmetic_v<T>)
                vetData->push_back(T(0));
            else
                vetData->push_back();

            if(vetInd->size() > 1)
            {
                int64_t size = vetInd->size();
                if((*vetInd)[size-2] > (*vetInd)[size-1])
                {
                    auto temp = std::views::zip(*vetInd, *vetData);
                    std::ranges::sort(temp);

                    id = binarySearch(*vetInd, j);
                    assertm(id < 0, "Id("<<j<<") should be find in line("<<i<<")\n");
                    return (*vetData)[id];

                }
            }

            return (*vetData)[vetData->size()-1];
        }

        T getVal(int64_t i, int64_t j) const
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i>=numLin, "i("<<i<<") is out of range in numLine("<<numLin<<")");
                assertm(j>=numCol, "j("<<j<<") is out of range in numCol("<<numCol<<")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j);
            if(id >= 0)
            {
                Vector<T>* vetData = vetVetData[i];
                return (*vetData)[(*vetInd)[id]];
            }

            return T(0);

        }

        // Create new lines in the matrix. Warning: EXPENSIVE, utilize vector emplace_back
        void addLines(int num)
        {
            assertm(num<=0, "num("<<num<<") should be > 0");

            for(int i=0; i < num; ++i)
            {
                vetVetData.emplace_back(new Vector<T>());
                vetVetInd.emplace_back(new Vector<int64_t>());
                vetNumElem.emplace_back(0);
            }

            numLin += num;
        }

        void addCollums(int num)
        {

            assertm(num<=0, "num("<<num<<") should be > 0");
            numCol += num;
        }

        // TODO
        SparseMatLin(const SparseMatLin<T> &mat)
        {
            assertm(1, "NAO IMPLEMENTADO");
        }

        // TODO
        SparseMatLin<T>& operator=(const SparseMatLin<T> &mat)
        {
            assertm(1, "NAO IMPLEMENTADO");
        }

        ~SparseMatLin()
        {

            for(int64_t i=0; i < numLin; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }
        }

        int64_t getNumLin()const{return numLin;}
        int64_t getNumCol()const{return numCol;}


    private:
        Vector<Vector<T>*> vetVetData;          // Have a ptr to a vector of T
        Vector<Vector<int64_t>*> vetVetInd;     // Have a ptr to a vector of indice
        Vector<int64_t> vetNumElem;             // Have the number of elements in each line

        int64_t numLin    = 0;
        int64_t numCol    = 0;


    };


    template<typename T>
    class SparseMatCol
    {
    public:

    };

    template<typename T>
    std::ostream& operator<<(std::ostream& os, const SparseMatLin<T>& sparseMat)
    {
        for(int64_t i=0; i < sparseMat.getNumLin(); ++i)
        {
            for(int64_t j=0; j < sparseMat.getNumCol(); ++j)
            {
                os<<sparseMat.getVal(i,j)<<" ";
            }

            os<<"\n";
        }

        return os;
    }



}


#endif //DW_SPARSEMAT_H
