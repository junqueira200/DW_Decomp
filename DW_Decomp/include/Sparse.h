//
// Created by igor on 12/10/24.
//

#ifndef DW_SPARSE_H
#define DW_SPARSE_H

#include "safe_vector.h"
#include "Aux.h"
#include <ranges>
#include <iostream>

namespace SparseNS
{

    template<typename T>
    int64_t binarySearch(Vector<T> &vet, T val, size_t size)
    {
        if(vet.empty())
            return -1;

        int64_t ini = 0;
        int64_t mid = size/2;
        int64_t end = size-1;

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

            capLin = numLin;

            vetVetData = Vector<Vector<T>*>(numLin);
            vetVetInd  = Vector<Vector<int64_t>*>(numLin);
            vetNumElem = Vector<int64_t>(numLin, 0);

            for(int64_t i=0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(0);
                vetVetInd[i]  = new Vector<int64_t>(0);
            }
        };

        // Warning: adds a new element if there isn't exist; Expensive if it called out of order in j!
        T& operator()(int64_t i, int64_t j)
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i>=numLin, "i("<<i<<") is out of range in numLine("<<numLin<<")");
                assertm(j>=numCol, "j("<<j<<") is out of range in numCol("<<numCol<<")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            if(id >= 0)
            {
                Vector<T>* vetData = vetVetData[i];
                return (*vetData)[id];
            }

            Vector<T>* vetData = vetVetData[i];

            if(vetInd->size() == vetNumElem[i])
            {
                vetInd->push_back(j);

                if constexpr(std::is_arithmetic_v<T>)
                    vetData->push_back(T(0));
                else
                    vetData->push_back();

                vetNumElem[i] += 1;
            }
            else
            {
                (*vetInd)[vetNumElem[i]] = j;
                if constexpr(std::is_arithmetic_v<T>)
                    (*vetData)[vetNumElem[i]] = T(0);

                vetNumElem[i] += 1;
                return (*vetData)[vetNumElem[i]];

            }


            if(vetInd->size() > 1)
            {
                int64_t size = vetNumElem[i];

                // Cheeks if vetInd is sorted
                if((*vetInd)[size-2] > (*vetInd)[size-1])
                {
                    auto temp = std::views::zip(*vetInd, *vetData);
                    std::ranges::sort(temp);

                    id = binarySearch(*vetInd, j, vetNumElem[i]);
                    assertm(id < 0, "Id("<<j<<") should be find in line("<<i<<")\n");
                    return (*vetData)[id];
                }
            }
            //std::cout<<"ret\n";
            return (*vetData)[vetNumElem[i]-1];
        }

        T getVal(int64_t i, int64_t j) const
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i>=numLin, "i("<<i<<") is out of range in numLine("<<numLin<<")");
                assertm(j>=numCol, "j("<<j<<") is out of range in numCol("<<numCol<<")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            //std::cout<<"\t\t\tvetInd: "<<*vetInd<<"\n";
            if(id >= 0)
            {
                Vector<T>* vetData = vetVetData[i];
                return (*vetData)[id];
            }

            return T(0);

        }


        void rmElement(int64_t i, int64_t j)
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i>=numLin, "i("<<i<<") is out of range in numLine("<<numLin<<")");
                assertm(j>=numCol, "j("<<j<<") is out of range in numCol("<<numCol<<")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            if(id < 0)
                return;

            (*vetInd)[id] = std::numeric_limits<int64_t>::infinity();
            Vector<T>* vetData = vetVetData[i];

            auto temp = std::views::zip(*vetInd, *vetData);
            std::ranges::sort(temp);
            vetNumElem[i] -= 1;

        }

        // Create new lines in the matrix. Warning: EXPENSIVE, utilize vector emplace_back
        void addLines(int num)
        {
            assertm(num<=0, "num("<<num<<") should be > 0");

            // Cheeks if the new number of lines is greater then capLin
            if((numLin+num) > capLin)
            {
                // For the lines that exist, put zero as the number of elements
                for(int i=0; i < (capLin-numLin); ++i)
                    vetNumElem[i+numLin] = 0;

                // Create new lines
                for(int i = 0; i < (num-(capLin-numLin)); ++i)
                {
                    vetVetData.emplace_back(new Vector<T>());
                    vetVetInd.emplace_back(new Vector<int64_t>());
                    vetNumElem.emplace_back(0);
                }

                numLin += num;
                capLin = numLin;
            }
            else
            {
                for(int i=numLin; i < (numLin+num); ++i)
                    vetNumElem[i] = 0;

                numLin += num;
            }
        }

        void rmLastLine()
        {
            numLin -= 1;
            if(numLin < 0)
                numLin = 0;
        }

        void addColums(int num)
        {

            assertm(num<=0, "num("<<num<<") should be positive and different from 0");
            numCol += num;
        }


        SparseMatLin(const SparseMatLin<T> &mat)
        {
            capLin = mat.capLin;
            numLin = mat.numLin;
            numCol = mat.numCol;

            // Create vectors
            vetVetData = Vector<Vector<T>*>(capLin);
            vetVetInd  = Vector<Vector<int64_t>*>(capLin);
            vetNumElem = Vector<int64_t>(capLin, 0);

            // Copy each line
            for(int64_t i=0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(mat.vetVetData[i]->size());
                vetVetInd[i]  = new Vector<int64_t>(mat.vetVetInd[i]->size());
                vetNumElem[i] = mat.vetNumElem[i];

                assertm(mat.vetVetData[i]->size()!=mat.vetVetInd[i]->size(), "size(vetData) != size(vetInd)!\n");

                // Copy each element
                for(int64_t t=0; t < mat.vetNumElem[i]; ++t)
                {
                    (*vetVetData[i])[t] = (*mat.vetVetData[i])[t];
                    (*vetVetInd[i])[t]  = (*mat.vetVetInd[i])[t];
                }
            }

        }

        // TODO
        SparseMatLin<T>& operator=(const SparseMatLin<T> &mat)
        {
            //std::cout<<"operator =\n";
            //assertm(1, "NAO IMPLEMENTADO");
            for(int64_t i=0; i < capLin; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }

            numLin = mat.numLin;
            numCol = mat.numCol;
            capLin = mat.capLin;

            vetVetData = Vector<Vector<T>*>(numLin);
            vetVetInd  = Vector<Vector<int64_t>*>(numLin);
            vetNumElem = Vector<int64_t>(mat.vetNumElem);

            for(int64_t i=0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(*mat.vetVetData[i]);
                vetVetInd[i]  = new Vector<int64_t>(*mat.vetVetInd[i]);
            }

            //std::cout<<"FIM operator =\n";

            return *this;

        }

        ~SparseMatLin()
        {

            for(int64_t i=0; i < capLin; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }
        }

        int64_t getNumLin()const{return numLin;}
        int64_t getCapLin()const{return capLin;}
        int64_t getNumCol()const{return numCol;}


//    private:
        Vector<Vector<T>*> vetVetData;          // Have a ptr to a vector of T
        Vector<Vector<int64_t>*> vetVetInd;     // Have a ptr to a vector of indice
        Vector<int64_t> vetNumElem;             // Have the number of elements in each line

        int64_t numLin    = 0;
        int64_t capLin    = 0;                  // capLin >= numLin
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


#endif //DW_SPARSE_H
