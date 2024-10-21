//
// Created by igor on 12/10/24.
//

#ifndef DW_SPARSE_H
#define DW_SPARSE_H

#include "safe_vector.h"
#include "Aux.h"
#include <utility>
#include <memory>
#include <ranges>
#include <iostream>

namespace SparseNS
{

    template<typename T>
    int64_t binarySearch(const Vector<T> &vet, T val, int64_t size)
    {
        //std::cout<<"binarySearch\n";
        if(size == 0)
            return -1;

        int64_t ini = 0;
        int64_t mid = size / 2;
        int64_t end = size - 1;

        while(vet[mid] != val)
        {
            if(val > vet[mid])
                ini = mid;
            else
                end = mid;
            int64_t temp = (end - ini + 1) / 2;
            mid = ini + temp;

            if(temp == 0)
                return -1;

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

        SparseMatLin(int64_t numLin_, int64_t numCol_) : numLin(numLin_), numCol(numCol_)
        {
            assertm(numLin == 0, "Number of lines have to be > 0");
            assertm(numCol == 0, "number of columns have to be > 0");

            capLin = numLin;

            vetVetData = Vector<Vector<T> *>(numLin);
            vetVetInd = Vector<Vector<int64_t> *>(numLin);
            vetNumElem = Vector<int64_t>(numLin, 0);

            for(int64_t i = 0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(0);
                vetVetInd[i] = new Vector<int64_t>(0);
            }
        };

        // Warning: adds a new element if there isn't exist; Expensive if it called out of order in j!
        T &operator()(int64_t i, int64_t j)
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
            assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
#endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            if(id >= 0)
            {
                Vector<T> *vetData = vetVetData[i];
                return (*vetData)[id];
            }

            Vector<T> *vetData = vetVetData[i];

            if(vetInd->size() == vetNumElem[i])
            {
                vetInd->push_back(j);

                if constexpr(std::is_arithmetic_v<T>)
                    vetData->push_back(T(0));
                else
                    vetData->push_back();

                vetNumElem[i] += 1;
            } else
            {
                (*vetInd)[vetNumElem[i]] = j;
                if constexpr(std::is_arithmetic_v<T>)
                    (*vetData)[vetNumElem[i]] = T(0);

                vetNumElem[i] += 1;
                //return (*vetData)[vetNumElem[i]];

            }


            if(vetInd->size() > 1)
            {
                int64_t size = vetNumElem[i];

                // Cheeks if vetInd is sorted
                if((*vetInd)[size - 2] > (*vetInd)[size - 1])
                {
                    auto temp = std::views::zip(*vetInd, *vetData);
                    std::ranges::sort(temp);

                    id = binarySearch(*vetInd, j, vetNumElem[i]);
                    assertm(id < 0, "Id(" << j << ") should be find in line(" << i << ")\n");
                    return (*vetData)[id];
                }
            }
            //std::cout<<"ret\n";
            return (*vetData)[vetNumElem[i] - 1];
        }

        T getVal(int64_t i, int64_t j) const
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
            assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
#endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            //std::cout<<"\t\t\tvetInd: "<<*vetInd<<"\n";
            if(id >= 0)
            {
                Vector<T> *vetData = vetVetData[i];
                return (*vetData)[id];
            }

            return T(0);

        }


        void rmElement(int64_t i, int64_t j)
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
            assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
#endif

            Vector<int64_t> *vetInd = vetVetInd[i];
            int64_t id = binarySearch(*vetInd, j, vetNumElem[i]);
            if(id < 0)
                return;

            (*vetInd)[id] = std::numeric_limits<int64_t>::infinity();
            Vector<T> *vetData = vetVetData[i];

            auto temp = std::views::zip(*vetInd, *vetData);
            std::ranges::sort(temp);
            vetNumElem[i] -= 1;

        }

        // Create new lines in the matrix. Warning: EXPENSIVE, utilize vector emplace_back
        void addLines(int num)
        {
            assertm(num <= 0, "num(" << num << ") should be > 0");

            // Cheeks if the new number of lines is greater then capLin
            if((numLin + num) > capLin)
            {
                // For the lines that exist, put zero as the number of elements
                for(int i = 0; i < (capLin - numLin); ++i)
                    vetNumElem[i + numLin] = 0;

                // Create new lines
                for(int i = 0; i < (num - (capLin - numLin)); ++i)
                {
                    vetVetData.emplace_back(new Vector<T>());
                    vetVetInd.emplace_back(new Vector<int64_t>());
                    vetNumElem.emplace_back(0);
                }

                numLin += num;
                capLin = numLin;
            } else
            {
                for(int i = numLin; i < (numLin + num); ++i)
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

            assertm(num <= 0, "num(" << num << ") should be positive and different from 0");
            numCol += num;
        }


        SparseMatLin(const SparseMatLin<T> &mat)
        {
            capLin = mat.capLin;
            numLin = mat.numLin;
            numCol = mat.numCol;

            // Create vectors
            vetVetData = Vector<Vector<T> *>(capLin);
            vetVetInd = Vector<Vector<int64_t> *>(capLin);
            vetNumElem = Vector<int64_t>(capLin, 0);

            // Copy each line
            for(int64_t i = 0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(mat.vetVetData[i]->size());
                vetVetInd[i] = new Vector<int64_t>(mat.vetVetInd[i]->size());
                vetNumElem[i] = mat.vetNumElem[i];

                assertm(mat.vetVetData[i]->size() != mat.vetVetInd[i]->size(), "size(vetData) != size(vetInd)!\n");

                // Copy each element
                for(int64_t t = 0; t < mat.vetNumElem[i]; ++t)
                {
                    (*vetVetData[i])[t] = (*mat.vetVetData[i])[t];
                    (*vetVetInd[i])[t] = (*mat.vetVetInd[i])[t];
                }
            }

        }


        SparseMatLin<T> &operator=(const SparseMatLin<T> &mat)
        {
            //std::cout<<"operator =\n";
            //assertm(1, "NAO IMPLEMENTADO");
            for(int64_t i = 0; i < capLin; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }

            numLin = mat.numLin;
            numCol = mat.numCol;
            capLin = mat.capLin;

            vetVetData = Vector<Vector<T> *>(numLin);
            vetVetInd = Vector<Vector<int64_t> *>(numLin);
            vetNumElem = Vector<int64_t>(mat.vetNumElem);

            for(int64_t i = 0; i < numLin; ++i)
            {
                vetVetData[i] = new Vector<T>(*mat.vetVetData[i]);
                vetVetInd[i] = new Vector<int64_t>(*mat.vetVetInd[i]);
            }

            //std::cout<<"FIM operator =\n";

            return *this;

        }

        ~SparseMatLin()
        {

            for(int64_t i = 0; i < capLin; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }
        }




        int64_t getNumLin() const
        { return numLin; }

        int64_t getCapLin() const
        { return capLin; }

        int64_t getNumCol() const
        { return numCol; }


//    private:
        Vector<Vector<T> *> vetVetData;          // Have a ptr to a vector of T
        Vector<Vector<int64_t> *> vetVetInd;     // Have a ptr to a vector of indice
        Vector<int64_t> vetNumElem;             // Have the number of elements in each line

        int64_t numLin = 0;
        int64_t capLin = 0;                  // capLin >= numLin
        int64_t numCol = 0;


    };


    template<typename T>
    class SparseMatCol
    {
    public:


        SparseMatCol() = default;

        SparseMatCol(int64_t numLin_, int64_t numCol_) : numLin(numLin_), numCol(numCol_)
        {
            assertm(numLin <= 0, "Number of lines have to be > 0");
            assertm(numCol <= 0, "number of columns have to be > 0");

            capCol = numCol;

            vetVetData = Vector<Vector<T> *>(numCol);
            vetVetInd = Vector<Vector<int64_t> *>(numCol);
            vetNumElem = Vector<int64_t>(numCol, 0);

            for(int64_t i = 0; i < numCol; ++i)
            {
                vetVetData[i] = new Vector<T>(0);
                vetVetInd[i] = new Vector<int64_t>(0);
            }

        };


        // Warning: adds a new element if there isn't exist; Expensive if it called out of order in i!
        T &operator()(int64_t i, int64_t j)
        {
            #ifdef VAR_VECTOR_SANITY_CHECK
                assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
                assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
            #endif

            Vector<int64_t> *vetInd = vetVetInd[j];
            int64_t id = binarySearch(*vetInd, i, vetNumElem[j]);

            if(id >= 0)
            {
                Vector<T> *vetData = vetVetData[j];
                return (*vetData)[id];
            }

            Vector<T> *vetData = vetVetData[j];


            if(vetInd->size() == vetNumElem[j])
            {
                vetInd->push_back(i);

                if constexpr(std::is_arithmetic_v<T>)
                    vetData->push_back(T(0));
                else
                    vetData->push_back();

                vetNumElem[j] += 1;
            } else
            {
                (*vetInd)[vetNumElem[j]] = i;
                if constexpr(std::is_arithmetic_v<T>)
                    (*vetData)[vetNumElem[j]] = T(0);

                vetNumElem[j] += 1;
                //return (*vetData)[vetNumElem[j]];

            }


            if(vetInd->size() > 1)
            {
                int64_t size = vetNumElem[j];

                // Cheeks if vetInd is sorted
                if((*vetInd)[size - 2] > (*vetInd)[size - 1])
                {
                    // sort vetInd and vetData by vetInd
                    auto temp = std::views::zip(*vetInd, *vetData);
                    std::ranges::sort(temp);

                    id = binarySearch(*vetInd, i, vetNumElem[j]);
                    assertm(id < 0, "Id(" << i << ") should be find in col(" << j << ")\n");
                    return (*vetData)[id];
                }
            }

            // Ponto do codigo nao eh atingido
            //std::cout<<"ret\n";
            return (*vetData)[vetNumElem[j] - 1];
        }

        T getVal(int64_t i, int64_t j) const
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
            assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
#endif

            Vector<int64_t> *vetInd = vetVetInd[j];
            int64_t id = binarySearch(*vetInd, i, vetNumElem[j]);
            //std::cout<<"\t\t\tvetInd: "<<*vetInd<<"\n";
            if(id >= 0)
            {
                Vector<T> *vetData = vetVetData[j];
                return (*vetData)[id];
            }

            return T(0);

        }


        void rmElement(int64_t i, int64_t j)
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= numLin, "i(" << i << ") is out of range in numLine(" << numLin << ")");
            assertm(j >= numCol, "j(" << j << ") is out of range in numCol(" << numCol << ")");
#endif

            Vector<int64_t> *vetInd = vetVetInd[j];
            int64_t id = binarySearch(*vetInd, i, vetNumElem[j]);
            if(id < 0)
                return;

            (*vetInd)[id] = std::numeric_limits<int64_t>::infinity();
            Vector<T> *vetData = vetVetData[j];

            auto temp = std::views::zip(*vetInd, *vetData);
            std::ranges::sort(temp);
            vetNumElem[j] -= 1;

        }


        void addLine(int num)
        {

            assertm(num <= 0, "num(" << num << ") should be positive and different from 0");
            numLin += num;
        }

        void addColumn(int num)
        {

            assertm(num <= 0, "num(" << num << ") should be > 0");

            // Cheeks if the new number of lines is greater then capLin
            if((numCol + num) > capCol)
            {
                // For the lines that exist, put zero as the number of elements
                for(int i = 0; i < (capCol - numCol); ++i)
                    vetNumElem[i + numCol] = 0;

                // Create new lines
                for(int i = 0; i < (num - (capCol - numCol)); ++i)
                {
                    vetVetData.emplace_back(new Vector<T>());
                    vetVetInd.emplace_back(new Vector<int64_t>());
                    vetNumElem.emplace_back(0);
                }

                numCol += num;
                capCol = numCol;
            } else
            {
                for(int i = numCol; i < (numCol + num); ++i)
                    vetNumElem[i] = 0;

                numCol += num;
            }
        }

        ~SparseMatCol()
        {

            for(int64_t i = 0; i < capCol; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }
        }


        SparseMatCol<T>& operator=(const SparseMatCol<T> &mat)
        {

            for(int64_t i = 0; i < capCol; ++i)
            {
                delete vetVetData[i];
                delete vetVetInd[i];
            }

            numLin = mat.numLin;
            numCol = mat.numCol;
            capCol = mat.capCol;

            vetNumElem = Vector<int64_t>(capCol);

            vetVetData = Vector<Vector<T>*>(capCol);
            vetVetInd = Vector<Vector<int64_t>*>(capCol);

            for(int64_t i=0; i < capCol; ++i)
            {

                vetVetData[i] = new Vector<T>(*mat.vetVetData[i]);
                vetVetInd[i]  = new Vector<int64_t>(*mat.vetVetInd[i]);
            }


            return *this;

        }




        int64_t getNumLin() const
        { return numLin; }

        int64_t getCapCol() const
        { return capCol; }

        int64_t getNumCol() const
        { return numCol; }


//    private:
        Vector<Vector<T> *> vetVetData;          // Have a ptr to a vector of T
        Vector<Vector<int64_t> *> vetVetInd;     // Have a ptr to a vector of index
        Vector<int64_t> vetNumElem;             // Have the number of elements in each line

        int64_t numLin = 0;
        int64_t capCol = 0;                  // capLin >= numCol
        int64_t numCol = 0;

    };

    template<typename T>
    class SparseVector
    {
    public:
        SparseVector() = default;

        SparseVector(int64_t n_)
        {
            n = n_;
            //cap = n;

            assertm(n <= 0, "Number of elements should be > 0");
        }

        // Warning: adds a new element if there isn't exist; Expensive if it called out of order in i!
        T &operator()(int64_t i)
        {
            //std::cout<<"i("<<i<<")\n";

#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= n, "i(" << i << ") is out of range in numElem(" << n << ")");
#endif
            int64_t id = binarySearch(vetInd, i, numElem);
            if(id >= 0)
            {
                return vetData[id];
            }

            if(vetInd.size() == numElem)
            {
                vetInd.push_back(i);

                if constexpr(std::is_arithmetic_v<T>)
                    vetData.push_back(T(0));
                else
                    vetData.push_back();

            } else
            {
                vetInd[n] = i;
                if constexpr(std::is_arithmetic_v<T>)
                    vetData[n] = T(0);

                //return (*vetData)[vetNumElem[j]];
            }


            numElem += 1;


            if(vetInd.size() > 1)
            {
                int64_t size = numElem;

                // Cheeks if vetInd is sorted
                if(vetInd[size - 2] > vetInd[size - 1])
                {
                    // sort vetInd and vetData by vetInd
                    auto temp = std::views::zip(vetInd, vetData);
                    std::ranges::sort(temp);

                    id = binarySearch(vetInd, i, numElem);
                    assertm(id < 0, "Id(" << i << ") should be find in vet\n");
                    return vetData[id];
                }
            }

            // Ponto do codigo nao eh atingido
            //std::cout<<"ret\n";
            return vetData[numElem - 1];
        }


        T getVal(int64_t i) const
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= n, "i(" << i << ") is out of range in numElem(" << n << ")");
#endif

            int64_t id = binarySearch(vetInd, i, numElem);
            //std::cout<<"\t\t\tvetInd: "<<*vetInd<<"\n";
            if(id >= 0)
            {
                return vetData[id];
            }

            return T(0);

        }


        void rmElement(int64_t i)
        {
#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(i >= n, "i(" << i << ") is out of range in numElem(" << n << ")");
#endif

            int64_t id = binarySearch(vetInd, i, numElem);
            if(id < 0)
                return;

            vetInd[id] = std::numeric_limits<int64_t>::infinity();
            auto temp = std::views::zip(vetInd, vetData);
            std::ranges::sort(temp);

            numElem -= 1;

        }

        void rmLastN_Elements(int64_t n_)
        {

#ifdef VAR_VECTOR_SANITY_CHECK
            assertm(n_ <= 0, "n(" << n << ") should be positive");
            assertm(n_ > numElem, "n(" << n << ") should be less then numElem(" << n << ")");
#endif

            for(int64_t i = (n_ - numElem); i < numElem; ++i)
                vetInd[i] = std::numeric_limits<int64_t>::infinity();

            numElem = numElem - n_;
        }

        void resetVector()
        {
            numElem = 0;
        }


        Vector<T> vetData;
        Vector<int64_t> vetInd;

        int64_t n = 0;
        int64_t numElem = 0;
        //int64_t cap     = 0;

    };

    template<typename T>
    std::ostream &operator<<(std::ostream &os, const SparseMatLin<T> &sparseMat)
    {
        for(int64_t i = 0; i < sparseMat.getNumLin(); ++i)
        {
            for(int64_t j = 0; j < sparseMat.getNumCol(); ++j)
            {
                os << sparseMat.getVal(i, j) << " ";
            }

            os << "\n";
        }

        return os;
    }


    template<typename T>
    std::ostream &operator<<(std::ostream &os, const SparseMatCol<T> &sparseMat)
    {
        for(int64_t i = 0; i < sparseMat.getNumLin(); ++i)
        {
            for(int64_t j = 0; j < sparseMat.getNumCol(); ++j)
            {
                os << sparseMat.getVal(i, j) << " ";
            }

            os << "\n";
        }

        return os;
    }


    template<typename T>
    std::ostream &operator<<(std::ostream &os, const SparseVector<T> &sparseVet)
    {
        for(int64_t i = 0; i < sparseVet.n; ++i)
        {
            os << sparseVet.getVal(i) << " ";
        }

        os << "\n";

        return os;
    }

    typedef SparseVector<double> SparseVectorD;
    typedef SparseMatLin<double> SparseMatLinD;
    typedef SparseMatCol<double> SparseMatColD;
}

#endif //DW_SPARSE_H
