/*
    Modificado por Igor de Andrade Junqueira,
    comentario original:

    Safe C++, Or How to Avoid Most Common Mistakes in C++ Code by Vladimir Kushnir, (OÕReilly).
    Copyright 2012 Vladimir Kushnir, ISBN 9781449320935.
    If you feel your use of code examples falls outside fair use or the
    permission given above, feel free to contact us at permissions@oreilly.com.
*/

#ifndef TTP_SAFE_MATRIX_H
#define TTP_SAFE_MATRIX_H

#include "iostream"
#include <vector>
#include <ostream>
#include <omp.h>
#include <cmath>
#include <cassert>


/*
class M_Coluna{M_Coluna()=default;};
class M_Linha{M_Linha()=default;};
class M_TrianngSup{M_TrianngSup()=default;};
class M_TrianngInf{M_TrianngInf()=default;};

*/

// Wrapper around std::vector, has temporary sanity checks in the operators [].
template <typename T, bool C=true>
class Matrix : public std::vector<T>
{
private:
    size_t n=0;
    size_t m=0;

    T get(size_t i)const{return std::vector<T>::operator[](i);}

public:

    explicit Matrix(size_t n_, size_t m_): std::vector<T>(n_*m_), n(n_), m(m_){}

    Matrix(size_t n_, size_t m_, const T& value): std::vector<T>(n_*m_, value), n(n_), m(m_){}

    Matrix(size_t n): std::vector<T>(n){}
    Matrix(){}

    Matrix(const Matrix<T> &temp):n(temp.n), m(temp.n), std::vector<T>(temp.n*temp.m)
    {
        for(int i=0; i < n*m; ++i)
            std::vector<T>::operator[](i) = temp.get(i);
    }


//    template <class InputIterator> Vector (InputIterator first, InputIterator last): std::vector<T>(first, last){}

    // Note: we do not provide a copy-ctor and assignment operator.
    // we rely on default versions of these methods generated by the compiler.

    T& operator[](size_t index) = delete;
    const T& operator[](size_t index) const = delete;

    // inline __attribute__((always_inline))
    virtual T operator()(size_t indexI, size_t indexJ) const
    {

        size_t index = indexI*m+indexJ;

        if constexpr(C)
        {

#if VAR_VECTOR_SANITY_CHECK
            if(indexI >= n)
            {
                std::cout << "Erro indice i: " << indexI << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }


            if(indexJ >= m)
            {
                std::cout << "Erro indice j: " << indexJ << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }

            if(index >= m * n)
            {
                std::cout << "Erro indice(" << index << ") >= " << m * n << "\n";
                std::cout << "i=" << indexI << "; j=" << indexJ << "\n";
                std::cout << "n=" << n << "\n";

                throw std::out_of_range("");
            }
#endif
        }

        return std::vector<T>::operator[](index);
    }

    virtual inline __attribute__((always_inline)) T& get(size_t indexI, size_t indexJ)
    {

        size_t index = indexI*m+indexJ;

        if constexpr(C)
        {
#if VAR_VECTOR_SANITY_CHECK
            if(indexI >= n)
            {
                std::cout << "Erro indice i: " << indexI << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }


            if(indexJ >= m)
            {
                std::cout << "Erro indice j: " << indexJ << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }


            if(index >= m * n)
            {
                std::cout << "Erro indice(" << index << ") >= " << m * n << "\n";
                std::cout << "i=" << indexI << "; j=" << indexJ << "\n";
                std::cout << "n=" << n << "\n";

                throw std::out_of_range("");
            }
#endif
        }

        return std::vector<T>::operator[](index);
    }

    [[nodiscard]] virtual size_t inline __attribute__((always_inline)) getNumLinhas() const
    {
        return n;
    }


    [[nodiscard]] virtual size_t inline __attribute__((always_inline)) getNumColunas() const
    {
        return m;
    }

    virtual void setVal(const T &val)
    {
        for(size_t i=0; i < m*n; ++i)
            std::vector<T>::operator[](i) = val;
    }
    virtual void setVal(const T &&val)
    {
        for(size_t i=0; i < m*n; ++i)
            std::vector<T>::operator[](i) = val;
    }

    inline __attribute__((always_inline))  auto getIterator(size_t indexI, size_t indexJ)
    {

        if constexpr(C)
        {
#if VAR_VECTOR_SANITY_CHECK
            if(indexI >= n)
            {
                std::cout << "Erro indice i: " << indexI << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }


            if(indexJ >= m)
            {
                std::cout << "Erro indice j: " << indexJ << " esta errado para matrix de tam " << n << " x " << m
                          << "\n";
                throw std::out_of_range("");
            }
#endif
        }
        return (std::vector<T>::begin()+(indexI*m+indexJ));
    }

    virtual void printVector()
    {
        for(size_t i=0; i < m*n; ++i)
            std::cout<<std::vector<T>::operator[](i)<<" ";
        std::cout<<"\n";
    }



    //Matrix(size_t n, const T &val): std::vector<T>(n, val){}




};


template <typename T, bool C>
inline __attribute__((always_inline)) std::ostream& operator << (std::ostream& os, Matrix<T,C>& mat)
{
    for(size_t i=0; i < mat.getNumLinhas(); ++i)
    {
        for(size_t j=0; j < mat.getNumColunas(); ++j)
        {
            const T val = mat(i, j);
            os<<val<<"\t";

        }

        os<<"\n";
    }

    return os;
}


template <typename T, bool C>
bool matrixsIguais(const Matrix<T, C> &mat0, const Matrix<T,C> &mat1)
{
    const T epsilon = 1E-5;

    if(mat0.getNumLinhas()  != mat1.getNumLinhas() || mat0.getNumColunas() != mat1.getNumColunas())
    {
        return false;
    }

    for(size_t i = 0; i < mat0.getNumLinhas(); ++i)
    {
            for(size_t j = 0; j < mat1.getNumColunas(); ++j)
            {
                if constexpr(std::is_same_v<T,double>)
                {
                    if(fabs(mat0(i,j)-mat1(i,j)) > epsilon)
                        return false;

                }
                else if constexpr(std::is_same_v<T,float>)
                {
                    if(fabsf(mat0(i,j)-mat1(i,j)) > epsilon)
                        return false;
                }
                else
                {
                    if(mat0(i, j) != mat1(i, j))
                    {
                        std::cout << mat0(i, j) << " != " << mat1(i, j) << "\n";
                        return false;
                    }
                }
            }
    }


    return true;
}


template <typename T, bool C>
void multMatrix(const Matrix<T,C> &mat0, const Matrix<T,C> &mat1, Matrix<T,C> &matResult)
{

    if(mat0.getNumLinhas() != mat1.getNumColunas() ||
       mat0.getNumColunas() != mat1.getNumLinhas())
    {
        std::cout<<"Erro, num linhas("<<mat0.getNumLinhas()<<") da mat0 != "<<
                                                             "num col("<<mat1.getNumColunas()<<") da mat1\n";
        throw "ERRO";
    }

    if(matResult.getNumLinhas() != mat0.getNumLinhas() || matResult.getNumColunas() != mat1.getNumColunas())
    {
#if VAR_VECTOR_SANITY_CHECK
        matResult = Matrix<T>(mat0.getNumLinhas(), mat1.getNumColunas(), 0.0);
#else
        matResult = Matrix<T>(mat0.getNumLinhas(), mat1.getNumColunas());
#endif
    }


     const size_t numLinhas = matResult.getNumLinhas();
     const size_t numCol    = matResult.getNumColunas();
     const size_t mat0NumCol = mat0.getNumColunas();


    {
        size_t i;
        #pragma omp parallel for simd firstprivate(numLinhas,numCol, mat0NumCol, mat0, mat1) shared(matResult) private(i)
        for(i = 0; i < numLinhas; ++i)
        {

            for(size_t j = 0; j < numCol; ++j)
            {
                T sum = T(0);


                for(size_t k = 0; k < mat0NumCol; ++k)
                {
                    sum += mat0(i, k) * mat1(k, j);
                    //std::cout<<"A("<<i<<","<<k<<") * B("<<k<<","<<j<<"\n";
                }

                //std::cout<<"Sum: "<<sum<<"\n";
                //exit(-1);

                //#pragma omp critical
                {
                    matResult(i, j) = sum;

                }
            }
        }

    }
}


#endif //TTP_SAFE_VECTOR_H
