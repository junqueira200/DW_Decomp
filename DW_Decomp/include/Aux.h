//
// Created by igor on 19/09/24.
//

#ifndef DW_DECOMP_AUX_H
#define DW_DECOMP_AUX_H
#include <Eigen/Eigen>

#define __PRETTYFILE__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define PRINT_DEBUG(inicio, texto) std::cout<<inicio<<"DEBUG: "<<texto<<"  FILE: "<<__PRETTYFILE__<<"  FUNC: "<<__PRETTY_FUNCTION__<<"  LINHA: "<<__LINE__<<"\n";

#define assertm(exp, msg) if(exp){std::cout<<msg<<"\n\nLINE: "<<__LINE__<<"\nFILE: "<<__PRETTYFILE__<<"\nFUNC: "<<__PRETTY_FUNCTION__<<"\n\n"; throw "ERROR";}


template<typename T, int option>
class TempSpMatPrint
{
public:

    Eigen::SparseMatrix<T,option> &m;
    TempSpMatPrint(Eigen::SparseMatrix<T,option> &m_):m(m_)
    {}

};


template<typename T, int option>
class TempSpVetPrint
{
public:

    Eigen::SparseVector<T,option> &m;
    TempSpVetPrint(Eigen::SparseVector<T,option> &m_):m(m_)
    {}

};


template<typename T, int option>
std::ostream & operator << (std::ostream & s, const TempSpMatPrint<T,option> &temp)
{
    s<<static_cast<const Eigen::SparseMatrixBase<Eigen::SparseMatrix<T,option>>&>(temp.m);
    return s;
}


template<typename T, int option>
std::ostream & operator << (std::ostream & s, const TempSpVetPrint<T,option> &temp)
{
    s<<static_cast<const Eigen::SparseMatrixBase<Eigen::SparseVector<T,option>>&>(temp.m);
    return s;
}

//typedef Eigen::Matrix<int, -1, -1, Eigen::RowMajor> EigenMatrixRow;
typedef Eigen::Matrix<int, -1, -1, Eigen::RowMajor> EigenMatrixRowI;
typedef Eigen::Matrix<double, -1, -1, Eigen::RowMajor> EigenMatrixRowD;


inline __attribute__((always_inline))
bool doubleLess(double a, double b, double ep=1E-3)
{
    return (a-b) < -ep;
}


inline __attribute__((always_inline))
bool doubleEqual(double a, double b, double ep=std::numeric_limits<double>::epsilon())
{
    //return std::fabs(a-b) < ep;
    constexpr double normal_min = std::numeric_limits<double>::min();

    if (!std::isfinite(a) || !std::isfinite(b)) // a = ±∞, NaN or b = ±∞, NaN
        return false;
    double diff = std::abs(a - b);
    // if "a" and "b" are near to zero, the relative error is less effective
    if (diff <= normal_min) // or also: user_epsilon * normal_min
        return true;

    double abs_a = std::abs(a);
    double abs_b = std::abs(b);
    return (diff / std::max(abs_a, abs_b)) <= ep;
}


inline __attribute__((always_inline))
bool doubleEqual(float a, float b, double ep=std::numeric_limits<double>::epsilon())
{
    //return std::fabs(a-b) < ep;
    constexpr float normal_min = std::numeric_limits<float>::min();

    if (!std::isfinite(a) || !std::isfinite(b)) // a = ±∞, NaN or b = ±∞, NaN
        return false;
    double diff = std::abs(a - b);
    // if "a" and "b" are near to zero, the relative error is less effective
    if (diff <= normal_min) // or also: user_epsilon * normal_min
        return true;

    float abs_a = std::abs(a);
    float abs_b = std::abs(b);
    return (diff / std::max(abs_a, abs_b)) <= ep;
}

#endif //DW_DECOMP_AUX_H
