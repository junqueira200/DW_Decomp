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
std::ostream & operator << (std::ostream & s, const TempSpMatPrint<T,option> &temp)
{
    s<<static_cast<const Eigen::SparseMatrixBase<Eigen::SparseMatrix<T,option>>&>(temp.m);
    return s;
}

#endif //DW_DECOMP_AUX_H
