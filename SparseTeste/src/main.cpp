#include <iostream>
#include "DW_Decomp.h"
#include "Sparse.h"
#include "Aux.h"
#include "SparseOp.h"


std::ostream& operator<<(std::ostream& os, const std::tuple<int, int, int> &tuple)
{
    os<<"("<<std::get<0>(tuple)<<"; "<<std::get<1>(tuple)<<"; "<<std::get<2>(tuple)<<") ";

    return os;
}


int main()
{
    SparseNS::SparseMatCol<double> sparseMatCol(5, 5);
    Vector<std::tuple<int, int, int>> vetTuple;

    for(int i=0; i < sparseMatCol.numLin; ++i)
    {
        for(int j=0; j < sparseMatCol.numCol; ++j)
        {
            vetTuple.push_back(std::make_tuple(i, j, i*sparseMatCol.numCol + j));
        }
    }


    for(int i=vetTuple.size()-1; i >= 0; --i)
    {
        auto tuple = vetTuple[i];
        int ii = std::get<0>(tuple);
        int jj = std::get<1>(tuple);
        int val = std::get<2>(tuple);

        sparseMatCol(ii, jj) = val;
    }

    std::cout<<sparseMatCol<<"\n\n";

    SparseNS::SparseVector<double> vet(5);
    SparseNS::SparseVector<double> vetResul(5);


    for(int i=0; i < vet.n; ++i)
    {
        vet(i) = i;
    }

    std::cout<<vet<<"\n";

    SparseOpNS::multSparseMatColTimesVetCol(sparseMatCol, vet, vetResul);

    std::cout<<vetResul;


    Eigen::SparseMatrix<double, Eigen::ColMajor> spMat(2, 2);
    spMat.coeffRef(0,0) = 2;

    Eigen::SparseVector<double, Eigen::ColMajor> spVet(2);
    Eigen::SparseVector<double, Eigen::RowMajor> spVet2(2);
    //Eigen::VectorXd spVet2(2);

    //spVet.coeffRef(0) = 1.0;
    //spVet.coeffRef(1) = 1.0;
    //spVet2 = spMat*spVet;

    Eigen::SparseMatrix<double, Eigen::ColMajor> r0(3, 1);
    r0 = (spMat*spVet).pruned();

    std::cout<<spMat*spVet<<"\n";
    std::cout<<spVet2*spMat<<"\n";

    //std::cout<<"r0: \n"<<r0<<"\n\n";
    std::cout<<"r0: \n"<<TempSpMatPrint(r0)<<"\n";//<<static_cast<const Eigen::SparseMatrixBase<Eigen::SparseMatrix<double>>&>(r0);
    //std::cout<<spVet<<"\n";


    return 0;
} // FIM main
