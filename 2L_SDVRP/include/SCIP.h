#ifndef SCIP_H
#define SCIP_H

#include "Instancia.h"
#include "Solucao.h"
#include <scippp/model.hpp>

namespace SCIP_NS
{

    typedef std::vector<scippp::Var> VectorVar;
    class MatrixVar
    {
    private:

        std::vector<scippp::Var> vetVar;
        size_t n, m;

    public:

        MatrixVar(scippp::Model& model, size_t n_, size_t m_, std::string&& prex, scippp::VarType type);
        MatrixVar()=default;

        INLINE
        scippp::Var& operator ()(size_t i, size_t j)
        {
            if(i > n)
                throw std::out_of_range(std::format("Error, i({}) > n({})", i, n));

            if(j > m)
                throw std::out_of_range(std::format("Error, i({}) > m({})", i, m));

            return vetVar[i*m+j];
        }

    };

    class Matrix3DVar
    {
    private:

        std::vector<scippp::Var> vetVar;
        size_t n, m, p;

    public:

        Matrix3DVar(scippp::Model& model, size_t n_, size_t m_, size_t p_, std::string&& prex, scippp::VarType type);
        Matrix3DVar()=default;

        INLINE
        scippp::Var& operator ()(size_t i, size_t j, size_t k)
        {
            if(i > n)
                throw std::out_of_range(std::format("Error, i({}) > n({})", i, n));

            if(j > m)
                throw std::out_of_range(std::format("Error, j({}) > m({})", j, m));

            if(k > p)
                throw std::out_of_range(std::format("Error, k({}) > p({})", k, p));

            return vetVar[(i*m*p)+(j*p)+k];
        }

    };

    class Scip3dPacking
    {
    public:

        const VectorI& vetItems;
        int numItems;
        const SolucaoNS::Rota& rota;

        scippp::Model model;
        SCIP *ptrScip = nullptr;

        VectorVar mStartPositionsX;
        VectorVar mEndPositionsX;
        VectorVar mStartPositionsY;
        VectorVar mEndPositionsY;
        VectorVar mStartPositionsZ;
        VectorVar mEndPositionsZ;

        Matrix3DVar mRelativeDirections;

        MatrixVar mItemsOverlapsXY; 	// mItemsOverlapsXY[i][j], bool, items i and j intersect in xy-plane
        MatrixVar mSupportXY; 			// mSupportXY[i][j], bool, 1, if item i is supported by item j xy-plane ? -> items
                                        // intersect AND item j is directly below item i
        MatrixVar mOverlapAreasXY; 		// mOverlapAreasXY[i][j], integer, size of intersection area in xy-plane of items i and j

        MatrixVar mWidths;
        MatrixVar mLengths;
        MatrixVar mHeights;

        Scip3dPacking(const VectorI& vetItems_, const int numItems_, const SolucaoNS::Rota& rota_);



    };


}

#endif // SCIP_H
