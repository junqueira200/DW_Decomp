/* ****************************************
 * ****************************************
 *  Date:    20/03/26
 *  File:    IBM_CpOptimizer.h
 *  Autor:   Igor de Andrade Junqueira
 *  Project: 2L-SDVRP
 * ****************************************
 * ****************************************/

#ifndef IBM_CPOPTIMIZER_H
#define IBM_CPOPTIMIZER_H

#include "Instancia.h"
#include "Solucao.h"
#include "InputOutput.h"
#include <ilcp/cp.h>
#include <ilconcert/ilomodel.h>
#include "AuxT.h"
#include "OPP_CP_3D.h"

namespace IBM_CpOptimizerNS
{
    using namespace ContainerLoading;
    /** ******************************************************************************************
     *  ******************************************************************************************
     * 	@brief Solve a single-bin packing problem using IBM ILOG CP Optimizer.
     *
     * 	This class models and solves a 3D bin packing (or container loading)
     * 	problem for a given set of items associated with a route. It attempts to
     * 	pack the specified items into a single bin while respecting all problem
     * 	constraints (e.g., dimensions, rotation, etc.).
     *
     *  ******************************************************************************************
     *  ******************************************************************************************
     */
    class CpOptimizer
    {
    public:

        static IloEnv* envPtr;
        IloModel model;

        const VectorI& vetItems;
        int numItems;
        const SolucaoNS::Rota& rota;

        IloIntVarArray mStartPositionsX;
        IloIntVarArray mEndPositionsX;

        IloIntVarArray mStartPositionsY;
        IloIntVarArray mEndPositionsY;

        IloIntVarArray mStartPositionsZ;
        IloIntVarArray mEndPositionsZ;

        IloIntVarArray mWidths;
        IloIntVarArray mLengths;
        IloIntVarArray mHeights;

        IloBoolVarArray mPlacedOnFloor;

        IloArray<IloArray<IloBoolVarArray>> mRelativeDirections; // mRelativeDirections[i][j][direction], bool, 1
                                                                 // if item i is placed relatively to item j in direction

        IloArray<IloBoolVarArray> mItemsOverlapsXY; // mItemsOverlapsXY[i][j], bool, items i and j intersect in xy-plane
        IloArray<IloBoolVarArray> mSupportXY;       // mSupportXY[i][j], bool, 1, if item i is supported by item j
                                                    // xy-plane ? -> items intersect AND item j is directly below item i
        IloArray<IloIntVarArray> mOverlapAreasXY;  // mOverlapAreasXY[i][j], integer, size of intersection area in
                                                    // xy-plane of items i and j
        IloArray<IloBoolVarArray> mOrientation;

        std::vector<Dimension> mDimensions = {{AxisY, RightY, LeftY}, {AxisX, InFrontX, BehindX}, {AxisZ, AboveZ, BelowZ}};


        CpOptimizer(const VectorI& vetItems, int numItems, const SolucaoNS::Rota& rota);
        bool solve(SolucaoNS::Bin& bin);
        void createVariables();
        void CreateNoOverlap();
        void CreateItemOrientations();
        void CreateSupportItem();
        void CreateSupportArea();
        void CreateXYIntersectionBool();
        void CreateXYIntersectionArea();
        void CreateLifo();
        void CreateOnFloorConstraints();
        void CreateAxleWeights();

        IloIntVar& getIntVars(DimensionType dimension, bool first, int i);

    };

    void inline
    NameVars(IloIntVarArray& a, const char * base) {
        for (IloInt i = 0; i < a.getSize(); i++) {
            char name[100];
            //sprintf(name, "%s[%ld]", base, (long)i);
            snprintf(name, 100, "%s[%ld]", base, (long)i);
            a[i].setName(name);
        }
    }

    void inline
    NameVars(IloIntervalVarArray& a, const char * base) {
        for (IloInt i = 0; i < a.getSize(); i++) {
            char name[100];
            //sprintf(name, "%s[%ld]", base, (long)i);
            snprintf(name, 100, "%s[%ld]", base, (long)i);
            a[i].setName(name);
        }
    }

    void inline
    NameVars(IloArray<IloIntVarArray>& a, const char * base) {
        for (IloInt i = 0; i < a.getSize(); i++) {
            for (IloInt j = 0; j < a[i].getSize(); j++) {
                char name[100];
                //sprintf(name, "%s[%ld][%ld]", base, (long)i, long(j));
                snprintf(name, 100, "%s[%ld][%ld]", base, (long)i, (long)j);
                a[i][j].setName(name);
            }
        }
    }


    void inline
    NameVars(IloArray<IloIntervalVarArray>& a, const char * base) {
        for (IloInt i = 0; i < a.getSize(); i++) {
            for (IloInt j = 0; j < a[i].getSize(); j++) {
                char name[100];
                //sprintf(name, "%s[%ld][%ld]", base, (long)i, long(j));
                snprintf(name, 100, "%s[%ld][%ld]", base, (long)i, (long)j);
                a[i][j].setName(name);
            }
        }
    }

    void testSCIP();

}

#endif // IBM_CPOPTIMIZER_H
