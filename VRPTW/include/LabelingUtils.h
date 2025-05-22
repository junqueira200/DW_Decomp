#ifndef LABELINGUTILS_H
#define LABELINGUTILS_H

#include "Label.h"
#include "DW_Decomp.h"

namespace LabelingAlgorithmNS
{


    void checkHeap(LabelHeap& heap, LabelingData& lData);
    void convertLabelBackwardToForward(Label* label, const ArrayResources &vetMaxResources, int numResources);
    void copyLabel(const Label& labelSrc, Label& labelDest, int numResorces);
    void whiteLabelIndex(Label* label);

    void removeCycles(Label &label, const int numCust);
    void removeCycles2(Label &label, const int numCust);
    void updateLabelCost(Label &label, const Vet3D_ResCost &vetMatResCost, FloatType labelStart);

    // Linear index for a nxn matrix
    inline __attribute__((always_inline))
    int getIndex(int i, int j, int numClie){return i*numClie+j;}

    bool checkDistance(const Eigen::Matrix<FloatType, -1, -1, Eigen::RowMajor> &matDist);
    bool containRoute(const Eigen::Array<Label*, 1, DW_DecompNS::NumMaxSolSubProb> &vetLabel, int numSol, Label* label);

    bool labelHaveRoute(std::vector<int> &vetRoute, Label *label);


}

#endif // LABELINGUTILS_H
