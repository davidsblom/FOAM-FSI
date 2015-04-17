
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFCoarsening_H
#define RBFCoarsening_H

#include "RBFInterpolation.H"
#include "fvCFD.H"
#include <fstream>

namespace rbf
{
    class RBFCoarsening
    {
public:

        RBFCoarsening( std::shared_ptr<RBFInterpolation> rbf );

        RBFCoarsening(
            std::shared_ptr<RBFInterpolation> rbf,
            bool enabled,
            bool livePointSelection,
            bool livePointSelectionSumValues,
            double tol,
            double tolLivePointSelection,
            int coarseningMinPoints,
            int coarseningMaxPoints,
            bool exportTxt
            );

        void greedySelection( const matrix & values );

        void compute(
            const matrix & positions,
            const matrix & positionsInterpolation
            );

        void interpolate(
            const matrix & values,
            matrix & valuesInterpolation
            );

        void setNbMovingAndStaticFaceCenters(
            int nbMovingFaceCenters,
            int nbStaticFaceCenters
            );

        std::shared_ptr<RBFInterpolation> rbf;
        std::shared_ptr<RBFInterpolation> rbfCoarse;
        bool enabled;
        bool livePointSelection;
        bool livePointSelectionSumValues;
        double tol;
        double tolLivePointSelection;
        int coarseningMinPoints;
        int coarseningMaxPoints;
        bool exportTxt;
        Eigen::VectorXi selectedPositions;
        int nbStaticFaceCentersRemove;
        matrix positions;
        matrix positionsInterpolation;
        matrix values;
        int nbMovingFaceCenters;
        int fileExportIndex;
    };
}

#endif
