
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "Coarsener.H"

namespace rbf
{
    class AdaptiveCoarsening : public Coarsener
    {
        public:
            AdaptiveCoarsening(
                double tol,
                double reselectionTol,
                int minPoints,
                int maxPoints
                );

            ~AdaptiveCoarsening();

            void compute(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<ElDistVector> positions,
                std::unique_ptr<ElDistVector> positionsInterpolation
                );

            void greedySelection( const std::unique_ptr<ElDistVector> & values );

            bool initialized();

            std::unique_ptr<ElDistVector> interpolate( const std::unique_ptr<ElDistVector> & values );

        private:
            std::pair<int, double> computeError( const std::unique_ptr<ElDistVector> & values );

            void selectData(
                const std::unique_ptr<ElDistVector> & data,
                std::unique_ptr<ElDistVector> & selection
                );

            const double tol;
            const double reselectionTol;
            const int minPoints;
            const int maxPoints;
            std::unique_ptr<ElRBFInterpolation> rbf;
            std::unique_ptr<ElRBFInterpolation> rbfCoarse;
            std::vector<size_t> selectedPositions;

            std::shared_ptr<RBFFunctionInterface> rbfFunction;
            std::unique_ptr<ElDistVector> positions;
            std::unique_ptr<ElDistVector> positionsInterpolation;
    };
}
