
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <memory>
#include <El.hpp>
#include "RBFFunctionInterface.H"

namespace rbf
{
    typedef El::DistMatrix<double, El::VR, El::STAR> ElDistVector;

    class ElRBFInterpolation
    {
        public:
            explicit ElRBFInterpolation();

            explicit ElRBFInterpolation(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<ElDistVector> positions,
                std::unique_ptr<ElDistVector> positionsInterpolation
                );

            ~ElRBFInterpolation();

            void compute(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<ElDistVector> positions,
                std::unique_ptr<ElDistVector> positionsInterpolation
                );

            bool initialized();

            std::unique_ptr<ElDistVector> interpolate( const std::unique_ptr<ElDistVector> & values );

        private:
            std::unique_ptr<El::DistMatrix<double> > H;
            std::unique_ptr<El::DistMatrix<double> > Phi;

            El::DistMatrix<double> HCopy;
            El::DistMatrix<double, El::MD, El::STAR> dSub;
            El::DistPermutation p;
    };
}
