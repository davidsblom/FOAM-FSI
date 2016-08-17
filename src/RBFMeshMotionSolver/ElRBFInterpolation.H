
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <memory>
#include "El.hpp"
#include "RBFFunctionInterface.H"

namespace rbf
{
    class ElRBFInterpolation
    {
        public:
            ElRBFInterpolation(
                std::unique_ptr<RBFFunctionInterface> & rbfFunction,
                std::unique_ptr<El::DistMultiVec<double> > & positions,
                std::unique_ptr<El::DistMultiVec<double> > & posInterpolation
                );

            ~ElRBFInterpolation();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > values );

        private:
            const std::unique_ptr<RBFFunctionInterface> rbfFunction;
            const std::unique_ptr<El::DistMultiVec<double> > positions;
            const std::unique_ptr<El::DistMultiVec<double> > positionsInterpolation;
            std::unique_ptr<El::DistMatrix<double> > H;
    };
}
