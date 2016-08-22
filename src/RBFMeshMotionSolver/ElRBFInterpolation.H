
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
            explicit ElRBFInterpolation(
                std::unique_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > posInterpolation
                );

            ~ElRBFInterpolation();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values );

        private:
            std::unique_ptr<El::DistMatrix<double> > H;
            std::unique_ptr<El::DistMatrix<double> > Phi;
    };
}
