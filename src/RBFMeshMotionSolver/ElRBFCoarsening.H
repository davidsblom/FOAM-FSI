
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "ElRBFInterpolation.H"
#include "Coarsener.H"

namespace rbf
{
    class ElRBFCoarsening
    {
        public:
            explicit ElRBFCoarsening( std::unique_ptr<Coarsener> coarsener );

            ~ElRBFCoarsening();

            void compute(
                std::unique_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
                );

            bool initialized();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values );

        private:
            std::unique_ptr<Coarsener> coarsener;
    };
}
