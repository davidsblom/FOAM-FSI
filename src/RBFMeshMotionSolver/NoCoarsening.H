
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "Coarsener.H"

namespace rbf
{
    class NoCoarsening : public Coarsener
    {
        public:
            explicit NoCoarsening( std::unique_ptr<ElRBFInterpolation> rbf );

            ~NoCoarsening();

            void compute(
                std::unique_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
                );

            bool initialized();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values );

        private:
            std::unique_ptr<ElRBFInterpolation> rbf;
    };
}
