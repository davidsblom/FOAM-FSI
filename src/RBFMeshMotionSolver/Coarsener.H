
#pragma once

#include "ElRBFInterpolation.H"

namespace rbf
{
    class Coarsener
    {
        public:
            virtual ~Coarsener(){}

            virtual void compute(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
                ) = 0;

            virtual bool initialized() = 0;

            virtual std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values ) = 0;
    };
}
