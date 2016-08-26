
#pragma once

#include "Coarsener.H"
#include "ElRBFInterpolation.H"

namespace rbf
{
    class UnitCoarsening : public Coarsener
    {
        public:
            UnitCoarsening(
                double tol,
                int minPoints,
                int maxPoints
                );

            ~UnitCoarsening();

            void compute(
                std::shared_ptr<RBFFunctionInterface> rbfFunction,
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
                );

            bool initialized();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values );

        private:
            void selectData(
                const std::unique_ptr<El::DistMatrix<double> > & data,
                std::unique_ptr<El::DistMatrix<double> > & selection
                );

            const double tol;
            const int minPoints;
            const int maxPoints;
            ElRBFInterpolation rbf;
            std::vector<size_t> selectedPositions;
    };
}
