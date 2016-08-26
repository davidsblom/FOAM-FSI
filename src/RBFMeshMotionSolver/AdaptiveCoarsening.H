
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
                std::unique_ptr<El::DistMatrix<double> > positions,
                std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
                );

            void greedySelection( const std::unique_ptr<El::DistMatrix<double> > & values );

            bool initialized();

            std::unique_ptr<El::DistMatrix<double> > interpolate( const std::unique_ptr<El::DistMatrix<double> > & values );

        private:
            std::pair<int, double> computeError( const std::unique_ptr<El::DistMatrix<double> > & values );

            void selectData(
                const std::unique_ptr<El::DistMatrix<double> > & data,
                std::unique_ptr<El::DistMatrix<double> > & selection
                );

            const double tol;
            const double reselectionTol;
            const int minPoints;
            const int maxPoints;
            std::unique_ptr<ElRBFInterpolation> rbf;
            std::unique_ptr<ElRBFInterpolation> rbfCoarse;
            std::vector<size_t> selectedPositions;

            std::shared_ptr<RBFFunctionInterface> rbfFunction;
            std::unique_ptr<El::DistMatrix<double> > positions;
            std::unique_ptr<El::DistMatrix<double> > positionsInterpolation;
    };

    template<typename Real>
    Real MaxAbs( const El::DistMatrix<Real> & A )
    {
        Real value = 0;

        if ( A.Participating() )
        {
            // Store the index/value of the local pivot candidate
            const El::Int mLocal = A.LocalHeight();
            const El::Int nLocal = A.LocalWidth();
            const Real * ABuf = A.LockedBuffer();
            const El::Int ALDim = A.LDim();

            for ( El::Int jLoc = 0; jLoc < nLocal; ++jLoc )
                for ( El::Int iLoc = 0; iLoc < mLocal; ++iLoc )
                    value = El::Max( std::abs( value ), ABuf[iLoc + jLoc * ALDim] );

            value = El::mpi::AllReduce( value, El::mpi::MAX, A.DistComm() );
        }

        El::mpi::Broadcast( value, A.Root(), A.CrossComm() );
        return value;
    }
}
