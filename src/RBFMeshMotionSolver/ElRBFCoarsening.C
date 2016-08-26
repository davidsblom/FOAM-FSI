
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFCoarsening.H"

namespace rbf
{
    ElRBFCoarsening::ElRBFCoarsening( std::unique_ptr<Coarsener> coarsener )
        :
        coarsener( std::move( coarsener ) )
    {}

    ElRBFCoarsening::~ElRBFCoarsening(){}

    void ElRBFCoarsening::compute(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<El::DistMatrix<double> > positions,
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation
        )
    {
        coarsener->compute( rbfFunction, std::move( positions ), std::move( positionsInterpolation ) );
    }

    bool ElRBFCoarsening::initialized()
    {
        return coarsener->initialized();
    }

    std::unique_ptr<El::DistMatrix<double> > ElRBFCoarsening::interpolate( const std::unique_ptr<El::DistMatrix<double> > & values )
    {
        return coarsener->interpolate( values );
    }
}
