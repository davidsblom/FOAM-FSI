
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "NoCoarsening.H"

namespace rbf
{
    NoCoarsening::NoCoarsening()
        :
        rbf( new ElRBFInterpolation() )
    {}

    NoCoarsening::~NoCoarsening(){}

    void NoCoarsening::compute(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        std::unique_ptr<ElDistVector> positions,
        std::unique_ptr<ElDistVector> positionsInterpolation
        )
    {
        rbf->compute( rbfFunction, std::move( positions ), std::move( positionsInterpolation ) );
    }

    bool NoCoarsening::initialized()
    {
        return rbf->initialized();
    }

    std::unique_ptr<ElDistVector> NoCoarsening::interpolate( const std::unique_ptr<ElDistVector> & values )
    {
        return rbf->interpolate( values );
    }
}
