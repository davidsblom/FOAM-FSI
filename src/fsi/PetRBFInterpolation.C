
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PetRBFInterpolation.H"
#include "PetscVector.H"
#include <vector>

namespace fsi
{
    PetRBFInterpolation::PetRBFInterpolation(
        const Eigen::MatrixXd & pos,
        const Eigen::MatrixXd & posInterp
        )
    {
        assert( pos.cols() == posInterp.cols() );
        assert( pos.cols() == 2 || pos.cols() == 3 );

        std::vector<PetscVector> positions, positionsInterpolation;

        const int dim = pos.cols();

        for ( int i = 0; i < dim; i++ )
        {
            positions.push_back( PetscVector( pos.rows(), false ) );
            positionsInterpolation.push_back( PetscVector( posInterp.rows(), false ) );
        }

        assert( positions[0].localSize() == pos.rows() );
        assert( positions[0].localSize() == positions[1].localSize() );
        assert( positionsInterpolation[0].localSize() == posInterp.rows() );

        for ( int i = 0; i < pos.rows(); i++ )
            for ( int j = 0; j < pos.cols(); j++ )
                positions[j].setLocal( i, pos( i, j ) );

        for ( int i = 0; i < posInterp.rows(); i++ )
            for ( int j = 0; j < posInterp.cols(); j++ )
                positionsInterpolation[j].setLocal( i, posInterp( i, j ) );

        for ( auto && vec : positions )
            vec.compress();

        for ( auto && vec : positionsInterpolation )
            vec.compress();

        assert( positions[0].localSize() == pos.rows() );
        assert( positions[0].localSize() == positions[1].localSize() );
        assert( positionsInterpolation[0].localSize() == posInterp.rows() );
    }

    PetRBFInterpolation::~PetRBFInterpolation()
    {}
}
