
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <cstdlib>
#include <ctime>
#include "gtest/gtest.h"
#include "PetRBFInterpolation.H"
#include "fvCFD.H"

TEST( RBFInterpolation, constructor )
{
    Eigen::MatrixXd pos, posInterp;
    pos.setRandom( 5, 3 );
    posInterp.setRandom( 23, pos.cols() );
    fsi::PetRBFInterpolation rbf( pos, posInterp );
}

TEST( RBFInterpolation, constructor_random_size )
{
    std::random_device rd;
    std::mt19937 gen( rd() );
    std::uniform_int_distribution<> dis( 10, 50 );

    Eigen::MatrixXd pos, posInterp;
    pos.setRandom( dis( gen ), 3 );
    posInterp.setRandom( dis( gen ), pos.cols() );
    fsi::PetRBFInterpolation rbf( pos, posInterp );
}

TEST( RBFInterpolation, constructor_zero )
{
    Eigen::MatrixXd pos, posInterp;
    pos.setRandom( 0, 3 );
    posInterp.setRandom( 0, pos.cols() );
    fsi::PetRBFInterpolation rbf( pos, posInterp );
}
