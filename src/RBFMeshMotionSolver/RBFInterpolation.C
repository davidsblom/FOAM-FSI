
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"

namespace rbf
{
  RBFInterpolation::RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction )
    :
    rbfFunction( rbfFunction ),
    computed( false ),
    n_A( 0 ),
    n_B( 0 ),
    dimGrid( 0 ),
    Hhat()
  {
    assert( rbfFunction );
  }

  void RBFInterpolation::compute(
    const matrix & positions,
    const matrix & positionsInterpolation
    )
  {
    // Verify input

    assert( positions.cols() == positionsInterpolation.cols() );
    assert( positions.rows() > 0 );
    assert( positions.cols() > 0 );
    assert( positionsInterpolation.rows() > 0 );

    n_A = positions.rows();
    n_B = positionsInterpolation.rows();
    Hhat.resize( n_B, n_A );
    dimGrid = positions.cols();

    // Radial basis function interpolation

    // Initialize variables

    matrix Q_A( n_A, dimGrid + 1 ), Q_B( n_B, dimGrid + 1 );

    // Initialize sparse matrices

    matrix H( n_A + dimGrid + 1, n_A + dimGrid + 1 ), Phi( n_B, n_A + dimGrid + 1 );
    H.setZero();

    // RBF function evaluation

    double r;

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = i; j < n_A; j++ )
      {
        r = ( positions.row( i ) - positions.row( j ) ).norm();
        H( j, i ) = rbfFunction->evaluate( r );
      }

      for ( int j = 0; j < i; j++ )
      {
        H( j, i ) = H( i, j );
      }
    }

    // Evaluate Q_A
    for ( int i = 0; i < n_A; i++ )
      Q_A( i, 0 ) = 1;

    Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

    // Build the matrix H

    H.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
    H.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

    // Evaluate Phi_BA which contains the evaluation of the radial basis function

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = 0; j < n_B; j++ )
      {
        r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
        Phi( j, i ) = rbfFunction->evaluate( r );
      }
    }

    // Evaluate Q_B

    for ( int i = 0; i < n_B; i++ )
      Q_B( i, 0 ) = 1;

    Q_B.block( 0, 1, n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

    // Determine the matrix Phi

    Phi.topRightCorner( Q_B.rows(), Q_B.cols() ) = Q_B;

    // Compute the LU decomposition of the matrix H

    Eigen::FullPivLU<matrix> lu( H );
    matrix Hinverse = lu.inverse();

    // Matrix matrix multiplication
    // Hhat = Phi * inv(H)

    Hhat.setZero();

    for ( int i = 0; i < Hhat.rows(); i++ )
      for ( int j = 0; j < Hhat.cols(); j++ )
        for ( int k = 0; k < Phi.cols(); k++ )
          Hhat( i, j ) += Phi( i, k ) * Hinverse( k, j );

    computed = true;
  }

  void RBFInterpolation::interpolate(
    const matrix & values,
    matrix & valuesInterpolation
    )
  {
    assert( computed );

    valuesInterpolation = Hhat * values;

    assert( valuesInterpolation.rows() == n_B );
    assert( values.cols() == valuesInterpolation.cols() );
  }

  /*
   * Compute interpolation matrix and directly interpolate the values.
   * The algorithms solves for the coefficients, and explicitly
   * uses the coefficients to interpolate the data to the new positions.
   */
  void RBFInterpolation::interpolate(
    const matrix & positions,
    const matrix & positionsInterpolation,
    const matrix & values,
    matrix & valuesInterpolation
    )
  {
    // Verify input

    assert( positions.cols() == positionsInterpolation.cols() );
    assert( positions.rows() > 0 );
    assert( positions.cols() > 0 );
    assert( positionsInterpolation.rows() > 0 );

    n_A = positions.rows();
    n_B = positionsInterpolation.rows();
    Hhat.resize( n_B, n_A );
    Hhat.setZero();
    dimGrid = positions.cols();

    // Radial basis function interpolation

    // Initialize matrices

    matrix H( n_A, n_A );
    H.setZero();
    Phi.conservativeResize( n_B, n_A );

    // RBF function evaluation

    double r;

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = i; j < n_A; j++ )
      {
        r = ( positions.row( i ) - positions.row( j ) ).norm();
        H( j, i ) = rbfFunction->evaluate( r );
      }

      for ( int j = 0; j < i; j++ )
      {
        H( j, i ) = H( i, j );
      }
    }

    // Calculate coefficients gamma and beta

    matrix B = H.fullPivLu().solve( values );

    // Evaluate Phi_BA which contains the evaluation of the radial basis function
    // This method is only used by the greedy algorithm, and the matrix Phi
    // is therefore enlarged at every greedy step.

    int i = Phi.cols() - 1;

    for ( int j = 0; j < n_B; j++ )
    {
      r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
      Phi( j, i ) = rbfFunction->evaluate( r );
    }

    valuesInterpolation = Phi * B;
  }
}
