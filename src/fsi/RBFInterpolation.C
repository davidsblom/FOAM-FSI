
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"

namespace fsi
{
  RBFInterpolation::RBFInterpolation()
    :
    computed( false ),
    n_A( 0 ),
    n_B( 0 ),
    dimGrid( 0 ),
    radius( 1 ),
    function( "thin_plate" ),
    Hhat()
  {}

  RBFInterpolation::RBFInterpolation(
    const std::string & function,
    const double radius
    )
    :
    computed( false ),
    n_A( 0 ),
    n_B( 0 ),
    dimGrid( 0 ),
    radius( radius ),
    function( function ),
    Hhat()
  {
    assert( function == "thin_plate" || function == "wendland_c2" );
    assert( radius > 0 );
  }

  RBFInterpolation::RBFInterpolation(
    const matrix & positions,
    const matrix & positionsInterpolation,
    const std::string & function,
    const double radius
    )
    :
    computed( false ),
    n_A( positions.rows() ),
    n_B( positionsInterpolation.rows() ),
    dimGrid( positions.cols() ),
    radius( radius ),
    function( function ),
    Hhat( n_B, n_A )
  {
    // Verify input
    assert( function == "thin_plate" || function == "wendland_c2" );
    assert( radius > 0 );

    compute( positions, positionsInterpolation );
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
    Hhat.setZero();
    dimGrid = positions.cols();

    // If position matrices are equal, the identity interpolation matrix
    // is used.
    if ( positions.rows() == positionsInterpolation.rows() )
    {
      bool equal = true;

      for ( int i = 0; i < positions.rows(); i++ )
      {
        if ( equal )
        {
          for ( int j = 0; j < positions.cols(); j++ )
          {
            if ( std::abs( positions( i, j ) - positionsInterpolation( i, j ) ) > 1.0e-15 )
            {
              equal = false;
              break;
            }
          }
        }
      }

      if ( equal )
      {
        Hhat = Eigen::MatrixXd::Identity( n_B, n_A );
        computed = true;
        return;
      }
    }

    // Radial basis function interpolation

    // Initialize variables

    matrix Q_A( n_A, dimGrid + 1 ), Q_B( n_B, dimGrid + 1 );
    Q_A.setZero();
    Q_B.setZero();

    // Initialize sparse matrices

    matrix H( n_A + dimGrid + 1, n_A + dimGrid + 1 ), Phi( n_B, n_A + dimGrid + 1 ), Phi_AA( n_A, n_A );
    H.setZero();

    // Evaluate Phi_AA which contains the evaluation of the radial basis function

    double r;

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = 0; j < n_A; j++ )
      {
        r = ( positions.row( i ) - positions.row( j ) ).norm();
        H( j, i ) = evaluateFunction( r );
      }
    }

    // Evaluate Q_A
    for ( int i = 0; i < n_A; i++ )
      Q_A( i, 0 ) = 1;

    Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

    // Build the matrix H

    H.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
    H.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

    // Compute the LU decomposition of the matrix H

    Eigen::FullPivLU<matrix> lu( H );

    // Evaluate Phi_BA which contains the evaluation of the radial basis function

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = 0; j < n_B; j++ )
      {
        r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
        Phi( j, i ) = evaluateFunction( r );
      }
    }

    // Evaluate Q_B

    for ( int i = 0; i < n_B; i++ )
      Q_B( i, 0 ) = 1;

    Q_B.block( 0, 1, n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

    // Determine the matrix Phi

    Phi.topRightCorner( Q_B.rows(), Q_B.cols() ) = Q_B;

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

  double RBFInterpolation::evaluateFunction( const double x )
  {
    if ( function == "thin_plate" )
      return TPS( x );

    // if ( function == "wendland_c2")
    return WendlandC2( x );
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

    assert( !computed );
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

    // Initialize variables

    matrix Q_A( n_A, dimGrid + 1 ), Q_B( n_B, dimGrid + 1 );
    Q_A.setZero();
    Q_B.setZero();

    // Initialize sparse matrices

    matrix H( n_A + dimGrid + 1, n_A + dimGrid + 1 ), Phi( n_B, n_A + dimGrid + 1 ), Phi_AA( n_A, n_A );
    H.setZero();

    /*
     * Evaluate Phi_AA which contains the evaluation of the radial basis function
     * Wendland's radial basis functions and the thin plate spline are positive definite.
     * Take advantage of that by reducing the number of radial basis evaluations by
     * a factor two.
     */

    double r;

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = i; j < n_A; j++ )
      {
        r = ( positions.row( i ) - positions.row( j ) ).norm();
        H( j, i ) = evaluateFunction( r );
      }

      for ( int j = 0; j < i; j++ )
        H( j, i ) = H( i, j );
    }

    // Evaluate Q_A
    for ( int i = 0; i < n_A; i++ )
      Q_A( i, 0 ) = 1;

    Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

    // Build the matrix H

    H.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
    H.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

    matrix W( values.rows() + values.cols() + 1, values.cols() );
    W.setZero();
    W.topLeftCorner( values.rows(), values.cols() ) = values;

    // Calculate coefficients gamma and beta

    if ( function == "wendland_c2" )
      B = H.ldlt().solve( W );
    else
      B = H.fullPivLu().solve( W );

    // Evaluate Phi_BA which contains the evaluation of the radial basis function

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = 0; j < n_B; j++ )
      {
        r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
        Phi( j, i ) = evaluateFunction( r );
      }
    }

    // Evaluate Q_B

    for ( int i = 0; i < n_B; i++ )
      Q_B( i, 0 ) = 1;

    Q_B.block( 0, 1, n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

    // Determine the matrix Phi

    Phi.topRightCorner( Q_B.rows(), Q_B.cols() ) = Q_B;

    valuesInterpolation = Phi * B;
  }

  void RBFInterpolation::interpolate(
    const matrix & positions,
    const matrix & positionsInterpolation,
    matrix & valuesInterpolation
    )
  {
    assert( positions.cols() == positionsInterpolation.cols() );
    assert( valuesInterpolation.rows() == positionsInterpolation.rows() );

    // Initialize variables

    n_A = positions.rows();
    n_B = positionsInterpolation.rows();
    dimGrid = positions.cols();

    matrix Phi( n_B, n_A + dimGrid + 1 ), Q_B( n_B, dimGrid + 1 );
    Q_B.setZero();

    // Evaluate Phi_BA which contains the evaluation of the radial basis function

    double r;

    for ( int i = 0; i < n_A; i++ )
    {
      for ( int j = 0; j < n_B; j++ )
      {
        r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
        Phi( j, i ) = evaluateFunction( r );
      }
    }

    // Evaluate Q_B

    for ( int i = 0; i < n_B; i++ )
      Q_B( i, 0 ) = 1;

    Q_B.block( 0, 1, n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

    // Determine the matrix Phi

    Phi.topRightCorner( Q_B.rows(), Q_B.cols() ) = Q_B;

    valuesInterpolation = Phi * B;
  }

  /*
   * Thin plate spline radial basis function
   */
  double RBFInterpolation::TPS( const double x )
  {
    if ( x > 0 )
      return std::log10( x ) * x * x;

    return 0;
  }

  /*
   * Wendland's C2. Wendland's radial basis functions
   * result in a positive definite matrix, which can be solved
   * twice as fast with a Cholesky decompositions instead of a full
   * LU decomposition for a non-symmetric matrix.
   */
  double RBFInterpolation::WendlandC2( const double x )
  {
    return std::pow( 1 - x, 4 ) * (4 * x + 1);
  }
}
