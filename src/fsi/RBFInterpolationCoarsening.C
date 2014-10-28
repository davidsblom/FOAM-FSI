
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolationCoarsening.H"

using std::deque;
using namespace fsi;

RBFInterpolationCoarsening::RBFInterpolationCoarsening( double tol )
  :
  tol( tol ),
  rbf( "wendland_c2", 1 ),
  positions( 0, 0 ),
  positionsInterpolation( 0, 0 )
{
  assert( tol > 1.0e-13 );
  assert( tol < 1 );
}

RBFInterpolationCoarsening::~RBFInterpolationCoarsening()
{}

/*
 * Use a greedy algorithm to select a subset of the positions in order
 * to decrease the problem size. The algorithm selects a set of points,
 * performs the interpolation and verifies the error at each positions.
 * The points with the greatest error are included, after which the
 * interpolation is performed again.
 * The algorithm keeps iterating until convergence.
 * The error is determined by looking at the distance between the interpolated
 * values and the actual values at the points not included in the interpolation.
 * The error is the Euclidean distance between the two values.
 */
void RBFInterpolationCoarsening::coarsen(
  const matrix & positions,
  const matrix & values,
  matrix & selectedPositions
  )
{
  assert( positions.rows() == values.rows() );

  // Initialize variables

  selectedPositions = positions.row( 0 );
  matrix selectedValues = values.row( 0 );
  matrix valuesInterpolation;
  fsi::vector error( positions.rows() );
  Eigen::VectorXi selectedPoints( 1 );
  selectedPoints << 0;

  int maxNbPoints = 2000;
  int minPoints = 100;
  maxNbPoints = std::min( 2000, static_cast<int>( positions.rows() ) );
  minPoints = std::min( 100, static_cast<int>( positions.rows() ) );

  assert( minPoints <= maxNbPoints );

  rbf.interpolate( selectedPositions, positions, selectedValues, valuesInterpolation );

  // Greedy algorithm

  for ( int i = 1; i < maxNbPoints; i++ )
  {
    // Evaluate the error

    for ( int j = 0; j < valuesInterpolation.rows(); j++ )
    {
      error( j ) = ( values.row( j ) - valuesInterpolation.row( j ) ).norm();
    }

    // Select the point with the largest error
    int index;
    error.maxCoeff( &index );

    double eps = std::pow( 2, -26 );
    double error = (values - valuesInterpolation).norm() / (values.norm() + eps);
    bool convergence = error < tol && i >= minPoints;

    if ( convergence )
      break;

    selectedPoints.conservativeResize( selectedPoints.rows() + 1 );
    selectedPoints( selectedPoints.rows() - 1 ) = index;

    selectedPositions.conservativeResize( selectedPositions.rows() + 1, selectedPositions.cols() );
    selectedPositions.row( selectedPositions.rows() - 1 ) = positions.row( index );

    selectedValues.conservativeResize( selectedValues.rows() + 1, selectedValues.cols() );
    selectedValues.row( selectedValues.rows() - 1 ) = values.row( index );

    rbf.interpolate( selectedPositions, positions, selectedValues, valuesInterpolation );
  }

  assert( selectedPositions.rows() <= positions.rows() );
}

/*
 * Perform the interpolation with radial basis functions.
 * First, select a subset of points with a greedy algorithm. Then,
 * perform the interpolation with the selected points.
 */
void RBFInterpolationCoarsening::interpolate(
  const matrix & positions,
  const matrix & positionsInterpolation,
  const matrix & values,
  matrix & valuesInterpolation
  )
{
  // Select points from positions matrix

  if ( this->positions.cols() == 0 && this->positions.rows() == 0 )
  {
    this->positions = positions;
    this->positionsInterpolation = positionsInterpolation;
  }

  matrix selectedPositions;

  coarsen( this->positions, values, selectedPositions );

  std::cout << "RBF interpolation with " << selectedPositions.rows()
            << " of " << positions.rows() << " points" << std::endl;

  // Perform interpolation with subset of positions

  rbf.interpolate( selectedPositions, this->positionsInterpolation, valuesInterpolation );
}
