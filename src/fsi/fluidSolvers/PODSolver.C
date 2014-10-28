
/*
 * Class
 *   FluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PODSolver.H"

PODSolver::PODSolver(
  shared_ptr<ImplicitMultiLevelFsiSolver> fsi,
  int participantId
  )
  :
  BaseMultiLevelSolver( 5, 3, 0 ),
  fsi( fsi ),
  participantId( participantId )
{
  assert( participantId == 0 || participantId == 1 );

  matrix readPositions;
  getReadPositions( readPositions );
  N = readPositions.rows();
  dim = readPositions.cols();
  data.resize( N, dim );
  data.setZero();
}

PODSolver::~PODSolver()
{}

void PODSolver::finalizeTimeStep()
{
  assert( init );

  init = false;
}

void PODSolver::getReadPositions( matrix & readPositions )
{
  if ( participantId == 0 )
    fsi->fsi->fluid->getReadPositions( readPositions );

  if ( participantId == 1 )
    fsi->fsi->solid->getReadPositions( readPositions );
}

void PODSolver::getWritePositions( matrix & writePositions )
{
  if ( participantId == 0 )
    fsi->fsi->fluid->getWritePositions( writePositions );

  if ( participantId == 1 )
    fsi->fsi->solid->getWritePositions( writePositions );
}

void PODSolver::initTimeStep()
{
  assert( !init );

  t = fsi->fsi->fluid->t;
  timeIndex = fsi->fsi->fluid->timeIndex;

  init = true;
}

bool PODSolver::isRunning()
{
  if ( participantId == 0 )
    return fsi->fsi->fluid->isRunning();

  if ( participantId == 1 )
    return fsi->fsi->solid->isRunning();

  return false;
}

void PODSolver::resetSolution()
{}

void PODSolver::solve(
  const matrix & input,
  matrix & output
  )
{
  // Count the number of snapshots

  int nbSnapshots = outputList.size();
  nbSnapshots = min( nbSnapshots, outputList.at( 0 ).rows() );

  assert( nbSnapshots >= 1 );

  // Build the snapshot matrix P

  matrix P( outputList.at( 0 ).rows(), nbSnapshots );
  matrix W( inputList.at(0).rows(), nbSnapshots);

  for ( unsigned i = 0; i < P.cols(); i++ )
  {
    P.col( i ) = outputList.at( outputList.size() - P.cols() + i );
    W.col( i ) = inputList.at( outputList.size() - P.cols() + i );
  }

  // Build first POD model using a singular value decomposition

  Eigen::JacobiSVD<matrix> svd_P( P, Eigen::ComputeThinU | Eigen::ComputeThinV );

  // POD modes

  matrix Phi = svd_P.matrixU();

  // Compute SVD decomposition of POD modes in order to solve for the
  // coefficient matrix T

  Eigen::JacobiSVD<matrix, Eigen::FullPivHouseholderQRPreconditioner> svd_Phi( Phi, Eigen::ComputeFullU | Eigen::ComputeFullV );

  // Solve for the coefficients T

  matrix T( nbSnapshots, Phi.cols() );

  for ( int i = 0; i < T.cols(); i++ )
    T.col( i ) = svd_Phi.solve( P.col( i ) );

  // Verify coefficient matrix T
  // Each column of P can be reconstructed, using the POD modes as a summation
  // of the coefficients T_i * Phi_i

  for ( int i = 0; i < P.cols(); i++ )
  {
    fsi::vector tmp( P.col( i ) );
    tmp.setZero();

    for ( int j = 0; j < T.rows(); j++ )
      tmp += T( j, i ) * Phi.col( j );

    assert( tmp.norm() - P.col( i ).norm() < 1.0e-13 );
  }

  // Build radial basis function network in order to interpolate the coefficients

  RBFInterpolation rbf;

  matrix values, valuesInterpolation;
  matrix positions, positionsInterpolation;

  positions = W.transpose();
  positionsInterpolation = input.transpose();

  values = T;
  valuesInterpolation = values;

  assert( positions.cols() == positionsInterpolation.cols() );

  rbf.compute( positions, positionsInterpolation );
  rbf.interpolate( values, valuesInterpolation );

  // Construct output from interpolated coefficients

  output.resize( P.rows(), 1 );
  output.setZero();

  for ( int i = 0; i < valuesInterpolation.rows(); i++ )
    output += valuesInterpolation( 0, i ) * Phi.col( i );
}
