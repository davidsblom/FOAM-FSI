
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AggressiveSpaceMapping.H"

using namespace fsi;

AggressiveSpaceMapping::AggressiveSpaceMapping(
  shared_ptr<SurrogateModel> fineModel,
  shared_ptr<SurrogateModel> coarseModel,
  int maxIter,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex,
  double singularityLimit
  )
  :
  SpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit )
{}

AggressiveSpaceMapping::~AggressiveSpaceMapping()
{}

void AggressiveSpaceMapping::performPostProcessing(
  const vector & y,
  const vector & x0,
  vector & xk,
  bool residualCriterium
  )
{
  assert( x0.rows() == xk.rows() );

  // Initialize variables

  int m = y.rows();
  vector yk( m ), output( m ), R( m ), zstar( m ), zk( m ), xkprev( m );
  xk = x0;
  xkprev = x0;
  coarseResiduals.resize( m, 0 );
  coarseResiduals.setZero();
  fineResiduals.resize( m, 0 );
  fineResiduals.setZero();

  // Determine optimum of coarse model zstar
  if ( residualCriterium )
  {
    assert( y.norm() < 1.0e-14 );
    coarseModel->optimize( x0, zstar );
  }

  if ( !residualCriterium )
    coarseModel->optimize( y, x0, zstar );

  if ( !coarseModel->allConverged() )
    Warning << "Surrogate model optimization process is not converged." << endl;

  xk = zstar;

  // Fine model evaluation

  fineModel->evaluate( xk, output, R );

  // Check convergence criteria
  if ( isConvergence( xk, xkprev, residualCriterium ) )
  {
    iterationsConverged();
    return;
  }

  // Parameter extraction

  coarseModel->optimize( R, zstar, zk );

  if ( !coarseModel->allConverged() )
    Warning << "Surrogate model optimization process is not converged." << endl;

  coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
  coarseResiduals.col( coarseResiduals.cols() - 1 ) = zk;

  fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
  fineResiduals.col( fineResiduals.cols() - 1 ) = xk;

  for ( int k = 0; k < maxIter - 1; k++ )
  {
    matrix J = Eigen::MatrixXd::Identity( m, m );

    // Determine the number of columns used to calculate the mapping matrix J

    int nbCols = k;
    int nbColsCurrentTimeStep = nbCols;
    int colIndex = 0;

    // Include information from previous time steps

    for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      nbCols += fineResidualsList.at( i ).cols() - 1;

    if ( nbCols > 0 )
    {
      Info << "Aggressive space mapping with ";
      Info << nbCols;
      Info << " cols for the Jacobian" << endl;

      // Include information from previous time steps

      for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      {
        assert( fineResidualsList.at( i ).cols() >= 2 );
        assert( coarseResidualsList.at( i ).cols() >= 2 );

        for ( int j = 0; j < fineResidualsList.at( i ).cols() - 1; j++ )
        {
          colIndex++;

          vector deltax, deltaz;

          deltax = fineResidualsList.at( i ).col( j + 1 ) - fineResidualsList.at( i ).col( j );
          deltaz = coarseResidualsList.at( i ).col( j + 1 ) - coarseResidualsList.at( i ).col( j );

          if ( deltax.norm() >= singularityLimit )
          {
            // Sherman–Morrison formula
            J += (deltax - J * deltaz) / (deltax.transpose() * J * deltaz) * (deltax.transpose() * J);
          }
        }
      }

      // Include information from previous iterations

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        colIndex++;

        vector deltax, deltaz;

        deltax = fineResiduals.col( i + 1 ) - fineResiduals.col( i );
        deltaz = coarseResiduals.col( i + 1 ) - coarseResiduals.col( i );

        if ( deltax.norm() >= singularityLimit )
        {
          // Sherman–Morrison formula
          J += (deltax - J * deltaz) / (deltax.transpose() * J * deltaz) * (deltax.transpose() * J);
        }
      }

      assert( colIndex == nbCols );
    }

    xkprev = xk;

    xk += J * (zstar - zk);

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
      break;

    // Parameter extraction

    coarseModel->optimize( R, zstar, zk );

    if ( !coarseModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
    coarseResiduals.col( coarseResiduals.cols() - 1 ) = zk;

    fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
    fineResiduals.col( fineResiduals.cols() - 1 ) = xk;
  }

  iterationsConverged();
}
