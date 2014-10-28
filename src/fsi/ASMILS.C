
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ASMILS.H"

using namespace fsi;

ASMILS::ASMILS(
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

ASMILS::~ASMILS()
{}

void ASMILS::fixedUnderRelaxation(
  vector & xk,
  vector & zstar,
  vector & zk
  )
{
  Info << "Fixed relaxation post processing" << endl;

  xk += zstar - zk;
}

void ASMILS::performPostProcessing(
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

  fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
  fineResiduals.col( fineResiduals.cols() - 1 ) = xk;

  // Parameter extraction

  coarseModel->optimize( R, zstar, zk );

  if ( !coarseModel->allConverged() )
    Warning << "Surrogate model optimization process is not converged." << endl;

  coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
  coarseResiduals.col( coarseResiduals.cols() - 1 ) = zk;

  for ( int k = 0; k < maxIter - 1; k++ )
  {
    xkprev = xk;

    // Determine the number of columns used to calculate the mapping matrix J

    int nbCols = fineResiduals.cols() - 1;
    nbCols = std::max( nbCols, 0 );

    // Include information from previous time steps

    for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      nbCols += fineResidualsList.at( i ).cols() - 1;

    nbCols = std::min( static_cast<int>( xk.rows() ), nbCols );

    assert( nbCols <= xk.rows() );

    // Use fixed under relaxation to startup the IQNILS algorithm
    if ( nbCols == 0 )
      fixedUnderRelaxation( xk, zstar, zk );

    if ( nbCols > 0 )
    {
      // IQNILS algorithm

      Info << "ASM-ILS post processing with " << nbCols << " cols for the Jacobian" << endl;

      // Construct the V and W matrices
      matrix V( xk.rows(), nbCols ), W( xk.rows(), nbCols );

      int nbColsCurrentTimeStep = std::max( static_cast<int>(fineResiduals.cols() - 1), 0 );
      nbColsCurrentTimeStep = std::min( nbColsCurrentTimeStep, nbCols );

      // Include information from previous iterations

      int colIndex = 0;

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        if ( colIndex >= V.cols() )
          continue;

        V.col( i ) = coarseResiduals.col( coarseResiduals.cols() - 2 - i ) - coarseResiduals.col( coarseResiduals.cols() - 1 );
        W.col( i ) = fineResiduals.col( fineResiduals.cols() - 2 - i ) - fineResiduals.col( fineResiduals.cols() - 1 );
        colIndex++;
      }

      // Include information from previous time steps

      for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      {
        assert( fineResidualsList.at( i ).cols() >= 2 );
        assert( coarseResidualsList.at( i ).cols() >= 2 );

        for ( int j = 0; j < fineResidualsList.at( i ).cols() - 1; j++ )
        {
          if ( colIndex >= V.cols() )
            continue;

          V.col( colIndex ) = coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 2 - j ) - coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 1 );
          W.col( colIndex ) = fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 2 - j ) - fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 1 );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      // Remove dependent columns of V and W

      int nbRemoveCols = 1;

      while ( nbRemoveCols > 0 )
      {
        nbRemoveCols = 0;

        Eigen::HouseholderQR<matrix> qr = V.householderQr();
        matrix QR_R = qr.matrixQR().triangularView <Eigen::Upper> ();
        Eigen::VectorXd diagonals = QR_R.diagonal();

        for ( int i = 0; i < diagonals.rows(); i++ )
        {
          if ( std::abs( diagonals[i] ) < singularityLimit )
          {
            // Remove the column from V and W
            removeColumnFromMatrix( V, i - nbRemoveCols );
            removeColumnFromMatrix( W, i - nbRemoveCols );

            nbRemoveCols++;
          }
        }

        if ( nbRemoveCols )
          Info << "ASM-ILS: remove " << nbRemoveCols << " columns from the Jacobian matrices" << endl;
      }

      assert( V.cols() == W.cols() );

      if ( V.cols() == 0 )
        fixedUnderRelaxation( xk, zstar, zk );

      if ( V.cols() > 0 )
      {
        // SVD decomposition
        Eigen::JacobiSVD<matrix, Eigen::FullPivHouseholderQRPreconditioner> svd( V, Eigen::ComputeFullU | Eigen::ComputeFullV );
        vector c = svd.solve( zstar - zk );

        // Update solution x
        xk += W * c + V * c + zk - zstar;
      }
    }

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );

    fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
    fineResiduals.col( fineResiduals.cols() - 1 ) = xk;

    // Parameter extraction

    coarseModel->optimize( R, zstar, zk );

    if ( !coarseModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
    coarseResiduals.col( coarseResiduals.cols() - 1 ) = zk;

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
      break;

    if ( k == maxIter - 2 )
      break;
  }

  iterationsConverged();
}

void ASMILS::removeColumnFromMatrix(
  matrix & A,
  int col
  )
{
  for ( int j = col; j < A.cols() - 1; j++ )
    A.col( j ) = A.col( j + 1 );

  A.conservativeResize( A.rows(), A.cols() - 1 );
}
