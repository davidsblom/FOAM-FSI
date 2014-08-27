
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "BroydenPostProcessing.H"

using namespace fsi;

BroydenPostProcessing::BroydenPostProcessing(
  shared_ptr<MultiLevelFsiSolver> fsi,
  double initialRelaxation,
  int maxIter,
  int maxUsedIterations,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex
  )
  :
  PostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex )
{}

BroydenPostProcessing::~BroydenPostProcessing()
{}

void BroydenPostProcessing::performPostProcessing(
  const vector & x0,
  vector & xk
  )
{
  vector y( x0.rows() );
  y.setZero();
  performPostProcessing( y, x0, xk, true );
}

void BroydenPostProcessing::performPostProcessing(
  const vector & y,
  const vector & x0,
  vector & xk
  )
{
  performPostProcessing( y, x0, xk, false );
}

void BroydenPostProcessing::performPostProcessing(
  const vector & y,
  const vector & x0,
  vector & xk,
  bool residualCriterium
  )
{
  assert( xk.rows() > 0 );
  assert( fsi->fluid->init );
  assert( fsi->solid->init );
  assert( y.rows() == x0.rows() );
  assert( y.rows() == xk.rows() );

  // Initialize variables
  vector xkprev = x0;
  xk = x0;
  xkprev.setZero();
  residuals.clear();
  sols.clear();
  vector yk = y;

  // Fsi evaluation
  vector output( xk.rows() ), R( xk.rows() );
  output.setZero();
  R.setZero();

  fsi->evaluate( x0, output, R );

  assert( x0.rows() == output.rows() );
  assert( x0.rows() == R.rows() );

  // Save output and residual
  residuals.push_back( R );
  sols.push_back( x0 );

  // Check convergence criteria
  if ( isConvergence( xk, xkprev, residualCriterium ) )
  {
    iterationsConverged();
    return;
  }

  for ( int iter = 0; iter < maxIter - 1; iter++ )
  {
    xkprev = xk;

    // Determine the number of columns of the V and W matrices
    int nbCols = residuals.size() - 1;
    nbCols = std::max( nbCols, 0 );

    // Include information from previous optimization solves
    for ( unsigned i = 0; i < solsList.size(); i++ )
      nbCols += solsList.at( i ).size() - 1;

    // Include information from previous time steps
    for ( unsigned i = 0; i < solsTimeList.size(); i++ )
      for ( unsigned j = 0; j < solsTimeList.at( i ).size(); j++ )
        nbCols += solsTimeList.at( i ).at( j ).size() - 1;

    int indexCounter = nbCols;
    nbCols = std::min( static_cast<int>( xk.rows() ), nbCols );
    nbCols = std::min( nbCols, maxUsedIterations );

    if ( nbCols == 0 )
    {
      // Use fixed under relaxation to startup the IQNILS algorithm

      Info << "Fixed relaxation post processing with factor " << initialRelaxation << endl;

      // Update solution x
      vector dx = initialRelaxation * (R - yk);

      xk += dx;
    }

    if ( nbCols > 0 )
    {
      Info << "Broyden post processing with " << nbCols << " cols for the Jacobian" << endl;

      matrix J = -Eigen::MatrixXd::Identity( R.rows(), R.rows() );

      int nbColsCurrentTimeStep = std::max( static_cast<int>(sols.size() - 1), 0 );

      int colIndex = 0;

      // Include information from previous time steps
      for ( unsigned i = residualsTimeList.size(); i-- > 0; )
      {
        for ( unsigned j = residualsTimeList.at( i ).size(); j-- > 0; )
        {
          for ( unsigned k = residualsTimeList.at( i ).at( j ).size() - 1; k-- > 0; )
          {
            assert( residualsTimeList.at( i ).at( j ).size() >= 2 );
            assert( residualsTimeList.at( i ).size() == solsTimeList.at( i ).size() );
            assert( residualsTimeList.at( i ).at( j ).size() == solsTimeList.at( i ).at( j ).size() );

            indexCounter--;

            if ( indexCounter >= nbCols )
              continue;

            fsi::vector dx = solsTimeList.at( i ).at( j ).at( k + 1 ) - solsTimeList.at( i ).at( j ).at( k );
            fsi::vector dR = residualsTimeList.at( i ).at( j ).at( k + 1 ) - residualsTimeList.at( i ).at( j ).at( k );

            // Sherman–Morrison formula
            J += (dx - J * dR) / (dx.transpose() * J * dR) * (dx.transpose() * J);

            colIndex++;
          }
        }
      }

      // Include information from previous optimization solves

      for ( unsigned i = residualsList.size(); i-- > 0; )
      {
        for ( unsigned j = residualsList.at( i ).size() - 1; j-- > 0; )
        {
          assert( residualsList.at( i ).size() >= 2 );
          assert( residualsList.at( i ).size() == solsList.at( i ).size() );

          indexCounter--;

          if ( indexCounter >= nbCols )
            continue;

          fsi::vector dx = solsList.at( i ).at( j + 1 ) - solsList.at( i ).at( j );
          fsi::vector dR = residualsList.at( i ).at( j + 1 ) - residualsList.at( i ).at( j );

          // Sherman–Morrison formula
          J += (dx - J * dR) / (dx.transpose() * J * dR) * (dx.transpose() * J);

          colIndex++;
        }
      }

      // Include information from previous iterations

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        assert( sols.size() >= 2 );

        indexCounter--;

        if ( indexCounter >= nbCols )
          continue;

        fsi::vector dx = sols.at( i + 1 ) - sols.at( i );
        fsi::vector dR = residuals.at( i + 1 ) - residuals.at( i );

        // Sherman–Morrison formula
        J += (dx - J * dR) / (dx.transpose() * J * dR) * (dx.transpose() * J);

        colIndex++;
      }

      assert( indexCounter == 0 );
      assert( colIndex == nbCols );

      xk += J * (y - R);
    }

    // Fsi evaluation
    fsi->evaluate( xk, output, R );

    assert( x0.rows() == output.rows() );
    assert( x0.rows() == R.rows() );

    // Save output and residual
    residuals.push_back( R );
    sols.push_back( xk );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
      break;

    assert( sols.size() == residuals.size() );
    assert( sols.at( 0 ).rows() == residuals.at( 0 ).rows() );
    assert( fsi->iter <= maxIter );
  }

  iterationsConverged();
}
