
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "OutputSpaceMapping.H"

using namespace fsi;

OutputSpaceMapping::OutputSpaceMapping(
  shared_ptr<SurrogateModel> fineModel,
  shared_ptr<SurrogateModel> surrogateModel,
  int maxIter,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex,
  double singularityLimit,
  int order
  )
  :
  SpaceMapping( fineModel, surrogateModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ),
  surrogateModel( surrogateModel ),
  order( order )
{
  assert( surrogateModel );
  assert( order == 0 || order == 1 || order == 2 );
}

OutputSpaceMapping::~OutputSpaceMapping()
{}

void OutputSpaceMapping::iterationsConverged()
{
  SpaceMapping::iterationsConverged();

  // Save input/output information for next time step
  if ( nbReuse > 0 && sols.size() >= 2 && timeIndex >= reuseInformationStartingFromTimeIndex )
    solsList.push_front( sols );

  // Remove the last item from the sols list and solutions
  // list in order to ensure that at maximum nbReuse time steps
  // are included.
  while ( static_cast<int>( solsList.size() ) > nbReuse )
    solsList.pop_back();

  sols.clear();

  assert( static_cast<int>( solsList.size() ) <= nbReuse );
  assert( solsList.size() == coarseResidualsList.size() );
  assert( solsList.size() == fineResidualsList.size() );
}

void OutputSpaceMapping::performPostProcessing(
  const vector & y,
  const vector & x0,
  vector & xk,
  bool residualCriterium
  )
{
  assert( x0.rows() == xk.rows() );

  // Initialize variables

  int m = y.rows();
  int n = x0.rows();
  vector yk( m ), output( m ), R( m );
  xk = x0;
  vector xkprev = x0;
  output.setZero();
  R.setZero();
  coarseResiduals.resize( m, 1 );
  coarseResiduals.setZero();
  fineResiduals.resize( m, 1 );
  fineResiduals.setZero();

  if ( timeIndex == 0 )
  {
    // Determine optimum of coarse model xstar
    if ( residualCriterium )
    {
      assert( y.norm() < 1.0e-14 );
      surrogateModel->optimize( x0, xk );
    }

    if ( !residualCriterium )
      surrogateModel->optimize( y, x0, xk );

    if ( !surrogateModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;
  }

  assert( xk.rows() == n );
  assert( x0.rows() == n );

  surrogateModel->setUseJacobian( false );

  // Initialize coarse model and fine model responses

  // Coarse model evaluation

  surrogateModel->evaluate( xk, output, R );
  coarseResiduals.col( 0 ) = R;
  assert( xk.rows() == n );
  assert( output.rows() == m );
  assert( R.rows() == m );

  // Fine model evaluation

  fineModel->evaluate( xk, output, R );
  fineResiduals.col( 0 ) = R;
  assert( output.rows() == m );
  assert( R.rows() == m );

  sols.push_back( xk );

  // Check convergence criteria
  if ( isConvergence( xk, xkprev, residualCriterium ) )
  {
    assert( fineModel->allConverged() );
    iterationsConverged();
    return;
  }

  assert( xk.rows() == n );

  for ( int k = 0; k < maxIter - 1; k++ )
  {
    // Determine the number of columns used to calculate the mapping matrix

    int nbCols = k;
    int nbColsCurrentTimeStep = nbCols;

    // Include information from previous time steps

    for ( unsigned i = 0; i < solsList.size(); i++ )
      nbCols += solsList.at( i ).size() - 1;

    // Update the design specification yk
    yk = coarseResiduals.col( k ) - (fineResiduals.col( k ) - y);

    if ( nbCols > 0 && order == 1 )
    {
      assert( fineResiduals.cols() == coarseResiduals.cols() );
      assert( static_cast<int>( sols.size() ) == fineResiduals.cols() );
      assert( solsList.size() == fineResidualsList.size() );
      assert( solsList.size() == coarseResidualsList.size() );

      // Initialize mapping matrix
      matrix J = -Eigen::MatrixXd::Identity( m, m );
      fsi::vector d, dprev, deltad, deltax;
      int colIndex = 0;

      Info << "Output space mapping with ";
      Info << nbCols;
      Info << " cols for the Jacobian" << endl;

      // Include information from previous time steps

      for ( unsigned i = solsList.size(); i-- > 0; )
      {
        for ( unsigned j = 0; j < solsList.at( i ).size() - 1; j++ )
        {
          colIndex++;

          d = fineResidualsList.at( i ).col( j + 1 ) - coarseResidualsList.at( i ).col( j + 1 );
          dprev = fineResidualsList.at( i ).col( j ) - coarseResidualsList.at( i ).col( j );
          deltad = d - dprev;
          deltax = solsList.at( i ).at( j + 1 ) - solsList.at( i ).at( j );

          // Broyden update for the Jacobian matrix

          if ( deltax.norm() < singularityLimit )
            continue;

          J += (deltad - J * deltax) / deltax.squaredNorm() * deltax.transpose();
        }
      }

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        colIndex++;

        d = fineResiduals.col( i + 1 ) - coarseResiduals.col( i + 1 );
        dprev = fineResiduals.col( i ) - coarseResiduals.col( i );
        deltad = d - dprev;
        deltax = sols.at( i + 1 ) - sols.at( i );

        // Broyden update for the Jacobian matrix

        if ( deltax.norm() < singularityLimit )
          continue;

        J += (deltad - J * deltax) / deltax.squaredNorm() * deltax.transpose();
      }

      assert( colIndex == nbCols );

      surrogateModel->setUseJacobian( true );
      surrogateModel->setSurrogateData( xk, J );
    }

    if ( nbCols > 0 && order == 2 )
    {
      nbCols = std::min( nbCols, n );

      matrix DeltaF( m, nbCols ), DeltaX( m, nbCols );

      int colIndex = 0;

      // Include information from current time step

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        if ( colIndex >= DeltaF.cols() )
          continue;

        DeltaF.col( colIndex ) = fineResiduals.col( k ) - coarseResiduals.col( k );
        DeltaF.col( colIndex ) -= fineResiduals.col( k - 1 - i ) - coarseResiduals.col( k - 1 - i );
        DeltaX.col( colIndex ) = sols.at( k ) - sols.at( k - 1 - i );
        colIndex++;
      }

      // Include information from previous time steps

      for ( unsigned i = 0; i < solsList.size(); i++ )
      {
        for ( unsigned j = 0; j < solsList.at( i ).size() - 1; j++ )
        {
          if ( colIndex >= DeltaF.cols() )
            continue;

          DeltaF.col( colIndex ) = fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 1 ) - coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 1 );

          DeltaF.col( colIndex ) -= fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 2 - j ) - coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 2 - j );

          DeltaX.col( colIndex ) = solsList.at( i ).at( solsList.at( i ).size() - 1 ) - solsList.at( i ).at( solsList.at( i ).size() - 2 - j );

          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      // Remove dependent columns of DeltaC and DeltaF

      int nbRemoveCols = 1;

      while ( nbRemoveCols > 0 )
      {
        nbRemoveCols = 0;

        if ( DeltaX.cols() == 1 )
          break;

        // Calculate singular value decomposition with Eigen

        Eigen::JacobiSVD<matrix> svd( DeltaX, Eigen::ComputeThinU | Eigen::ComputeThinV );

        vector singularValues = svd.singularValues();

        for ( int i = 0; i < singularValues.rows(); i++ )
        {
          if ( std::abs( singularValues( i ) ) <= singularityLimit && DeltaX.cols() > 1 )
          {
            // Remove the column from DeltaC and DeltaF
            removeColumnFromMatrix( DeltaX, i - nbRemoveCols );
            removeColumnFromMatrix( DeltaF, i - nbRemoveCols );

            nbRemoveCols++;
          }
        }

        if ( nbRemoveCols )
          Info << "Output space mapping: remove " << nbRemoveCols << " columns from the Jacobian matrices" << endl;
      }

      matrix I = Eigen::MatrixXd::Identity( m, m );

      Eigen::JacobiSVD<matrix> svd( DeltaX, Eigen::ComputeThinU | Eigen::ComputeThinV );

      matrix pseudoDeltaX = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();

      matrix J = -I + (DeltaF + DeltaX) * pseudoDeltaX;

      surrogateModel->setUseJacobian( true );
      surrogateModel->setSurrogateData( xk, J );
    }

    xkprev = xk;

    // Update the fine model optimum
    output.resize( n );
    surrogateModel->optimize( yk, xk, output );
    assert( output.rows() == n );
    assert( yk.rows() == m );
    assert( xk.rows() == n );

    if ( !surrogateModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    xk = output;

    surrogateModel->setUseJacobian( false );

    // Coarse model evaluation
    output.resize( m );
    surrogateModel->evaluate( xk, output, R );
    assert( xk.rows() == n );
    assert( output.rows() == m );
    assert( R.rows() == m );
    coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
    coarseResiduals.col( coarseResiduals.cols() - 1 ) = R;

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );
    assert( xk.rows() == n );
    assert( output.rows() == m );
    assert( R.rows() == m );
    fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
    fineResiduals.col( fineResiduals.cols() - 1 ) = R;

    sols.push_back( xk );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
    {
      assert( fineModel->allConverged() );
      break;
    }
  }

  iterationsConverged();
}

void OutputSpaceMapping::removeColumnFromMatrix(
  matrix & A,
  int col
  )
{
  for ( int j = col; j < A.cols() - 1; j++ )
    A.col( j ) = A.col( j + 1 );

  A.conservativeResize( A.rows(), A.cols() - 1 );
}
