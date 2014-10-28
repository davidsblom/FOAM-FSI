
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ManifoldMapping.H"

using namespace fsi;

ManifoldMapping::ManifoldMapping(
  shared_ptr<SurrogateModel> fineModel,
  shared_ptr<SurrogateModel> coarseModel,
  int maxIter,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex,
  double singularityLimit
  )
  :
  SpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ),
  iter( 0 )
{}

ManifoldMapping::~ManifoldMapping()
{}

void ManifoldMapping::performPostProcessing(
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
  iter = 0;

  // Determine optimum of coarse model xstar
  if ( residualCriterium )
  {
    assert( y.norm() < 1.0e-14 );
    coarseModel->optimize( x0, xk );
  }

  if ( !residualCriterium )
    coarseModel->optimize( y, x0, xk );

  if ( !coarseModel->allConverged() )
    Warning << "Surrogate model optimization process is not converged." << endl;

  assert( xk.rows() == n );
  assert( x0.rows() == n );

  // Initialize coarse model and fine model responses

  // Coarse model evaluation

  coarseModel->evaluate( xk, output, R );
  coarseResiduals.col( 0 ) = R;
  assert( xk.rows() == n );
  assert( output.rows() == m );
  assert( R.rows() == m );

  // Fine model evaluation

  fineModel->evaluate( xk, output, R );
  fineResiduals.col( 0 ) = R;
  assert( output.rows() == m );
  assert( R.rows() == m );

  iter++;

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
    xkprev = xk;

    // Determine the number of columns used to calculate the mapping matrix

    int nbCols = std::min( k, n );
    int nbColsCurrentTimeStep = nbCols;

    // Include information from previous time steps

    for ( std::deque<matrix>::iterator it = fineResidualsList.begin(); it != fineResidualsList.end(); ++it )
      nbCols += it->cols() - 1;

    nbCols = std::min( nbCols, n );

    if ( nbCols == 0 )
      // Update the design specification yk
      yk = coarseResiduals.col( k ) - (fineResiduals.col( k ) - y);

    if ( nbCols >= 1 )
    {
      // Update the mapping matrix

      matrix DeltaF( m, nbCols ), DeltaC( m, nbCols );

      int colIndex = 0;

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        DeltaF.col( i ) = fineResiduals.col( k ) - fineResiduals.col( k - 1 - i );
        DeltaC.col( i ) = coarseResiduals.col( k ) - coarseResiduals.col( k - 1 - i );
        colIndex++;
      }

      // Include information from previous time steps

      int tmpIndex = colIndex;

      for ( std::deque<matrix>::iterator it = fineResidualsList.begin(); it != fineResidualsList.end(); ++it )
      {
        assert( it->cols() >= 2 );

        for ( int i = 0; i < it->cols() - 1; i++ )
        {
          if ( colIndex >= DeltaF.cols() )
            continue;

          DeltaF.col( colIndex ) = it->col( it->cols() - 1 ) - it->col( it->cols() - 2 - i );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      colIndex = tmpIndex;

      for ( std::deque<matrix>::iterator it = coarseResidualsList.begin(); it != coarseResidualsList.end(); ++it )
      {
        assert( it->cols() >= 2 );

        for ( int i = 0; i < it->cols() - 1; i++ )
        {
          if ( colIndex >= DeltaC.cols() )
            continue;

          DeltaC.col( colIndex ) = it->col( it->cols() - 1 ) - it->col( it->cols() - 2 - i );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      Info << "Manifold mapping with ";
      Info << nbCols;
      Info << " cols for the Jacobian" << endl;

      // Initialize variables for singular value decomposition

      vector S_F, S_C;
      matrix V_F, U_F, V_C, U_C, Sigma_F, pseudoSigma_F;

      // Remove dependent columns of DeltaC and DeltaF

      int nbRemoveCols = 1;

      while ( nbRemoveCols > 0 )
      {
        nbRemoveCols = 0;

        if ( DeltaF.cols() == 0 )
          break;

        // Calculate singular value decomposition with Eigen

        Eigen::JacobiSVD<matrix> svd( DeltaF, Eigen::ComputeThinU | Eigen::ComputeThinV );

        vector singularValues = svd.singularValues();

        for ( int i = 0; i < singularValues.rows(); i++ )
        {
          if ( std::abs( singularValues( i ) ) <= singularityLimit )
          {
            // Remove the column from DeltaC and DeltaF
            removeColumnFromMatrix( DeltaC, i - nbRemoveCols );
            removeColumnFromMatrix( DeltaF, i - nbRemoveCols );

            nbRemoveCols++;
          }
        }

        if ( nbRemoveCols )
          Info << "Manifold mapping: remove " << nbRemoveCols << " columns from the Jacobian matrices" << endl;
      }

      assert( DeltaC.cols() == DeltaF.cols() );

      if ( DeltaF.cols() > 0 )
      {
        // Calculate singular value decomposition with Eigen

        Eigen::JacobiSVD<matrix> svd_F( DeltaF, Eigen::ComputeThinU | Eigen::ComputeThinV );
        Eigen::JacobiSVD<matrix> svd_C( DeltaC, Eigen::ComputeThinU | Eigen::ComputeThinV );

        matrix pseudoSigma_F = svd_F.singularValues().asDiagonal();

        for ( int i = 0; i < pseudoSigma_F.cols(); i++ )
          pseudoSigma_F( i, i ) = 1.0 / pseudoSigma_F( i, i );

        matrix pseudoDeltaF = svd_F.matrixV() * pseudoSigma_F * svd_F.matrixU().transpose();

        // Update the design specification yk

        U_F = svd_F.matrixU();
        U_C = svd_C.matrixU();

        // matrix I = Eigen::MatrixXd::Identity( m, m );
        // matrix Tk = DeltaC * pseudoDeltaF + ( I - U_C * U_C.transpose() ) * ( I - U_F * U_F.transpose() );
        // yk = coarseResiduals.col( k ) - Tk * ( fineResiduals.col(k) - y );

        yk = coarseResiduals.col( k );
        yk -= DeltaC * ( pseudoDeltaF * (fineResiduals.col( k ) - y) );
        yk += U_F * ( U_F.transpose() * (fineResiduals.col( k ) - y) );
        yk += U_C * ( U_C.transpose() * (fineResiduals.col( k ) - y) );
        yk -= fineResiduals.col( k ) - y;
        yk -= U_C * ( U_C.transpose() * ( U_F * ( U_F.transpose() * (fineResiduals.col( k ) - y) ) ) );
      }

      if ( DeltaF.cols() == 0 )
      {
        // Update the design specification yk
        yk = coarseResiduals.col( k ) - (fineResiduals.col( k ) - y);
      }
    }

    // Update the fine model optimum
    output.resize( n );
    coarseModel->optimize( yk, xk, output );
    assert( output.rows() == n );
    assert( yk.rows() == m );
    assert( xk.rows() == n );

    if ( !coarseModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    xk = output;

    // Coarse model evaluation
    output.resize( m );
    coarseModel->evaluate( xk, output, R );
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

    iter++;

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
    {
      assert( fineModel->allConverged() );
      break;
    }
  }

  iterationsConverged();
}

void ManifoldMapping::removeColumnFromMatrix(
  matrix & A,
  int col
  )
{
  for ( int j = col; j < A.cols() - 1; j++ )
    A.col( j ) = A.col( j + 1 );

  A.conservativeResize( A.rows(), A.cols() - 1 );
}
