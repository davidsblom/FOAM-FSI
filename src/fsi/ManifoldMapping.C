
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
  double singularityLimit,
  bool updateJacobian
  )
  :
  SpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ),
  updateJacobian( updateJacobian ),
  initialSolutionCoarseModel( false ),
  scaling( false ),
  iter( 0 ),
  scalingFactors( Eigen::VectorXd::Ones( 2 ) ),
  sizeVar0( 0 ),
  sizeVar1( 0 )
{}

ManifoldMapping::ManifoldMapping(
  shared_ptr<SurrogateModel> fineModel,
  shared_ptr<SurrogateModel> coarseModel,
  int maxIter,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex,
  double singularityLimit,
  bool updateJacobian,
  bool initialSolutionCoarseModel
  )
  :
  SpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ),
  updateJacobian( updateJacobian ),
  initialSolutionCoarseModel( initialSolutionCoarseModel ),
  scaling( false ),
  iter( 0 ),
  scalingFactors( Eigen::VectorXd::Ones( 2 ) ),
  sizeVar0( 0 ),
  sizeVar1( 0 )
{}

ManifoldMapping::~ManifoldMapping()
{}

void ManifoldMapping::applyScaling( vector & vec )
{
  vec.head( sizeVar0 ) /= scalingFactors( 0 );
  vec.tail( sizeVar1 ) /= scalingFactors( 1 );
}

void ManifoldMapping::applyScaling( matrix & mat )
{
  mat.topLeftCorner( sizeVar0, mat.cols() ).array() /= scalingFactors( 0 );
  mat.bottomLeftCorner( sizeVar1, mat.cols() ).array() /= scalingFactors( 1 );
}

void ManifoldMapping::determineScalingFactors( const vector & output )
{
  std::shared_ptr<ImplicitMultiLevelFsiSolver> fsi;
  fsi = std::dynamic_pointer_cast<ImplicitMultiLevelFsiSolver>( fineModel );

  if ( fsi )
    scaling = fsi->fsi->parallel;

  if ( scaling && timeIndex <= reuseInformationStartingFromTimeIndex )
  {
    // Use scaling if the fluid and solid are coupled in parallel
    sizeVar0 = fsi->fsi->solidSolver->couplingGridSize * fsi->fsi->solid->dim;
    sizeVar1 = fsi->fsi->fluidSolver->couplingGridSize * fsi->fsi->fluid->dim;

    assert( output.rows() == sizeVar0 + sizeVar1 );

    scalingFactors( 0 ) = output.head( sizeVar0 ).norm();
    scalingFactors( 1 ) = output.tail( sizeVar1 ).norm();

    for ( int i = 0; i < scalingFactors.rows(); i++ )
      if ( std::abs( scalingFactors( i ) ) < 1.0e-13 )
        scalingFactors( i ) = 1;

    Info << "Parallel coupling of fluid and solid solvers with scaling factors ";
    Info << scalingFactors( 0 ) << " and " << scalingFactors( 1 ) << endl;

    // Reset the mapping matrix Tprev since a different scaling factor is used
    Tkprev.resize( 0, 0 );
  }
}

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
  matrix Tk;

  if ( timeIndex == 0 || initialSolutionCoarseModel )
  {
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
  }

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

  determineScalingFactors( output );

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

    for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      nbCols += fineResidualsList.at( i ).cols() - 1;

    nbCols = std::min( nbCols, n );

    // Update the design specification yk

    vector alpha = fineResiduals.col( k ) - y;
    yk = coarseResiduals.col( k );

    // Apply scaling
    if ( scaling )
    {
      applyScaling( alpha );
      applyScaling( yk );
    }

    if ( nbCols == 0 )
    {
      // Update the design specification yk

      if ( updateJacobian && Tkprev.rows() == m )
        yk -= Tkprev * alpha;
      else
        yk -= alpha;
    }

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

      for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      {
        assert( fineResidualsList.at( i ).cols() >= 2 );

        for ( unsigned j = 0; j < fineResidualsList.at( i ).cols() - 1; j++ )
        {
          if ( colIndex >= DeltaF.cols() )
            continue;

          DeltaF.col( colIndex ) = fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 1 ) - fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 2 - j );
          DeltaC.col( colIndex ) = coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 1 ) - coarseResidualsList.at( i ).col( coarseResidualsList.at( i ).cols() - 2 - j );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      // Apply scaling
      if ( scaling )
      {
        applyScaling( DeltaF );
        applyScaling( DeltaC );
      }

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

        Eigen::JacobiSVD<matrix> svd_C( DeltaC, Eigen::ComputeThinU | Eigen::ComputeThinV );
        Eigen::JacobiSVD<matrix> svd_F( DeltaF, Eigen::ComputeThinU | Eigen::ComputeThinV );

        matrix pseudoSigma_F = svd_F.singularValues().asDiagonal();

        for ( int i = 0; i < pseudoSigma_F.cols(); i++ )
          pseudoSigma_F( i, i ) = 1.0 / pseudoSigma_F( i, i );

        matrix pseudoDeltaF = svd_F.matrixV() * pseudoSigma_F * svd_F.matrixU().transpose();

        U_F = svd_F.matrixU();
        U_C = svd_C.matrixU();

        if ( !updateJacobian )
        {
          vector beta = U_F * (U_F.transpose() * alpha);

          yk -= alpha;
          yk -= DeltaC * (pseudoDeltaF * alpha);
          yk += U_C * ( U_C.transpose() * (alpha - beta) );
          yk += beta;
        }

        if ( updateJacobian )
        {
          matrix I = Eigen::MatrixXd::Identity( DeltaF.rows(), DeltaF.rows() );

          if ( Tkprev.rows() == DeltaF.rows() )
            Tk = Tkprev + (DeltaC - Tkprev * DeltaF) * pseudoDeltaF;
          else
            Tk = DeltaC * pseudoDeltaF + ( I - U_C * U_C.transpose() ) * ( I - U_F * U_F.transpose() );

          yk -= Tk * alpha;
        }
      }

      if ( DeltaF.cols() == 0 )
      {
        // Update the design specification yk

        if ( updateJacobian && Tkprev.rows() == m )
          yk -= Tkprev * alpha;
        else
          yk -= alpha;
      }
    }

    // Apply scaling
    if ( scaling )
    {
      yk.head( sizeVar0 ) *= scalingFactors( 0 );
      yk.tail( sizeVar1 ) *= scalingFactors( 1 );
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

    determineScalingFactors( output );

    iter++;

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
    {
      assert( fineModel->allConverged() );

      if ( updateJacobian && Tk.cols() > 0 && timeIndex >= reuseInformationStartingFromTimeIndex )
        Tkprev = Tk;

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
