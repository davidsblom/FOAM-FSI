
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ResponseParameterMapping.H"

using namespace fsi;

ResponseParameterMapping::ResponseParameterMapping(
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

ResponseParameterMapping::~ResponseParameterMapping()
{}

void ResponseParameterMapping::iterationsConverged()
{
  SpaceMapping::iterationsConverged();

  assert( xktildeList.size() == xkList.size() );

  // Save input/output information for next time step
  if ( nbReuse > 0 && xktildeList.size() >= 2 && timeIndex >= reuseInformationStartingFromTimeIndex )
  {
    xktildePreviousTimeStepsList.push_front( xktildeList );
    xkPreviousTimeStepsList.push_front( xkList );
  }

  // Remove the last item from the sols list and solutions
  // list in order to ensure that at maximum nbReuse time steps
  // are included.
  while ( static_cast<int>( xktildePreviousTimeStepsList.size() ) > nbReuse )
  {
    xktildePreviousTimeStepsList.pop_back();
    xkPreviousTimeStepsList.pop_back();
  }

  xktildeList.clear();
  xkList.clear();

  assert( static_cast<int>( xktildePreviousTimeStepsList.size() ) <= nbReuse );
  assert( xkPreviousTimeStepsList.size() == coarseResidualsList.size() );
  assert( xktildePreviousTimeStepsList.size() == fineResidualsList.size() );
  assert( xktildePreviousTimeStepsList.size() == xkPreviousTimeStepsList.size() );
}

void ResponseParameterMapping::performPostProcessing(
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
  matrix B( x0.rows(), x0.rows() );
  B.setZero();
  vector yk( m ), output( m ), R( m );
  vector xktilde = x0;
  xk = x0;
  vector xkprev = x0;
  output.setZero();
  R.setZero();
  coarseResiduals.resize( m, 1 );
  coarseResiduals.setZero();
  fineResiduals.resize( m, 1 );
  fineResiduals.setZero();
  iter = 0;
  xktildeList.clear();
  xkList.clear();
  mappedCoarseList.clear();

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

  // Initialize coarse model and fine model responses

  // Fine model evaluation

  fineModel->evaluate( xk, output, R );
  fineResiduals.col( 0 ) = R;

  iter++;

  // Parameter extraction
  Info << nl << "Response and parameter mapping: parameter extraction process" << endl;
  coarseModel->optimize( R, xk, xktilde );
  xktildeList.push_back( xktilde );
  xkList.push_back( xk );

  // Coarse model evaluation

  coarseModel->evaluate( xktilde, output, R );
  coarseResiduals.col( 0 ) = R;

  // Check convergence criteria
  if ( isConvergence( xk, xkprev, residualCriterium ) )
  {
    iterationsConverged();
    return;
  }

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
      // Update the mapping matrix B

      matrix DeltaXTilde( m, nbCols ), DeltaX( m, nbCols );

      int colIndex = 0;

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        DeltaXTilde.col( i ) = xktildeList.at( k - 1 - i ) - xktildeList.at( k );
        DeltaX.col( i ) = xkList.at( k - 1 - i ) - xkList.at( k );
        colIndex++;
      }

      // Include information from previous time steps for the mapping matrix B

      assert( xkPreviousTimeStepsList.size() == xktildePreviousTimeStepsList.size() );
      assert( static_cast<int> ( fineResidualsList.size() ) == static_cast<int> ( xkPreviousTimeStepsList.size() ) );

      for ( unsigned i = 0; i < xkPreviousTimeStepsList.size(); i++ )
      {
        assert( xkPreviousTimeStepsList.at( i ).size() == xktildePreviousTimeStepsList.at( i ).size() );
        assert( static_cast<int>( xkPreviousTimeStepsList.at( i ).size() ) == fineResidualsList.at( i ).cols() );

        for ( unsigned j = 0; j < xkPreviousTimeStepsList.at( i ).size() - 1; j++ )
        {
          if ( colIndex >= DeltaXTilde.cols() )
            continue;

          DeltaXTilde.col( colIndex ) = xktildePreviousTimeStepsList.at( i ).at( xktildePreviousTimeStepsList.at( i ).size() - 2 - j ) - xktildePreviousTimeStepsList.at( i ).at( xktildePreviousTimeStepsList.at( i ).size() - 1 );
          DeltaX.col( colIndex ) = xkPreviousTimeStepsList.at( i ).at( xkPreviousTimeStepsList.at( i ).size() - 2 - j ) - xkPreviousTimeStepsList.at( i ).at( xkPreviousTimeStepsList.at( i ).size() - 1 );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      // Remove dependent columns of DeltaX and DeltaXTilde

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
          if ( std::abs( singularValues( i ) ) <= singularityLimit )
          {
            // Remove the column from DeltaC and DeltaF
            removeColumnFromMatrix( DeltaX, i - nbRemoveCols );
            removeColumnFromMatrix( DeltaXTilde, i - nbRemoveCols );

            nbRemoveCols++;
          }
        }

        if ( nbRemoveCols )
          Info << "Response and parameter mapping: remove " << nbRemoveCols << " columns from the DeltaX and DeltaXTilde matrices" << endl;
      }

      assert( DeltaX.cols() == DeltaXTilde.cols() );
      assert( DeltaX.cols() >= 1 );

      // Calculate singular value decomposition with Eigen

      Eigen::JacobiSVD<matrix> svd_X( DeltaX, Eigen::ComputeThinU | Eigen::ComputeThinV );
      Eigen::JacobiSVD<matrix> svd_XTilde( DeltaXTilde, Eigen::ComputeThinU | Eigen::ComputeThinV );

      matrix pseudoSigma_X = svd_X.singularValues().asDiagonal();

      for ( int i = 0; i < pseudoSigma_X.cols(); i++ )
        pseudoSigma_X( i, i ) = 1.0 / pseudoSigma_X( i, i );

      matrix pseudoDeltaX = svd_X.matrixV() * pseudoSigma_X * svd_X.matrixU().transpose();
      matrix I = Eigen::MatrixXd::Identity( m, m );

      // Calculate the new mapping matrix B

      B = DeltaXTilde * pseudoDeltaX + ( I - svd_XTilde.matrixU() * svd_XTilde.matrixU().transpose() ) * ( I - svd_X.matrixU() * svd_X.matrixU().transpose() );

      // Evaluate the coarse model for each xk to fill the mapped coarse response list mappedCoarseList
      mappedCoarseList.clear();
      vector xkTmp, outputTmp( m ), Rtmp( m );

      for ( unsigned i = 0; i < xkList.size(); i++ )
      {
        xkTmp = xkList.at( i );
        coarseModel->evaluate( xktilde + B * (xkTmp - xk), outputTmp, Rtmp );
        mappedCoarseList.push_back( Rtmp );
      }

      assert( mappedCoarseList.size() == xkList.size() );

      // Evaluate the coarse model for each xk for the previous time steps

      mappedCoarsePreviousTimeStepsList.clear();

      for ( unsigned i = 0; i < xkPreviousTimeStepsList.size(); i++ )
      {
        std::deque<fsi::vector> mappedCoarsePreviousTimeStepList;

        for ( unsigned j = 0; j < xkPreviousTimeStepsList.at( i ).size(); j++ )
        {
          xkTmp = xkPreviousTimeStepsList.at( i ).at( j );
          coarseModel->evaluate( xktilde + B * (xkTmp - xk), outputTmp, Rtmp );
          mappedCoarsePreviousTimeStepList.push_back( Rtmp );
        }

        mappedCoarsePreviousTimeStepsList.push_back( mappedCoarsePreviousTimeStepList );

        assert( mappedCoarsePreviousTimeStepList.size() == xkPreviousTimeStepsList.at( i ).size() );
      }

      assert( mappedCoarsePreviousTimeStepsList.size() == xkPreviousTimeStepsList.size() );

      // Update the mapping matrix T

      matrix DeltaF( m, nbCols ), DeltaC( m, nbCols );

      colIndex = 0;

      for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
      {
        DeltaF.col( i ) = fineResiduals.col( k - 1 - i ) - fineResiduals.col( k );
        DeltaC.col( i ) = mappedCoarseList.at( k - 1 - i ) - mappedCoarseList.at( k );
        colIndex++;
      }

      // Include information from previous time steps for the mapping matrix T

      for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
      {
        for ( unsigned j = 0; j < fineResidualsList.at( i ).cols() - 1; j++ )
        {
          if ( colIndex >= DeltaF.cols() )
            continue;

          DeltaF.col( colIndex ) = fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 2 - j ) - fineResidualsList.at( i ).col( fineResidualsList.at( i ).cols() - 1 );
          DeltaC.col( colIndex ) = mappedCoarsePreviousTimeStepsList.at( i ).at( mappedCoarsePreviousTimeStepsList.at( i ).size() - 2 - j ) - mappedCoarsePreviousTimeStepsList.at( i ).at( mappedCoarsePreviousTimeStepsList.at( i ).size() - 1 );
          colIndex++;
        }
      }

      assert( colIndex == nbCols );

      // Remove dependent columns of DeltaC and DeltaF

      nbRemoveCols = 1;

      while ( nbRemoveCols > 0 )
      {
        nbRemoveCols = 0;

        if ( DeltaF.cols() == 1 )
          break;

        // Calculate singular value decomposition with Eigen

        Eigen::JacobiSVD<matrix> svd_F( DeltaF, Eigen::ComputeThinU | Eigen::ComputeThinV );

        vector singularValues = svd_F.singularValues();

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
          Info << "Response and parameter mapping: remove " << nbRemoveCols << " columns from the DeltaC and DeltaF matrices" << endl;
      }

      assert( DeltaC.cols() == DeltaF.cols() );
      assert( DeltaC.cols() >= 1 );

      // Calculate singular value decomposition with Eigen

      Eigen::JacobiSVD<matrix> svd_F( DeltaF, Eigen::ComputeThinU | Eigen::ComputeThinV );
      Eigen::JacobiSVD<matrix> svd_C( DeltaC, Eigen::ComputeThinU | Eigen::ComputeThinV );

      matrix pseudoSigma_F = svd_F.singularValues().asDiagonal();

      for ( int i = 0; i < pseudoSigma_F.cols(); i++ )
        pseudoSigma_F( i, i ) = 1.0 / pseudoSigma_F( i, i );

      matrix pseudoDeltaF = svd_F.matrixV() * pseudoSigma_F * svd_F.matrixU().transpose();

      // Update the design specification yk

      matrix Tk = DeltaC * pseudoDeltaF + ( I - svd_C.matrixU() * svd_C.matrixU().transpose() ) * ( I - svd_F.matrixU() * svd_F.matrixU().transpose() );
      yk = coarseResiduals.col( k ) - Tk * (fineResiduals.col( k ) - y);
    }

    // Update the fine model optimum

    Info << nl << "Response and parameter mapping: update the fine model optimum by minimizing the coarse model" << endl;

    output.resize( n );

    if ( nbCols == 0 )
      coarseModel->optimize( yk, xk, output );

    if ( nbCols >= 1 )
      coarseModel->optimize( yk, xk, output, B, xktilde, xk );

    if ( !coarseModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    xk = output;

    output.resize( m );

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );
    fineResiduals.conservativeResize( fineResiduals.rows(), fineResiduals.cols() + 1 );
    fineResiduals.col( fineResiduals.cols() - 1 ) = R;

    iter++;

    // Parameter extraction
    Info << nl << "Response and parameter mapping: parameter extraction process" << endl;
    coarseModel->optimize( R, xk, xktilde );
    xktildeList.push_back( xktilde );
    xkList.push_back( xk );

    if ( !coarseModel->allConverged() )
      Warning << "Surrogate model optimization process is not converged." << endl;

    // Coarse model evaluation

    coarseModel->evaluate( xktilde, output, R );
    coarseResiduals.conservativeResize( coarseResiduals.rows(), coarseResiduals.cols() + 1 );
    coarseResiduals.col( coarseResiduals.cols() - 1 ) = R;

    assert( coarseResiduals.cols() == fineResiduals.cols() );
    assert( static_cast<int>( xktildeList.size() ) == coarseResiduals.cols() );
    assert( static_cast<int>( xkList.size() ) == coarseResiduals.cols() );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
      break;
  }

  iterationsConverged();
}

void ResponseParameterMapping::removeColumnFromMatrix(
  matrix & A,
  int col
  )
{
  for ( int j = col; j < A.cols() - 1; j++ )
    A.col( j ) = A.col( j + 1 );

  A.conservativeResize( A.rows(), A.cols() - 1 );
}
