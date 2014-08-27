
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "IQNILSPostProcessing.H"

namespace fsi
{
  IQNILSPostProcessing::IQNILSPostProcessing(
    shared_ptr<MultiLevelFsiSolver> fsi,
    int maxIter,
    double initialRelaxation,
    int maxUsedIterations,
    int nbReuse,
    double singularityLimit,
    int reuseInformationStartingFromTimeIndex,
    bool scaling
    )
    :
    PostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ),
    scaling( scaling ),
    singularityLimit( singularityLimit ),
    scalingFactors( Eigen::VectorXd::Ones( 2 ) )
  {
    assert( fsi );
    assert( singularityLimit > 0 );
    assert( singularityLimit < 1 );

    if ( scaling )
      assert( fsi->parallel );
  }

  void IQNILSPostProcessing::fixedUnderRelaxation(
    vector & xk,
    vector & R,
    vector & yk
    )
  {
    Info << "Fixed relaxation post processing with factor " << initialRelaxation << endl;

    // Use scaling if the fluid and solid are coupled in parallel
    int sizeVar0 = fsi->solidSolver->couplingGridSize * fsi->solid->dim;
    int sizeVar1 = fsi->fluidSolver->couplingGridSize * fsi->fluid->dim;

    // Update solution x
    fsi::vector dx = initialRelaxation * (R - yk);

    if ( scaling )
    {
      dx.head( sizeVar0 ).array() *= scalingFactors( 0 );
      dx.tail( sizeVar1 ).array() *= scalingFactors( 1 );
    }

    xk += dx;
  }

  void IQNILSPostProcessing::performPostProcessing(
    const vector & x0,
    vector & xk
    )
  {
    vector y( x0.rows() );
    y.setZero();
    performPostProcessing( y, x0, xk, true );
  }

  void IQNILSPostProcessing::performPostProcessing(
    const vector & y,
    const vector & x0,
    vector & xk
    )
  {
    performPostProcessing( y, x0, xk, false );
  }

  /*
   * Minimize the FSI residual.
   * Two different convergence criteria can be used: based on the sequence of
   * iterations, or on the FSI residual
   * bool residualCriterium = false: criterium based on the sequence of iterations
   * residualCriterium = true: based on the input/output information of the
   * FSI residual.
   */
  void IQNILSPostProcessing::performPostProcessing(
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
    xkprev.setZero();
    xk = x0;
    residuals.clear();
    sols.clear();
    vector yk = y;

    // Fsi evaluation
    vector output( xk.rows() ), R( xk.rows() );
    output.setZero();
    R.setZero();

    fsi->evaluate( x0, output, R );

    // Use scaling if the fluid and solid are coupled in parallel
    int sizeVar0 = fsi->solidSolver->couplingGridSize * fsi->solid->dim;
    int sizeVar1 = fsi->fluidSolver->couplingGridSize * fsi->fluid->dim;

    if ( scaling && timeIndex < 10 )
    {
      assert( output.rows() == sizeVar0 + sizeVar1 );

      scalingFactors( 0 ) = output.head( sizeVar0 ).norm();
      scalingFactors( 1 ) = output.tail( sizeVar1 ).norm();

      for ( int i = 0; i < scalingFactors.rows(); i++ )
        if ( std::abs( scalingFactors( i ) ) < 1.0e-13 )
          scalingFactors( i ) = 1;

      Info << "Parallel coupling of fluid and solid solvers with scaling factors ";
      Info << scalingFactors( 0 ) << " and " << scalingFactors( 1 ) << endl;
    }

    assert( x0.rows() == output.rows() );
    assert( x0.rows() == R.rows() );

    // Save output and residual
    residuals.push_front( R );
    sols.push_front( output );

    if ( scaling )
    {
      R.head( sizeVar0 ).array() /= scalingFactors( 0 );
      R.tail( sizeVar1 ).array() /= scalingFactors( 1 );
      yk.head( sizeVar0 ).array() /= scalingFactors( 0 );
      yk.tail( sizeVar1 ).array() /= scalingFactors( 1 );
    }

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

      nbCols = std::min( static_cast<int>( xk.rows() ), nbCols );
      nbCols = std::min( nbCols, maxUsedIterations );

      assert( nbCols <= xk.rows() );

      if ( nbCols == 0 )
      {
        // Use fixed under relaxation to startup the IQNILS algorithm

        fixedUnderRelaxation( xk, R, yk );
      }

      if ( nbCols > 0 )
      {
        // IQNILS algorithm

        Info << "IQN-ILS post processing with " << nbCols << " cols for the Jacobian" << endl;

        // Construct the V and W matrices
        matrix V( xk.rows(), nbCols ), W( xk.rows(), nbCols );

        int nbColsCurrentTimeStep = std::max( static_cast<int>(sols.size() - 1), 0 );
        nbColsCurrentTimeStep = std::min( nbColsCurrentTimeStep, nbCols );

        // Include information from previous iterations

        int colIndex = 0;

        for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
        {
          if ( colIndex >= V.cols() )
            continue;

          V.col( i ) = residuals.at( i ) - residuals.at( i + 1 );
          W.col( i ) = sols.at( i ) - sols.at( i + 1 );
          colIndex++;
        }

        // Include information from previous optimization solves
        for ( unsigned i = 0; i < residualsList.size(); i++ )
        {
          assert( residualsList.at( i ).size() >= 2 );

          for ( unsigned j = 0; j < residualsList.at( i ).size() - 1; j++ )
          {
            if ( colIndex >= V.cols() )
              continue;

            V.col( colIndex ) = residualsList.at( i ).at( j ) - residualsList.at( i ).at( j + 1 );
            W.col( colIndex ) = solsList.at( i ).at( j ) - solsList.at( i ).at( j + 1 );
            colIndex++;
          }
        }

        // Include information from previous time steps

        for ( unsigned i = 0; i < residualsTimeList.size(); i++ )
        {
          for ( unsigned j = 0; j < residualsTimeList.at( i ).size(); j++ )
          {
            assert( residualsTimeList.at( i ).at( j ).size() >= 2 );

            for ( unsigned k = 0; k < residualsTimeList.at( i ).at( j ).size() - 1; k++ )
            {
              if ( colIndex >= V.cols() )
                continue;

              V.col( colIndex ) = residualsTimeList.at( i ).at( j ).at( k ) - residualsTimeList.at( i ).at( j ).at( k + 1 );
              W.col( colIndex ) = solsTimeList.at( i ).at( j ).at( k ) - solsTimeList.at( i ).at( j ).at( k + 1 );
              colIndex++;
            }
          }
        }

        assert( colIndex == nbCols );

        // Apply scaling
        if ( scaling )
        {
          V.topLeftCorner( sizeVar0, V.cols() ).array() /= scalingFactors( 0 );
          W.topLeftCorner( sizeVar0, W.cols() ).array() /= scalingFactors( 0 );
          V.bottomLeftCorner( sizeVar1, V.cols() ).array() /= scalingFactors( 1 );
          W.bottomLeftCorner( sizeVar1, W.cols() ).array() /= scalingFactors( 1 );
        }

        // Remove dependent columns of V and W

        int nbRemoveCols = 1;

        while ( nbRemoveCols > 0 )
        {
          nbRemoveCols = 0;

          Eigen::HouseholderQR<matrix> qr = V.householderQr();
          matrix R = qr.matrixQR().triangularView <Eigen::Upper> ();
          Eigen::VectorXd diagonals = R.diagonal();

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
            Info << "IQN-ILS: remove " << nbRemoveCols << " columns from the Jacobian matrices" << endl;
        }

        assert( V.cols() == W.cols() );

        if ( V.cols() == 0 )
          fixedUnderRelaxation( xk, R, yk );

        else
        {
          // SVD decomposition

          Eigen::JacobiSVD<matrix, Eigen::FullPivHouseholderQRPreconditioner> svd( V, Eigen::ComputeFullU | Eigen::ComputeFullV );
          vector c = svd.solve( yk - R );

          // Update solution x
          vector dx = W * c + R - yk;

          if ( scaling )
          {
            dx.head( sizeVar0 ).array() *= scalingFactors( 0 );
            dx.tail( sizeVar1 ).array() *= scalingFactors( 1 );
          }

          xk += dx;
        }
      }

      // Fsi evaluation
      fsi->evaluate( xk, output, R );

      assert( x0.rows() == output.rows() );
      assert( x0.rows() == R.rows() );

      // Save output and residual
      residuals.push_front( R );
      sols.push_front( output );

      if ( scaling )
      {
        R.head( sizeVar0 ).array() /= scalingFactors( 0 );
        R.tail( sizeVar1 ).array() /= scalingFactors( 1 );
      }

      // Check convergence criteria
      if ( isConvergence( xk, xkprev, residualCriterium ) )
        break;

      assert( sols.size() == residuals.size() );
      assert( sols.at( 0 ).rows() == residuals.at( 0 ).rows() );
      assert( fsi->iter <= maxIter );
    }

    iterationsConverged();
  }

  void IQNILSPostProcessing::removeColumnFromMatrix(
    matrix & A,
    int col
    )
  {
    for ( int j = col; j < A.cols() - 1; j++ )
      A.col( j ) = A.col( j + 1 );

    A.conservativeResize( A.rows(), A.cols() - 1 );
  }
}
