
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AndersonPostProcessing.H"

namespace fsi
{
    AndersonPostProcessing::AndersonPostProcessing(
        shared_ptr<MultiLevelFsiSolver> fsi,
        int maxIter,
        scalar initialRelaxation,
        int maxUsedIterations,
        int nbReuse,
        scalar singularityLimit,
        int reuseInformationStartingFromTimeIndex,
        bool scaling,
        scalar beta,
        bool updateJacobian
        )
        :
        PostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ),
        scaling( scaling ),
        beta( beta ),
        singularityLimit( singularityLimit ),
        updateJacobian( updateJacobian ),
        scalingFactors( fsi::vector::Ones( 2 ) ),
        sizeVar0( 0 ),
        sizeVar1( 0 )
    {
        assert( fsi );
        assert( singularityLimit > 0 );
        assert( singularityLimit < 1 );
        assert( beta > 0 );

        if ( scaling )
            assert( fsi->parallel );
    }

    void AndersonPostProcessing::applyScaling( vector & vec )
    {
        vec.head( sizeVar0 ) /= scalingFactors( 0 );
        vec.tail( sizeVar1 ) /= scalingFactors( 1 );
    }

    void AndersonPostProcessing::applyScaling( matrix & mat )
    {
        mat.topLeftCorner( sizeVar0, mat.cols() ).array() /= scalingFactors( 0 );
        mat.bottomLeftCorner( sizeVar1, mat.cols() ).array() /= scalingFactors( 1 );
    }

    void AndersonPostProcessing::determineScalingFactors( const vector & output )
    {
        if ( scaling && timeIndex <= reuseInformationStartingFromTimeIndex )
        {
            // Use scaling if the fluid and solid are coupled in parallel
            sizeVar0 = fsi->solidSolver->couplingGridSize * fsi->solid->dim;
            sizeVar1 = fsi->fluidSolver->couplingGridSize * fsi->fluid->dim;

            assert( output.rows() == sizeVar0 + sizeVar1 );

            scalingFactors( 0 ) = output.head( sizeVar0 ).norm();
            scalingFactors( 1 ) = output.tail( sizeVar1 ).norm();

            for ( int i = 0; i < scalingFactors.rows(); i++ )
                if ( std::abs( scalingFactors( i ) ) < 1.0e-13 )
                    scalingFactors( i ) = 1;

            Info << "Parallel coupling of fluid and solid solvers with scaling factors ";
            Info << scalingFactors( 0 ) << " and " << scalingFactors( 1 ) << endl;

            // Reset jacobian matrix J since a difference scaling factor is used
            Jprev.resize( 0, 0 );
        }
    }

    void AndersonPostProcessing::fixedUnderRelaxation(
        vector & xk,
        vector & R,
        vector & yk
        )
    {
        // Update solution x
        fsi::vector dx;

        if ( updateJacobian && Jprev.rows() == yk.rows() )
        {
            Info << "Anderson mixing method: reuse Jacobian of previous time step or optimization" << endl;
            dx = Jprev * (yk - R);
        }
        else
        {
            Info << "Fixed relaxation post processing with factor " << initialRelaxation << endl;
            dx = initialRelaxation * (R - yk);
        }

        if ( scaling )
        {
            dx.head( sizeVar0 ).array() *= scalingFactors( 0 );
            dx.tail( sizeVar1 ).array() *= scalingFactors( 1 );
        }

        xk += dx;
    }

    void AndersonPostProcessing::performPostProcessing(
        const vector & x0,
        vector & xk
        )
    {
        vector y( x0.rows() );
        y.setZero();
        performPostProcessing( y, x0, xk, true );
    }

    void AndersonPostProcessing::performPostProcessing(
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
    void AndersonPostProcessing::performPostProcessing(
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
        matrix J;

        // Fsi evaluation
        vector output( xk.rows() ), R( xk.rows() );
        output.setZero();
        R.setZero();

        fsi->evaluate( x0, output, R );

        determineScalingFactors( output );

        assert( x0.rows() == output.rows() );
        assert( x0.rows() == R.rows() );

        // Save output and residual
        residuals.push_front( R );
        sols.push_front( x0 );

        // Check convergence criteria
        if ( isConvergence( output, output + y - R, residualCriterium ) )
        {
            bool keepIterations = residualCriterium || solsList.size() == 0;
            iterationsConverged( keepIterations );
            return;
        }

        if ( scaling )
        {
            applyScaling( R );
            applyScaling( yk );
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

            // Use fixed under relaxation to startup the IQNILS algorithm
            if ( nbCols == 0 )
                fixedUnderRelaxation( xk, R, yk );

            if ( nbCols > 0 )
            {
                // Anderson mixing method

                Info << "Anderson mixing method: post processing with " << nbCols << " cols for the Jacobian" << endl;

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
                    applyScaling( V );
                    applyScaling( W );
                }

                // Remove dependent columns of V and W

                int nbRemoveCols = 1;

                while ( nbRemoveCols > 0 )
                {
                    nbRemoveCols = 0;

                    Eigen::HouseholderQR<matrix> qr = V.householderQr();
                    matrix R = qr.matrixQR().triangularView <Eigen::Upper> ();
                    fsi::vector diagonals = R.diagonal();

                    for ( int i = 0; i < diagonals.rows(); i++ )
                    {
                        if ( std::abs( diagonals[i] ) < singularityLimit && V.cols() > 5 )
                        {
                            // Remove the column from V and W
                            removeColumnFromMatrix( V, i - nbRemoveCols );
                            removeColumnFromMatrix( W, i - nbRemoveCols );

                            nbRemoveCols++;
                        }
                    }

                    if ( nbRemoveCols )
                        Info << "Anderson mixing method: remove " << nbRemoveCols << " columns from the Jacobian matrices" << endl;
                }

                assert( V.cols() == W.cols() );

                if ( V.cols() == 0 )
                    fixedUnderRelaxation( xk, R, yk );

                else
                {
                    vector dx;

                    if ( updateJacobian )
                    {
                        Eigen::JacobiSVD<matrix> svd( V, Eigen::ComputeThinU | Eigen::ComputeThinV );

                        matrix Vinverse = svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose();

                        matrix I = fsi::matrix::Identity( V.rows(), V.rows() );

                        if ( Jprev.cols() == R.rows() )
                        {
                            Info << "Anderson mixing method: reuse Jacobian of previous time step or optimization" << endl;
                            J = Jprev + (W - Jprev * V) * Vinverse;
                        }
                        else
                            J = (V + W) * Vinverse - I;

                        dx = J * (yk - R);
                    }

                    if ( !updateJacobian )
                    {
                        vector c = V.fullPivHouseholderQr().solve( yk - R );
                        dx = beta * (R - yk) + W * c + beta * V * c;
                    }

                    // Update solution x

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
            sols.push_front( xk );

            // Check convergence criteria
            if ( isConvergence( output, output + y - R, residualCriterium ) )
            {
                bool keepIterations = residualCriterium || solsList.size() == 0;
                iterationsConverged( keepIterations );

                if ( updateJacobian && J.cols() > 0 && timeIndex >= reuseInformationStartingFromTimeIndex )
                    Jprev = J;

                break;
            }

            determineScalingFactors( output );

            if ( scaling )
            {
                yk = y;
                applyScaling( R );
                applyScaling( yk );
            }

            assert( sols.size() == residuals.size() );
            assert( sols.at( 0 ).rows() == residuals.at( 0 ).rows() );
            assert( fsi->iter <= maxIter );
        }
    }
}

void AndersonPostProcessing::removeColumnFromMatrix(
    matrix & A,
    int col
    )
{
    for ( int j = col; j < A.cols() - 1; j++ )
        A.col( j ) = A.col( j + 1 );

    A.conservativeResize( A.rows(), A.cols() - 1 );
}
