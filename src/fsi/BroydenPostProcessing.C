
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "BroydenPostProcessing.H"

using namespace fsi;

BroydenPostProcessing::BroydenPostProcessing(
    shared_ptr<MultiLevelFsiSolver> fsi,
    int maxIter,
    scalar initialRelaxation,
    int maxUsedIterations,
    int nbReuse,
    scalar singularityLimit,
    int reuseInformationStartingFromTimeIndex
    )
    :
    PostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ),
    J(),
    singularityLimit( singularityLimit )
{
    assert( singularityLimit > 0 );
    assert( singularityLimit < 1 );
}

BroydenPostProcessing::~BroydenPostProcessing()
{}

void BroydenPostProcessing::finalizeTimeStep()
{
    assert( sols.size() == 0 );

    PostProcessing::finalizeTimeStep();

    assert( residualsList.size() == 0 );

    Info << "Broyden post processing: rebuild Jacobian with information from " << nbReuse << " previous time steps" << endl;

    // Rebuild the Jacobian at the end of each time step

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

    J = -fsi::matrix::Identity( J.rows(), J.cols() );

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

                colIndex++;

                fsi::vector dx = solsTimeList.at( i ).at( j ).at( k + 1 ) - solsTimeList.at( i ).at( j ).at( k );
                fsi::vector dR = residualsTimeList.at( i ).at( j ).at( k + 1 ) - residualsTimeList.at( i ).at( j ).at( k );

                if ( dx.norm() < singularityLimit )
                    continue;

                // Sherman–Morrison formula
                J += (dx - J * dR) / (dx.transpose() * J * dR) * (dx.transpose() * J);
            }
        }
    }

    assert( colIndex == nbCols );
}

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
    if ( isConvergence( output, output + y - R, residualCriterium ) )
    {
        bool keepIterations = residualCriterium;
        iterationsConverged( keepIterations );
        return;
    }

    if ( J.cols() != R.rows() || J.rows() != R.rows() )
        J = -fsi::matrix::Identity( R.rows(), R.rows() );

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

        int nbColsCurrentTimeStep = std::max( static_cast<int>(sols.size() - 1), 0 );

        assert( nbCols >= 0 );

        if ( nbCols == 0 )
        {
            // Use fixed under relaxation to startup the IQNILS algorithm

            Info << "Fixed relaxation post processing with factor " << initialRelaxation << endl;

            // Update solution x
            xk += initialRelaxation * (R - y);
        }

        if ( nbCols > 0 )
        {
            Info << "Broyden post processing with " << nbCols << " cols for the Jacobian" << endl;

            if ( nbColsCurrentTimeStep )
            {
                fsi::vector dx, dR;
                dx = sols.at( sols.size() - 1 ) - sols.at( sols.size() - 2 );
                dR = residuals.at( residuals.size() - 1 ) - residuals.at( residuals.size() - 2 );

                if ( dx.norm() >= singularityLimit )
                    J += (dx - J * dR) / (dx.transpose() * J * dR) * (dx.transpose() * J);
            }

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
        if ( isConvergence( output, output + y - R, residualCriterium ) )
        {
            bool keepIterations = residualCriterium;
            iterationsConverged( keepIterations );
            break;
        }

        assert( sols.size() == residuals.size() );
        assert( sols.at( 0 ).rows() == residuals.at( 0 ).rows() );
        assert( fsi->iter <= maxIter );
    }

    // Do not save the iteration vectors since the scheme did not converge
    // A large number of those vectors do not contain any information which
    // can speedup subsequent optimizations
    bool keepIterations = false;
    iterationsConverged( keepIterations );
}
