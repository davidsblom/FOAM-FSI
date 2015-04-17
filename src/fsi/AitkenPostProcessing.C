
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <math.h>

#include "AitkenPostProcessing.H"

using namespace fsi;

AitkenPostProcessing::AitkenPostProcessing(
    shared_ptr<MultiLevelFsiSolver> fsi,
    double initialRelaxation,
    int maxIter,
    int maxUsedIterations,
    int nbReuse,
    int reuseInformationStartingFromTimeIndex
    )
    :
    PostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ),
    aitkenFactor( initialRelaxation )
{
    assert( aitkenFactor > 0 );
    assert( aitkenFactor < 1 );
}

AitkenPostProcessing::~AitkenPostProcessing()
{}

void AitkenPostProcessing::performPostProcessing(
    const vector & x0,
    vector & xk
    )
{
    vector y( x0.rows() );
    y.setZero();
    performPostProcessing( y, x0, xk, true );
}

void AitkenPostProcessing::performPostProcessing(
    const vector & y,
    const vector & x0,
    vector & xk
    )
{
    performPostProcessing( y, x0, xk, false );
}

void AitkenPostProcessing::performPostProcessing(
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
    residuals.clear();
    sols.clear();
    vector yk = y;
    aitkenFactor = initialRelaxation;

    // Fsi evaluation
    vector output( xk.rows() ), R( xk.rows() ), Rprev;
    output.setZero();
    R.setZero();

    fsi->evaluate( x0, output, R );

    assert( x0.rows() == output.rows() );
    assert( x0.rows() == R.rows() );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
    {
        bool keepIterations = residualCriterium;
        iterationsConverged( keepIterations );
        return;
    }

    for ( int iter = 0; iter < maxIter - 1; iter++ )
    {
        if ( iter == 0 )
            xk += initialRelaxation * (R - yk);

        if ( iter > 0 )
        {
            vector deltaR = R - Rprev;

            aitkenFactor *= -deltaR.dot( Rprev ) / deltaR.squaredNorm();

            // Safe guards for the aitken relaxation factor
            aitkenFactor = std::max( aitkenFactor, initialRelaxation );
            aitkenFactor = std::min( aitkenFactor, 2.0 );

            assert( !std::isnan( aitkenFactor ) );

            xk += aitkenFactor * (R - yk);
        }

        Info << "Aitken post processing with factor " << aitkenFactor << endl;

        Rprev = R;

        // Fsi evaluation
        fsi->evaluate( xk, output, R );

        assert( x0.rows() == output.rows() );
        assert( x0.rows() == R.rows() );

        // Check convergence criteria
        if ( isConvergence( xk, xkprev, residualCriterium ) )
        {
            bool keepIterations = residualCriterium;
            iterationsConverged( keepIterations );
            break;
        }

        assert( fsi->iter <= maxIter );
    }
}
