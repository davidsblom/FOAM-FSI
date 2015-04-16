
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ImplicitMultiLevelFsiSolver.H"

using namespace fsi;

ImplicitMultiLevelFsiSolver::ImplicitMultiLevelFsiSolver(
    shared_ptr<MultiLevelFsiSolver> fsi,
    shared_ptr<PostProcessing> postProcessing
    )
    :
    fsi( fsi ),
    postProcessing( postProcessing ),
    init( false )
{
    assert( fsi );
    assert( postProcessing );
    assert( fsi->fluid );
    assert( fsi->solid );
    assert( fsi->solid->data.rows() > 0 );
    assert( fsi->solid->data.cols() > 0 );
    assert( fsi->fluid->data.rows() > 0 );
    assert( fsi->fluid->data.cols() > 0 );
}

ImplicitMultiLevelFsiSolver::~ImplicitMultiLevelFsiSolver()
{}

bool ImplicitMultiLevelFsiSolver::allConverged()
{
    return fsi->allConverged;
}

void ImplicitMultiLevelFsiSolver::evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    )
{
    assert( init );

    assert( x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( output.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || output.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( R.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || R.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( x.rows() == output.rows() );
    assert( x.rows() == R.rows() );

    fsi->evaluate( x, output, R );
}

void ImplicitMultiLevelFsiSolver::finalizeTimeStep()
{
    assert( init );

    fsi->finalizeTimeStep();

    postProcessing->finalizeTimeStep();

    init = false;
}

void ImplicitMultiLevelFsiSolver::initTimeStep()
{
    assert( !init );
    assert( !fsi->init );
    assert( !fsi->fluid->init );
    assert( !fsi->solid->init );

    fsi->initTimeStep();

    init = true;

    assert( fsi->init );
    assert( fsi->fluid->init );
    assert( fsi->solid->init );
    assert( fsi->fluid->timeIndex > 0 );
}

bool ImplicitMultiLevelFsiSolver::isConvergence()
{
    return fsi->isConvergence();
}

bool ImplicitMultiLevelFsiSolver::isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    )
{
    assert( init );

    assert( x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( xprev.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || xprev.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( x.rows() == xprev.rows() );

    return fsi->isConvergence( x, xprev );
}

void ImplicitMultiLevelFsiSolver::optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    )
{
    assert( init );

    assert( x0.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || x0.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( xk.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || xk.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( x0.rows() == xk.rows() );

    fsi->newMeasurementSeries();

    postProcessing->performPostProcessing( x0, xk );
}

void ImplicitMultiLevelFsiSolver::optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    )
{
    assert( init );

    assert( y.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || y.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( x0.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || x0.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( xk.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || xk.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    assert( x0.rows() == xk.rows() );
    assert( y.rows() == xk.rows() );

    fsi->newMeasurementSeries();

    postProcessing->performPostProcessing( y, x0, xk );
}

void ImplicitMultiLevelFsiSolver::run()
{
    assert( !init );

    time = std::clock();

    while ( fsi->isRunning() )
        solveTimeStep();
}

void ImplicitMultiLevelFsiSolver::setSurrogateData(
    fsi::vector & xf,
    matrix & J
    )
{
    fsi->setSurrogateData( xf, J );
}

void ImplicitMultiLevelFsiSolver::setUseJacobian( bool useJacobian )
{
    fsi->setUseJacobian( useJacobian );
}

void ImplicitMultiLevelFsiSolver::solve()
{
    assert( init );

    assert( fsi->x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        || fsi->x.rows() == fsi->solidSolver->couplingGridSize * fsi->solid->dim
        + fsi->fluidSolver->couplingGridSize * fsi->fluid->dim );

    fsi->newMeasurementSeries();

    // Initial solution
    fsi::vector x0 = fsi->x;

    postProcessing->performPostProcessing( x0, fsi->x );
}

void ImplicitMultiLevelFsiSolver::solveTimeStep()
{
    assert( !init );

    initTimeStep();
    solve();
    finalizeTimeStep();
}
