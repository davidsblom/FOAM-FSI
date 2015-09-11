
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SpaceMappingSolver.H"

using namespace fsi;

SpaceMappingSolver::SpaceMappingSolver(
    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel,
    shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel,
    shared_ptr<SpaceMapping> spaceMapping
    )
    :
    ImplicitMultiLevelFsiSolver( fineModel->fsi, fineModel->postProcessing ),
    fineModel( fineModel ),
    coarseModel( coarseModel ),
    spaceMapping( spaceMapping )
{
    assert( fineModel );
    assert( coarseModel );
    assert( spaceMapping );
    assert( fineModel->fsi->fluidSolver->couplingGridSize > 0 );
    assert( coarseModel->fsi->fluidSolver->couplingGridSize > 0 );
    assert( fineModel->fsi->fluidSolver->level == fineModel->fsi->solidSolver->level );
    assert( coarseModel->fsi->fluidSolver->level == coarseModel->fsi->solidSolver->level );
    assert( fineModel->fsi->fluidSolver->level > coarseModel->fsi->fluidSolver->level );
    assert( fineModel->fsi->N <= fineModel->fsi->fluidSolver->couplingGridSize );
}

SpaceMappingSolver::~SpaceMappingSolver()
{}

bool SpaceMappingSolver::allConverged()
{
    return fineModel->fsi->allConverged;
}

void SpaceMappingSolver::evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    )
{
    assert( init );

    assert( x.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || x.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( output.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || output.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( R.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || R.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( x.rows() == output.rows() );
    assert( x.rows() == R.rows() );

    fineModel->evaluate( x, output, R );
}

void SpaceMappingSolver::finalizeTimeStep()
{
    assert( init );

    fineModel->finalizeTimeStep();
    coarseModel->finalizeTimeStep();

    spaceMapping->finalizeTimeStep();

    init = false;
}

void SpaceMappingSolver::initTimeStep()
{
    assert( !init );
    assert( !fineModel->init );
    assert( !coarseModel->init );

    fineModel->initTimeStep();
    coarseModel->initTimeStep();

    init = true;

    assert( fineModel->init );
    assert( fineModel->fsi->init );
    assert( fineModel->fsi->fluid->init );
    assert( fineModel->fsi->fluid->timeIndex > 0 );
    assert( fineModel->fsi->fluid->t > 0 );
    assert( coarseModel->init );
}

bool SpaceMappingSolver::isConvergence()
{
    return fineModel->isConvergence();
}

bool SpaceMappingSolver::isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    )
{
    return fineModel->isConvergence( x, xprev );
}

bool SpaceMappingSolver::isRunning()
{
    assert( !init );

    return fineModel->fsi->isRunning();
}

void SpaceMappingSolver::optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    )
{
    fsi::vector y( xk.rows() );
    y.setZero();

    optimize( y, x0, xk );
}

void SpaceMappingSolver::optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    )
{
    assert( init );

    assert( y.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || y.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( x0.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || x0.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( xk.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || xk.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    assert( y.rows() == x0.rows() );
    assert( y.rows() == xk.rows() );

    fineModel->fsi->newMeasurementSeries();

    spaceMapping->performPostProcessing( y, x0, xk );
}

void SpaceMappingSolver::run()
{
    assert( !init );

    time = std::clock();

    while ( isRunning() )
        solveTimeStep();
}

void SpaceMappingSolver::setSurrogateData(
    fsi::vector & xf,
    matrix & J
    )
{
    fineModel->setSurrogateData( xf, J );
}

void SpaceMappingSolver::setUseJacobian( bool useJacobian )
{
    fineModel->setUseJacobian( useJacobian );
}

void SpaceMappingSolver::solve()
{
    assert( init );
    assert( fineModel->fsi->x.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        || fineModel->fsi->x.rows() == fineModel->fsi->solidSolver->couplingGridSize * fineModel->fsi->solid->dim
        + fineModel->fsi->fluidSolver->couplingGridSize * fineModel->fsi->fluid->dim );

    // Initialize solution and design specification
    fsi::vector x0( fineModel->fsi->x.rows() ), xk( fineModel->fsi->x.rows() );

    x0 = fineModel->fsi->x;

    fineModel->fsi->newMeasurementSeries();

    spaceMapping->performPostProcessing( x0, xk );

    fineModel->fsi->x = xk;
}

void SpaceMappingSolver::solveTimeStep()
{
    assert( !init );

    initTimeStep();
    solve();
    finalizeTimeStep();
}
