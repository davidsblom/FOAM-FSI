
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCFsiSolver.H"

SDCFsiSolver::SDCFsiSolver(
    std::shared_ptr<sdc::SDCFsiSolverInterface> fluid,
    std::shared_ptr<sdc::SDCFsiSolverInterface> solid,
    std::shared_ptr<PostProcessing> postProcessing
    )
    :
    fluid( fluid ),
    solid( solid ),
    postProcessing( postProcessing ),
    dofFluid( fluid->getDOF() ),
    dofSolid( solid->getDOF() ),
    k( 0 ),
    xStages()
{
    assert( fluid );
    assert( solid );
    assert( postProcessing );
    assert( dofFluid > 0 );
    assert( dofSolid > 0 );
}

SDCFsiSolver::~SDCFsiSolver()
{}

void SDCFsiSolver::evaluateFunction(
    const int k,
    const fsi::vector & q,
    const scalar t,
    fsi::vector & f
    )
{
    fsi::vector fFluid( dofFluid ), fSolid( dofSolid );

    Eigen::Map<const fsi::vector> qFluid( q.data(), dofFluid );
    Eigen::Map<const fsi::vector> qSolid( q.data() + dofFluid, dofSolid );

    fluid->evaluateFunction( k, qFluid, t, fFluid );
    solid->evaluateFunction( k, qSolid, t, fSolid );

    f.head( dofFluid ) = fFluid;
    f.tail( dofSolid ) = fSolid;
}

void SDCFsiSolver::finalizeTimeStep()
{
    postProcessing->fsi->finalizeTimeStep();
}

int SDCFsiSolver::getDOF()
{
    return dofFluid + dofSolid;
}

scalar SDCFsiSolver::getScalingFactor()
{
    return 1;
}

void SDCFsiSolver::getSolution(
    fsi::vector & solution,
    fsi::vector & f
    )
{
    assert( dofFluid > 0 );
    assert( dofSolid > 0 );
    assert( dofFluid + dofSolid == solution.rows() );
    assert( solution.rows() == f.rows() );

    fsi::vector solFluid( dofFluid ), solSolid( dofSolid ), fFluid( dofFluid ), fSolid( dofSolid );
    fluid->getSolution( solFluid, fFluid );
    solid->getSolution( solSolid, fSolid );

    solution.head( dofFluid ) = solFluid;
    solution.tail( dofSolid ) = solSolid;
    f.head( dofFluid ) = fFluid;
    f.tail( dofSolid ) = fSolid;
}

void SDCFsiSolver::setSolution(
    const fsi::vector & solution,
    const fsi::vector & f
    )
{
    Eigen::Map<const fsi::vector> solFluid( solution.data(), dofFluid );
    Eigen::Map<const fsi::vector> solSolid( solution.data() + dofFluid, dofSolid );

    Eigen::Map<const fsi::vector> fFluid( f.data(), dofFluid );
    Eigen::Map<const fsi::vector> fSolid( f.data() + dofFluid, dofSolid );

    fluid->setSolution( solFluid, fFluid );
    solid->setSolution( solSolid, fSolid );
}

scalar SDCFsiSolver::getEndTime()
{
    return fluid->getEndTime();
}

scalar SDCFsiSolver::getTimeStep()
{
    return fluid->getTimeStep();
}

void SDCFsiSolver::nextTimeStep()
{
    fluid->nextTimeStep();
    solid->nextTimeStep();

    for ( int i = 0; i < k; i++ )
        xStages.at( i ) = postProcessing->fsi->x;
}

void SDCFsiSolver::initTimeStep()
{
    postProcessing->fsi->initTimeStep();
}

void SDCFsiSolver::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    fluid->setNumberOfImplicitStages( k );
    solid->setNumberOfImplicitStages( k );
    postProcessing->setNumberOfImplicitStages( k );

    xStages.clear();

    for ( int i = 0; i < k + 1; i++ )
        xStages.push_back( postProcessing->fsi->x );
}

void SDCFsiSolver::implicitSolve(
    bool corrector,
    const int k,
    const int kold,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    )
{
    Eigen::Map<const fsi::vector> qoldFluid( qold.data(), dofFluid );
    Eigen::Map<const fsi::vector> qoldSolid( qold.data() + dofFluid, dofSolid );

    Eigen::Map<const fsi::vector> rhsFluid( rhs.data(), dofFluid );
    Eigen::Map<const fsi::vector> rhsSolid( rhs.data() + dofFluid, dofSolid );

    fluid->prepareImplicitSolve( corrector, k, kold, t, dt, qoldFluid, rhsFluid );
    solid->prepareImplicitSolve( corrector, k, kold, t, dt, qoldSolid, rhsSolid );

    // Perform FSI iterations to solve the coupled problem

    postProcessing->fsi->newMeasurementSeries();

    // Initial solution
    fsi::vector x0;

    if ( corrector )
        x0 = xStages.at( k + 1 );
    else
        x0 = xStages.at( k );

    postProcessing->initStage( k );
    postProcessing->performPostProcessing( x0, postProcessing->fsi->x );
    postProcessing->finalizeStage();

    getSolution( result, f );

    xStages.at( k + 1 ) = postProcessing->fsi->x;
}

scalar SDCFsiSolver::getStartTime()
{
    return fluid->getStartTime();
}

void SDCFsiSolver::getVariablesInfo(
    std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    )
{
    std::deque<int> dofFluid, dofSolid;
    std::deque<bool> enabledFluid, enabledSolid;
    std::deque<std::string> namesFluid, namesSolid;
    fluid->getVariablesInfo( dofFluid, enabledFluid, namesFluid );
    solid->getVariablesInfo( dofSolid, enabledSolid, namesSolid );

    assert( dofFluid.size() == enabledFluid.size() );
    assert( dofFluid.size() == namesFluid.size() );
    assert( dofSolid.size() == enabledSolid.size() );
    assert( dofSolid.size() == namesSolid.size() );

    for ( unsigned int i = 0; i < dofFluid.size(); i++ )
    {
        dof.push_back( dofFluid.at( i ) );
        enabled.push_back( enabledFluid.at( i ) );
        names.push_back( namesFluid.at( i ) );
    }

    for ( unsigned int i = 0; i < dofSolid.size(); i++ )
    {
        dof.push_back( dofSolid.at( i ) );
        enabled.push_back( enabledSolid.at( i ) );
        names.push_back( namesSolid.at( i ) );
    }
}

bool SDCFsiSolver::isConverged()
{
    return fluid->isConverged() && solid->isConverged();
}
