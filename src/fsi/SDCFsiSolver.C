
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCFsiSolver.H"

SDCFsiSolver::SDCFsiSolver(
    std::shared_ptr<SDCSolver> fluid,
    std::shared_ptr<SDCSolver> solid,
    std::shared_ptr<PostProcessing> postProcessing
    )
    :
    fluid( fluid ),
    solid( solid ),
    postProcessing( postProcessing ),
    dofFluid( fluid->getDOF() ),
    dofSolid( solid->getDOF() )
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
    fsi::vector qFluid( dofFluid ), qSolid( dofSolid ), fFluid( dofFluid ), fSolid( dofSolid );
    qFluid = q.head( dofFluid );
    qSolid = q.tail( dofSolid );
    fluid->evaluateFunction( k, qFluid, t, fFluid );
    solid->evaluateFunction( k, qSolid, t, fSolid );
    f.head( dofFluid ) = fFluid;
    f.tail( dofSolid ) = fSolid;
}

void SDCFsiSolver::finalizeTimeStep()
{
    fluid->finalizeTimeStep();
    solid->finalizeTimeStep();
}

int SDCFsiSolver::getDOF()
{
    return dofFluid + dofSolid;
}

scalar SDCFsiSolver::getScalingFactor()
{
    return 1;
}

void SDCFsiSolver::getSolution( fsi::vector & solution )
{
    fsi::vector solFluid( dofFluid ), solSolid( dofSolid );
    fluid->getSolution( solFluid );
    solid->getSolution( solSolid );
    solution.head( dofFluid ) = solFluid;
    solution.tail( dofSolid ) = solSolid;
}

void SDCFsiSolver::setSolution(
    const fsi::vector & solution,
    const fsi::vector & f
    )
{
    fsi::vector solFluid( dofFluid ), solSolid( dofSolid ), fFluid( dofFluid ), fSolid( fSolid );
    solFluid = solution.head( dofFluid );
    solSolid = solution.tail( dofSolid );
    fFluid = f.head( dofFluid );
    fSolid = f.tail( dofSolid );
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
}

void SDCFsiSolver::initTimeStep()
{
    fluid->initTimeStep();
    solid->initTimeStep();
}

void SDCFsiSolver::setNumberOfImplicitStages( int k )
{
    fluid->setNumberOfImplicitStages( k );
    solid->setNumberOfImplicitStages( k );
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
{}

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
