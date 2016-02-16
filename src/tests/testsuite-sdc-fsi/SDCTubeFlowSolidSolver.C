
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCTubeFlowSolidSolver.H"

namespace tubeflow
{
    SDCTubeFlowSolidSolver::SDCTubeFlowSolidSolver(
        scalar a0,
        scalar cmk,
        scalar p0,
        scalar rho,
        scalar L,
        int N
        )
        :
        TubeFlowSolidSolver( a0, cmk, p0, rho, L, N )
    {}

    SDCTubeFlowSolidSolver::~SDCTubeFlowSolidSolver()
    {}

    void SDCTubeFlowSolidSolver::evaluateFunction(
        const int /*k*/,
        const fsi::vector & /*q*/,
        const scalar /*t*/,
        fsi::vector & /*f*/
        )
    {}

    void SDCTubeFlowSolidSolver::finalizeTimeStep()
    {
        TubeFlowSolidSolver::finalizeTimeStep();
    }

    int SDCTubeFlowSolidSolver::getDOF()
    {
        return 0;
    }

    void SDCTubeFlowSolidSolver::getSolution(
        fsi::vector & /*solution*/,
        fsi::vector & /*f*/
        )
    {}

    void SDCTubeFlowSolidSolver::setSolution(
        const fsi::vector & /*solution*/,
        const fsi::vector & /*f*/
        )
    {
        assert( false );
    }

    scalar SDCTubeFlowSolidSolver::getEndTime()
    {
        assert( false );
        return 0;
    }

    scalar SDCTubeFlowSolidSolver::getTimeStep()
    {
        assert( false );
        return 0;
    }

    void SDCTubeFlowSolidSolver::nextTimeStep()
    {}

    void SDCTubeFlowSolidSolver::initTimeStep()
    {
        TubeFlowSolidSolver::initTimeStep();
    }

    void SDCTubeFlowSolidSolver::setNumberOfImplicitStages( int /*k*/ )
    {}

    void SDCTubeFlowSolidSolver::implicitSolve(
        bool corrector,
        const int k,
        const int kold,
        const scalar t,
        const scalar dt,
        const fsi::vector & qold,
        const fsi::vector & rhs,
        fsi::vector & /*f*/,
        fsi::vector & /*result*/
        )
    {
        prepareImplicitSolve( corrector, k, kold, t, dt, qold, rhs );
    }

    void SDCTubeFlowSolidSolver::prepareImplicitSolve(
        bool /*corrector*/,
        const int /*k*/,
        const int /*kold*/,
        const scalar /*t*/,
        const scalar /*dt*/,
        const fsi::vector & /*qold*/,
        const fsi::vector & /*rhs*/
        )
    {}

    void SDCTubeFlowSolidSolver::getVariablesInfo(
        std::deque<int> & dof,
        std::deque<bool> & enabled,
        std::deque<std::string> & names
        )
    {
        dof.push_back( 0 );
        enabled.push_back( false );
        names.push_back( "solid" );
    }
}
