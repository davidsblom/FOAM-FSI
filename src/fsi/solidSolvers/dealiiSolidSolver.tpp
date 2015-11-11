
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

template <int dimension>
dealiiSolidSolver<dimension>::dealiiSolidSolver(
    double time_step,
    double final_time,
    double theta,
    double degree,
    double gravity,
    double distributed_load,
    double rho,
    double E,
    double nu,
    unsigned int n_global_refines
    )
    :
    Step23::LinearElasticity<dimension>( time_step, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines ),
    BaseMultiLevelSolver( 5, dimension, 0 )
{
    matrix tmp;
    getWritePositions( tmp );

    BaseMultiLevelSolver::N = tmp.rows();
    BaseMultiLevelSolver::data.resize( N, dimension );
    BaseMultiLevelSolver::data.setZero();
}

template <int dimension>
dealiiSolidSolver<dimension>::~dealiiSolidSolver()
{}

template <int dimension>
void dealiiSolidSolver<dimension>::finalizeTimeStep()
{
    assert( BaseMultiLevelSolver::init );

    Step23::LinearElasticity<dimension>::finalizeTimeStep();

    BaseMultiLevelSolver::init = false;
}

template <int dimension>
void dealiiSolidSolver<dimension>::getReadPositions( matrix & readPositions )
{
    Step23::LinearElasticity<dimension>::getReadPositions( readPositions );
}

template <int dimension>
void dealiiSolidSolver<dimension>::getWritePositions( matrix & writePositions )
{
    Step23::LinearElasticity<dimension>::getWritePositions( writePositions );
}

template <int dimension>
void dealiiSolidSolver<dimension>::initTimeStep()
{
    assert( !BaseMultiLevelSolver::init );

    Step23::LinearElasticity<dimension>::initTimeStep();

    BaseMultiLevelSolver::init = true;
}

template <int dimension>
bool dealiiSolidSolver<dimension>::isRunning()
{
    return Step23::LinearElasticity<dimension>::isRunning();
}

template <int dimension>
void dealiiSolidSolver<dimension>::resetSolution()
{}

template <int dimension>
void dealiiSolidSolver<dimension>::solve(
    const matrix & input,
    matrix & output
    )
{
    Info << "Solve solid domain with deal.II" << endl;

    Step23::LinearElasticity<dimension>::setTraction( input );

    Step23::LinearElasticity<dimension>::solve();

    Step23::LinearElasticity<dimension>::getDisplacement( output );

    BaseMultiLevelSolver::data = output;
}
