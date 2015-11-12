
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

template <int dimension>
dealiiSolidSolver<dimension>::dealiiSolidSolver( dealiifsi::DataStorage & data )
    :
    dealiifsi::LinearElasticity<dimension>( data ),
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

    dealiifsi::LinearElasticity<dimension>::finalizeTimeStep();

    BaseMultiLevelSolver::init = false;
}

template <int dimension>
void dealiiSolidSolver<dimension>::getReadPositions( matrix & readPositions )
{
    dealiifsi::LinearElasticity<dimension>::getReadPositions( readPositions );
}

template <int dimension>
void dealiiSolidSolver<dimension>::getWritePositions( matrix & writePositions )
{
    dealiifsi::LinearElasticity<dimension>::getWritePositions( writePositions );
}

template <int dimension>
void dealiiSolidSolver<dimension>::initTimeStep()
{
    assert( !BaseMultiLevelSolver::init );

    dealiifsi::LinearElasticity<dimension>::initTimeStep();

    BaseMultiLevelSolver::init = true;
}

template <int dimension>
bool dealiiSolidSolver<dimension>::isRunning()
{
    return dealiifsi::LinearElasticity<dimension>::isRunning();
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

    dealiifsi::LinearElasticity<dimension>::setTraction( input );

    dealiifsi::LinearElasticity<dimension>::solve();

    dealiifsi::LinearElasticity<dimension>::getDisplacement( output );

    BaseMultiLevelSolver::data = output;
}
