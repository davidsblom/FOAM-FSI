
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

template <int dimension>
dealiiSolidSolver<dimension>::dealiiSolidSolver( dealiifsi::DataStorage & data )
    :
    dealiifsi::LinearElasticity<dimension>( data ),
    BaseMultiLevelSolver( 5, dimension, 0 ),
    k( 0 ),
    UStages(),
    VStages()
{
    matrix tmp;
    getWritePositions( tmp );

    BaseMultiLevelSolver::N = tmp.rows();
    BaseMultiLevelSolver::data.resize( N, dimension );
    BaseMultiLevelSolver::data.setZero();
}

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
    dealiifsi::LinearElasticity<dimension>( time_step, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines ),
    BaseMultiLevelSolver( 5, dimension, 0 ),
    k( 0 ),
    UStages(),
    VStages()
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

    if ( UStages.size() > 0 )
    {
        UStages.at( kindex + 1 ) = dealiifsi::LinearElasticity<dimension>::solution_u;
        VStages.at( kindex + 1 ) = dealiifsi::LinearElasticity<dimension>::solution_v;
    }
}

template <int dimension>
void dealiiSolidSolver<dimension>::evaluateFunction(
    const int k,
    const fsi::vector & q,
    const scalar t,
    fsi::vector & f
    )
{}

template <int dimension>
int dealiiSolidSolver<dimension>::getDOF()
{
    return dealiifsi::LinearElasticity<dimension>::solution_u.size() + dealiifsi::LinearElasticity<dimension>::solution_v.size();
}

template <int dimension>
void dealiiSolidSolver<dimension>::getSolution(
    fsi::vector & solution,
    fsi::vector & f
    )
{
    copy( dealiifsi::LinearElasticity<dimension>::solution_u, solution, 0 );
    copy( dealiifsi::LinearElasticity<dimension>::solution_v, solution, dealiifsi::LinearElasticity<dimension>::solution_u.size() );
    copy( dealiifsi::LinearElasticity<dimension>::u_f, f, 0 );
    copy( dealiifsi::LinearElasticity<dimension>::v_f, f, dealiifsi::LinearElasticity<dimension>::u_f.size() );
}

template <int dimension>
void dealiiSolidSolver<dimension>::setSolution(
    const fsi::vector & solution,
    const fsi::vector & f
    )
{
    copy( solution, dealiifsi::LinearElasticity<dimension>::solution_u, 0 );
    copy( solution, dealiifsi::LinearElasticity<dimension>::solution_v, dealiifsi::LinearElasticity<dimension>::solution_u.size() );
    copy( f, dealiifsi::LinearElasticity<dimension>::u_f, 0 );
    copy( f, dealiifsi::LinearElasticity<dimension>::v_f, dealiifsi::LinearElasticity<dimension>::u_f.size() );
}

template <int dimension>
scalar dealiiSolidSolver<dimension>::getEndTime()
{
    return dealiifsi::LinearElasticity<dimension>::final_time;
}

template <int dimension>
scalar dealiiSolidSolver<dimension>::getTimeStep()
{
    return dealiifsi::LinearElasticity<dimension>::time_step;
}

template <int dimension>
void dealiiSolidSolver<dimension>::nextTimeStep()
{
    for ( int i = 0; i < k; i++ )
    {
        UStages.at( i ) = dealiifsi::LinearElasticity<dimension>::solution_u;
        VStages.at( i ) = dealiifsi::LinearElasticity<dimension>::solution_v;
    }
}

template <int dimension>
void dealiiSolidSolver<dimension>::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    UStages.clear();
    VStages.clear();

    for ( int i = 0; i < k + 1; i++ )
    {
        UStages.push_back( dealiifsi::LinearElasticity<dimension>::solution_u );
        VStages.push_back( dealiifsi::LinearElasticity<dimension>::solution_v );
    }
}

template <int dimension>
void dealiiSolidSolver<dimension>::prepareImplicitSolve(
    bool corrector,
    const int k,
    const int kold,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs
    )
{
    dealiifsi::LinearElasticity<dimension>::time_step = dt;
    dealiifsi::LinearElasticity<dimension>::time = t;

    if ( corrector )
    {
        dealiifsi::LinearElasticity<dimension>::solution_u = UStages.at( k + 1 );
        dealiifsi::LinearElasticity<dimension>::solution_v = VStages.at( k + 1 );
    }

    copy( qold, dealiifsi::LinearElasticity<dimension>::old_solution_u, 0 );
    copy( qold, dealiifsi::LinearElasticity<dimension>::old_solution_v, dealiifsi::LinearElasticity<dimension>::old_solution_u.size() );
    copy( rhs, dealiifsi::LinearElasticity<dimension>::u_rhs, 0 );
    copy( rhs, dealiifsi::LinearElasticity<dimension>::v_rhs, dealiifsi::LinearElasticity<dimension>::u_rhs.size() );
}

template <int dimension>
void dealiiSolidSolver<dimension>::implicitSolve(
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
    kindex = k;
    prepareImplicitSolve( corrector, k, kold, t, dt, qold, rhs );

    dealiifsi::LinearElasticity<dimension>::solve();

    UStages.at( k + 1 ) = dealiifsi::LinearElasticity<dimension>::solution_u;
    VStages.at( k + 1 ) = dealiifsi::LinearElasticity<dimension>::solution_v;

    getSolution( result, f );
}

template <int dimension>
void dealiiSolidSolver<dimension>::getVariablesInfo(
    std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    )
{
    dof.push_back( dealiifsi::LinearElasticity<dimension>::solution_u.size() );
    dof.push_back( dealiifsi::LinearElasticity<dimension>::solution_v.size() );
    enabled.push_back( true );
    enabled.push_back( true );
    names.push_back( "solid U" );
    names.push_back( "solid V" );
}

void copy(
    const dealii::Vector<scalar> & source,
    fsi::vector & target,
    unsigned int targetOffset
    )
{
    for ( unsigned int i = 0; i < source.size(); ++i )
        target( i + targetOffset ) = source[i];
}

void copy(
    const fsi::vector & source,
    dealii::Vector<scalar> & target,
    unsigned int sourceOffset
    )
{
    for ( unsigned int i = 0; i < target.size(); ++i )
        target[i] = source( i + sourceOffset );
}
