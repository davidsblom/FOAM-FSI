
/*
 * Class
 *   SDCLaplacianSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCLaplacianSolver.H"

SDCLaplacianSolver::SDCLaplacianSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    )
    :
    init( false ),
    name( name ),
    args( args ),
    runTime( runTime ),
    mesh
    (
    Foam::IOobject
    (
        name,
        runTime->timeName(),
        *runTime,
        Foam::IOobject::MUST_READ
    )
    ),
    T
    (
    IOobject
    (
        "T",
        runTime->timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    ),
    rhsT
    (
    IOobject
    (
        "rhsT",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "rhsT", dimTemperature / dimTime, 0 )
    ),
    transportProperties
    (
    IOobject
    (
        "transportProperties",
        runTime->constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
    ),
    DT
    (
    transportProperties.lookup( "DT" )
    ),
    nNonOrthCorr( mesh.solutionDict().subDict( "SIMPLE" ).lookupOrDefault<int>("nNonOrthogonalCorrectors", 0) )
{
    assert( args );
    assert( runTime );
}

SDCLaplacianSolver::~SDCLaplacianSolver()
{}

void SDCLaplacianSolver::evaluateFunction(
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    )
{
    for ( int i = 0; i < T.size(); i++ )
        T[i] = q( i );

    volScalarField F = fvc::laplacian( DT, T );

    for ( int i = 0; i < F.size(); i++ )
        f( i ) = F[i];
}

void SDCLaplacianSolver::finalizeTimeStep()
{
    assert( init );

    init = false;
}

int SDCLaplacianSolver::getDOF()
{
    return T.size();
}

int SDCLaplacianSolver::getNbCells()
{
    return mesh.globalData().nTotalCells();
}

double SDCLaplacianSolver::getScalingFactor()
{
    return 1;
}

void SDCLaplacianSolver::getSolution( Eigen::VectorXd & solution )
{
    for ( int i = 0; i < T.size(); i++ )
        solution( i ) = T[i];
}

void SDCLaplacianSolver::setSolution(
    const Eigen::VectorXd & solution,
    const Eigen::VectorXd & f
    )
{
    for ( int i = 0; i < T.size(); i++ )
        T[i] = solution( i );
}

double SDCLaplacianSolver::getEndTime()
{
    return runTime->endTime().value() - runTime->startTime().value();
}

double SDCLaplacianSolver::getTimeStep()
{
    return runTime->deltaT().value();
}

void SDCLaplacianSolver::nextTimeStep()
{
    if ( TStages.size() == static_cast<unsigned>(k) )
        for ( int i = 0; i < k; i++ )
            TStages.at( i ) = T;
}

void SDCLaplacianSolver::initTimeStep()
{
    assert( !init );

    init = true;
}

bool SDCLaplacianSolver::isRunning()
{
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void SDCLaplacianSolver::setDeltaT( double dt )
{
    runTime->setDeltaT( dt );
}

void SDCLaplacianSolver::setNumberOfStages( int k )
{
    this->k = k;

    for ( int i = 0; i < k; i++ )
        TStages.push_back( volScalarField( T ) );
}

void SDCLaplacianSolver::implicitSolve(
    bool corrector,
    const int k,
    const double t,
    const double dt,
    const Eigen::VectorXd & qold,
    const Eigen::VectorXd & rhs,
    Eigen::VectorXd & f,
    Eigen::VectorXd & result
    )
{
    assert( init );

    runTime->setDeltaT( dt );

    if ( corrector )
        T = TStages.at( k + 1 );
    else
        // predictor
        T = TStages.at( k );

    for ( int i = 0; i < T.size(); i++ )
        T.oldTime()[i] = qold( i );

    for ( int i = 0; i < T.size(); i++ )
        rhsT[i] = rhs( i );

    for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
    {
        solve
        (
            fvm::ddt( T ) - fvm::laplacian( DT, T ) - rhsT / dt
        );
    }

    TStages.at( k + 1 ) = T;

    for ( int i = 0; i < T.size(); i++ )
        result( i ) = T[i];

    volScalarField F = fvc::laplacian( DT, T );

    for ( int i = 0; i < F.size(); i++ )
        f( i ) = F[i];
}
