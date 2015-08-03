
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
    nNonOrthCorr( mesh.solutionDict().subDict( "SIMPLE" ).lookupOrDefault<int>("nNonOrthogonalCorrectors", 0) ),
    k( 0 )
{
    assert( args );
    assert( runTime );
}

SDCLaplacianSolver::~SDCLaplacianSolver()
{}

void SDCLaplacianSolver::evaluateFunction(
    const int k,
    const fsi::vector & q,
    const scalar t,
    fsi::vector & f
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

scalar SDCLaplacianSolver::getScalingFactor()
{
    return 1;
}

void SDCLaplacianSolver::getSolution( fsi::vector & solution )
{
    for ( int i = 0; i < T.size(); i++ )
        solution( i ) = T[i];
}

void SDCLaplacianSolver::setSolution(
    const fsi::vector & solution,
    const fsi::vector & f
    )
{
    for ( int i = 0; i < T.size(); i++ )
        T[i] = solution( i );
}

scalar SDCLaplacianSolver::getEndTime()
{
    return runTime->endTime().value() - runTime->startTime().value();
}

scalar SDCLaplacianSolver::getTimeStep()
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

void SDCLaplacianSolver::setDeltaT( scalar dt )
{
    runTime->setDeltaT( dt );
}

void SDCLaplacianSolver::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    for ( int i = 0; i < k + 1; i++ )
        TStages.push_back( volScalarField( T ) );
}

void SDCLaplacianSolver::implicitSolve(
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
