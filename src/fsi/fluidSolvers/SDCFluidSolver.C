
/*
 * Class
 *   SDCFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCFluidSolver.H"

SDCFluidSolver::SDCFluidSolver(
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    foamFluidSolver( name, args, runTime ),
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
    pimple( mesh.solutionDict().subDict( "PIMPLE" ) ),
    nu( transportProperties.lookup( "nu" ) ),
    rho( transportProperties.lookup( "rho" ) ),
    p
    (
    IOobject
    (
        "p",
        runTime->timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    ),
    U
    (
    IOobject
    (
        "U",
        runTime->timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    ),
    phi
    (
    IOobject
    (
        "phi",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate( U ) & mesh.Sf()
    ),
    Uf
    (
    IOobject
    (
        "Uf",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate( U )
    ),
    AU
    (
    IOobject
    (
        "AU",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1.0 / runTime->deltaT(),
    zeroGradientFvPatchScalarField::typeName
    ),
    HU
    (
    IOobject
    (
        "HU",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    U.dimensions() / runTime->deltaT().dimensions(),
    zeroGradientFvPatchVectorField::typeName
    ),
    rhsU
    (
    IOobject
    (
        "rhsU",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector( "rhsU", dimVelocity / dimTime, Foam::vector::zero ),
    zeroGradientFvPatchVectorField::typeName
    ),
    rhsUf
    (
    IOobject
    (
        "rhsUf",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate( rhsU )
    ),
    sumLocalContErr( 0 ),
    globalContErr( 0 ),
    cumulativeContErr( 0 ),
    convergenceTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "convergenceTolerance" ) ) ),
    k( 0 ),
    pStages(),
    phiStages(),
    UStages(),
    UfStages(),
    UF
    (
    IOobject
    (
        "UF",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "UF", dimVelocity / dimTime, Foam::vector::zero )
    ),
    UfF
    (
    IOobject
    (
        "UfF",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector( "UfF", dimVelocity / dimTime, Foam::vector::zero )
    )
{
    // Ensure that the absolute tolerance of the linear solver is less than the
    // used convergence tolerance for the non-linear system.
    scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
    assert( absTolerance < convergenceTolerance );

    absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "p" ).lookup( "tolerance" ) );
    assert( absTolerance < convergenceTolerance );

    initialize();
}

SDCFluidSolver::~SDCFluidSolver()
{}

void SDCFluidSolver::continuityErrs()
{
    volScalarField contErr = fvc::div( phi );

    sumLocalContErr = runTime->deltaT().value() *
        mag( contErr ) ().weightedAverage( mesh.V() ).value();

    globalContErr = runTime->deltaT().value() *
        contErr.weightedAverage( mesh.V() ).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr
         << endl;
}

void SDCFluidSolver::courantNo()
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;

    if ( mesh.nInternalFaces() )
    {
        surfaceScalarField magPhi = mag( phi );

        surfaceScalarField SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs() * magPhi;

        const scalar deltaT = runTime->deltaT().value();

        CoNum = max( SfUfbyDelta / mesh.magSf() ).value() * deltaT;

        meanCoNum = ( sum( SfUfbyDelta ) / sum( mesh.magSf() ) ).value() * deltaT;

        velMag = max( magPhi / mesh.magSf() ).value();
    }

    Info << "Courant Number mean: " << meanCoNum
         << " max: " << CoNum
         << " velocity magnitude: " << velMag
         << endl;
}

void SDCFluidSolver::createFields()
{
    U.oldTime();
    Uf.oldTime();
    phi.oldTime();

    // Read pressure properties and create turbulence model
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

double SDCFluidSolver::evaluateMomentumResidual( double dt )
{
    volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) - fvc::laplacian( nu, U ) + fvc::grad( p ) - rhsU / dt;

    scalarField magResU = mag( residual.internalField() );
    scalar momentumResidual = std::sqrt( gSumSqr( magResU ) / mesh.globalData().nTotalCells() );
    scalar rmsU = std::sqrt( gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells() );
    rmsU /= runTime->deltaT().value();

    // Scale the residual by the root mean square of the velocity field
    momentumResidual /= rmsU;

    return momentumResidual;
}

void SDCFluidSolver::getAcousticsDensityLocal( matrix & data )
{
    assert( false );
}

void SDCFluidSolver::getAcousticsVelocityLocal( matrix & data )
{
    assert( false );
}

void SDCFluidSolver::getAcousticsPressureLocal( matrix & data )
{
    assert( false );
}

void SDCFluidSolver::getTractionLocal( matrix & traction )
{
    assert( false );
}

void SDCFluidSolver::getWritePositionsLocalAcoustics( matrix & writePositions )
{
    assert( false );
}

void SDCFluidSolver::initialize()
{
    assert( !init );

    readPIMPLEControls();

    createFields();
}

void SDCFluidSolver::initTimeStep()
{
    assert( !init );

    readPIMPLEControls();

    init = true;
}

bool SDCFluidSolver::isRunning()
{
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void SDCFluidSolver::readPIMPLEControls()
{
    nOuterCorr = readInt( pimple.lookup( "nOuterCorrectors" ) );
    nCorr = readInt( pimple.lookup( "nCorrectors" ) );

    nNonOrthCorr = pimple.lookupOrDefault<int>( "nNonOrthogonalCorrectors", 0 );
}

void SDCFluidSolver::resetSolution()
{
    assert( false );
}

void SDCFluidSolver::setDeltaT( double dt )
{
    runTime->setDeltaT( dt );
}

void SDCFluidSolver::setNumberOfStages( int k )
{
    this->k = k;

    for ( int i = 0; i < k; i++ )
    {
        pStages.push_back( volScalarField( p ) );
        phiStages.push_back( surfaceScalarField( phi ) );
        UStages.push_back( volVectorField( U ) );
        UfStages.push_back( surfaceVectorField( Uf ) );
    }
}

void SDCFluidSolver::nextTimeStep()
{
    timeIndex++;

    if ( pStages.size() == static_cast<unsigned>(k) )
    {
        for ( int i = 0; i < k; i++ )
        {
            pStages.at( i ) = p;
            phiStages.at( i ) = phi;
            UStages.at( i ) = U;
            UfStages.at( i ) = Uf;
        }
    }
}

void SDCFluidSolver::solve()
{
    assert( false );
}

void SDCFluidSolver::finalizeTimeStep()
{
    foamFluidSolver::finalizeTimeStep();
}

int SDCFluidSolver::getDOF()
{
    return U.size() * 3 + Uf.size() * 3;
}

void SDCFluidSolver::getSolution( Eigen::VectorXd & solution )
{
    for ( int i = 0; i < U.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            solution( i * 3 + j ) = U[i][j];

    for ( int i = 0; i < Uf.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            solution( i * 3 + j + U.size() * 3 ) = Uf[i][j];
}

double SDCFluidSolver::getTimeStep()
{
    return runTime->deltaT().value();
}

void SDCFluidSolver::evaluateFunction(
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    )
{
    for ( int i = 0; i < UF.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            f( i * 3 + j ) = UF[i][j];

    for ( int i = 0; i < UfF.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            f( i * 3 + j + UF.size() * 3 ) = UfF[i][j];
}

void SDCFluidSolver::implicitSolve(
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
    runTime->setDeltaT( dt );

    if ( corrector )
    {
        p = pStages.at( k + 1 );
        phi = phiStages.at( k + 1 );
        U = UStages.at( k + 1 );
        Uf = UfStages.at( k + 1 );
    }
    else
    {
        // predictor
        p = pStages.at( k );
        phi = phiStages.at( k );
        U = UStages.at( k );
        Uf = UfStages.at( k );
    }

    for ( int i = 0; i < U.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            U.oldTime()[i][j] = qold( i * 3 + j );

    for ( int i = 0; i < Uf.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            Uf.oldTime()[i][j] = qold( i * 3 + j + U.size() * 3 );

    for ( int i = 0; i < U.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            rhsU[i][j] = rhs( i * 3 + j );

    for ( int i = 0; i < Uf.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            rhsUf[i][j] = rhs( i * 3 + j + U.size() * 3 );

    courantNo();

    if ( runTime->timeIndex() > 1 || k > 0 || corrector )
    {
        double initMomentumResidual = evaluateMomentumResidual( dt );
        convergenceTolerance = std::max( 1.0e-2 * initMomentumResidual, 1.0e-15 );

        Info << "root mean square residual norm = " << initMomentumResidual;
        Info << ", tolerance = " << convergenceTolerance;
        Info << ", iteration = 0, convergence = false" << endl;
    }

    // PIMPLE algorithm

    for ( label oCorr = 0; oCorr < nOuterCorr; oCorr++ )
    {
        // for computing Hp
        fvVectorMatrix UEqn
        (
            fvm::div( phi, U )
            - fvm::laplacian( nu, U )
        );

        // for computing \tilde{a}_{P}
        fvVectorMatrix UEqnt
        (
            fvm::ddt( U )
            + fvm::div( phi, U )
            - fvm::laplacian( nu, U )
        );

        {
            // To ensure S0 and B0 are thrown out of memory
            // Source and boundaryCoeffs need to be saved when relexation is applied
            // to still obtain time consistent behavior.
            // Only source is affected by relaxation, boundaryCoeffs is not relaxation
            // dependent.
            // BoundaryCoeffs needs to be saved to generate the correct UEqn after
            // solving. Explicit terms (depending on U(n)) need to remain depending
            // on U(n) and not on new solution)
            vectorField S0 = UEqnt.source();
            FieldField<Field, Foam::vector> B0 = UEqnt.boundaryCoeffs();

            UEqnt.relax();

            // Eqn 25
            Foam::solve( UEqnt == -fvc::grad( p ) + rhsU / dt );

            // Reset equation to ensure relaxation parameter is not causing problems for time order
            UEqnt =
                (
                fvm::ddt( U )
                + fvm::div( phi, U )
                - fvm::laplacian( nu, U )
                );
            UEqnt.source() = S0;
            UEqnt.boundaryCoeffs() = B0;
        }

        // --- PISO loop

        // Relative convergence measure for the PISO loop:
        // Perform at maximum nCorr PISO corrections.
        // If the relative residual with respect to the initial
        // residual is decreased by factor tol: assume convergence.

        double initResidual = 1;
        double currResidual = 1;
        double pressureResidual = 1;
        double tol = 1.0e-2;

        for ( label corr = 0; corr < nCorr; corr++ )
        {
            p.storePrevIter();

            volVectorField HU = UEqn.H(); // Eqn 6
            volScalarField AU = UEqnt.A(); // Eqn 14

            dimensionedScalar rDeltaT = 1.0 / runTime->deltaT();

            // Eqn 15
            U = HU / AU + rDeltaT / AU * U.oldTime() + 1 / AU * rhsU / dt;

            // Eqn 16
            Uf = fvc::interpolate( HU ) / fvc::interpolate( AU ) + rDeltaT / fvc::interpolate( AU ) * Uf.oldTime() + 1 / fvc::interpolate( AU ) * rhsUf / dt;

            {
                forAll( Uf.boundaryField(), patchI )
                {
                    if ( !Uf.boundaryField()[patchI].coupled() )
                    {
                        Uf.boundaryField()[patchI] =
                            (
                            U.boundaryField()[patchI]
                            );
                    }
                }

                phi = Uf & mesh.Sf();

                forAll( phi.boundaryField(), patchI )
                {
                    if ( !phi.boundaryField()[patchI].coupled() )
                    {
                        phi.boundaryField()[patchI] =
                            (
                            U.boundaryField()[patchI]
                            & mesh.Sf().boundaryField()[patchI]
                            );
                    }
                }
            }

            for ( label nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        1.0 / fvc::interpolate( AU ), p,
                        "laplacian((1|A(U)),p)"
                    )
                    == fvc::div( phi )
                );

                pEqn.setReference( pRefCell, pRefValue );
                pressureResidual = pEqn.solve().initialResidual();

                if ( corr == 0 && nonOrth == 0 )
                    initResidual = pressureResidual;
                else
                if ( nonOrth == 0 )
                    currResidual = pressureResidual;

                if ( nonOrth == nNonOrthCorr )
                    phi -= pEqn.flux();
            }

            p.relax();

            // Eqn 30
            U -= fvc::grad( p ) / AU;

            U.correctBoundaryConditions();

            if ( currResidual < std::max( tol * initResidual, 1.0e-15 ) )
                break;
        }

        surfaceVectorField nf = mesh.Sf() / mesh.magSf();
        surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
        surfaceVectorField Unor = phi / mesh.magSf() * nf;

        Uf = Utang + Unor;

        scalar momentumResidual = evaluateMomentumResidual( dt );

        bool convergence = momentumResidual <= convergenceTolerance;

        Info << "root mean square residual norm = " << momentumResidual;
        Info << ", tolerance = " << convergenceTolerance;
        Info << ", iteration = " << oCorr + 1;
        Info << ", convergence = ";

        if ( convergence )
            Info << "true";
        else
            Info << "false";

        Info << endl;

        if ( convergence )
            break;
    }

    continuityErrs();

    pStages.at( k + 1 ) = p;
    phiStages.at( k + 1 ) = phi;
    UStages.at( k + 1 ) = U;
    UfStages.at( k + 1 ) = Uf;

    for ( int i = 0; i < U.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            result( i * 3 + j ) = U[i][j];

    for ( int i = 0; i < Uf.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            result( i * 3 + j + U.size() * 3 ) = Uf[i][j];

    dimensionedScalar rDeltaT = 1.0 / runTime->deltaT();
    UF = rDeltaT * ( U - U.oldTime() ) - rhsU / dt;
    UfF = rDeltaT * ( Uf - Uf.oldTime() ) - rhsUf / dt;

    for ( int i = 0; i < UF.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            f( i * 3 + j ) = UF[i][j];

    for ( int i = 0; i < UfF.size(); i++ )
        for ( int j = 0; j < 3; j++ )
            f( i * 3 + j + UF.size() * 3 ) = UfF[i][j];
}

int SDCFluidSolver::getNbCells()
{
    return mesh.globalData().nTotalCells();
}

double SDCFluidSolver::getScalingFactor()
{
    scalar rmsU = gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells();
    rmsU += gSumSqr( mag( Uf.internalField() ) ) / mesh.globalData().nTotalFaces();
    rmsU = std::sqrt( rmsU );
    rmsU /= runTime->deltaT().value();

    return rmsU;
}
