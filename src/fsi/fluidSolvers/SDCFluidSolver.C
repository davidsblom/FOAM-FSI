
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
    Uf & mesh.Sf()
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
    dimensionedVector( "rhsU", dimVelocity, Foam::vector::zero )
    ),
    rhsPhi
    (
    IOobject
    (
        "rhsPhi",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate( rhsU ) & mesh.Sf()
    ),
    nCorr( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "nCorrectors" ) ) ),
    nNonOrthCorr( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "nNonOrthogonalCorrectors" ) ) ),
    minIter( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "minIter" ) ) ),
    maxIter( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "maxIter" ) ) ),
    absoluteTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "tolerance" ) ) ),
    relativeTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "relTol" ) ) ),
    sumLocalContErr( 0 ),
    globalContErr( 0 ),
    cumulativeContErr( 0 ),
    laminarTransport( U, phi ),
    turbulence( autoPtr<incompressible::turbulenceModel>
    (
        incompressible::turbulenceModel::New( U, phi, laminarTransport )
    ) ),
    k( 0 ),
    pStages(),
    phiStages(),
    UStages(),
    UfStages(),
    UFHeader
    (
    "UF",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
    phiFHeader
    (
    "phiF",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
    UF
    (
    UFHeader,
    mesh,
    dimensionedVector( "UF", dimVelocity / dimTime, Foam::vector::zero )
    ),
    phiF
    (
    phiFHeader,
    fvc::interpolate( UF ) & mesh.Sf()
    ),
    turbulenceSwitch( true ),
    explicitFirstStage( false )
{
    if ( UFHeader.headerOk() && phiFHeader.headerOk() )
        explicitFirstStage = true;

    // Ensure that the absolute tolerance of the linear solver is less than the
    // used convergence tolerance for the non-linear system.
    scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "p" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    {
        IOdictionary dict
        (
            IOobject
            (
                "turbulenceProperties",
                runTime->constant(),
                *runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if ( word( dict.lookup( "simulationType" ) ) == "laminar" )
            turbulenceSwitch = false;

        if ( word( dict.lookup( "simulationType" ) ) == "RASModel" )
        {
            IOdictionary dict
            (
                IOobject
                (
                    "RASProperties",
                    runTime->constant(),
                    *runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            if ( word( dict.lookup( "RASModel" ) ) == "laminar" )
                turbulenceSwitch = false;
        }
    }

    Info << "Turbulence ";

    if ( turbulenceSwitch )
        Info << "enabled";
    else
        Info << "disabled";

    Info << endl;

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

    surfaceVectorField nf = mesh.Sf() / mesh.magSf();
    surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
    surfaceVectorField Unor = phi / mesh.magSf() * nf;

    Uf = Utang + Unor;

    // Read pressure properties and create turbulence model
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

double SDCFluidSolver::evaluateMomentumResidual()
{
    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();
    volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + fvc::grad( p ) - rDeltaT * rhsU;

    if ( turbulenceSwitch )
        residual += turbulence->divDevReff( U ) & U;
    else
        residual += -fvc::laplacian( nu, U );

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

    createFields();
}

void SDCFluidSolver::initTimeStep()
{
    assert( !init );

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
    int index = 0;
    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                index++;
            }
        }
    }

    forAll( phi.internalField(), i )
    {
        index++;
    }

    forAll( phi.boundaryField(), patchI )
    {
        forAll( phi.boundaryField()[patchI], i )
        {
            index++;
        }
    }

    return index;
}

void SDCFluidSolver::getSolution( Eigen::VectorXd & solution )
{
    int index = 0;

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            solution( index ) = U.internalField()[i][j];
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                solution( index ) = U.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    forAll( phi.internalField(), i )
    {
        solution( index ) = phi[i];
        index++;
    }

    forAll( phi.boundaryField(), patchI )
    {
        forAll( phi.boundaryField()[patchI], i )
        {
            solution( index ) = phi.boundaryField()[patchI][i];
            index++;
        }
    }

    assert( index == solution.rows() );
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
    if ( explicitFirstStage )
    {
        UF = -fvc::div( phi, U ) - fvc::grad( p ) + fvc::laplacian( nu, U );
        phiF = fvc::interpolate( UF ) & mesh.Sf();
        explicitFirstStage = false;
    }

    int index = 0;

    forAll( UF.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            f( index ) = UF.internalField()[i][j];
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
    {
        forAll( UF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                f( index ) = UF.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    forAll( phiF.internalField(), i )
    {
        f( index ) = phiF[i];
        index++;
    }

    forAll( phiF.boundaryField(), patchI )
    {
        forAll( phiF.boundaryField()[patchI], i )
        {
            f( index ) = phiF.boundaryField()[patchI][i];
            index++;
        }
    }

    assert( index == f.rows() );
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
    bool convergence = false;
    runTime->setDeltaT( dt );

    double told = runTime->value();
    runTime->setTime( t, runTime->timeIndex() );

    if ( corrector )
    {
        p = pStages.at( k + 1 );
        phi = phiStages.at( k + 1 );
        U = UStages.at( k + 1 );
    }

    int index = 0;

    forAll( U.oldTime().internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            U.oldTime().internalField()[i][j] = qold( index );
            index++;
        }
    }

    forAll( U.oldTime().boundaryField(), patchI )
    {
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                U.oldTime().boundaryField()[patchI][i][j] = qold( index );
                index++;
            }
        }
    }

    forAll( phi.oldTime().internalField(), i )
    {
        phi.oldTime().internalField()[i] = qold( index );
        index++;
    }

    forAll( phi.oldTime().boundaryField(), patchI )
    {
        forAll( phi.boundaryField()[patchI], i )
        {
            phi.oldTime().boundaryField()[patchI][i] = qold( index );
            index++;
        }
    }

    index = 0;

    forAll( rhsU.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            rhsU.internalField()[i][j] = rhs( index );
            index++;
        }
    }

    forAll( rhsU.boundaryField(), patchI )
    {
        forAll( rhsU.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                rhsU.boundaryField()[patchI][i][j] = rhs( index );
                index++;
            }
        }
    }

    forAll( rhsPhi.internalField(), i )
    {
        rhsPhi.internalField()[i] = rhs( index );
        index++;
    }

    forAll( rhsPhi.boundaryField(), patchI )
    {
        forAll( rhsPhi.boundaryField()[patchI], i )
        {
            rhsPhi.boundaryField()[patchI][i] = rhs( index );
            index++;
        }
    }

    // -------------------------------------------------------------------------

    scalar convergenceTolerance = absoluteTolerance;

    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();

    // --- PIMPLE loop
    for ( label oCorr = 0; oCorr < maxIter; oCorr++ )
    {
        U.storePrevIter();

        fvVectorMatrix UEqn
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
            vectorField S0 = UEqn.source();
            FieldField<Field, Foam::vector> B0 = UEqn.boundaryCoeffs();

            UEqn.relax();

            Foam::solve( UEqn == -fvc::grad( p ) + rDeltaT * rhsU );

            // Reset equation to ensure relaxation parameter is not causing problems for time order
            UEqn =
                (
                fvm::ddt( U )
                + fvm::div( phi, U )
                - fvm::laplacian( nu, U )
                );

            UEqn.source() = S0;
            UEqn.boundaryCoeffs() = B0;
        }

        // Relative convergence measure for the PISO loop:
        // Perform at maximum nCorr PISO corrections.
        // If the relative residual with respect to the initial
        // residual is decreased by factor tol: assume convergence.

        double initResidual = 1;
        double currResidual = 1;
        double pressureResidual = 1;
        double tol = 1.0e-2;

        // --- PISO loop
        for ( int corr = 0; corr < nCorr; corr++ )
        {
            p.storePrevIter();

            HU = UEqn.H();
            AU = UEqn.A();
            U = HU / AU + 1 / AU * rDeltaT * rhsU;

            {
                phi = ( fvc::interpolate( HU ) / fvc::interpolate( AU ) ) & mesh.Sf();

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

                // === Set boundaries correct of U === //
                surfaceScalarField ddtPhiCoeff
                (
                    IOobject
                    (
                        "ddtPhiCoeff",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensioned<scalar>( "1", dimless, 1.0 )
                );

                forAll( U.boundaryField(), patchI )
                {
                    if ( U.boundaryField()[patchI].fixesValue() )
                    {
                        ddtPhiCoeff.boundaryField()[patchI] = 0.0;
                    }
                }

                surfaceScalarField phi0 = phi.oldTime() - ( fvc::interpolate( U.oldTime() ) & mesh.Sf() ) + rhsPhi;

                phi += rDeltaT * ddtPhiCoeff * phi0 / fvc::interpolate( AU );
            }

            for ( int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++ )
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian( 1.0 / fvc::interpolate( AU ), p, "laplacian((1|A(U)),p)" ) == fvc::div( phi )
                );

                pEqn.setReference( pRefCell, pRefValue );

                pressureResidual = pEqn.solve().initialResidual();

                if ( corr == 0 && nonOrth == 0 )
                    initResidual = pressureResidual;
                else
                if ( nonOrth == 0 )
                    currResidual = pressureResidual;

                if ( nonOrth == nNonOrthCorr )
                {
                    phi -= pEqn.flux();
                }
            }

            p.relax();

            U -= (1.0 / AU) * fvc::grad( p );
            U.correctBoundaryConditions();

            if ( currResidual < std::max( tol * initResidual, 1.0e-15 ) )
                break;
        }

        scalar momentumResidual = evaluateMomentumResidual();

        if ( oCorr == 0 )
            convergenceTolerance = std::max( relativeTolerance * momentumResidual, absoluteTolerance );

        convergence = momentumResidual <= convergenceTolerance && oCorr >= minIter - 1;

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

    // -------------------------------------------------------------------------

    if ( static_cast<int>( pStages.size() ) >= k + 2 )
    {
        pStages.at( k + 1 ) = p;
        phiStages.at( k + 1 ) = phi;
        UStages.at( k + 1 ) = U;
    }

    UF = rDeltaT * (U - U.oldTime() - rhsU);
    phiF = rDeltaT * (phi - phi.oldTime() - rhsPhi);

    getSolution( result );
    evaluateFunction( k + 1, qold, t, f );

    runTime->setTime( told, runTime->timeIndex() );
}

double SDCFluidSolver::getScalingFactor()
{
    return 1;
    scalar rmsU = gSumSqr( mag( U.internalField() ) ) / mesh.globalData().nTotalCells();
    rmsU += gSumSqr( mag( phi.internalField() ) ) / mesh.globalData().nTotalFaces();
    rmsU = std::sqrt( rmsU );
    rmsU /= runTime->deltaT().value();

    return rmsU;
}
