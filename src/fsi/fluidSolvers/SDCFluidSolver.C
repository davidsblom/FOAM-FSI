
/*
 * Class
 *   SDCFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCFluidSolver.H"

SDCFluidSolver::SDCFluidSolver(
    std::string name,
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
    pRefCell( 0 ),
    pRefValue( 0.0 ),
    laminarTransport( U, phi ),
    turbulence( autoPtr<incompressible::turbulenceModel>
        (
            incompressible::turbulenceModel::New( U, phi, laminarTransport )
        ) ),
    k( 0 ),
    pStages(),
    phiStages(),
    UStages(),
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
    explicitFirstStage( true ),
    writeInterval( readLabel( runTime->controlDict().lookup( "writeInterval" ) ) )
{
    if ( UFHeader.headerOk() && phiFHeader.headerOk() )
        explicitFirstStage = false;

    // Ensure that the absolute tolerance of the linear solver is less than the
    // used convergence tolerance for the non-linear system.
    scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    if ( absTolerance >= absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver of the momentum equation (U) should be smaller than PIMPLE::tolerance in order to reach convergence of the non-linear system" );

    absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "p" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    if ( absTolerance >= absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver of the Poisson equation (p) should be smaller than PIMPLE::tolerance in order to reach convergence of the non-linear system" );

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

    const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
    const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );
    word ddtScheme;

    if ( ddtSchemes.found( "ddt(U)" ) )
        ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
    else
        ddtScheme = word( ddtSchemes.lookup( "default" ) );

    assert( ddtScheme == "bdf1" );
    assert( word( runTime->controlDict().lookup( "writeControl" ) ) == word( "timeStep" ) );

    initialize();
}

SDCFluidSolver::~SDCFluidSolver()
{}

void SDCFluidSolver::continuityErrs()
{
    tmp<volScalarField> contErr = fvc::div( phi );

    sumLocalContErr = runTime->deltaT().value() *
        mag( contErr() ) ().weightedAverage( mesh.V() ).value();

    globalContErr = runTime->deltaT().value() *
        contErr->weightedAverage( mesh.V() ).value();

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
        tmp<surfaceScalarField> magPhi = mag( phi );

        surfaceScalarField SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs() * magPhi();

        const scalar deltaT = runTime->deltaT().value();

        CoNum = max( SfUfbyDelta / mesh.magSf() ).value() * deltaT;

        meanCoNum = ( sum( SfUfbyDelta ) / sum( mesh.magSf() ) ).value() * deltaT;

        velMag = max( magPhi() / mesh.magSf() ).value();
    }

    Info << "Courant Number mean: " << meanCoNum
         << " max: " << CoNum
         << " velocity magnitude: " << velMag
         << endl;
}

void SDCFluidSolver::createFields()
{
    U.oldTime();
    phi.oldTime();

    // Read pressure properties and create turbulence model
    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

scalar SDCFluidSolver::evaluateMomentumResidual()
{
    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();
    tmp<volVectorField> residual = fvc::ddt( U ) + fvc::div( phi, U ) + fvc::grad( p ) - rDeltaT * rhsU;

    if ( turbulenceSwitch )
        residual() += turbulence->divDevReff() & U;
    else
        residual() += -fvc::laplacian( nu, U );

    scalarField magResU = mag( residual->internalField() );
    return std::sqrt( gSumSqr( magResU ) / ( mesh.globalData().nTotalCells() * mesh.nGeometricD() ) );
}

void SDCFluidSolver::getAcousticsDensityLocal( matrix & )
{
    assert( false );
}

void SDCFluidSolver::getAcousticsVelocityLocal( matrix & )
{
    assert( false );
}

void SDCFluidSolver::getAcousticsPressureLocal( matrix & )
{
    assert( false );
}

void SDCFluidSolver::getTractionLocal( matrix & )
{
    assert( false );
}

void SDCFluidSolver::getWritePositionsLocalAcoustics( matrix & )
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

    t = runTime->time().value();

    init = true;
}

bool SDCFluidSolver::isRunning()
{
    return runTime->run();
}

void SDCFluidSolver::resetSolution()
{
    assert( false );
}

void SDCFluidSolver::setDeltaT( scalar dt )
{
    runTime->setDeltaT( dt );
}

void SDCFluidSolver::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    pStages.clear();
    phiStages.clear();
    UStages.clear();

    for ( int i = 0; i < k + 1; i++ )
    {
        pStages.push_back( volScalarField( p ) );
        phiStages.push_back( surfaceScalarField( phi ) );
        UStages.push_back( volVectorField( U ) );
    }
}

void SDCFluidSolver::nextTimeStep()
{
    timeIndex++;
    (*runTime)++;

    for ( int i = 0; i < k; i++ )
    {
        pStages.at( i ) = p;
        phiStages.at( i ) = phi;
        UStages.at( i ) = U;
    }
}

void SDCFluidSolver::solve()
{
    assert( false );
}

void SDCFluidSolver::finalizeTimeStep()
{
    assert( init );

    if ( runTime->timeIndex() % writeInterval == 0 )
        runTime->writeNow();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    init = false;
}

int SDCFluidSolver::getDOF()
{
    int index = 0;
    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    index++;
                }
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

void SDCFluidSolver::getSolution(
    fsi::vector & solution,
    fsi::vector & f
    )
{
    int index = 0;

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            solution( index ) = U.internalField()[i][j];
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    solution( index ) = U.boundaryField()[patchI][i][j];
                    index++;
                }
            }
        }
    }

    forAll( phi.internalField(), i )
    {
        solution( index ) = phi.internalField()[i];
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

    index = 0;

    forAll( UF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = UF.internalField()[i][j];
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( UF.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    f( index ) = UF.boundaryField()[patchI][i][j];
                    index++;
                }
            }
        }
    }

    forAll( phiF.internalField(), i )
    {
        f( index ) = phiF.internalField()[i];
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

void SDCFluidSolver::setSolution(
    const fsi::vector & solution,
    const fsi::vector & f
    )
{
    p = pStages.at( 0 );

    int index = 0;

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            U.internalField()[i][j] = solution( index );
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    U.boundaryField()[patchI][i][j] = solution( index );
                    index++;
                }
            }
        }
    }

    forAll( phi.internalField(), i )
    {
        phi.internalField()[i] = solution( index );
        index++;
    }

    forAll( phi.boundaryField(), patchI )
    {
        forAll( phi.boundaryField()[patchI], i )
        {
            phi.boundaryField()[patchI][i] = solution( index );
            index++;
        }
    }

    assert( index == solution.rows() );

    index = 0;

    forAll( UF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            UF.internalField()[i][j] = f( index );
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( UF.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    UF.boundaryField()[patchI][i][j] = f( index );
                    index++;
                }
            }
        }
    }

    forAll( phiF.internalField(), i )
    {
        phiF.internalField()[i] = f( index );
        index++;
    }

    forAll( phiF.boundaryField(), patchI )
    {
        forAll( phiF.boundaryField()[patchI], i )
        {
            phiF.boundaryField()[patchI][i] = f( index );
            index++;
        }
    }

    assert( index == f.rows() );
}

scalar SDCFluidSolver::getEndTime()
{
    return runTime->endTime().value();
}

scalar SDCFluidSolver::getStartTime()
{
    return runTime->startTime().value();
}

scalar SDCFluidSolver::getTimeStep()
{
    return runTime->deltaT().value();
}

void SDCFluidSolver::evaluateFunction(
    const int,
    const fsi::vector &,
    const scalar,
    fsi::vector & f
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
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = UF.internalField()[i][j];
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( UF.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    f( index ) = UF.boundaryField()[patchI][i][j];
                    index++;
                }
            }
        }
    }

    forAll( phiF.internalField(), i )
    {
        f( index ) = phiF.internalField()[i];
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
    const int,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    )
{
    runTime->setDeltaT( dt );
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
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            U.oldTime().internalField()[i][j] = qold( index );
            index++;
        }
    }

    forAll( U.oldTime().boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    U.oldTime().boundaryField()[patchI][i][j] = qold( index );
                    index++;
                }
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

    assert( index == qold.rows() );

    index = 0;

    forAll( rhsU.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            rhsU.internalField()[i][j] = rhs( index );
            index++;
        }
    }

    forAll( rhsU.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" )
        {
            forAll( rhsU.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    rhsU.boundaryField()[patchI][i][j] = rhs( index );
                    index++;
                }
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

    assert( index == rhs.rows() );

    // -------------------------------------------------------------------------

    courantNo();

    scalar convergenceTolerance = absoluteTolerance;

    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();

    // --- PIMPLE loop
    for ( label oCorr = 0; oCorr < maxIter; oCorr++ )
    {
        U.storePrevIter();

        fvVectorMatrix UEqn
        (
            fvm::div( phi, U )
            - fvm::laplacian( nu, U )
        );

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

            Foam::solve( UEqnt == -fvc::grad( p ) + rDeltaT * rhsU );

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

        // Relative convergence measure for the PISO loop:
        // Perform at maximum nCorr PISO corrections.
        // If the relative residual with respect to the initial
        // residual is decreased by factor tol: assume convergence.

        scalar initResidual = 1;
        scalar currResidual = 1;
        scalar pressureResidual = 1;
        scalar tol = 1.0e-2;

        // --- PISO loop
        for ( int corr = 0; corr < nCorr; corr++ )
        {
            p.storePrevIter();

            HU = UEqn.H();
            AU = UEqnt.A();

            U = (UEqnt.H() + rDeltaT * rhsU) / AU;

            phi = fvc::interpolate( HU ) / fvc::interpolate( AU ) & mesh.Sf();
            phi += rDeltaT * phi.oldTime() / fvc::interpolate( AU );
            phi += rDeltaT * rhsPhi / fvc::interpolate( AU );

            forAll( phi.boundaryField(), patchI )
            {
                if ( !phi.boundaryField()[patchI].coupled() )
                {
                    phi.boundaryField()[patchI] = U.boundaryField()[patchI] & mesh.Sf().boundaryField()[patchI];
                }
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

            if ( currResidual < std::max( tol * initResidual, scalar( 1.0e-15 ) ) )
                break;
        }

        scalar momentumResidual = evaluateMomentumResidual();

        if ( oCorr == 0 )
            convergenceTolerance = std::max( relativeTolerance * momentumResidual, absoluteTolerance );

        bool convergence = momentumResidual <= convergenceTolerance && oCorr >= minIter - 1;

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

    pStages.at( k + 1 ) = p;
    phiStages.at( k + 1 ) = phi;
    UStages.at( k + 1 ) = U;

    UF = rDeltaT * (U - U.oldTime() - rhsU);
    phiF = rDeltaT * (phi - phi.oldTime() - rhsPhi);

    getSolution( result, f );
}

void SDCFluidSolver::getVariablesInfo(
    std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    )
{
    dof.push_back( 0 );
    enabled.push_back( true );
    names.push_back( "U" );

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            dof.back() += 1;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "oscillatingCavityFixedValue" && U.boundaryField().types()[patchI] != "empty" )
        {
            dof.push_back( 0 );
            enabled.push_back( true );
            names.push_back( "patch " + mesh.boundaryMesh()[patchI].name() );

            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    dof.back() += 1;
                }
            }
        }
    }

    dof.push_back( 0 );
    enabled.push_back( true );
    names.push_back( "phi" );

    forAll( phi.internalField(), i )
    {
        dof.back() += 1;
    }

    forAll( phi.boundaryField(), patchI )
    {
        forAll( phi.boundaryField()[patchI], i )
        {
            dof.back() += 1;
        }
    }
}

void SDCFluidSolver::prepareImplicitSolve(
    bool,
    const int,
    const int,
    const scalar,
    const scalar,
    const fsi::vector &,
    const fsi::vector &
    )
{
    assert( false );
}
