
/*
 * Class
 *   SDCDynamicMeshFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCDynamicMeshFluidSolver.H"

SDCDynamicMeshFluidSolver::SDCDynamicMeshFluidSolver(
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
    rhsMeshPhi
    (
    IOobject
    (
        "rhsMeshPhi",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "0", dimVolume / dimTime, 0.0 )
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
    pointsStages(),
    volumeStages(),
    interpolateVolumeStages(),
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
    "UfF",
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
    UfF
    (
    phiFHeader,
    fvc::interpolate( UF )
    ),
    meshPhiF
    (
    IOobject
    (
        "rhsMeshPhi",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar( "0", dimVolume / dimTime, 0.0 )
    ),
    explicitFirstStage( false )
{
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

    const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
    const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );
    word ddtScheme;

    if ( ddtSchemes.found( "ddt(U)" ) )
        ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
    else
        ddtScheme = word( ddtSchemes.lookup( "default" ) );

    assert( ddtScheme == "bdf1" );

    initialize();
}

SDCDynamicMeshFluidSolver::~SDCDynamicMeshFluidSolver()
{}

void SDCDynamicMeshFluidSolver::continuityErrs()
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

void SDCDynamicMeshFluidSolver::createFields()
{
    surfaceVectorField nf = mesh.Sf() / mesh.magSf();
    surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
    surfaceVectorField Unor = phi / mesh.magSf() * nf;

    Uf = Utang + Unor;

    U.oldTime();
    Uf.oldTime();
    phi.oldTime();

    // Read pressure properties and create turbulence model
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

double SDCDynamicMeshFluidSolver::evaluateMomentumResidual()
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

void SDCDynamicMeshFluidSolver::getAcousticsDensityLocal( matrix & data )
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::getAcousticsVelocityLocal( matrix & data )
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::getAcousticsPressureLocal( matrix & data )
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::getTractionLocal( matrix & traction )
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::getWritePositionsLocalAcoustics( matrix & writePositions )
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::initialize()
{
    assert( !init );

    createFields();
}

void SDCDynamicMeshFluidSolver::initTimeStep()
{}

bool SDCDynamicMeshFluidSolver::isRunning()
{
    return runTime->run();
}

void SDCDynamicMeshFluidSolver::resetSolution()
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::setDeltaT( double dt )
{
    runTime->setDeltaT( dt );
}

void SDCDynamicMeshFluidSolver::setNumberOfStages( int k )
{
    this->k = k;

    volScalarField V
    (
        IOobject
        (
            "V",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );

    V.internalField() = mesh.V();
    V.correctBoundaryConditions();

    for ( int i = 0; i < k; i++ )
    {
        pStages.push_back( volScalarField( p ) );
        phiStages.push_back( surfaceScalarField( phi ) );
        UStages.push_back( volVectorField( U ) );
        UfStages.push_back( surfaceVectorField( Uf ) );
        pointsStages.push_back( mesh.points() );
        volumeStages.push_back( mesh.V() );
        interpolateVolumeStages.push_back( surfaceScalarField( fvc::interpolate( V ) ) );
    }
}

void SDCDynamicMeshFluidSolver::nextTimeStep()
{
    timeIndex++;
    (*runTime)++;

    if ( pStages.size() == static_cast<unsigned>(k) )
    {
        volScalarField V
        (
            IOobject
            (
                "V",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless,
            zeroGradientFvPatchScalarField::typeName
        );

        V.internalField() = mesh.V();
        V.correctBoundaryConditions();

        for ( int i = 0; i < k; i++ )
        {
            pStages.at( i ) = p;
            phiStages.at( i ) = phi;
            UStages.at( i ) = U;
            UfStages.at( i ) = Uf;
            pointsStages.at( i ) = mesh.points();
            volumeStages.at( i ) = mesh.V();
            interpolateVolumeStages.at( i ) = fvc::interpolate( V );
        }
    }
}

void SDCDynamicMeshFluidSolver::solve()
{
    assert( false );
}

void SDCDynamicMeshFluidSolver::finalizeTimeStep()
{
    runTime->writeNow();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;
}

int SDCDynamicMeshFluidSolver::getDOF()
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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    index++;
                }
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
            index++;
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
                index++;
        }
    }

    forAll( mesh.phi().internalField(), i )
    {
        index++;
    }

    forAll( mesh.phi().boundaryField(), patchI )
    {
        forAll( mesh.phi().boundaryField()[patchI], i )
        {
            index++;
        }
    }

    return index;
}

void SDCDynamicMeshFluidSolver::getSolution( Eigen::VectorXd & solution )
{
    volScalarField V
    (
        IOobject
        (
            "V",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );

    V.internalField() = mesh.V();
    V.correctBoundaryConditions();
    surfaceScalarField interpolateV = fvc::interpolate( V );

    int index = 0;

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            solution( index ) = U.internalField()[i][j] * V.internalField()[i];
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    solution( index ) = U.boundaryField()[patchI][i][j] * V.boundaryField()[patchI][i];
                    index++;
                }
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            solution( index ) = Uf.internalField()[i][j] * interpolateV.internalField()[i];
            index++;
        }
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                solution( index ) = Uf.boundaryField()[patchI][i][j] * interpolateV.boundaryField()[patchI][i];
                index++;
            }
        }
    }

    // Unknown how to compute the actual value of the swept volumes
    forAll( mesh.phi().internalField(), i )
    {
        solution( index ) = 0.0;
        index++;
    }

    forAll( mesh.phi().boundaryField(), patchI )
    {
        forAll( mesh.phi().boundaryField()[patchI], i )
        {
            solution( index ) = 0.0;
            index++;
        }
    }

    assert( index == solution.rows() );
}

void SDCDynamicMeshFluidSolver::setSolution(
    const Eigen::VectorXd & solution,
    const Eigen::VectorXd & f
    )
{
    p = pStages.at( 0 );

    int index = 0;

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            U.internalField()[i][j] = solution( index );
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                U.boundaryField()[patchI][i][j] = solution( index );
                index++;
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            Uf.internalField()[i][j] = solution( index );
            index++;
        }
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                Uf.boundaryField()[patchI][i][j] = solution( index );
                index++;
            }
        }
    }

    assert( index == solution.rows() );

    index = 0;

    forAll( UF.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            UF.internalField()[i][j] = f( index );
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
    {
        forAll( UF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                UF.boundaryField()[patchI][i][j] = f( index );
                index++;
            }
        }
    }

    forAll( UfF.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            UfF.internalField()[i][j] = f( index );
            index++;
        }
    }

    forAll( UfF.boundaryField(), patchI )
    {
        forAll( UfF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                UfF.boundaryField()[patchI][i][j] = f( index );
                index++;
            }
        }
    }

    forAll( meshPhiF.internalField(), i )
    {
        meshPhiF.internalField()[i] = f( index );
        index++;
    }

    forAll( meshPhiF.boundaryField(), patchI )
    {
        forAll( meshPhiF.boundaryField()[patchI], i )
        {
            meshPhiF.boundaryField()[patchI][i] = f( index );
            index++;
        }
    }

    assert( index == f.rows() );
}

double SDCDynamicMeshFluidSolver::getEndTime()
{
    return runTime->endTime().value();
}

double SDCDynamicMeshFluidSolver::getStartTime()
{
    return runTime->startTime().value();
}

double SDCDynamicMeshFluidSolver::getTimeStep()
{
    return runTime->deltaT().value();
}

void SDCDynamicMeshFluidSolver::evaluateFunction(
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    )
{
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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity" )
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
    }

    forAll( UfF.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            f( index ) = UfF.internalField()[i][j];
            index++;
        }
    }

    forAll( UfF.boundaryField(), patchI )
    {
        forAll( UfF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                f( index ) = UfF.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    forAll( meshPhiF.internalField(), i )
    {
        f( index ) = meshPhiF.internalField()[i];
        index++;
    }

    forAll( meshPhiF.boundaryField(), patchI )
    {
        forAll( meshPhiF.boundaryField()[patchI], i )
        {
            f( index ) = meshPhiF.boundaryField()[patchI][i];
            index++;
        }
    }

    assert( index == f.rows() );
}

void SDCDynamicMeshFluidSolver::implicitSolve(
    bool corrector,
    const int k,
    const int kold,
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
    runTime->setTime( t, runTime->timeIndex() );

    if ( corrector )
    {
        p = pStages.at( k + 1 );
        phi = phiStages.at( k + 1 );
        U = UStages.at( k + 1 );
        Uf = UfStages.at( k + 1 );
    }

    Uf.oldTime() = UfStages.at( kold );
    U.oldTime() = UStages.at( kold );
    phi.oldTime() = phiStages.at( kold );

    // Moving wall boundary condition

    forAll( U.boundaryField().types(), i )
    {
        assert( U.boundaryField().types()[i] != "movingWallVelocity" );

        if ( U.boundaryField().types()[i] == "SDCMovingWallVelocity" )
        {
            SDCMovingWallVelocityFvPatchVectorField & bc = static_cast<SDCMovingWallVelocityFvPatchVectorField &>(U.boundaryField()[i]);
            bc.setSDCInfo( corrector, k );
        }
    }

    int index = 0;

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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity" )
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
    }

    forAll( rhsUf.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            rhsUf.internalField()[i][j] = rhs( index );
            index++;
        }
    }

    forAll( rhsUf.boundaryField(), patchI )
    {
        forAll( rhsUf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
            {
                rhsUf.boundaryField()[patchI][i][j] = rhs( index );
                index++;
            }
        }
    }

    forAll( rhsMeshPhi.internalField(), i )
    {
        rhsMeshPhi.internalField()[i] = rhs( index );
        index++;
    }

    forAll( rhsMeshPhi.boundaryField(), patchI )
    {
        forAll( rhsMeshPhi.boundaryField()[patchI], i )
        {
            rhsMeshPhi.boundaryField()[patchI][i] = rhs( index );
            index++;
        }
    }

    assert( index == rhs.rows() );

    // Update mesh.phi()
    {
        // Reset the mesh point locations to the old stage
        pointField pointsOld = pointsStages.at( kold );
        tmp<scalarField> sweptVols = mesh.movePoints( pointsOld );

        mesh.setOldPoints( pointsOld );
        mesh.setV0() = volumeStages.at( kold );

        mesh.update();

        scalar rDeltaT = 1.0 / runTime->deltaT().value();
        mesh.setPhi() -= rDeltaT * rhsMeshPhi;
    }

    volScalarField V
    (
        IOobject
        (
            "V",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );

    V.internalField() = mesh.V();
    V.correctBoundaryConditions();
    volScalarField V0 = V;
    V0.internalField() = mesh.V0();
    V0.correctBoundaryConditions();

    rhsU /= V;
    rhsUf /= fvc::interpolate( V );

    // -------------------------------------------------------------------------

    scalar convergenceTolerance = absoluteTolerance;

    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();

    // --- PIMPLE loop
    for ( label oCorr = 0; oCorr < maxIter; oCorr++ )
    {
        // Make the fluxes relative to the mesh motion
        // The absolute fluxes are saved in memory between stages
        fvc::makeRelative( phi, U );

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

        Foam::solve( UEqnt == -fvc::grad( p ) + rDeltaT * rhsU );

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
            AU = UEqnt.A();

            U = (UEqnt.H() + rDeltaT * rhsU) / AU;
            Uf = (fvc::interpolate( HU ) + rDeltaT * Uf.oldTime() * interpolateVolumeStages.at( kold ) / fvc::interpolate( V ) + rDeltaT * rhsUf) / fvc::interpolate( AU );

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

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative( phi, U );

            U -= fvc::grad( p ) / AU;
            U.correctBoundaryConditions();

            if ( currResidual < std::max( tol * initResidual, 1.0e-15 ) )
                break;
        }

        // Update the face velocities
        fvc::makeAbsolute( phi, U );
        {
            surfaceVectorField nf = mesh.Sf() / mesh.magSf();
            surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
            surfaceVectorField Unor = phi / mesh.magSf() * nf;

            Uf = Utang + Unor;
        }
        fvc::makeRelative( phi, U );

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

        // Make the fluxes absolute to the mesh motion
        fvc::makeAbsolute( phi, U );

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
        UfStages.at( k + 1 ) = Uf;
        pointsStages.at( k + 1 ) = mesh.points();
        volumeStages.at( k + 1 ) = mesh.V();
        interpolateVolumeStages.at( k + 1 ) = fvc::interpolate( V );
    }

    UF = rDeltaT * (U * V - U.oldTime() * V0 - rhsU * V);
    UfF = rDeltaT * ( Uf * fvc::interpolate( V ) - Uf.oldTime() * interpolateVolumeStages.at( kold ) - rhsUf * fvc::interpolate( V ) );
    meshPhiF = mesh.phi();

    getSolution( result );
    evaluateFunction( k + 1, qold, t, f );
}

double SDCDynamicMeshFluidSolver::getScalingFactor()
{
    return 1;
}

void SDCDynamicMeshFluidSolver::getVariablesInfo(
    std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    )
{
    dof.push_back( 0 );
    dof.push_back( 0 );
    dof.push_back( 0 );
    enabled.push_back( true );
    enabled.push_back( true );
    enabled.push_back( false );
    names.push_back( "U" );
    names.push_back( "Uf" );
    names.push_back( "meshPhi" );

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
        {
            dof.at( 0 ) += 1;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    dof.at( 0 ) += 1;
                }
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < 3; j++ )
            dof.at( 1 ) += 1;
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < 3; j++ )
                dof.at( 1 ) += 1;
        }
    }

    forAll( mesh.phi().internalField(), i )
    {
        dof.at( 2 ) += 1;
    }

    forAll( mesh.phi().boundaryField(), patchI )
    {
        forAll( mesh.phi().boundaryField()[patchI], i )
        {
            dof.at( 2 ) += 1;
        }
    }
}

bool SDCDynamicMeshFluidSolver::isConverged()
{
    bool convergence = true;

    forAll( U.boundaryField().types(), i )
    {
        assert( U.boundaryField().types()[i] != "movingWallVelocity" );

        if ( U.boundaryField().types()[i] == "SDCMovingWallVelocity" )
        {
            SDCMovingWallVelocityFvPatchVectorField & bc = static_cast<SDCMovingWallVelocityFvPatchVectorField &>(U.boundaryField()[i]);
            bc.sdc->outputResidual();

            if ( not bc.sdc->isConverged() )
                convergence = false;
        }
    }

    return convergence;
}
