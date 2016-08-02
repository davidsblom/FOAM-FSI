
/*
 * Class
 *   FluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "FluidSolver.H"
#include <stdexcept>

FluidSolver::FluidSolver(
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
    UfHeader
    (
        "Uf",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Uf
    (
        UfHeader,
        mesh,
        dimensionedVector( "0", U.dimensions(), Foam::vector::zero )
    ),
    AU
    (
        IOobject
        (
            "AU",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        U.dimensions() / runTime->deltaT().dimensions(),
        zeroGradientFvPatchVectorField::typeName
    ),
    nCorr( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "nCorrectors" ) ) ),
    nNonOrthCorr( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "nNonOrthogonalCorrectors" ) ) ),
    minIter( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "minIter" ) ) ),
    maxIter( readInt( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "maxIter" ) ) ),
    absoluteTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "tolerance" ) ) ),
    relativeTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "relTol" ) ) ),
    pisoTolerance( readScalar( mesh.solutionDict().subDict( "PIMPLE" ).lookup( "pisoTol" ) ) ),
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
    CoNum( 0 ),
    meanCoNum( 0 ),
    velMag( 0 ),
    turbulenceSwitch( true )
{
    assert( absoluteTolerance < 1 );
    assert( absoluteTolerance > 0 );
    assert( nCorr > 0 );
    assert( maxIter >= 1 );
    assert( nNonOrthCorr >= 0 );
    assert( relativeTolerance < 1 );
    assert( minIter <= maxIter );
    assert( minIter >= 0 );
    assert( pisoTolerance < 1 );
    assert( pisoTolerance > 0 );

    // Ensure that the absolute tolerance of the linear solver is less than the
    // used convergence tolerance for the non-linear system.
    scalar absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    if ( absTolerance >= absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver of the momentum equation (U) should be smaller than PIMPLE::absoluteTolerance in order to reach convergence of the non-linear system" );

    absTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "p" ).lookup( "tolerance" ) );
    assert( absTolerance < absoluteTolerance );

    if ( absTolerance >= absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver of the Poisson equation (p) should be smaller than PIMPLE::absoluteTolerance in order to reach convergence of the non-linear system" );

    checkTimeDiscretisationScheme();

    U.oldTime();
    Uf.oldTime();
    phi.oldTime();

    if ( !UfHeader.headerOk() )
    {
        tmp<surfaceVectorField> nf( mesh.Sf() / mesh.magSf() );
        tmp<surfaceVectorField> Unor( phi / mesh.magSf() * nf() );
        tmp<surfaceVectorField> Utang( fvc::interpolate( U ) - nf * ( fvc::interpolate( U ) & nf() ) );

        Uf = Utang() + Unor();
    }

    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );

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
}

FluidSolver::~FluidSolver()
{}

void FluidSolver::checkTimeDiscretisationScheme()
{
    const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
    const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );
    word ddtScheme;

    if ( ddtSchemes.found( "ddt(U)" ) )
        ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
    else
        ddtScheme = word( ddtSchemes.lookup( "default" ) );

    if ( ddtScheme != word( "bdf1" ) && ddtScheme != word( "bdf2" ) && ddtScheme != word( "bdf3" ) )
    {
        FatalErrorIn( "checkTimeDiscretisationScheme: " )
            << " ddt(U) scheme is currently set to " << ddtScheme
            << ". This solver only works with ddt(U) scheme = bdf1, bdf2 or bdf3."
            << abort( FatalError );
    }

    forAll( U.boundaryField().types(), i )
    {
        assert( U.boundaryField().types()[i] != "movingWallVelocity" );
        assert( U.boundaryField().types()[i] != "SDCMovingWallVelocity" );
    }
}

void FluidSolver::continuityErrs()
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

void FluidSolver::courantNo()
{
    if ( mesh.nInternalFaces() )
    {
        tmp<surfaceScalarField> magPhi = mag( phi );

        tmp<surfaceScalarField> SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs() * magPhi();

        CoNum = max( SfUfbyDelta() / mesh.magSf() )
            .value() * runTime->deltaT().value();

        meanCoNum = ( sum( SfUfbyDelta() ) / sum( mesh.magSf() ) )
            .value() * runTime->deltaT().value();

        velMag = max( magPhi() / mesh.magSf() ).value();
    }

    Info << "Courant Number mean: " << meanCoNum
         << " max: " << CoNum
         << " velocity magnitude: " << velMag
         << endl;
}

scalar FluidSolver::evaluateMomentumResidual()
{
    tmp<volVectorField> residual = fvc::ddt( U ) + fvc::div( phi, U ) + fvc::grad( p );

    if ( turbulenceSwitch )
        residual() += turbulence->divDevReff( U ) & U;
    else
        residual() += -fvc::laplacian( nu, U );

    tmp<scalarField> magResU = mag( residual->internalField() );
    return std::sqrt( gSumSqr( magResU ) / ( mesh.globalData().nTotalCells() * mesh.nGeometricD() ) );
}

void FluidSolver::getAcousticsDensityLocal( matrix & )
{
    assert( false );
}

void FluidSolver::getAcousticsVelocityLocal( matrix & )
{
    assert( false );
}

void FluidSolver::getAcousticsPressureLocal( matrix & )
{
    assert( false );
}

void FluidSolver::getTractionLocal( matrix & traction )
{
    int size = 0;

    forAll( movingPatchIDs, patchI )
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();
    }

    vectorField tractionField( size, Foam::vector::zero );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        tmp<vectorField> tractionFieldPatchI = -rho.value() * nu.value()
            * U.boundaryField()[movingPatchIDs[patchI]].snGrad()
            + rho.value() * p.boundaryField()[movingPatchIDs[patchI]]
            * mesh.boundary()[movingPatchIDs[patchI]].nf();

        forAll( tractionFieldPatchI(), i )
        {
            tractionField[i + offset] = tractionFieldPatchI()[i];
        }

        offset += size;
    }

    traction.resize( tractionField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < traction.rows(); i++ )
        for ( int j = 0; j < traction.cols(); j++ )
            traction( i, j ) = tractionField[i][j];
}

void FluidSolver::getWritePositionsLocalAcoustics( matrix & )
{
    assert( false );
}

void FluidSolver::initTimeStep()
{
    assert( !init );

    timeIndex++;
    t = runTime->time().value();

    Info << "\nTime = " << runTime->value() << endl;

    courantNo();

    init = true;
}

void FluidSolver::finalizeTimeStep()
{
    foamFluidSolver::finalizeTimeStep();
    runTime->write();
}

bool FluidSolver::isRunning()
{
    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void FluidSolver::resetSolution()
{}

void FluidSolver::solve()
{
    Info << "Solve fluid domain" << endl;

    mesh.update();

    scalar convergenceTolerance = absoluteTolerance;

    // --- PIMPLE loop
    for ( label oCorr = 0; oCorr < maxIter; oCorr++ )
    {
        // Make the fluxes relative to the mesh motion
        fvc::makeRelative( phi, U );

        U.storePrevIter();

        fvVectorMatrix UEqn
        (
            fvm::ddt( U )
            + fvm::div( phi, U )
        );

        if ( turbulenceSwitch )
            UEqn += turbulence->divDevReff( U );
        else
            UEqn += -fvm::laplacian( nu, U );

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

            Foam::solve( UEqn == -fvc::grad( p ) );

            // Reset equation to ensure relaxation parameter is not causing problems for time order
            UEqn =
                (
                fvm::ddt( U )
                + fvm::div( phi, U )
                );

            if ( turbulenceSwitch )
                UEqn += turbulence->divDevReff( U );
            else
                UEqn += -fvm::laplacian( nu, U );

            UEqn.source() = S0;
            UEqn.boundaryCoeffs() = B0;
        }

        // Relative convergence measure for the PISO loop:
        // Perform at maximum nCorr PISO corrections.
        // If the relative residual with respect to the initial
        // residual is decreased by factor tol: assume convergence.

        scalar initResidual = 1;
        scalar currResidual = 1;
        scalar pressureResidual = 1;

        // --- PISO loop
        for ( int corr = 0; corr < nCorr; corr++ )
        {
            p.storePrevIter();

            HU = UEqn.H();
            AU = UEqn.A();
            U = HU / AU;

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

                phi += fvc::ddtPhiCorr( 1.0 / AU, U, phi );
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

            U -= (1.0 / AU) * fvc::grad( p );
            U.correctBoundaryConditions();

            if ( currResidual < std::max( pisoTolerance * initResidual, scalar( 1.0e-15 ) ) )
                break;
        }

        if ( turbulenceSwitch )
            turbulence->correct();

        // Update the face velocities
        fvc::makeAbsolute( phi, U );
        {
            tmp<surfaceVectorField> nf( mesh.Sf() / mesh.magSf() );
            tmp<surfaceVectorField> Unor( phi / mesh.magSf() * nf() );
            tmp<surfaceVectorField> Utang( fvc::interpolate( U ) - nf * ( fvc::interpolate( U ) & nf() ) );

            Uf = Utang() + Unor();
        }
        fvc::makeRelative( phi, U );

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

        // Make the fluxes absolute to the mesh motion
        fvc::makeAbsolute( phi, U );

        if ( convergence )
            break;
    }

    continuityErrs();
}
