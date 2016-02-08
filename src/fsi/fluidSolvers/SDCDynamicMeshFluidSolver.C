
/*
 * Class
 *   SDCDynamicMeshFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCDynamicMeshFluidSolver.H"

SDCDynamicMeshFluidSolver::SDCDynamicMeshFluidSolver(
    std::string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    SDCFluidSolver( name, args, runTime ),
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
    fvc::interpolate( U )
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
    UfStages(),
    pointsStages(),
    volumeStages(),
    interpolateVolumeStages(),
    UfFHeader
    (
    "UfF",
    runTime->timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
    ),
    UfF
    (
    UfFHeader,
    fvc::interpolate( UF )
    ),
    meshPhiF
    (
    IOobject
    (
        "meshPhiF",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar( "0", dimVolume / dimTime, 0.0 )
    ),
    kold( 0 ),
    indexk( 0 )
{
    initialize();
}

SDCDynamicMeshFluidSolver::~SDCDynamicMeshFluidSolver()
{}

void SDCDynamicMeshFluidSolver::createFields()
{
    surfaceVectorField nf = mesh.Sf() / mesh.magSf();
    surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
    surfaceVectorField Unor = phi / mesh.magSf() * nf;

    Uf = Utang + Unor;

    U.oldTime();
    Uf.oldTime();

    // Read pressure properties and create turbulence model
    pRefCell = 0;
    pRefValue = 0.0;
    setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );
}

void SDCDynamicMeshFluidSolver::initialize()
{
    assert( !init );

    createFields();
}

void SDCDynamicMeshFluidSolver::resetSolution()
{}

void SDCDynamicMeshFluidSolver::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

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

    for ( int i = 0; i < k + 1; i++ )
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

    assert( timeIndex == runTime->timeIndex() );

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

int SDCDynamicMeshFluidSolver::getDOF()
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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
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

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
            index++;
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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

void SDCDynamicMeshFluidSolver::getSolution(
    fsi::vector & solution,
    fsi::vector & f
    )
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
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            solution( index ) = U.internalField()[i][j] * V.internalField()[i];
            index++;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    solution( index ) = U.boundaryField()[patchI][i][j] * V.boundaryField()[patchI][i];
                    index++;
                }
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            solution( index ) = Uf.internalField()[i][j] * interpolateV.internalField()[i];
            index++;
        }
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
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

    forAll( UfF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = UfF.internalField()[i][j];
            index++;
        }
    }

    forAll( UfF.boundaryField(), patchI )
    {
        forAll( UfF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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

void SDCDynamicMeshFluidSolver::setSolution(
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
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                U.boundaryField()[patchI][i][j] = solution( index );
                index++;
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            Uf.internalField()[i][j] = solution( index );
            index++;
        }
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            UF.internalField()[i][j] = f( index );
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
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

    forAll( UfF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            UfF.internalField()[i][j] = f( index );
            index++;
        }
    }

    forAll( UfF.boundaryField(), patchI )
    {
        forAll( UfF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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

void SDCDynamicMeshFluidSolver::evaluateFunction(
    const int,
    const fsi::vector &,
    const scalar,
    fsi::vector & f
    )
{
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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
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

    forAll( UfF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = UfF.internalField()[i][j];
            index++;
        }
    }

    forAll( UfF.boundaryField(), patchI )
    {
        forAll( UfF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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

void SDCDynamicMeshFluidSolver::prepareImplicitSolve(
    bool corrector,
    const int k,
    const int kold,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs
    )
{
    assert( k < this->k );
    assert( kold < this->k );
    assert( kold <= k );
    assert( qold.rows() == rhs.rows() );
    this->kold = kold;
    this->indexk = k;

    runTime->setDeltaT( dt );
    runTime->setTime( t, runTime->timeIndex() );
    this->t = runTime->time().value();

    if ( corrector )
    {
        p = pStages.at( k + 1 );
        phi = phiStages.at( k + 1 );
        U = UStages.at( k + 1 );
        Uf = UfStages.at( k + 1 );
    }

    Uf.oldTime() = UfStages.at( kold );
    U.oldTime() = UStages.at( kold );

    // Moving wall boundary condition

    forAll( U.boundaryField().types(), i )
    {
        assert( U.boundaryField().types()[i] != "movingWallVelocity" );
        assert( U.boundaryField().types()[i] != "myMovingWallVelocity" );

        if ( U.boundaryField().types()[i] == "SDCMovingWallVelocity" )
        {
            SDCMovingWallVelocityFvPatchVectorField & bc = static_cast<SDCMovingWallVelocityFvPatchVectorField &>(U.boundaryField()[i]);
            bc.setSDCInfo( corrector, k );
        }
    }

    int index = 0;

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
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
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

    forAll( rhsUf.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            rhsUf.internalField()[i][j] = rhs( index );
            index++;
        }
    }

    forAll( rhsUf.boundaryField(), patchI )
    {
        forAll( rhsUf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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
        mesh.setV0() = volumeStages.at( kold );

        RBFMeshMotionSolver & motionSolver =
            const_cast<RBFMeshMotionSolver &>
            (
            mesh.lookupObject<RBFMeshMotionSolver>( "dynamicMeshDict" )
            );
        motionSolver.corrector = corrector;
        motionSolver.k = k;
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
}

void SDCDynamicMeshFluidSolver::implicitSolve(
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
    prepareImplicitSolve( corrector, k, kold, t, dt, qold, rhs );

    solve();

    getSolution( result, f );
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
    names.push_back( "fluid U" );
    names.push_back( "fluid Uf" );
    names.push_back( "fluid meshPhi" );

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            dof.at( 0 ) += 1;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        if ( U.boundaryField().types()[patchI] != "SDCMovingWallVelocity"
            && U.boundaryField().types()[patchI] != "transitionalParabolicVelocity" )
        {
            forAll( U.boundaryField()[patchI], i )
            {
                for ( int j = 0; j < mesh.nGeometricD(); j++ )
                {
                    dof.at( 0 ) += 1;
                }
            }
        }
    }

    forAll( Uf.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
            dof.at( 1 ) += 1;
    }

    forAll( Uf.boundaryField(), patchI )
    {
        forAll( Uf.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
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
        assert( U.boundaryField().types()[i] != "myMovingWallVelocity" );

        if ( U.boundaryField().types()[i] == "SDCMovingWallVelocity" )
        {
            SDCMovingWallVelocityFvPatchVectorField & bc = static_cast<SDCMovingWallVelocityFvPatchVectorField &>(U.boundaryField()[i]);
            bc.timeIntegrationScheme->outputResidual( "moving wall velocity" );

            if ( not bc.timeIntegrationScheme->isConverged() )
                convergence = false;
        }
    }

    RBFMeshMotionSolver & motionSolver =
        const_cast<RBFMeshMotionSolver &>
        (
        mesh.lookupObject<RBFMeshMotionSolver>( "dynamicMeshDict" )
        );

    if ( motionSolver.timeIntegrationScheme )
    {
        motionSolver.timeIntegrationScheme->outputResidual( "mesh motion" );

        if ( not motionSolver.timeIntegrationScheme->isConverged() )
            convergence = false;
    }

    return convergence;
}

void SDCDynamicMeshFluidSolver::solve()
{
    Info << "Solve fluid domain" << endl;

    // Update mesh.phi()
    {
        mesh.update();

        scalar rDeltaT = 1.0 / runTime->deltaT().value();

        // Create swept volumes
        const faceList & f = mesh.faces();

        scalarField sweptVols( f.size() );

        forAll( f, faceI )
        {
            sweptVols[faceI] = f[faceI].sweptVol( pointsStages.at( kold ), mesh.points() );
        }

        mesh.setPhi().internalField() = scalarField::subField( sweptVols, mesh.nInternalFaces() );
        mesh.setPhi().internalField() *= rDeltaT;

        const fvPatchList & patches = mesh.boundary();

        forAll( patches, patchI )
        {
            mesh.setPhi().boundaryField()[patchI] = patches[patchI].patchSlice( sweptVols );
            mesh.setPhi().boundaryField()[patchI] *= rDeltaT;
        }

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

    // -------------------------------------------------------------------------

    courantNo();

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
                    Uf -= pEqn.flux() * mesh.Sf() / ( mesh.magSf() * mesh.magSf() );
                    phi = Uf & mesh.Sf();
                }
            }

            p.relax();

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative( phi, U );

            U -= fvc::grad( p ) / AU;
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

        // Make the fluxes absolute to the mesh motion
        fvc::makeAbsolute( phi, U );

        if ( convergence )
            break;
    }

    continuityErrs();

    // -------------------------------------------------------------------------

    pStages.at( indexk + 1 ) = p;
    phiStages.at( indexk + 1 ) = phi;
    UStages.at( indexk + 1 ) = U;
    UfStages.at( indexk + 1 ) = Uf;
    pointsStages.at( indexk + 1 ) = mesh.points();
    volumeStages.at( indexk + 1 ) = mesh.V();
    interpolateVolumeStages.at( indexk + 1 ) = fvc::interpolate( V );

    UF = rDeltaT * (U * V - U.oldTime() * V0 - rhsU * V);
    UfF = rDeltaT * ( Uf * fvc::interpolate( V ) - Uf.oldTime() * interpolateVolumeStages.at( kold ) - rhsUf * fvc::interpolate( V ) );
    meshPhiF = mesh.phi();
}

void SDCDynamicMeshFluidSolver::getTractionLocal( matrix & traction )
{
    vectorField tractionField( getInterfaceSizeLocal(), Foam::vector::zero );

    int offset = 0;

    forAll( movingPatchIDs, patchI )
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        vectorField tractionFieldPatchI = -rho.value() * nu.value()
            * U.boundaryField()[movingPatchIDs[patchI]].snGrad()
            + rho.value() * p.boundaryField()[movingPatchIDs[patchI]]
            * mesh.boundary()[movingPatchIDs[patchI]].nf();

        forAll( tractionFieldPatchI, i )
        {
            tractionField[i + offset] = tractionFieldPatchI[i];
        }

        offset += size;
    }

    assert( tractionField.size() == nGlobalCenters[Pstream::myProcNo()] );

    traction.resize( tractionField.size(), mesh.nGeometricD() );

    for ( int i = 0; i < traction.rows(); i++ )
        for ( int j = 0; j < traction.cols(); j++ )
            traction( i, j ) = tractionField[i][j];
}
