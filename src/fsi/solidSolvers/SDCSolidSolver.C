
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDCSolidSolver.H"

SDCSolidSolver::SDCSolidSolver (
    std::string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    foamSolidSolver( name, args, runTime ),
    gradU( fvc::grad( U ) ),
    epsilon
    (
        IOobject
        (
            "epsilonGreen",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor( "zero", dimless, symmTensor::zero )
    ),
    sigma
    (
        IOobject
        (
            "sigma",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor( "zero", dimForce / dimArea, symmTensor::zero )
    ),
    rheology( sigma, U ),
    rho( rheology.rho() ),
    mu( rheology.mu() ),
    lambda( rheology.lambda() ),
    muf( fvc::interpolate( mu, "mu" ) ),
    lambdaf( fvc::interpolate( lambda, "lambda" ) ),
    n( mesh.Sf() / mesh.magSf() ),
    V
    (
        IOobject
        (
            "V",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt( U )
    ),
    minIter( readInt( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "minIter" ) ) ),
    maxIter( readInt( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "maxIter" ) ) ),
    absoluteTolerance( readScalar( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "tolerance" ) ) ),
    relativeTolerance( readScalar( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "relTol" ) ) ),
    interpolator( nullptr ),
    k( 0 ),
    indexk( 0 ),
    UStages(),
    VStages(),
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
        dimensionedVector( "rhsU", dimLength, Foam::vector::zero )
    ),
    rhsV
    (
        IOobject
        (
            "rhsV",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector( "rhsV", dimVelocity, Foam::vector::zero )
    ),
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
        dimensionedVector( "UF", dimLength / dimTime, Foam::vector::zero )
    ),
    VF
    (
        IOobject
        (
            "VF",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector( "VF", dimVelocity / dimTime, Foam::vector::zero )
    )
{
    const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
    const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );
    word ddtScheme;

    if ( ddtSchemes.found( "ddt(U)" ) )
        ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
    else
        ddtScheme = word( ddtSchemes.lookup( "default" ) );

    assert( ddtScheme == "bdf1" );
}

SDCSolidSolver::SDCSolidSolver (
    std::string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime,
    std::shared_ptr<rbf::RBFCoarsening> interpolator
    )
    :
    foamSolidSolver( name, args, runTime ),
    gradU( fvc::grad( U ) ),
    epsilon
    (
        IOobject
        (
            "epsilonGreen",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor( "zero", dimless, symmTensor::zero )
    ),
    sigma
    (
        IOobject
        (
            "sigma",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor( "zero", dimForce / dimArea, symmTensor::zero )
    ),
    rheology( sigma, U ),
    rho( rheology.rho() ),
    mu( rheology.mu() ),
    lambda( rheology.lambda() ),
    muf( fvc::interpolate( mu, "mu" ) ),
    lambdaf( fvc::interpolate( lambda, "lambda" ) ),
    n( mesh.Sf() / mesh.magSf() ),
    V
    (
        IOobject
        (
            "V",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt( U )
    ),
    minIter( readInt( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "minIter" ) ) ),
    maxIter( readInt( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "maxIter" ) ) ),
    absoluteTolerance( readScalar( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "tolerance" ) ) ),
    relativeTolerance( readScalar( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "relTol" ) ) ),
    interpolator( interpolator ),
    k( 0 ),
    indexk( 0 ),
    UStages(),
    VStages(),
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
        dimensionedVector( "rhsU", dimLength, Foam::vector::zero )
    ),
    rhsV
    (
        IOobject
        (
            "rhsV",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector( "rhsV", dimVelocity, Foam::vector::zero )
    ),
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
        dimensionedVector( "UF", dimLength / dimTime, Foam::vector::zero )
    ),
    VF
    (
        IOobject
        (
            "VF",
            runTime->timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector( "VF", dimVelocity / dimTime, Foam::vector::zero )
    )
{
    const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
    const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );
    word ddtScheme;

    if ( ddtSchemes.found( "ddt(U)" ) )
        ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
    else
        ddtScheme = word( ddtSchemes.lookup( "default" ) );

    assert( ddtScheme == "bdf1" );
}

SDCSolidSolver::~SDCSolidSolver()
{}

void SDCSolidSolver::calculateEpsilonSigma()
{
    // - Green finite strain tensor
    epsilon = symm( gradU ) + 0.5 * symm( gradU & gradU.T() );

    // - second Piola-Kirchhoff stress tensor
    sigma = 2 * mu * epsilon + lambda * ( I * tr( epsilon ) );
}

void SDCSolidSolver::initialize()
{
    readCouplingProperties();
}

void SDCSolidSolver::initTimeStep()
{
    assert( !init );

    readSolidMechanicsControls();

    init = true;
}

bool SDCSolidSolver::interpolateVolField( std::shared_ptr<BaseMultiLevelSolver> solver )
{
    std::shared_ptr<SDCSolidSolver> fineModel;
    fineModel = std::dynamic_pointer_cast<SDCSolidSolver>( solver );

    assert( fineModel );

    // Interpolate the displacement of the fineMesh onto the own
    // mesh

    Info << "Mesh to mesh volume interpolation of field U for the solid domain" << endl;

    const volVectorField & fieldSource = fineModel->mesh.lookupObject<volVectorField>( "U" );

    // Gather all the cell centers of the source

    labelList cellCentresSourceSize( Pstream::nProcs(), 0 );

    const Foam::vectorField & cellCentresSource = fineModel->mesh.cellCentres();

    cellCentresSourceSize[Pstream::myProcNo()] = cellCentresSource.size();

    // Reduce (gather, scatter) over procs

    reduce( cellCentresSourceSize, sumOp<labelList>() );

    assert( interpolator );

    if ( !interpolator->rbf->computed )
    {
        // Gather all the cell centers of the source

        const Foam::vectorField & cellCentresTarget = mesh.cellCentres();

        vectorField globalControlPoints( sum( cellCentresSourceSize ), Foam::vector::zero );

        label startIndex = 0;

        for ( int i = 0; i < Pstream::myProcNo(); i++ )
            startIndex += cellCentresSourceSize[i];

        for ( int i = 0; i < cellCentresSourceSize[Pstream::myProcNo()]; i++ )
            globalControlPoints[startIndex + i] = cellCentresSource[i];

        reduce( globalControlPoints, sumOp<vectorField>() );

        matrix positions( sum( cellCentresSourceSize ), 3 );
        matrix positionsInterpolation( cellCentresTarget.size(), 3 );

        for ( int i = 0; i < positions.rows(); i++ )
            for ( int j = 0; j < positions.cols(); j++ )
                positions( i, j ) = globalControlPoints[i][j];

        for ( int i = 0; i < positionsInterpolation.rows(); i++ )
            for ( int j = 0; j < positionsInterpolation.cols(); j++ )
                positionsInterpolation( i, j ) = cellCentresTarget[i][j];

        interpolator->compute( positions, positionsInterpolation );
    }

    vectorField globalFieldSource( fineModel->mesh.globalData().nTotalCells(), Foam::vector::zero );

    label startIndex = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        startIndex += cellCentresSourceSize[i];

    for ( int i = 0; i < cellCentresSourceSize[Pstream::myProcNo()]; i++ )
        globalFieldSource[startIndex + i] = fieldSource[i];

    reduce( globalFieldSource, sumOp<vectorField>() );

    // Initialize variables for interpolation

    matrix values( globalFieldSource.size(), 3 );
    matrix valuesInterpolation;

    for ( int i = 0; i < values.rows(); i++ )
        for ( int j = 0; j < values.cols(); j++ )
            values( i, j ) = globalFieldSource[i][j];

    interpolator->interpolate( values, valuesInterpolation );

    for ( int i = 0; i < valuesInterpolation.rows(); i++ )
        for ( int j = 0; j < valuesInterpolation.cols(); j++ )
            U[i][j] = valuesInterpolation( i, j );

    gradU = fvc::grad( U );

    calculateEpsilonSigma();

    return true;
}

bool SDCSolidSolver::isRunning()
{
    return runTime->run();
}

void SDCSolidSolver::readSolidMechanicsControls()
{
    // Ensure that the absolute tolerance of the linear solver is less
    // than the used convergence tolerance for the non-linear system.
    scalar linearTolerance = readScalar( mesh.solutionDict().subDict( "solvers" ).subDict( "U" ).lookup( "tolerance" ) );
    assert( linearTolerance < absoluteTolerance );
    assert( relativeTolerance < 1 );
    assert( absoluteTolerance > 0 );
    assert( absoluteTolerance < 1 );
    assert( minIter < maxIter );
    assert( maxIter > 0 );
    assert( minIter >= 0 );

    if ( not linearTolerance < absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver (U) should be smaller than solidMechanics::absoluteTolerance in order to reach convergence of the non-linear system" );
}

void SDCSolidSolver::resetSolution()
{
    U == U.oldTime();
    gradU = fvc::grad( U );
}

void SDCSolidSolver::solve()
{
    Info << "Solve solid domain" << endl;

    scalar iCorr = 0;
    scalar displacementResidual = 1;
    scalar velocityResidual = 1;
    scalar residual = 1;
    scalar initialResidual = 1;
    lduMatrix::solverPerformance solverPerf;
    lduMatrix::debug = 0;
    scalar convergenceTolerance = absoluteTolerance;

    gradU = fvc::grad( U );

    calculateEpsilonSigma();

    dimensionedVector gravity( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "gravity" ) );

    dimensionedScalar deltaT = runTime->deltaT();

    for ( iCorr = 0; iCorr < maxIter; iCorr++ )
    {
        U.storePrevIter();
        V.storePrevIter();

        surfaceTensorField shearGradU = ( (I - n * n) & fvc::interpolate( gradU ) );

        fvVectorMatrix UEqn
        (
            fvm::ddt( U )
            ==
            deltaT / rho * (
                fvm::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
                + fvc::div(
                    mesh.magSf()
                    * (
                        -(muf + lambdaf) * ( fvc::snGrad( U ) & (I - n * n) )
                        + lambdaf * tr( shearGradU & (I - n * n) ) * n
                        + muf * (shearGradU & n)
                        + muf * ( n & fvc::interpolate( gradU & gradU.T() ) )
                        + 0.5 * lambdaf
                        * ( n * tr( fvc::interpolate( gradU & gradU.T() ) ) )
                        + ( n & fvc::interpolate( sigma & gradU ) )
                        )
                    )
                )
            + V.oldTime()
            + rhsV
            + rhsU / deltaT
        );

        // Add gravity

        UEqn -= deltaT * gravity;

        solverPerf = UEqn.solve();

        gradU = fvc::grad( U );
        shearGradU = ( (I - n * n) & fvc::interpolate( gradU ) );
        calculateEpsilonSigma();

        V = deltaT / rho * (
            fvc::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
            + fvc::div(
                mesh.magSf()
                * (
                    -(muf + lambdaf) * ( fvc::snGrad( U ) & (I - n * n) )
                    + lambdaf * tr( shearGradU & (I - n * n) ) * n
                    + muf * (shearGradU & n)
                    + muf * ( n & fvc::interpolate( gradU & gradU.T() ) )
                    + 0.5 * lambdaf
                    * ( n * tr( fvc::interpolate( gradU & gradU.T() ) ) )
                    + ( n & fvc::interpolate( sigma & gradU ) )
                    )
                )
            )
            + V.oldTime()
            + rhsV
            + deltaT * gravity
        ;

        displacementResidual = gSumMag( U.internalField() - U.prevIter().internalField() ) / (gSumMag( U.internalField() - U.oldTime().internalField() ) + SMALL);
        velocityResidual = gSumMag( V.internalField() - V.prevIter().internalField() ) / (gSumMag( V.internalField() - V.oldTime().internalField() ) + SMALL);
        residual = max( displacementResidual, solverPerf.initialResidual() );
        residual = max( residual, velocityResidual );

        if ( iCorr == 0 )
        {
            initialResidual = residual;
            convergenceTolerance = std::max( relativeTolerance * residual, absoluteTolerance );
            assert( convergenceTolerance > 0 );
            assert( convergenceTolerance < 1 );
        }

        bool convergence = residual <= convergenceTolerance && iCorr >= minIter - 1;

        if ( convergence )
            break;
    }

    lduMatrix::debug = 1;

    Info << "Solving for " << U.name();
    Info << ", Initial residual = " << initialResidual;
    Info << ", Final residual velocity = " << velocityResidual;
    Info << ", Final residual displacement = " << displacementResidual;
    Info << ", No outer iterations " << iCorr << endl;

    // -------------------------------------------------------------------------

    UStages.at( indexk + 1 ) = U;
    VStages.at( indexk + 1 ) = V;

    dimensionedScalar rDeltaT = 1.0 / mesh.time().deltaT();
    UF = rDeltaT * (U - U.oldTime() - rhsU);
    VF = rDeltaT * (V - V.oldTime() - rhsV);
}

void SDCSolidSolver::evaluateFunction(
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
        forAll( UF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                f( index ) = UF.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    forAll( VF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = VF.internalField()[i][j];
            index++;
        }
    }

    forAll( VF.boundaryField(), patchI )
    {
        forAll( VF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                f( index ) = VF.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    assert( index == f.rows() );
}

void SDCSolidSolver::finalizeTimeStep()
{
    foamSolidSolver::finalizeTimeStep();

    runTime->writeNow();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;
}

int SDCSolidSolver::getDOF()
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
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                index++;
            }
        }
    }

    forAll( V.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            index++;
        }
    }

    forAll( V.boundaryField(), patchI )
    {
        forAll( V.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                index++;
            }
        }
    }

    return index;
}

void SDCSolidSolver::getSolution(
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
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                solution( index ) = U.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    forAll( V.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            solution( index ) = V.internalField()[i][j];
            index++;
        }
    }

    forAll( V.boundaryField(), patchI )
    {
        forAll( V.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                solution( index ) = V.boundaryField()[patchI][i][j];
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
            f( index ) = UF.internalField()[i][j];
            index++;
        }
    }

    forAll( UF.boundaryField(), patchI )
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

    forAll( VF.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            f( index ) = VF.internalField()[i][j];
            index++;
        }
    }

    forAll( VF.boundaryField(), patchI )
    {
        forAll( VF.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                f( index ) = VF.boundaryField()[patchI][i][j];
                index++;
            }
        }
    }

    assert( index == f.rows() );
}

void SDCSolidSolver::setSolution(
    const fsi::vector & solution,
    const fsi::vector &
    )
{
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

    forAll( V.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            V.internalField()[i][j] = solution( index );
            index++;
        }
    }

    forAll( V.boundaryField(), patchI )
    {
        forAll( V.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                V.boundaryField()[patchI][i][j] = solution( index );
                index++;
            }
        }
    }

    assert( index == solution.rows() );
}

scalar SDCSolidSolver::getEndTime()
{
    return runTime->endTime().value();
}

scalar SDCSolidSolver::getTimeStep()
{
    return runTime->deltaT().value();
}

void SDCSolidSolver::nextTimeStep()
{
    timeIndex++;

    if ( runTime->timeIndex() != timeIndex )
        (*runTime)++;

    assert( runTime->timeIndex() == timeIndex );

    for ( int i = 0; i < k; i++ )
    {
        UStages.at( i ) = U;
        VStages.at( i ) = V;
    }
}

void SDCSolidSolver::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    UStages.clear();
    VStages.clear();

    for ( int i = 0; i < k + 1; i++ )
    {
        UStages.push_back( volVectorField( U ) );
        VStages.push_back( volVectorField( V ) );
    }
}

void SDCSolidSolver::prepareImplicitSolve(
    bool corrector,
    const int k,
    const int,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs
    )
{
    runTime->setDeltaT( dt );
    runTime->setTime( t, runTime->timeIndex() );
    indexk = k;

    if ( corrector )
    {
        U = UStages.at( k + 1 );
        V = VStages.at( k + 1 );
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
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                U.oldTime().boundaryField()[patchI][i][j] = qold( index );
                index++;
            }
        }
    }

    forAll( V.oldTime().internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            V.oldTime().internalField()[i][j] = qold( index );
            index++;
        }
    }

    forAll( V.oldTime().boundaryField(), patchI )
    {
        forAll( V.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                V.oldTime().boundaryField()[patchI][i][j] = qold( index );
                index++;
            }
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
        forAll( rhsU.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                rhsU.boundaryField()[patchI][i][j] = rhs( index );
                index++;
            }
        }
    }

    forAll( rhsV.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            rhsV.internalField()[i][j] = rhs( index );
            index++;
        }
    }

    forAll( rhsV.boundaryField(), patchI )
    {
        forAll( rhsV.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                rhsV.boundaryField()[patchI][i][j] = rhs( index );
                index++;
            }
        }
    }

    assert( index == rhs.rows() );
}

void SDCSolidSolver::implicitSolve(
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

scalar SDCSolidSolver::getStartTime()
{
    return runTime->startTime().value();
}

void SDCSolidSolver::getVariablesInfo(
    std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    )
{
    dof.push_back( 0 );
    dof.push_back( 0 );
    enabled.push_back( true );
    enabled.push_back( true );
    names.push_back( "solid U" );
    names.push_back( "solid V" );

    forAll( U.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            dof.at( 0 ) += 1;
        }
    }

    forAll( U.boundaryField(), patchI )
    {
        forAll( U.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                dof.at( 0 ) += 1;
            }
        }
    }

    forAll( V.internalField(), i )
    {
        for ( int j = 0; j < mesh.nGeometricD(); j++ )
        {
            dof.at( 1 ) += 1;
        }
    }

    forAll( V.boundaryField(), patchI )
    {
        forAll( V.boundaryField()[patchI], i )
        {
            for ( int j = 0; j < mesh.nGeometricD(); j++ )
            {
                dof.at( 1 ) += 1;
            }
        }
    }
}
