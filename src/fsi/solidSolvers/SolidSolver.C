
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SolidSolver.H"

SolidSolver::SolidSolver (
    string name,
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
    minIter( 0 ),
    maxIter( 0 ),
    absoluteTolerance( 0 ),
    relativeTolerance( 0 ),
    interpolator( nullptr )
{}

SolidSolver::SolidSolver (
    string name,
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
    minIter( 0 ),
    maxIter( 0 ),
    absoluteTolerance( 0 ),
    relativeTolerance( 0 ),
    interpolator( interpolator )
{}


SolidSolver::~SolidSolver()
{}

void SolidSolver::calculateEpsilonSigma()
{
    // - Green finite strain tensor
    epsilon = symm( gradU ) + 0.5 * symm( gradU & gradU.T() );

    // - second Piola-Kirchhoff stress tensor
    sigma = 2 * mu * epsilon + lambda * ( I * tr( epsilon ) );
}

void SolidSolver::initialize()
{
    readCouplingProperties();
}

void SolidSolver::initTimeStep()
{
    assert( !init );

    timeIndex++;
    t = runTime->time().value();

    readSolidMechanicsControls();

    init = true;
}

bool SolidSolver::interpolateVolField( std::shared_ptr<BaseMultiLevelSolver> solver )
{
    std::shared_ptr<SolidSolver> fineModel;
    fineModel = std::dynamic_pointer_cast<SolidSolver>( solver );

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

void SolidSolver::finalizeTimeStep()
{
    foamSolidSolver::finalizeTimeStep();

    if ( runTime->outputTime() )
        runTime->write();
}

bool SolidSolver::isRunning()
{
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void SolidSolver::readSolidMechanicsControls()
{
    const dictionary & stressControl =
        mesh.solutionDict().subDict( "solidMechanics" );

    minIter = readInt( stressControl.lookup( "minIter" ) );
    maxIter = readInt( stressControl.lookup( "maxIter" ) );
    absoluteTolerance = readScalar( stressControl.lookup( "tolerance" ) );
    relativeTolerance = readScalar( stressControl.lookup( "relTol" ) );

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

    if ( linearTolerance >= absoluteTolerance )
        throw std::runtime_error( "The absolute tolerance for the linear solver (U) should be smaller than solidMechanics::absoluteTolerance in order to reach convergence of the non-linear system" );
}

void SolidSolver::resetSolution()
{}

void SolidSolver::solve()
{
    Info << "Solve solid domain" << endl;

    scalar iCorr = 0;
    scalar displacementResidual = 1;
    scalar initialResidual = 1;
    lduMatrix::solverPerformance solverPerf;
    lduMatrix::debug = 0;
    scalar convergenceTolerance = absoluteTolerance;

    gradU = fvc::grad( U );

    calculateEpsilonSigma();

    dimensionedVector gravity( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "gravity" ) );

    for ( iCorr = 0; iCorr < maxIter; iCorr++ )
    {
        U.storePrevIter();

        tmp<surfaceTensorField> shearGradU =
            ( (I - n * n) & fvc::interpolate( gradU ) );

        fvVectorMatrix UEqn
        (
            rho * fvm::d2dt2( U )
            ==
            fvm::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
            + fvc::div(
                mesh.magSf()
                * (
                    -(muf + lambdaf) * ( fvc::snGrad( U ) & (I - n * n) )
                    + lambdaf * tr( shearGradU() & (I - n * n) ) * n
                    + muf * (shearGradU() & n)
                    + muf * ( n & fvc::interpolate( gradU & gradU.T() ) )
                    + 0.5 * lambdaf
                    * ( n * tr( fvc::interpolate( gradU & gradU.T() ) ) )
                    + ( n & fvc::interpolate( sigma & gradU ) )
                    )
                )
        );

        // Add gravity

        UEqn -= rho * gravity;

        solverPerf = UEqn.solve();

        U.relax();

        gradU = fvc::grad( U );

        calculateEpsilonSigma();

        displacementResidual = gSumMag( U.internalField() - U.prevIter().internalField() ) / (gSumMag( U.internalField() ) + SMALL);
        displacementResidual = max( displacementResidual, solverPerf.initialResidual() );

        if ( iCorr == 0 )
        {
            initialResidual = displacementResidual;
            convergenceTolerance = std::max( relativeTolerance * displacementResidual, absoluteTolerance );
            assert( convergenceTolerance > 0 );
        }

        bool convergence = displacementResidual <= convergenceTolerance && iCorr >= minIter - 1;

        if ( convergence )
            break;
    }

    lduMatrix::debug = 1;

    Info << "Solving for " << U.name();
    Info << ", Initial residual = " << initialResidual;
    Info << ", Final residual = " << displacementResidual;
    Info << ", No outer iterations " << iCorr << endl;
}
