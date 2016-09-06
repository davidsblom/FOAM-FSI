
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElasticSolidSolver.H"

ElasticSolidSolver::ElasticSolidSolver (
    const string & name,
    const std::shared_ptr<argList> & args,
    const std::shared_ptr<Time> & runTime
    )
    :
    foamSolidSolver( name, args, runTime ),
    gradU
    (
        IOobject
        (
            "grad(U)",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor( "zero", dimless, tensor::zero )
    ),
    snGradU
    (
        IOobject
        (
            "snGrad(U)",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector( "zero", dimless, Foam::vector::zero )
    ),
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
    gradV( fvc::ddt( gradU ) ),
    snGradV( ( snGradU - snGradU.oldTime() ) / runTime->deltaT() ),
    epsilon
    (
        IOobject
        (
            "epsilon",
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
    divSigmaExp
    (
        IOobject
        (
            "divSigmaExp",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector( "zero", dimForce / dimVolume, Foam::vector::zero )
    ),
    rheology( sigma, U ),
    rho( rheology.rho() ),
    mu( rheology.mu() ),
    lambda( rheology.lambda() ),
    muf( rheology.muf() ),
    lambdaf( rheology.lambdaf() ),
    n( mesh.Sf() / mesh.magSf() ),
    nCorr( readInt( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "nCorrectors" ) ) ),
    convergenceTolerance( readScalar( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "U" ) ) ),
    divSigmaExpMethod( mesh.solutionDict().subDict( "solidMechanics" ).lookup( "divSigmaExp" ) )
{}

ElasticSolidSolver::~ElasticSolidSolver()
{}

void ElasticSolidSolver::finalizeTimeStep()
{
    assert( init );

    V = fvc::ddt( U );
    gradV = fvc::ddt( gradU );
    snGradV = ( snGradU - snGradU.oldTime() ) / runTime->deltaT();

    init = false;
}

void ElasticSolidSolver::initTimeStep()
{
    assert( !init );

    timeIndex++;
    t = runTime->time().value();

    U += V * runTime->deltaT();
    gradU += gradV * runTime->deltaT();
    snGradU += snGradV * runTime->deltaT();

    init = true;
}

bool ElasticSolidSolver::isRunning()
{
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void ElasticSolidSolver::resetSolution()
{}

void ElasticSolidSolver::solve()
{
    Info << "Solve solid domain" << endl;

    int iCorr = 0;
    lduMatrix::solverPerformance solverPerf;
    scalar initialResidual = 1.0;
    lduMatrix::debug = 0;

    do
    {
        U.storePrevIter();

        if ( divSigmaExpMethod == "standard" )
        {
            divSigmaExp = fvc::div
                (
                mu * gradU.T() + lambda * ( I * tr( gradU ) ) - (mu + lambda) * gradU,
                "div(sigma)"
                );
        }
        else
        if ( divSigmaExpMethod == "surface" )
        {
            divSigmaExp = fvc::div
                (
                muf * ( mesh.Sf() & fvc::interpolate( gradU.T() ) )
                + lambdaf * ( mesh.Sf() & I * fvc::interpolate( tr( gradU ) ) )
                - (muf + lambdaf) * ( mesh.Sf() & fvc::interpolate( gradU ) )
                );
        }
        else
        if ( divSigmaExpMethod == "decompose" )
        {
            snGradU = fvc::snGrad( U );

            surfaceTensorField shearGradU = ( (I - n * n) & fvc::interpolate( gradU ) );

            divSigmaExp = fvc::div
                (
                mesh.magSf() *
                (
                    -(muf + lambdaf) * ( snGradU & (I - n * n) )
                    + lambdaf * tr( shearGradU & (I - n * n) ) * n
                    + muf * (shearGradU & n)
                )
                );
        }
        else
        if ( divSigmaExpMethod == "expLaplacian" )
        {
            divSigmaExp =
                -fvc::laplacian( mu + lambda, U, "laplacian(DU,U)" )
                + fvc::div( mu * gradU.T() + lambda * ( I * tr( gradU ) ), "div(sigma)" );
        }
        else
        {
            FatalErrorIn( args->executable() )
                << "divSigmaExp method " << divSigmaExpMethod << " not found!" << endl;
        }

        // linear momentum equation
        fvVectorMatrix UEqn
        (
            rho * fvm::d2dt2( U )
            ==
            fvm::laplacian( 2 * muf + lambdaf, U, "laplacian(DU,U)" )
            + divSigmaExp
        );

        solverPerf = UEqn.solve();

        if ( iCorr == 0 )
            initialResidual = solverPerf.initialResidual();

        U.relax();

        gradU = fvc::grad( U );
    }
    while
    (
        iCorr++ == 0
        ||
        (
            solverPerf.initialResidual() > convergenceTolerance
            && iCorr < nCorr
        )
    );

    lduMatrix::debug = 1;

    Info << "Solving for " << U.name();
    Info << ", Initial residual = " << initialResidual;
    Info << ", Final residual = " << solverPerf.initialResidual();
    Info << ", No outer iterations " << iCorr << endl;
}
