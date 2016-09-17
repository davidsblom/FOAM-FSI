
/*
 * Class
 *   CoupledFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "CoupledFluidSolver.H"

CoupledFluidSolver::CoupledFluidSolver(const string & name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    )
    :
    foamFluidSolver(name, args, runTime),
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
            IOobject::NO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    ),
    Up
    (
        IOobject
        (
            "Up",
            runTime->timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector4("zero", dimless, vector4::zero)
    ),
    laminarTransport(U, phi),
    turbulence(incompressible::turbulenceModel::New(U, phi, laminarTransport)),
    sumLocalContErr(0),
    globalContErr(0),
    cumulativeContErr(0),
    pRefCell(0),
    pRefValue(0),
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
    nu(transportProperties.lookup("nu")),
    rho(transportProperties.lookup("rho")),
    convergenceTolerance(readScalar(mesh.solutionDict().subDict("blockSolver").lookup("convergenceTolerance"))),
    nOuterCorr(readLabel(mesh.solutionDict().subDict("blockSolver").lookup("nOuterCorrectors"))),
    totalVolume(sum(mesh.V()).value()),
    CoNum(0),
    meanCoNum(0),
    velMag(0) {
    // Ensure that the absolute tolerance of the linear solver is less than the
    // used convergence tolerance for the non-linear system.
    scalar absTolerance = readScalar(mesh.solutionDict().subDict("solvers").subDict("Up").lookup("tolerance"));
    assert(absTolerance < convergenceTolerance);

    if (absTolerance >= convergenceTolerance)
        throw std::runtime_error("The absolute tolerance for the linear solver Up should be smaller than blockSolver::convergenceTolerance in order to reach convergence of the non-linear system");

    readBlockSolverControls();
}

CoupledFluidSolver::~CoupledFluidSolver() {}

void CoupledFluidSolver::checkTotalVolume() {
    scalar newTotalVolume = gSum(mesh.cellVolumes());
    scalar totalVolRatio = newTotalVolume / totalVolume;

    Info << "Volume: new = " << newTotalVolume << " old = " << totalVolume
         << " change = " << Foam::mag(newTotalVolume - totalVolume)
         << " ratio = " << totalVolRatio - 1 << endl;

    totalVolume = newTotalVolume;
}

void CoupledFluidSolver::continuityErrs() {
    tmp<volScalarField> contErr = fvc::div(phi);

    sumLocalContErr = runTime->deltaT().value() *
        mag(contErr()) ().weightedAverage(mesh.V()).value();

    globalContErr = runTime->deltaT().value() *
        contErr->weightedAverage(mesh.V()).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr
         << endl;
}

void CoupledFluidSolver::courantNo() {
    CoNum = 0.0;
    meanCoNum = 0.0;
    velMag = 0.0;

    if (mesh.nInternalFaces()) {
        tmp<surfaceScalarField> magPhi = mag(phi);

        tmp<surfaceScalarField> SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs() * magPhi();

        CoNum = max(SfUfbyDelta() / mesh.magSf())
            .value() * runTime->deltaT().value();

        meanCoNum = (sum(SfUfbyDelta()) / sum(mesh.magSf()))
            .value() * runTime->deltaT().value();

        velMag = max(magPhi() / mesh.magSf()).value();
    }

    Info << "Courant Number mean: " << meanCoNum
         << " max: " << CoNum
         << " velocity magnitude: " << velMag
         << endl;
}

void CoupledFluidSolver::getAcousticsDensityLocal(matrix &) {
    assert(false);
}

void CoupledFluidSolver::getAcousticsVelocityLocal(matrix &) {
    assert(false);
}

void CoupledFluidSolver::getAcousticsPressureLocal(matrix &) {
    assert(false);
}

void CoupledFluidSolver::getTractionLocal(matrix & traction) {
    int size = 0;

    forAll(movingPatchIDs, patchI)
    {
        size += mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();
    }

    vectorField tractionField(size, Foam::vector::zero);

    int offset = 0;

    forAll(movingPatchIDs, patchI)
    {
        int size = mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres().size();

        tmp<vectorField> tractionFieldPatchI = -rho.value() * nu.value()
            * U.boundaryField()[movingPatchIDs[patchI]].snGrad()
            + rho.value() * p.boundaryField()[movingPatchIDs[patchI]]
            * mesh.boundary()[movingPatchIDs[patchI]].nf();

        forAll(tractionFieldPatchI(), i)
        {
            tractionField[i + offset] = tractionFieldPatchI()[i];
        }

        offset += size;
    }

    assert(tractionField.size() == nGlobalCenters[Pstream::myProcNo()]);

    traction.resize(tractionField.size(), mesh.nGeometricD());

    for (int i = 0; i < traction.rows(); i++)
        for (int j = 0; j < traction.cols(); j++)
            traction(i, j) = tractionField[i][j];
}

void CoupledFluidSolver::getWritePositionsLocalAcoustics(matrix &) {
    assert(false);
}

void CoupledFluidSolver::initTimeStep() {
    assert(!init);

    timeIndex++;
    t = timeIndex * runTime->deltaT().value();

    Info << "\nTime = " << runTime->timeName() << endl;

    checkTotalVolume();
    courantNo();

    init = true;
}

bool CoupledFluidSolver::isRunning() {
    runTime->write();

    Info << "ExecutionTime = " << runTime->elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime->elapsedClockTime() << " s"
         << endl << endl;

    return runTime->loop();
}

void CoupledFluidSolver::readBlockSolverControls() {
    setRefCell(p, mesh.solutionDict().subDict("blockSolver"), pRefCell, pRefValue);
}

void CoupledFluidSolver::resetSolution() {
    U == U.oldTime();
    p == p.oldTime();
    phi == phi.oldTime();
}

void CoupledFluidSolver::solve() {
    Info << "Solve fluid domain" << endl;

    mesh.update();

    int oCorr;
    bool underrelaxation = true;

    // Outer correction loop to solve the non-linear system
    for (oCorr = 0; oCorr < nOuterCorr; oCorr++) {
        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        p.storePrevIter();
        U.storePrevIter();

        // Initialize the Up block system (matrix, source and reference to Up)
        fvBlockMatrix<vector4> UpEqn(Up);

        // Momentum equation
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
            + fvm::div(phi, U)
            + turbulence->divDevReff(U)
        );

        if (underrelaxation)
            UEqn.relax();

        UpEqn.insertEquation(0, UEqn);

        // Assemble and insert pressure equation

        // Pressure parts of the continuity equation
        surfaceScalarField rUAf
        (
            "rUAf",
            fvc::interpolate(1.0 / UEqn.A())
        );

        surfaceScalarField presSource
        (
            "presSource",
            rUAf * (fvc::interpolate(fvc::grad(p, "grad(pSource)")) & mesh.Sf())
        );

        fvScalarMatrix pEqn
        (
            -fvm::laplacian(rUAf, p)
            ==
            -fvc::div(presSource)
        );

        pEqn.setReference(pRefCell, pRefValue);

        UpEqn.insertEquation(3, pEqn);

        // Calculate grad p coupling matrix. Needs to be here if one uses
        // gradient schemes with limiters.  VV, 9/June/2014
        BlockLduSystem<Foam::vector, Foam::vector> pInU(fvm::grad(p));

        // Calculate div U coupling.  Could be calculated only once since
        // it is only geometry dependent.  VV, 9/June/2014
        BlockLduSystem<Foam::vector, scalar> UInp(fvm::UDiv(U));

        // Last argument in insertBlockCoupling says if the column direction
        // should be incremented. This is needed for arbitrary positioning
        // of U and p in the system. This could be better. VV, 30/April/2014
        UpEqn.insertBlockCoupling(0, 3, pInU, true);
        UpEqn.insertBlockCoupling(3, 0, UInp, false);

        // Solve the block matrix
        UpEqn.solve();

        // Retrieve solution
        UpEqn.retrieveSolution(0, U.internalField());
        UpEqn.retrieveSolution(3, p.internalField());

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

        if (underrelaxation)
            p.relax();

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        turbulence->correct();

        tmp<volVectorField> residual = fvc::ddt(U) + fvc::div(phi, U) + (turbulence->divDevReff(U) & U) + fvc::grad(p);

        tmp<scalarField> magResU = mag(residual->internalField());
        scalar momentumResidual = std::sqrt(gSumSqr(magResU) / mesh.globalData().nTotalCells());
        scalar rmsU = std::sqrt(gSumSqr(mag(U.internalField())) / mesh.globalData().nTotalCells());
        rmsU /= runTime->deltaT().value();

        // Scale the residual by the root mean square of the velocity field
        momentumResidual /= rmsU;

        int minIter = 2;
        bool convergence = momentumResidual <= convergenceTolerance && oCorr >= minIter - 1;

        Info << "root mean square residual norm = " << momentumResidual;
        Info << ", tolerance = " << convergenceTolerance;
        Info << ", iteration = " << oCorr + 1;
        Info << ", convergence = ";

        if (convergence)
            Info << "true";
        else
            Info << "false";

        Info << endl;

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        if (convergence)
            break;
    }

    continuityErrs();
}
