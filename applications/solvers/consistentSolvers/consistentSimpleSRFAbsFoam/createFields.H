
Info << "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << "Reading/calculating face flux field phi\n" << endl;
#include "createPhi.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell( p, mesh.solutionDict().subDict( "SIMPLE" ), pRefCell, pRefValue );

singlePhaseTransportModel laminarTransport( U, phi );

autoPtr<incompressible::RASModel> turbulence
(
    incompressible::RASModel::New( U, phi, laminarTransport )
);

Info << "Creating absolute SRF model\n" << endl;
autoPtr<SRF::SRFModelAbs> SRF
(
    SRF::SRFModelAbs::New( U )
);
SRF->correctBoundaryVelocity( U );

/*volVectorField Fcoriolis
 * (
 *  IOobject
 *  (
 *      "Fcoriolis",
 *      runTime.timeName(),
 *      mesh,
 *      IOobject::NO_READ,
 *      IOobject::AUTO_WRITE
 *  ),
 *  mesh,
 *  SRF->Fcoriolis()().dimensions(),
 *  fixedValueFvPatchVectorField::typeName
 * );
 * Fcoriolis.internalField() = SRF->Fcoriolis();*/

// Make phi relative
SRF->relativeFlux( phi );
