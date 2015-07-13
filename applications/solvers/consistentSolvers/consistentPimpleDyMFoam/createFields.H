
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
U.oldTime();

#include "createPhi.H"
phi.oldTime();

#include "createUf.H"
#include "createV000.H"

// Read pressure properties and create turbulence model
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell( p, mesh.solutionDict().subDict( "PIMPLE" ), pRefCell, pRefValue );

singlePhaseTransportModel laminarTransport( U, phi );

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New( U, phi, laminarTransport )
);

// Read or create fields for poison equation for the pressure
Info << "Reading field AU if present\n" << endl;
volScalarField AU
(
    IOobject
    (
        "AU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    1.0 / runTime.deltaT(),
    zeroGradientFvPatchScalarField::typeName
);

Info << "Reading field HU if present\n" << endl;
volVectorField HU
(
    IOobject
    (
        "HU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    U.dimensions() / runTime.deltaT().dimensions(),
    zeroGradientFvPatchVectorField::typeName
);
