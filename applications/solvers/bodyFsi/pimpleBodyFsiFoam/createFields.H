Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);

dimensionedScalar rhoFluid //Needed to calculate pressure at structure (p is actually p/rho)
(
    transportProperties.lookup("rho")
);


Info<< "Reading field p\n" << endl;
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

Info<< "Reading field U\n" << endl;
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

#   include "createPhi.H"
phi.oldTime();

//Creating field to store old face areas needed voor ddtPhiCorr
#   include "createUf.H"

//Read pressure properties and create turbulence model
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

Info<< "Reading field AU if present\n" << endl;
volScalarField AU
(
    IOobject
    (
        "AU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1.0/runTime.deltaT(),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Reading field HU if present\n" << endl;
volVectorField HU
(
    IOobject
    (
        "HU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    U.dimensions()/runTime.deltaT().dimensions(),
    zeroGradientFvPatchVectorField::typeName
);
