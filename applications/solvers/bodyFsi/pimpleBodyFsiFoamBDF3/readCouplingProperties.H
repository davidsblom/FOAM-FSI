    Info << "Reading coupling properties" << endl;

    IOdictionary couplingProperties
    (
        IOobject
        (
            "couplingProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

	// === Setting up interface properties === //
	fsiInterface fsiInter(mesh,couplingProperties,bodyCol.getMovingPatchIDs());
	fsiInter.initializeForce(rhoFluid,nu,U,p);