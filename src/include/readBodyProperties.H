    Info<< "Reading solid body properties\n" << endl;

    IOdictionary bodyMotionDict
    (
        IOobject
        (
            "bodyMotionDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    bodyCollector bodyCol(bodyMotionDict,mesh);

    // === Setting up fluid mesh motionSolver === //
    RBFMotionSolverExt& ms =
        const_cast<RBFMotionSolverExt&>
        (
            mesh.lookupObject<RBFMotionSolverExt>("dynamicMeshDict")
        );