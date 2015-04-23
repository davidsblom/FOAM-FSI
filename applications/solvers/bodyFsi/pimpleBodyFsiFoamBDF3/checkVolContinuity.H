if (mesh.moving())
{
    // The ddt term constructed by hand because it must be Euler
    const scalar dt = runTime.deltaT().value();

    const scalarField& V = mesh.V().field();
    const scalarField& V0 = mesh.V0().field();
    scalar sumV = gSum(V);

    /*scalarField volChange = (1.0 - V0/V)/dt;
    scalarField divFlux = - fvc::div(mesh.phi())().internalField();
    scalarField conserve = (volChange + divFlux);
    scalar sumLocalContErr = dt*gSum(mag(conserve*mesh.V()))/sumV;
    scalar globalContErr = dt*gSum(conserve*mesh.V())/sumV;*/
    
    scalarField volChange = (V - V0)/dt;
    scalarField divFlux = - V*fvc::div(mesh.phi())().internalField();
    scalarField conserve = (volChange + divFlux);
    scalar sumLocalContErr = dt*gSum(mag(conserve))/sumV;
    scalar globalContErr = dt*gSum(conserve)/sumV;

    label nMotionErrors = 0;
    scalar maxMotionError = 0;

    forAll (conserve, cellI)
    {
        maxMotionError = Foam::max(maxMotionError, mag(conserve[cellI]));

        /*if (mag(conserve[cellI]) > max(1e-12/dt,1e-8))
        {
            Info<< "Motion conservation error in cell " << cellI << ": "
                << conserve[cellI] << " > " << max(1e-12/dt,1e-8)
                << " V: " << V[cellI] << " V0: " << V0[cellI]
                << " volume change: " << volChange[cellI]
                << " div flux: " << divFlux[cellI] << endl;

            nMotionErrors++;
        }*/
    }

    if (nMotionErrors > 0)
    {
        Info<< "Number of MOTION ERRORS: " << nMotionErrors
            << " out of nCells = " << mesh.nCells() << endl;
    }

    Info<< "volume continuity errors : "
        << "volume = " << sumV
        << ", max error = " << maxMotionError
        << ", sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}
