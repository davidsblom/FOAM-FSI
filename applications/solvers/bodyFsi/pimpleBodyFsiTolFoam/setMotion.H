	//Get motion based on forces from fsiInterface. fsiInter.getPatchForces() has forces in order of movingPatchIDs_ from bodyCollector
    //solve solid
    ms.setMotion(bodyCol.getPatchMotion(fsiInter.getPatchForces()));
    t = std::clock() - t;
    runTimeSolid += static_cast<float>(t) / CLOCKS_PER_SEC;
    t = std::clock();

    //move mesh
    bool meshChanged = mesh.update();
    t = std::clock() - t;
    runTimeMesh += static_cast<float>(t) / CLOCKS_PER_SEC;
    t = std::clock();

#	include "checkVolContinuity.H"//Check moving mesh continuity errors
