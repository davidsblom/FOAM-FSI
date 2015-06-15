{
    const labelList& movingPatchIDs = bodyCol.getMovingPatchIDs();
    const scalar& pi = mathematicalConstant::pi;

    scalar nMovingFaces = 0;
    forAll(movingPatchIDs,id)
    {
        nMovingFaces += mesh.boundaryMesh()[movingPatchIDs[id]].size();
    }

    //Initialize force field
    Field<vectorField> patchForces(mesh.boundaryMesh().size(),vectorField(0,vector::zero));
    forAll(patchForces,ipatch){
    	patchForces[ipatch]=vectorField(mesh.boundaryMesh()[ipatch].size(),vector::zero);
    }

    //Set forces if forcing is requisted
    if(forcing)
    {
        forAll(movingPatchIDs,id)
        {
            if(prescribedFrequency>SMALL)
            {
                patchForces[movingPatchIDs[id]] = prescribedForces/nMovingFaces*0.5*(1.0-Foam::cos(2*pi*prescribedFrequency*runTime.value()));
            }else{
                patchForces[movingPatchIDs[id]] = prescribedForces/nMovingFaces;
            }
        }
    }

    bodyCol.getPatchMotion(patchForces);
}
