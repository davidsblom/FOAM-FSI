// for all internal faces check theat the d dot S product is positive
const vectorField& centres = mesh.cellCentres();
const vectorField& areas = mesh.faceAreas();

const labelList& own = mesh.faceOwner();
const labelList& nei = mesh.faceNeighbour();

// Severe nonorthogonality threshold
//const scalar severeNonorthogonalityThreshold = ::cos(nonOrthThreshold_()/180.0*mathematicalConstant::pi);
const scalar nonOrthThreshold = 70;
const scalar severeNonorthogonalityThreshold = ::cos(nonOrthThreshold/180.0*mathematicalConstant::pi);

scalar minDDotS = GREAT;
scalar sumDDotS = 0;
label severeNonOrth = 0;
label errorNonOrth = 0;

Info << "Calculating internal field" << endl;
forAll (nei, faceI)
{
    vector d = centres[nei[faceI]] - centres[own[faceI]];
    const vector& s = areas[faceI];

    scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

    //Ensure round of errors do not cause value larger than 1
    dDotS = min(dDotS,1);

    //Save to internalField
    nonOrthoFaces.internalField()[faceI] = ::acos(dDotS)/mathematicalConstant::pi*180.0;

    if (dDotS < severeNonorthogonalityThreshold)
    {
        if (dDotS > SMALL)
        {
            severeNonOrth++;
        }
        else
        {
            errorNonOrth++;
        }
    }

    if (dDotS < minDDotS)
    {
        minDDotS = dDotS;
    }

    sumDDotS += dDotS;
}

reduce(minDDotS, minOp<scalar>());
reduce(sumDDotS, sumOp<scalar>());
reduce(severeNonOrth, sumOp<label>());
reduce(errorNonOrth, sumOp<label>());

label neiSize = nei.size();
reduce(neiSize, sumOp<label>());

if (neiSize > 0)
{
    Info<< "    Internal non-orthogonality Max: "
        << ::acos(minDDotS)/mathematicalConstant::pi*180.0
        << " average: " <<
        ::acos(sumDDotS/neiSize)/mathematicalConstant::pi*180.0
        << " Threshold = " << nonOrthThreshold
        << endl;
}

if (severeNonOrth > 0)
{
    Info<< "   *Number of severely non-orthogonal faces: "
        << severeNonOrth << "." << endl;
}

/*if (errorNonOrth > 0)
{
    Info<< " ***Number of non-orthogonality errors: "
        << errorNonOrth << "." << endl;
}
else
{

    Info<< "    Non-orthogonality check OK." << endl;
}*/

//Checking boundary
Info << "Check boundaries" << endl;
scalarList minDDotSPatches(nonOrthoFaces.boundaryField().size(),GREAT);
scalarList sumDDotSPatches(nonOrthoFaces.boundaryField().size(),0.0);
labelList psize(nonOrthoFaces.boundaryField().size(),0.0);

label pcounter = 0;
forAll(nonOrthoFaces.boundaryField(), iPatch)
{
    scalar pminDDotS = GREAT;
    scalar psumDDotS = 0;
    if(mesh.boundaryMesh()[iPatch].type()!="empty")
    {
        Info << "    Boundary patch: " << mesh.boundaryMesh()[iPatch].name() << endl;
        const polyPatch& patch = mesh.boundaryMesh()[iPatch];
        const vectorField& fcentres = patch.faceCentres();
        const labelList& fcells = patch.faceCells();

        forAll(fcentres,iface)
        {
            vector d = fcentres[iface] - centres[fcells[iface]];
            const vector& s = patch.faceAreas()[iface];

            scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

            //Ensure round of errors do not cause value larger than 1
            dDotS = min(dDotS,1);

            //Save to internalField
            nonOrthoFaces.boundaryField()[iPatch][iface] = ::acos(dDotS)/mathematicalConstant::pi*180.0;

            if (dDotS < pminDDotS)
            {
                pminDDotS = dDotS;
            }
            psumDDotS += dDotS;
        }

        if(fcentres.size() > 0)
        {
            Info<< "        Patch non-orthogonality Max: " << ::acos(pminDDotS)/mathematicalConstant::pi*180.0
                << " average: " << ::acos(psumDDotS/fcentres.size())/mathematicalConstant::pi*180.0 << endl;
        }

        minDDotSPatches[pcounter] = pminDDotS;
        sumDDotSPatches[pcounter] = psumDDotS;
        psize[pcounter] = fcentres.size();
        pcounter++;
    }
}

//Set size correctly and reduce for parallel processing
minDDotSPatches.setSize(pcounter);
sumDDotSPatches.setSize(pcounter);
psize.setSize(pcounter);
reduce(minDDotSPatches, minOp<scalarList>());
reduce(sumDDotSPatches, sumOp<scalarList>());
reduce(psize, sumOp<labelList>());

//Average non-orthogonality based on internalfield and boundary field
scalar avNonOrtho = sumDDotS;
label avSize = neiSize;
forAll(psize,ipatch)
{
    avNonOrtho += sumDDotSPatches[ipatch];
    avSize += psize[ipatch];
}
avNonOrtho /= avSize;
avNonOrtho = ::acos(avNonOrtho)/mathematicalConstant::pi*180.0;

//Max non-orthogonality based on internalfield and boundary field
scalar maxNonOrtho = ::acos(min(minDDotS, min(minDDotSPatches)))/mathematicalConstant::pi*180.0;

Info<< "Mesh non-orthogonality Max: " << maxNonOrtho << " average: " << avNonOrtho << endl;
