
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "CoupledSolidSolver.H"

using namespace fsi;

CoupledSolidSolver::CoupledSolidSolver(
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    )
    :
    foamSolidSolver( name, args, runTime ),
    solidMesh( mesh ),
    pMesh( mesh ),
    volToPointInterp( mesh ),
    gradU
    (
    IOobject
    (
        "grad(U)",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimless, tensor::zero )
    ),
    gradUf
    (
    IOobject
    (
        "gradUf",
        runTime->timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor( "zero", dimless, tensor::zero )
    ),
    epsilon
    (
    IOobject
    (
        "epsilon",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor( "zero", dimless, symmTensor::zero )
    ),
    sigma
    (
    IOobject
    (
        "sigma",
        runTime->timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor( "zero", dimForce / dimArea, symmTensor::zero )
    ),
    n( mesh.Sf() / mesh.magSf() ),
    rheology( sigma ),
    rho( rheology.rho() ),
    mu( rheology.mu() ),
    lambda( rheology.lambda() ),
    muf( fvc::interpolate( mu, "mu" ) ),
    lambdaf( fvc::interpolate( lambda, "lambda" ) ),
    stressControl( mesh.solutionDict().subDict( "solidMechanics" ) ),
    nCorr( stressControl.lookupOrDefault<int>("nCorrectors", 1000) )
{
    createNonLinearGeometry();

    BlockLduMatrix<Foam::vector>::debug = 5;
}

CoupledSolidSolver::~CoupledSolidSolver()
{}

void CoupledSolidSolver::createNonLinearGeometry()
{
    /*
     * Add nonlinear geometry entry to solution dict
     * This can then be accessed by boundary conditions, constitutive model, etc to
     * know what type of solver is being used.
     */
    Info << "Setting nonLinear to off" << nl << endl;
    mesh.solutionDict().subDict( "stressedFoam" ).set
    (
        "nonLinear",
        nonLinearGeometry::nonLinearNames_
        [
            nonLinearGeometry::OFF
        ]
    );
}

void CoupledSolidSolver::solve()
{
    Info << "Solve solid domain" << endl;

    BlockLduMatrix<Foam::vector> blockM( solidMesh );

    // Grab block diagonal and set it to zero

    Field<tensor> & d = blockM.diag().asSquare();
    d = tensor::zero;

    // Grab linear off-diagonal and set it to zero

    Field<tensor> & l = blockM.lower().asSquare();
    Field<tensor> & u = blockM.upper().asSquare();
    u = tensor::zero;
    l = tensor::zero;

    // Solution vector

    vectorField & blockX = U.internalField();

    // Source vector

    vectorField blockB( mesh.nCells(), Foam::vector::zero );

    // Discretise by hand

    Info << "Calculating required addressing" << endl;
    const unallocLabelList & fvOwner = mesh.owner();
    const unallocLabelList & fvNeighbour = mesh.neighbour();
    const scalarField & deltaCoeffsI = mesh.deltaCoeffs().internalField();
    const scalarField & magSfI = mesh.magSf().internalField();
    const scalarField & mufI = muf.internalField();
    const scalarField & lambdafI = lambdaf.internalField();
    const vectorField & points = mesh.points();
    const vectorField & nI = n.internalField();
    const faceList & faces = mesh.faces();
    const cellList & cells = mesh.cells();
    const labelListList & pointCells = mesh.pointCells();
    const labelListList & edgeCells = mesh.edgeCells();
    const vectorField CI = mesh.C().internalField();
    const labelListList & faceEdges = mesh.faceEdges();

    const bool nonOrthogonalMesh = !mesh.orthogonal();
    const vectorField * corrVecIPtr = NULL;

    if ( nonOrthogonalMesh )
    {
        corrVecIPtr = &mesh.correctionVectors().internalField();
    }

    const labelList & fvMap = solidMesh.fvMeshAddressingMap();
    const edgeList & imagEdges = solidMesh.imagEdges();
    const scalarListList & pointCellWeights = solidMesh.pointWeights();
    const scalarListList & edgeCellWeights = solidMesh.edgeWeights();
    const labelListList & cellImagEdges = solidMesh.cellImagEdges();


    Info << "Adding normal derivatve terms to matrix" << endl;
    forAll( fvOwner, faceI )
    {
        const label own = fvOwner[faceI];
        const label nei = fvNeighbour[faceI];

        const Foam::vector & faceN = nI[faceI];
        const scalar faceMu = mufI[faceI];
        const scalar faceLambda = lambdafI[faceI];
        const scalar faceMagSf = magSfI[faceI];
        const scalar faceDeltaCoeff = deltaCoeffsI[faceI];

        // Normal derivative terms
        const tensor coeff =
            sqr( faceN ) * (2 * faceMu + faceLambda) * faceMagSf * faceDeltaCoeff
            + ( I - sqr( faceN ) ) * faceMu * faceMagSf * faceDeltaCoeff;

        d[own] -= coeff;
        d[nei] -= coeff;

        u[fvMap[faceI]] += coeff;
        l[fvMap[faceI]] += coeff;
    }


    // Tangential derivative terms
    // forAll edges of this face
    // add contribution to pointCells of edge start and end
    // contribution for each pointCellI of start/end is:
    // a += 0.5*mu*w*edgeN*faceN + 0.5*lambda*w*faceN*edgeN
    // where w is the pointCell weight
    // Note that magSf does not feature as it is cancelled out for
    // finite volume Gauss (multiply by magSf) and finite area Gauss
    // (divide by magSf)
    // challenges: addressing and efficiency

    // Method
    // for all cells
    // for all faces of cellI
    // for all edges of faceI
    // for all pointCells of edge.start() and edge.end()
    // add constributions to d, u and l
    // end
    // end
    // end
    // end

    // Create bool list with points on boundary
    boolList pointsOnBoundary( mesh.allPoints().size(), false );
    forAll( mesh.boundaryMesh(), patchI )
    {
        const labelList & meshPoints = mesh.boundaryMesh()[patchI].meshPoints();

        forAll( meshPoints, pointI )
        {
            pointsOnBoundary[meshPoints[pointI]] = true;
        }
    }
    boolList edgesOnBoundary( mesh.nEdges(), false );
    forAll( mesh.boundaryMesh(), patchI )
    {
        const labelList & meshEdges = mesh.boundaryMesh()[patchI].meshEdges();

        forAll( meshEdges, edgeI )
        {
            edgesOnBoundary[meshEdges[edgeI]] = true;
        }
    }


    Info << "Adding tangential derivative terms to the matrix" << endl;
    forAll( cells, cellI )
    {
        // Info<< "cell " << cellI << endl;

        // const cell& curCell = cells[cellI];
        const labelList & curCellFaces = cells[cellI];

        forAll( curCellFaces, faceI )
        {
            const label curFaceID = curCellFaces[faceI];

            if ( mesh.isInternalFace( curFaceID ) )
            {
                // Info<< "    face " << curFaceID << " n " << nI[curFaceID] << endl;
                const face & curFace = faces[curFaceID];
                const edgeList curFaceEdges = curFace.edges();

                const label own = fvOwner[curFaceID];

                // const label nei = fvNeighbour[curFaceID];

                const scalar faceMu = mufI[curFaceID];
                const scalar faceLambda = lambdafI[curFaceID];

                // Remember if cellI is not the owner of curFace
                // Flip face normal if cellI is not the owner
                const Foam::vector & faceN = nI[curFaceID];
                bool flipFace = false;

                if ( cellI != own )
                {
                    flipFace = true;
                }

                // Non-orthogonal correction vector
                const Foam::vector * faceKPtr = NULL;

                if ( nonOrthogonalMesh )
                {
                    faceKPtr = &( (*corrVecIPtr)[curFaceID] );
                }

                // Store (2*mu + lambda) as it is used for each edge
                const scalar faceTwoMuLambda = 2 * faceMu + faceLambda;

                forAll( curFaceEdges, edgeI )
                {
                    const edge & curEdge = curFaceEdges[edgeI];

                    // Calculate Le*faceN which features in all the coefficients

                    // Projected edge vector
                    Foam::vector e = curEdge.vec( points );
                    e -= faceN * (faceN & e);

                    // Edge length vector
                    Foam::vector Le = e ^ faceN;
                    Le *= curFace.edgeDirection( curEdge );

                    // LeFaceN term is the tensor component of the coeff
                    // Change sign if cellI is not the owner as faceN is opposite
                    tensor LeFaceN = Le * faceN;

                    if ( flipFace )
                    {
                        LeFaceN *= -1;
                    }

                    // Non-orthogonal correction component
                    scalar * kDotLePtr = NULL;

                    if ( nonOrthogonalMesh )
                    {
                        kDotLePtr = new scalar( (*faceKPtr) & Le );

                        if ( flipFace )
                        {
                            (*kDotLePtr) *= -1;
                        }
                    }

                    // There are three types of edges:
                    // - internal edges
                    // no points on the boundary
                    // - near boundary edges
                    // internal edge with one or two points on the
                    // boundary
                    // - boundary edges
                    // edge lies along the boundary
                    // We must perform the discretisation separately on
                    // edge type
                    const label curEdgeID = faceEdges[curFaceID][edgeI];
                    const bool boundaryEdge = edgesOnBoundary[curEdgeID];
                    const label startPointID = curEdge.start();
                    const label endPointID = curEdge.end();
                    const bool startOnBoundary = pointsOnBoundary[startPointID];
                    const bool endOnBoundary = pointsOnBoundary[endPointID];
                    const bool nearBoundaryEdge =
                        !boundaryEdge && (startOnBoundary || endOnBoundary);

                    // const bool bothPointsOnBoundary =
                    // startOnBoundary && endOnBoundary;

                    // Info<< "        edge " << curEdgeID
                    // << " s " << startPointID << " e " << endPointID << endl;

                    if ( nearBoundaryEdge )
                    {
                        // If the start or end cell is on the boundary, then
                        // its value is calculated by extrapolating from the
                        // edge centre:
                        // Umiddle = 0.5*(Ustart + Uend)
                        // Therefore:
                        // Uend = 2*Umiddle - Ustart
                        // where here we assume that the end is the boundary
                        // Ustart is implicitly interpolated from pointCells
                        // using inverse distance interpolation, and Umiddle
                        // is implicitly interpolated from edgeCells using
                        // inverse distance interpolation.
                        // This is equivalent to interpolating directly to the edge
                        // centre from the edgeCells using inverse distance.
                        // So we only need to add coeff contributions to the edge
                        // cells.

                        const scalarList & curEdgeCellWeights =
                            edgeCellWeights[curEdgeID];
                        const labelList & curEdgeCells = edgeCells[curEdgeID];

                        forAll( curEdgeCells, edgeCellI )
                        {
                            const scalar wCell = curEdgeCellWeights[edgeCellI];

                            // Info<< "        edgeCell " << curEdgeCells[edgeCellI]
                            // << " w " << wCell << endl;

                            // There is no (1/2) factor as we directly interpolate
                            // edge centre value i.e. coeff is twice that of
                            // internal edge
                            tensor coeff =
                                faceMu * wCell * LeFaceN + faceLambda * wCell * LeFaceN.T();

                            // Add non-orthogonal correction term
                            if ( nonOrthogonalMesh )
                            {
                                coeff +=
                                    tensor
                                    (
                                    sqr( faceN ) *
                                    faceTwoMuLambda * wCell * (*kDotLePtr)
                                    + ( I - sqr( faceN ) )
                                    * faceMu * wCell * (*kDotLePtr)
                                    );
                            }

                            // Add coeff to the matrix
                            if ( curEdgeCells[edgeCellI] == cellI )
                            {
                                // Contributions to diag
                                // Info<< "        nbe d[" << cellI
                                // << "] += " << coeff << endl;
                                d[cellI] += coeff;
                            }
                            else
                            {
                                // Find which imaginary edge connects cellI to
                                // pointCellI
                                const label curImEdgeID =
                                    solidMesh.findImaginaryEdge
                                    (
                                    cellI,
                                    curEdgeCells[edgeCellI]
                                    );

                                if ( curImEdgeID == -1 )
                                {
                                    FatalErrorIn
                                    (
                                        "insertTangential.H"
                                    ) << "curImEdge not found"
                                      << abort( FatalError );
                                }

                                // Contributions to upper/lower
                                if ( cellI > curEdgeCells[edgeCellI] )
                                {
                                    l[curImEdgeID] += coeff;
                                }
                                else
                                {
                                    u[curImEdgeID] += coeff;
                                }
                            }
                        } // forAll edgeCells fo current edge
                    } // if near-boundary edge
                    else
                    if ( boundaryEdge )
                    {
                        // Edge lies along boundary
                        // Point value are extrapolated along the adjoining
                        // edges within the current face

                        // startEdge is the edge connected to the start of the
                        // curEdge, and endEdge is the edge connected to the end
                        // of the curEdge

                        // We add 2*coeff to startEdge and endEdge edgeCells, and
                        // -coeff to startEdge and endEdge otherVertex pointCells

                        // Find startEdge and endEdge
                        label startEdgeLocalID = -1;
                        label endEdgeLocalID = -1;
                        forAll( curFaceEdges, eI )
                        {
                            if ( eI != edgeI )
                            {
                                const edge & otherEdge = curFaceEdges[eI];

                                if
                                (
                                    otherEdge.start() == curEdge.start()
                                    || otherEdge.end() == curEdge.start()
                                )
                                {
                                    startEdgeLocalID = eI;
                                }
                                else
                                if
                                (
                                    otherEdge.start() == curEdge.end()
                                    || otherEdge.end() == curEdge.end()
                                )
                                {
                                    endEdgeLocalID = eI;
                                }
                            }
                        }

                        if ( startEdgeLocalID == -1 || endEdgeLocalID == -1 )
                        {
                            FatalErrorIn
                            (
                                "insertTangential.H"
                            ) << "startEdge and/or endEdge not found"
                              << abort( FatalError );
                        }

                        // Add 2*coeff to edgeCells of startEdge and endEdge
                        const label startEdgeID =
                            faceEdges[curFaceID][startEdgeLocalID];
                        const label endEdgeID =
                            faceEdges[curFaceID][endEdgeLocalID];

                        labelList seEdgeIDs( 2, -1 );
                        seEdgeIDs[0] = startEdgeID;
                        seEdgeIDs[1] = endEdgeID;

                        forAll( seEdgeIDs, seEdgeI )
                        {
                            const label curSeEdgeID = seEdgeIDs[seEdgeI];

                            const scalarList & seEdgeCellWeights =
                                edgeCellWeights[curSeEdgeID];

                            const labelList & seEdgeCells = edgeCells[curSeEdgeID];

                            forAll( seEdgeCells, edgeCellI )
                            {
                                const scalar wCell = seEdgeCellWeights[edgeCellI];

                                // Note: coeff is 2*coeff of internal edge
                                // pointCells
                                tensor coeff =
                                    faceMu * wCell * LeFaceN
                                    + faceLambda * wCell * LeFaceN.T();

                                // Add non-orthogonal correction term
                                // Note: term is 2*coeff of internal edge
                                if ( nonOrthogonalMesh )
                                {
                                    coeff +=
                                        tensor
                                        (
                                        sqr( faceN )
                                        * faceTwoMuLambda * wCell * (*kDotLePtr)
                                        + ( I - sqr( faceN ) )
                                        * faceMu * wCell * (*kDotLePtr)
                                        );
                                }

                                // Add coeff to the matrix
                                if ( seEdgeCells[edgeCellI] == cellI )
                                {
                                    // Contributions to diag
                                    // Info<< "        be d[" << cellI
                                    // << "] += " << coeff << endl;
                                    d[cellI] += coeff;
                                }
                                else
                                {
                                    // Find which imaginary edge connects cellI to
                                    // pointCellI
                                    const label curImEdgeID =
                                        solidMesh.findImaginaryEdge
                                        (
                                        cellI,
                                        seEdgeCells[edgeCellI]
                                        );

                                    if ( curImEdgeID == -1 )
                                    {
                                        FatalErrorIn
                                        (
                                            "insertTangential.H"
                                        ) << "curImEdge not found"
                                          << abort( FatalError );
                                    }

                                    // Contributions to upper/lower
                                    if ( cellI > seEdgeCells[edgeCellI] )
                                    {
                                        l[curImEdgeID] += coeff;
                                    }
                                    else
                                    {
                                        u[curImEdgeID] += coeff;
                                    }
                                }
                            } // forAll edgeCells of seEdge
                        } // for startEdge and endEdge


                        // Now add -coeff to pointCells of far vertex of startEdge
                        // and endEdge

                        // Find start and end other vertices
                        const edge & startEdge = curFaceEdges[startEdgeLocalID];
                        const edge & endEdge = curFaceEdges[endEdgeLocalID];
                        const label startOtherVertexID =
                            startEdge.otherVertex( startPointID );
                        const label endOtherVertexID =
                            endEdge.otherVertex( endPointID );

                        labelList sePointIDs( 2, -1 );
                        sePointIDs[0] = startOtherVertexID;
                        sePointIDs[1] = endOtherVertexID;

                        forAll( sePointIDs, pointI )
                        {
                            // otherVertex of startEdge or endEdge
                            const label sePointID = sePointIDs[pointI];

                            const scalarList & seCellWeights =
                                pointCellWeights[sePointID];
                            const labelList & sePointCells = pointCells[sePointID];

                            forAll( sePointCells, pointCellI )
                            {
                                const scalar wCell = seCellWeights[pointCellI];

                                // Note: coeff is negative for otherVertex
                                tensor coeff =
                                    -0.5 * faceMu * wCell * LeFaceN
                                    - 0.5 * faceLambda * wCell * LeFaceN.T();

                                // Add non-orthogonal correction term
                                if ( nonOrthogonalMesh )
                                {
                                    coeff -=
                                        tensor
                                        (
                                        sqr( faceN ) *
                                        0.5 * faceTwoMuLambda * wCell * (*kDotLePtr)
                                        + ( I - sqr( faceN ) ) *
                                        0.5 * faceMu * wCell * (*kDotLePtr)
                                        );
                                }

                                if ( sePointCells[pointCellI] == cellI )
                                {
                                    // Contributions to diag
                                    // Info<< "        -be d[" << cellI
                                    // << "] += " << coeff << endl;
                                    d[cellI] += coeff;
                                }
                                else
                                {
                                    // Find which imaginary edge connects cellI to
                                    // pointCellI
                                    const label curImEdgeID =
                                        solidMesh.findImaginaryEdge
                                        (
                                        cellI,
                                        sePointCells[pointCellI]
                                        );

                                    if ( curImEdgeID == -1 )
                                    {
                                        FatalErrorIn
                                        (
                                            "insertTangential.H"
                                        ) << "curImEdge not found"
                                          << abort( FatalError );
                                    }

                                    // Contributions to upper/lower
                                    if ( cellI > sePointCells[pointCellI] )
                                    {
                                        l[curImEdgeID] += coeff;
                                    }
                                    else
                                    {
                                        u[curImEdgeID] += coeff;
                                    }
                                }
                            } // forAll point cells
                        } // for start and end other vertices
                    }
                    else // internal edge
                    {
                        // Put start and end ID in a list so we will perform same
                        // actions to both start and end
                        labelList sePointIDs( 2, -1 );
                        sePointIDs[0] = startPointID;
                        sePointIDs[1] = endPointID;

                        // For internal edges, we add coeff contributions to the
                        // pointCells of the start and end of the edge
                        forAll( sePointIDs, pointI )
                        {
                            // start or end point of edge
                            const label sePointID = sePointIDs[pointI];

                            const scalarList & seCellWeights =
                                pointCellWeights[sePointID];
                            const labelList & sePointCells = pointCells[sePointID];

                            forAll( sePointCells, pointCellI )
                            {
                                const scalar wCell = seCellWeights[pointCellI];

                                tensor coeff =
                                    0.5 * faceMu * wCell * LeFaceN
                                    + 0.5 * faceLambda * wCell * LeFaceN.T();

                                // Add non-orthogonal correction term
                                if ( nonOrthogonalMesh )
                                {
                                    coeff +=
                                        tensor
                                        (
                                        sqr( faceN ) *
                                        0.5 * faceTwoMuLambda * wCell * (*kDotLePtr)
                                        + ( I - sqr( faceN ) ) *
                                        0.5 * faceMu * wCell * (*kDotLePtr)
                                        );
                                }

                                if ( sePointCells[pointCellI] == cellI )
                                {
                                    // Contributions to diag
                                    d[cellI] += coeff;
                                }
                                else
                                {
                                    // Find which imaginary edge connects cellI to
                                    // pointCellI
                                    const label curImEdgeID =
                                        solidMesh.findImaginaryEdge
                                        (
                                        cellI,
                                        sePointCells[pointCellI]
                                        );

                                    if ( curImEdgeID == -1 )
                                    {
                                        FatalErrorIn
                                        (
                                            "insertTangential.H"
                                        ) << "curImEdge not found"
                                          << abort( FatalError );
                                    }

                                    // Contributions to upper/lower
                                    if ( cellI > sePointCells[pointCellI] )
                                    {
                                        l[curImEdgeID] += coeff;
                                    }
                                    else
                                    {
                                        u[curImEdgeID] += coeff;
                                    }
                                }
                            } // forAll sePointCells
                        } // for start and end of edge
                    } // internal edge
                } // forAll edges of current face
            } // if internal face
        } // forAll faces of current cell
    } // forAll cells


    // Boundary conditions
    // Note: currently there is no non-orthogonal correction on fixedValue or
    // symmetryPlane; solidTraction does not need any.
    Info << "Adding boundary condition contributions" << endl;
    forAll( mesh.boundary(), patchI )
    {
        if
        (
            U.boundaryField()[patchI].type()
            == solidTractionFvPatchVectorField::typeName
        )
        {
            const solidTractionFvPatchVectorField & tracPatch =
                refCast<solidTractionFvPatchVectorField>
                (
                U.boundaryField()[patchI]
                );

            const unallocLabelList & faceCells =
                mesh.boundary()[patchI].faceCells();

            const vectorField pForce =
                (
                tracPatch.traction()
                - n.boundaryField()[patchI] * tracPatch.pressure()
                ) * mesh.magSf().boundaryField()[patchI];

            forAll( mesh.boundary()[patchI], faceI )
            {
                const label own = faceCells[faceI];

                // Add force to source
                blockB[own] -= pForce[faceI];
            }
        }
        else
        if
        (
            U.boundaryField()[patchI].type()
            == fixedValueFvPatchVectorField::typeName
        )
        {
            const vectorField & pU = U.boundaryField()[patchI];
            const vectorField & pN = n.boundaryField()[patchI];
            const scalarField & pMu = muf.boundaryField()[patchI];
            const scalarField & pLambda = lambdaf.boundaryField()[patchI];
            const scalarField & pDeltaCoeffs =
                mesh.deltaCoeffs().boundaryField()[patchI];
            const scalarField & pMagSf =
                mesh.magSf().boundaryField()[patchI];
            const unallocLabelList & faceCells =
                mesh.boundary()[patchI].faceCells();

            forAll( mesh.boundary()[patchI], faceI )
            {
                const label own = faceCells[faceI];

                const Foam::vector & faceN = pN[faceI];
                const scalar faceMu = pMu[faceI];
                const scalar faceLambda = pLambda[faceI];
                const scalar faceMagSf = pMagSf[faceI];
                const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

                const tensor coeff =
                    sqr( faceN ) *
                    (2 * faceMu + faceLambda) * faceMagSf * faceDeltaCoeff
                    + ( I - sqr( faceN ) ) * faceMu * faceMagSf * faceDeltaCoeff;

                // Diag contribution
                d[own] -= coeff;

                // Source contribution
                blockB[own] -= coeff & pU[faceI];
            }
        }
        else
        if ( U.boundaryField()[patchI].type() == "symmetryPlane" )
        {
            const vectorField & pN = n.boundaryField()[patchI];
            const scalarField & pMu = muf.boundaryField()[patchI];
            const scalarField & pLambda = lambdaf.boundaryField()[patchI];
            const scalarField & pDeltaCoeffs =
                mesh.deltaCoeffs().boundaryField()[patchI];
            const scalarField & pMagSf =
                mesh.magSf().boundaryField()[patchI];
            const unallocLabelList & faceCells =
                mesh.boundary()[patchI].faceCells();

            forAll( mesh.boundary()[patchI], faceI )
            {
                const label own = faceCells[faceI];

                const Foam::vector & faceN = pN[faceI];
                const scalar faceMu = pMu[faceI];
                const scalar faceLambda = pLambda[faceI];
                const scalar faceMagSf = pMagSf[faceI];
                const scalar faceDeltaCoeff = pDeltaCoeffs[faceI];

                // Normal term
                d[own] -=
                    tensor
                    (
                    sqr( faceN ) * (2 * faceMu + faceLambda) * faceMagSf * faceDeltaCoeff
                    );

                // Tangential term
                const face & curFace =
                    faces[mesh.boundaryMesh()[patchI].start() + faceI];
                const edgeList curFaceEdges = curFace.edges();

                forAll( curFaceEdges, edgeI )
                {
                    const edge & curEdge = curFaceEdges[edgeI];

                    // Projected edge vector
                    Foam::vector e = curEdge.vec( points );
                    e -= faceN * (faceN & e);

                    // Edge length vector
                    Foam::vector Le = e ^ faceN;
                    Le *= curFace.edgeDirection( curEdge );

                    // Weights for pointCells
                    const label startPointID = curEdge.start();
                    const label endPointID = curEdge.end();
                    const scalarList & startCellWeights =
                        pointCellWeights[startPointID];
                    const scalarList & endCellWeights =
                        pointCellWeights[endPointID];

                    tensor LeFaceN = Le * faceN;

                    const labelList & startPointCells = pointCells[startPointID];
                    const labelList & endPointCells = pointCells[endPointID];

                    const labelList & curCellImagEdges = cellImagEdges[own];

                    // Start point cell contributions
                    forAll( startPointCells, pointCellI )
                    {
                        const scalar wCell = startCellWeights[pointCellI];

                        // Shear term
                        tensor coeff = 0.5 * faceLambda * wCell * LeFaceN.T();

                        if ( startPointCells[pointCellI] == own )
                        {
                            // Contributions to diag from start of edge
                            d[own] += coeff;
                        }
                        else // contribution to lower or upper
                        {
                            // Find which imaginary edge connects own to this
                            // pointCell
                            label curImEdgeLocalID = -1;
                            label imEdgeOwn = own;
                            label imEdgeNei = startPointCells[pointCellI];
                            bool flipImEdge = false;

                            if ( imEdgeOwn > imEdgeNei )
                            {
                                flipImEdge = true;
                            }

                            forAll( curCellImagEdges, imEdgeI )
                            {
                                const edge curImEdge =
                                    imagEdges[curCellImagEdges[imEdgeI]];

                                if
                                (
                                    (curImEdge.end() == imEdgeNei
                                        && !flipImEdge)
                                    || (curImEdge.start() == imEdgeNei
                                        && flipImEdge)
                                )
                                {
                                    curImEdgeLocalID = imEdgeI;
                                    break;
                                }
                            }

                            if ( curImEdgeLocalID == -1 )
                            {
                                FatalError
                                << "curImEdge not found"
                                << abort( FatalError );
                            }

                            const label curImEdgeID =
                                curCellImagEdges[curImEdgeLocalID];

                            // Contributions to upper/lower from start of edge
                            if ( flipImEdge )
                            {
                                l[curImEdgeID] += coeff;
                            }
                            else
                            {
                                u[curImEdgeID] += coeff;
                            }
                        }
                    } // forAll startPointCells

                    // End point cell contributions
                    forAll( endPointCells, pointCellI )
                    {
                        const scalar wCell = endCellWeights[pointCellI];

                        tensor coeff = 0.5 * faceLambda * wCell * LeFaceN.T();

                        if ( endPointCells[pointCellI] == own )
                        {
                            // Contributions to diag from end of edge
                            d[own] += coeff;
                        }
                        else
                        {
                            // Find which imag edge connects own to this
                            // pointCell
                            label curImEdgeLocalID = -1;
                            label imEdgeOwn = own;
                            label imEdgeNei = endPointCells[pointCellI];
                            bool flipImEdge = false;

                            if ( imEdgeOwn > imEdgeNei )
                            {
                                flipImEdge = true;
                            }

                            forAll( curCellImagEdges, imEdgeI )
                            {
                                const edge curImEdge =
                                    imagEdges[curCellImagEdges[imEdgeI]];

                                if
                                (
                                    (curImEdge.end() == imEdgeNei
                                        && !flipImEdge)
                                    || (curImEdge.start() == imEdgeNei
                                        && flipImEdge)
                                )
                                {
                                    curImEdgeLocalID = imEdgeI;
                                    break;
                                }
                            }

                            if ( curImEdgeLocalID == -1 )
                            {
                                FatalError
                                << "curImEdge not found"
                                << abort( FatalError );
                            }

                            const label curImEdgeID =
                                curCellImagEdges[curImEdgeLocalID];

                            // Contributions to upper/lower from end of edge
                            if ( flipImEdge )
                            {
                                l[curImEdgeID] += coeff;
                            }
                            else
                            {
                                u[curImEdgeID] += coeff;
                            }
                        }
                    } // forAll endPointCells
                } // forAll edges of current face
            } // forAll faces
        }
        else
        if
        (
            U.boundaryField()[patchI].type()
            == processorFvPatchVectorField::typeName
        )
        {
            Pout << "Patch " << patchI << " is a processor boundary" << endl;
        }
        else
        if
        (
            U.boundaryField()[patchI].type() != emptyFvPatchVectorField::typeName
        )
        {
            FatalError
            << "Boundary condition " << U.boundaryField()[patchI].type()
            << " not allowed" << nl
            << "Currently implemented boundary conditions:" << nl
            << fixedValueFvPatchVectorField::typeName << nl
            << solidTractionFvPatchVectorField::typeName << nl
            << "symmetryPlane" << nl
            << abort( FatalError );
        }
    }


    // - Block coupled solver call
    bool converged = false;
    int iCorr = 0;

    /*
     * Curently the GMRES solver is the only implemented linear
     * solver for non-diagonally dominant matrices
     * But GMRES convergence stalls unless it is restarted, so
     * we use a loop here to restart GMRES until convergence is
     * achieved
     */
    do
    {
        BlockSolverPerformance<Foam::vector> solverPerf =
            BlockLduSolver<Foam::vector>::New
            (
            U.name(),
            blockM,
            mesh.solutionDict().solver( "blockU" )
            )->solve( blockX, blockB );

        converged = solverPerf.converged();

        solverPerf.print();
    }
    while ( !converged && ++iCorr < nCorr );

    U.correctBoundaryConditions();
}
