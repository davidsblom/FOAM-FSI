/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Class applies a two-dimensional correction to mesh motion point field.
    The correction guarantees that the mesh does not get twisted during motion
    and thus introduce a third dimension into a 2-D problem.

\*---------------------------------------------------------------------------*/

#include "twoDPointCorrectorRBF.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Set marker to -1 if is on the "back" side
void twoDPointCorrectorRBF::setMarker()
{
    if (!required()) return;

    // Change size of useablePointIDs_ and shadowPointIDs_
    useablePointIDs_.setSize(mesh_.nPoints()/2);
    shadowPointIDs_.setSize(mesh_.nPoints()/2);

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();
    const pointField& points(mesh_.points());

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    forAll (neIndices, edgeI)
    {
        const label& pStartInd = meshEdges[neIndices[edgeI]].start();
        const label& pEndInd = meshEdges[neIndices[edgeI]].end();
        const point& pStart = points[pStartInd];
        const point& pEnd = points[pEndInd];

        // calculate average point position
        const point A = 0.5*(pStart + pEnd);

        // Calculate inner product with plane normal
        scalar pStartInner = (pn & (pStart - A));
        scalar pEndInner = (pn & (pEnd - A));
        if(pStartInner > 0 && pEndInner < 0){
            pointMarker_[pEndInd] = -1;
            useablePointIDs_[edgeI] = pStartInd;
            shadowPointIDs_[edgeI] = pEndInd;
        }else if(pEndInner > 0 && pStartInner < 0)
        {
            pointMarker_[pStartInd] = -1;
            useablePointIDs_[edgeI] = pEndInd;
            shadowPointIDs_[edgeI] = pStartInd;
        }else
        {
            FatalErrorIn("void twoDPointCorrectorRBF::setMarker()")
            << "Both points give back a negative value with the inner product. Programming error?"
            << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDPointCorrectorRBF::twoDPointCorrectorRBF(const polyMesh& mesh)
:
    twoDPointCorrector(mesh),
    mesh_(mesh),
    pointMarker_(mesh.nPoints(),0),
    useablePointIDs_(mesh.nPoints(),0),
    shadowPointIDs_(mesh.nPoints(),0)
{
    setMarker();
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoDPointCorrectorRBF::~twoDPointCorrectorRBF()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& twoDPointCorrectorRBF::marker() const{
    return pointMarker_;
}

void twoDPointCorrectorRBF::setShadowSide(vectorField& newpoints) const{
    forAll(useablePointIDs_,ipoint){
        newpoints[shadowPointIDs_[ipoint]] = newpoints[useablePointIDs_[ipoint]];
    }
}

/*void twoDPointCorrectorRBF::correctPoints(pointField& p) const
{
    if (!required_) return;

    // Algorithm:
    // Loop through all edges. Calculate the average point position A for
    // the front and the back. Correct the position of point P (in two planes)
    // such that vectors AP and planeNormal are parallel

    // Get reference to edges
    const edgeList&  meshEdges = mesh_.edges();

    const labelList& neIndices = normalEdgeIndices();
    const vector& pn = planeNormal();

    forAll (neIndices, edgeI)
    {
        point& pStart = p[meshEdges[neIndices[edgeI]].start()];

        point& pEnd = p[meshEdges[neIndices[edgeI]].end()];

        // calculate average point position
        const point A = 0.5*(pStart + pEnd);

        // correct point locations
        pStart = A + pn*(pn & (pStart - A));
        pEnd = A + pn*(pn & (pEnd - A));
    }
}*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
