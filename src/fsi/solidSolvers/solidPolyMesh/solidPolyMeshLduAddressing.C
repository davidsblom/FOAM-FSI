/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    solidPolyMeshLduAddressing

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "solidPolyMeshLduAddressing.H"
#include "globalMeshData.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPolyMeshLduAddressing::solidPolyMeshLduAddressing
(
    const solidPolyMesh& mesh
)
:
    lduAddressing(mesh.nCells()),
    lowerAddr_(),
    upperAddr_(),
    patchAddr_(0), //mesh.boundary().size()),
    patchSchedule_(mesh.globalData().patchSchedule())
{
    // Set upper and lower addressing using imaginary edges connectivity

    const edgeList& imEdges = mesh.imagEdges();

    upperAddr_.setSize(imEdges.size());
    lowerAddr_.setSize(imEdges.size());

    forAll(lowerAddr_, imEdgeI)
    {
        lowerAddr_[imEdgeI] = imEdges[imEdgeI].start();
        upperAddr_[imEdgeI] = imEdges[imEdgeI].end();
    }

    // Add dummy boundary addressing
    forAll (patchAddr_, patchI)
    {
        patchAddr_[patchI].setSize(0);
    }
}


// ************************************************************************* //
