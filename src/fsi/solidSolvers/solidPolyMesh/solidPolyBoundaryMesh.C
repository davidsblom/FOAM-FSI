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

\*---------------------------------------------------------------------------*/

#include "solidPolyBoundaryMesh.H"
#include "polyBoundaryMesh.H"
//#include "faceTetPolyPatch.H"
//#include "globalTetPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyBoundaryMesh
solidPolyBoundaryMesh::solidPolyBoundaryMesh
(
    const solidPolyMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    //solidPolyPatchList(basicBdry.size()),
    mesh_(m),
    polyBdry_(basicBdry)
{
    if (Pstream::parRun())
    {
        //FatalErrorIn
        WarningIn
        (
            "solidPolyBoundaryMesh::solidPolyBoundaryMesh"
        )   << " lduInterface not implemented yet for coupled solid solver"
            //<< abort(FatalError);
            << endl;
    }

    // PC: Do we need to create special polyPatches...? or can we just use
    // polyBoundaryMesh functionality...

    // Set boundary patches
    //solidPolyPatchList& Patches = *this;

    // forAll(Patches, patchI)
    // {
    //     Patches.set
    //     (
    //         patchI,
    //         faceTetPolyPatch::New(basicBdry[patchI], *this).ptr()
    //     );
    // }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

lduInterfacePtrsList solidPolyBoundaryMesh::interfaces() const
{
    //lduInterfacePtrsList interfaces(size());
    //lduInterfacePtrsList interfaces(0);
    lduInterfacePtrsList interfaces(polyBdry_.size());

    Warning
        << "solidPolyBoundaryMesh::solidPolyBoundaryMesh"
        << " lduInterface size is 0 for now" << endl;

    forAll (interfaces, patchi)
    {
        //if (isA<lduInterface>(this->operator[](patchi)))
        if (isA<lduInterface>(polyBdry_[patchi]))
        {
            interfaces.set
            (
                patchi,
                //&refCast<const lduInterface>(this->operator[](patchi))
                &refCast<const lduInterface>(polyBdry_[patchi])
            );
        }
    }

    return interfaces;
}


// const globalTetPolyPatch&
// solidPolyBoundaryMesh::globalPatch() const
// {
//     const solidPolyPatchList& patches = *this;

//     forAll (patches, patchI)
//     {
//         if (isA<globalTetPolyPatch>(patches[patchI]))
//         {
//             return refCast<const globalTetPolyPatch>
//             (
//                 patches[patchI]
//             );
//         }
//     }

//     FatalErrorIn
//     (
//         "const globalTetPolyPatch&"
//         "solidPolyBoundaryMesh::globalPatch() const"
//     )   << "patch not found.  Is this case running in parallel?"
//         << abort(FatalError);

//     // Dummy return
//     return refCast<const globalTetPolyPatch>(patches[0]);
// }


// faceListList solidPolyBoundaryMesh::boundaryTriFaces() const
// {
//     faceListList result(size());

//     forAll (result, patchI)
//     {
//         result[patchI] = operator[](patchI).triFaces();
//     }

//     return result;
// }


// void solidPolyBoundaryMesh::updateMesh()
// {
//     solidPolyPatchList& Patches = *this;

//     forAll(Patches, patchI)
//     {
//         Patches[patchI].updateMesh();
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
