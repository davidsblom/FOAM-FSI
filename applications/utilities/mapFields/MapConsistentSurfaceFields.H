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

#ifndef MapConsistentSurfaceFields_H
#define MapConsistentSurfaceFields_H

#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void MapConsistentSurfaceFields
(
    const IOobjectList& objects,
    const fvMesh& meshSource,
    const fvMesh& meshTarget
)
{
    word fieldClassName
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    for
    (
        IOobjectList::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        Info<< "    interpolating " << fieldIter()->name()
            << endl;

        // Read field
        GeometricField<Type, fvsPatchField, surfaceMesh> fieldSource
        (
            *fieldIter(),
            meshSource
        );

        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.headerOk())
        {
            // Read fieldTarget
            GeometricField<Type, fvsPatchField, surfaceMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            // Copy Field
            fieldTarget.internalField() = fieldSource.internalField();
            forAll(fieldTarget.boundaryField(),iBoundary){
                fieldTarget.boundaryField()[iBoundary] = fieldSource.boundaryField()[iBoundary];
            }
            /*meshToMeshInterp.interpolate
            (
                fieldTarget,
                fieldSource,
                meshToMesh::INTERPOLATE
            );*/



            // Write field
            fieldTarget.write();
        }
        else
        {
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;

            // Interpolate field
            GeometricField<Type, fvsPatchField, surfaceMesh> fieldTarget
            (
                fieldTargetIOobject,
                fieldSource
            );

            // Write field
            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
