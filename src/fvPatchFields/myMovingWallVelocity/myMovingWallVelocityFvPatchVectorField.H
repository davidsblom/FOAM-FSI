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

Class
    Foam::myMovingWallVelocityFvPatchVectorField

Description
    Foam::myMovingWallVelocityFvPatchVectorField

SourceFiles
    myMovingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myMovingWallVelocityFvPatchVectorField_H
#define myMovingWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class myMovingWallVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class myMovingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    //To keep track when to copy patch face centres
    label myTimeIndex_;

    //Old face centres of patch. Needed for (second order) time integration
    vectorField Fc_;
    vectorField oldFc_;
    vectorField oldoldFc_;


public:

    //- Runtime type information
    TypeName("myMovingWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        myMovingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myMovingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given myMovingWallVelocityFvPatchVectorField
        //  onto a new patch
        myMovingWallVelocityFvPatchVectorField
        (
            const myMovingWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myMovingWallVelocityFvPatchVectorField
        (
            const myMovingWallVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myMovingWallVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myMovingWallVelocityFvPatchVectorField
        (
            const myMovingWallVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new myMovingWallVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
