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
    Foam::myMovingWallVelocity3FvPatchVectorField

Description
    Foam::myMovingWallVelocity3FvPatchVectorField

SourceFiles
    myMovingWallVelocity3FvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myMovingWallVelocity3FvPatchVectorField_H
#define myMovingWallVelocity3FvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class myMovingWallVelocity3FvPatch Declaration
\*---------------------------------------------------------------------------*/

class myMovingWallVelocity3FvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    //To keep track when to copy patch face centres
    label myTimeIndex_;

    //Old face centres of patch. Needed for (second order) time integration
    vectorField Fc_;
    vectorField Fc0_;
    vectorField Fc00_;
    vectorField Fc000_;

    // Private members

    scalarField bdf3Coeffs(const volVectorField& vf) const;


public:

    //- Runtime type information
    TypeName("myMovingWallVelocity3");


    // Constructors

        //- Construct from patch and internal field
        myMovingWallVelocity3FvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myMovingWallVelocity3FvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given myMovingWallVelocity3FvPatchVectorField
        //  onto a new patch
        myMovingWallVelocity3FvPatchVectorField
        (
            const myMovingWallVelocity3FvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myMovingWallVelocity3FvPatchVectorField
        (
            const myMovingWallVelocity3FvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myMovingWallVelocity3FvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myMovingWallVelocity3FvPatchVectorField
        (
            const myMovingWallVelocity3FvPatchVectorField&,
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
                new myMovingWallVelocity3FvPatchVectorField(*this, iF)
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
