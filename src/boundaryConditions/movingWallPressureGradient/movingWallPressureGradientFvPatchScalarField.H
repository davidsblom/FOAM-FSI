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
    Foam::movingWallPressureGradientFvPatchScalarField

Description
    Foam::movingWallPressureGradientFvPatchScalarField

SourceFiles
    movingWallPressureGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallPressureGradientFvPatchScalarFields_H
#define movingWallPressureGradientFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class movingWallPressureGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

class movingWallPressureGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        scalar deltaT_() const;

        scalar deltaT0_() const;

        scalar deltaT0_( const volVectorField& vf ) const;

        vectorField ddtU(const volVectorField& U) const;

        vectorField ddtU(const volVectorField& U, const surfaceVectorField& Uf) const;

public:

    //- Runtime type information
    TypeName("movingWallPressureGradient");


    // Constructors

        //- Construct from patch and internal field
        movingWallPressureGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallPressureGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallPressureGradientFvPatchScalarField onto
        //  a new patch
        movingWallPressureGradientFvPatchScalarField
        (
            const movingWallPressureGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        movingWallPressureGradientFvPatchScalarField
        (
            const movingWallPressureGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        movingWallPressureGradientFvPatchScalarField
        (
            const movingWallPressureGradientFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new movingWallPressureGradientFvPatchScalarField(*this, iF)
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
