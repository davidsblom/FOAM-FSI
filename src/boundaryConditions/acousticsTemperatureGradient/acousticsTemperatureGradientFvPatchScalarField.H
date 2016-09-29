/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::acousticsTemperatureGradientFvPatchScalarField

Description
    Foam::acousticsTemperatureGradientFvPatchScalarField

SourceFiles
    acousticsTemperatureGradientFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef acousticsTemperatureGradientFvPatchScalarFields_H
#define acousticsTemperatureGradientFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class acousticsTemperatureGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

class acousticsTemperatureGradientFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data


public:

    //- Runtime type information
    TypeName("acousticsTemperatureGradient");


    // Constructors

        //- Construct from patch and internal field
        acousticsTemperatureGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        acousticsTemperatureGradientFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given acousticsTemperatureGradientFvPatchScalarField onto
        //  a new patch
        acousticsTemperatureGradientFvPatchScalarField
        (
            const acousticsTemperatureGradientFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        acousticsTemperatureGradientFvPatchScalarField
        (
            const acousticsTemperatureGradientFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new acousticsTemperatureGradientFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        acousticsTemperatureGradientFvPatchScalarField
        (
            const acousticsTemperatureGradientFvPatchScalarField&,
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
                new acousticsTemperatureGradientFvPatchScalarField(*this, iF)
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
