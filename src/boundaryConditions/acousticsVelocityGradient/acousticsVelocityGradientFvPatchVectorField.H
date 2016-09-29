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
    Foam::acousticsVelocityGradientFvPatchVectorField

Description
    Foam::acousticsVelocityGradientFvPatchVectorField

SourceFiles
    acousticsVelocityGradientFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef acousticsVelocityGradientFvPatchVectorFields_H
#define acousticsVelocityGradientFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef IOList<vector> vectorIOList;

/*---------------------------------------------------------------------------*\
                   Class acousticsVelocityGradientFvPatch Declaration
\*---------------------------------------------------------------------------*/

class acousticsVelocityGradientFvPatchVectorField
:
    public fixedGradientFvPatchField<vector>
{
    // Private data


public:

    //- Runtime type information
    TypeName("acousticsVelocityGradient");


    // Constructors

        //- Construct from patch and internal field
        acousticsVelocityGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        acousticsVelocityGradientFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given acousticsVelocityGradientFvPatchVectorField onto
        //  a new patch
        acousticsVelocityGradientFvPatchVectorField
        (
            const acousticsVelocityGradientFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        acousticsVelocityGradientFvPatchVectorField
        (
            const acousticsVelocityGradientFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new acousticsVelocityGradientFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        acousticsVelocityGradientFvPatchVectorField
        (
            const acousticsVelocityGradientFvPatchVectorField&,
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
                new acousticsVelocityGradientFvPatchVectorField(*this, iF)
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
