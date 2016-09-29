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

\*---------------------------------------------------------------------------*/

#include "acousticsVelocityGradientFvPatchVectorField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::acousticsVelocityGradientFvPatchVectorField::acousticsVelocityGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF)
{}


Foam::acousticsVelocityGradientFvPatchVectorField::acousticsVelocityGradientFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary&
)
:
    fixedGradientFvPatchVectorField(p, iF)
{

        fvPatchField<vector>::operator=(patchInternalField());
        gradient() = Foam::vector::zero;

}


Foam::acousticsVelocityGradientFvPatchVectorField::acousticsVelocityGradientFvPatchVectorField
(
    const acousticsVelocityGradientFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::acousticsVelocityGradientFvPatchVectorField::acousticsVelocityGradientFvPatchVectorField
(
    const acousticsVelocityGradientFvPatchVectorField& ptf
)
:
    fixedGradientFvPatchVectorField(ptf)
{}


Foam::acousticsVelocityGradientFvPatchVectorField::acousticsVelocityGradientFvPatchVectorField
(
    const acousticsVelocityGradientFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::acousticsVelocityGradientFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const vectorIOList& velocity = db().time().lookupObject<vectorIOList>("velocityAcoustics");

    gradient() = velocity;

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void Foam::acousticsVelocityGradientFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        acousticsVelocityGradientFvPatchVectorField
    );
}

// ************************************************************************* //
