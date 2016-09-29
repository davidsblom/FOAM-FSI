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

#include "acousticsTemperatureGradientFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::acousticsTemperatureGradientFvPatchScalarField::acousticsTemperatureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::acousticsTemperatureGradientFvPatchScalarField::acousticsTemperatureGradientFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary&
)
:
    fixedGradientFvPatchScalarField(p, iF)
{

        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;

}


Foam::acousticsTemperatureGradientFvPatchScalarField::acousticsTemperatureGradientFvPatchScalarField
(
    const acousticsTemperatureGradientFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::acousticsTemperatureGradientFvPatchScalarField::acousticsTemperatureGradientFvPatchScalarField
(
    const acousticsTemperatureGradientFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf)
{}


Foam::acousticsTemperatureGradientFvPatchScalarField::acousticsTemperatureGradientFvPatchScalarField
(
    const acousticsTemperatureGradientFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::acousticsTemperatureGradientFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarIOList& temperature = db().time().lookupObject<scalarIOList>("temperatureAcoustics");

    gradient() = temperature;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::acousticsTemperatureGradientFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        acousticsTemperatureGradientFvPatchScalarField
    );
}

// ************************************************************************* //
