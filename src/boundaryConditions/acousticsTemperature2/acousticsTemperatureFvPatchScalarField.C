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

#include "acousticsTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::acousticsTemperatureFvPatchScalarField::acousticsTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::acousticsTemperatureFvPatchScalarField::acousticsTemperatureFvPatchScalarField
(
    const acousticsTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::acousticsTemperatureFvPatchScalarField::acousticsTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::acousticsTemperatureFvPatchScalarField::acousticsTemperatureFvPatchScalarField
(
    const acousticsTemperatureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::acousticsTemperatureFvPatchScalarField::acousticsTemperatureFvPatchScalarField
(
    const acousticsTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::acousticsTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::acousticsTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::acousticsTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarIOList& temperature = db().time().lookupObject<scalarIOList>("temperatureAcoustics");

    scalarField T ( patch().Cf().size(), scalar(0) );
    T = temperature;

    operator==
    (
        this->patchInternalField() + T/this->patch().deltaCoeffs()
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::acousticsTemperatureFvPatchScalarField::snGrad() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


Foam::tmp<Foam::scalarField>
Foam::acousticsTemperatureFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


Foam::tmp<Foam::scalarField>
Foam::acousticsTemperatureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


void Foam::acousticsTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        acousticsTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
