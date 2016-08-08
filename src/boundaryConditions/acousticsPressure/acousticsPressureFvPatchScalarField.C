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

#include "acousticsPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::acousticsPressureFvPatchScalarField::acousticsPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    psiName_("psi"),
    gamma_(0.0),
    T0_(p.size(), 0.0)
{}


Foam::acousticsPressureFvPatchScalarField::acousticsPressureFvPatchScalarField
(
    const acousticsPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    T0_(ptf.T0_, mapper)
{}


Foam::acousticsPressureFvPatchScalarField::acousticsPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    gamma_(readScalar(dict.lookup("gamma"))),
    T0_("T0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(T0_);
    }
}


Foam::acousticsPressureFvPatchScalarField::acousticsPressureFvPatchScalarField
(
    const acousticsPressureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_)
{}


Foam::acousticsPressureFvPatchScalarField::acousticsPressureFvPatchScalarField
(
    const acousticsPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::acousticsPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    T0_.autoMap(m);
}


void Foam::acousticsPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const acousticsPressureFvPatchScalarField& tiptf =
        refCast<const acousticsPressureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::acousticsPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarIOList& pressure = db().lookupObject<scalarIOList>("pressureAcoustics");

    scalarField p ( patch().Cf().size(), scalar(0) );
    p = pressure;

    operator==
    (
        p
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::tmp<Foam::scalarField>
Foam::acousticsPressureFvPatchScalarField::snGrad() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


Foam::tmp<Foam::scalarField>
Foam::acousticsPressureFvPatchScalarField::gradientInternalCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


Foam::tmp<Foam::scalarField>
Foam::acousticsPressureFvPatchScalarField::gradientBoundaryCoeffs() const
{
    return tmp<scalarField>
    (
        new scalarField(this->size(), 0.0)
    );
}


void Foam::acousticsPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    if (UName_ != "U")
    {
        os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    }
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (psiName_ != "psi")
    {
        os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    T0_.writeEntry("T0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        acousticsPressureFvPatchScalarField
    );
}

// ************************************************************************* //
