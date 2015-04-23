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

#include "rotatingWallVelocitySmoothFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Private Members  * * * * * * * * * * * * * * //

scalar rotatingWallVelocitySmoothFvPatchVectorField::getOmega() const
{
    const scalar& time = dimensionedInternalField().mesh().time().value();
    if(T0_ > time && T0_>SMALL){
        scalar smoothFactor = 0.5*(1-Foam::cos(mathematicalConstant::pi/T0_*time));
        return smoothFactor*omega_;
    }else{
        return omega_;
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatingWallVelocitySmoothFvPatchVectorField::rotatingWallVelocitySmoothFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(vector::zero),
    axis_(vector::zero),
    omega_(0),
    T0_(0)
{}


rotatingWallVelocitySmoothFvPatchVectorField::rotatingWallVelocitySmoothFvPatchVectorField
(
    const rotatingWallVelocitySmoothFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_),
    T0_(ptf.T0_)
{}


rotatingWallVelocitySmoothFvPatchVectorField::rotatingWallVelocitySmoothFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(readScalar(dict.lookup("omega"))),
    T0_(readScalar(dict.lookup("T0")))
{
    // Evaluate the wall velocity
    updateCoeffs();
}


rotatingWallVelocitySmoothFvPatchVectorField::rotatingWallVelocitySmoothFvPatchVectorField
(
    const rotatingWallVelocitySmoothFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchField<vector>(pivpvf),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_),
    T0_(pivpvf.T0_)
{}


rotatingWallVelocitySmoothFvPatchVectorField::rotatingWallVelocitySmoothFvPatchVectorField
(
    const rotatingWallVelocitySmoothFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(pivpvf, iF),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_),
    T0_(pivpvf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rotatingWallVelocitySmoothFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar omega = getOmega();

    // Calculate the rotating wall velocity from the specification of the motion
    vectorField Up = (-omega)*((patch().Cf() - origin_) ^ (axis_/mag(axis_)));

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    vectorField n = patch().nf();
    vectorField::operator=(Up - n*(n & Up));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void rotatingWallVelocitySmoothFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    os.writeKeyword("omega") << omega_ << token::END_STATEMENT << nl;
    os.writeKeyword("T0") << T0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    rotatingWallVelocitySmoothFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
