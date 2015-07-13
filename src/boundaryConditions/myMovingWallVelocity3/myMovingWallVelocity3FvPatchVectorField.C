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

\*---------------------------------------------------------------------------*/

#include "myMovingWallVelocity3FvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * * * //

scalarField myMovingWallVelocity3FvPatchVectorField::bdf3Coeffs(const volVectorField& vf) const
{
    scalarField coeffs(4,0.0);

    //if first time step
    if (vf.oldTime().timeIndex() == vf.oldTime().oldTime().timeIndex() || vf.oldTime().oldTime().timeIndex() < 0)
    {
        coeffs[0]=1.0;
        coeffs[1]=-1.0;
    }
    else if(vf.oldTime().oldTime().timeIndex() == vf.oldTime().oldTime().oldTime().timeIndex() || vf.oldTime().oldTime().oldTime().timeIndex() < 0)
    {
        coeffs[0]=3.0/2.0;
        coeffs[1]=-2.0;
        coeffs[2]=1.0/2.0;
    }
    else
    {
        coeffs[0]=11.0/6.0;
        coeffs[1]=-3.0;
        coeffs[2]=3.0/2.0;
        coeffs[3]=-1.0/3.0;
    }

    return coeffs;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myMovingWallVelocity3FvPatchVectorField::myMovingWallVelocity3FvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    Fc0_(p.patch().size(),vector::zero),
    Fc00_(p.patch().size(),vector::zero),
    Fc000_(p.patch().size(),vector::zero)
{}


myMovingWallVelocity3FvPatchVectorField::myMovingWallVelocity3FvPatchVectorField
(
    const myMovingWallVelocity3FvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    myTimeIndex_(ptf.myTimeIndex_),
    Fc_(p.patch().size(),vector::zero),
    Fc0_(p.patch().size(),vector::zero),
    Fc00_(p.patch().size(),vector::zero),
    Fc000_(p.patch().size(),vector::zero)
{}


myMovingWallVelocity3FvPatchVectorField::myMovingWallVelocity3FvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    Fc0_(p.patch().size(),vector::zero),
    Fc00_(p.patch().size(),vector::zero),
    Fc000_(p.patch().size(),vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    Fc_ = p.patch().faceCentres();
    Fc0_ = p.patch().faceCentres();
    Fc00_ = p.patch().faceCentres();
    Fc000_ = p.patch().faceCentres();
}


myMovingWallVelocity3FvPatchVectorField::myMovingWallVelocity3FvPatchVectorField
(
    const myMovingWallVelocity3FvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc_),
    Fc0_(pivpvf.Fc0_),
    Fc00_(pivpvf.Fc00_),
    Fc000_(pivpvf.Fc000_)
{}


myMovingWallVelocity3FvPatchVectorField::myMovingWallVelocity3FvPatchVectorField
(
    const myMovingWallVelocity3FvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc0_),
    Fc0_(pivpvf.Fc0_),
    Fc00_(pivpvf.Fc00_),
    Fc000_(pivpvf.Fc000_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myMovingWallVelocity3FvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Get current mesh data
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();

    //Update the face center locations in time
    if(myTimeIndex_< mesh.time().timeIndex())
    {
        Fc000_ = Fc00_;
        Fc00_ = Fc0_;
        Fc0_ = Fc_;
        Fc_ = pp.faceCentres();
        myTimeIndex_ = mesh.time().timeIndex();
    }

    const volVectorField& U = mesh.lookupObject<volVectorField>(dimensionedInternalField().name());

    scalarField coeffs = bdf3Coeffs(U);

    vectorField Up = (coeffs[0]*Fc_ + coeffs[1]*Fc0_ + coeffs[2]*Fc00_ + coeffs[3]*Fc000_)/mesh.time().deltaT().value();

    scalarField phip = p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    vectorField::operator=(Up + n*(Un - (n & Up)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void myMovingWallVelocity3FvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    myMovingWallVelocity3FvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
