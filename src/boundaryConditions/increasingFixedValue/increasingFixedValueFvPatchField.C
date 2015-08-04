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

#include "increasingFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar increasingFixedValueFvPatchField<Type>::currentScale() const
{
	scalar currentScale = 0.0;
    if(this->db().time().value()<time0_)
    {
        currentScale = 0.0;
    }
	else if(this->db().time().value()<(time0_+dT_))
	{
		currentScale = 	0.5*(1.0 - cos(2*mathematicalConstant::pi*(1/dT_/2)*(this->db().time().value()-time0_) ));
	}
    else
    {
        currentScale = 1.0;
    }
    return currentScale;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
increasingFixedValueFvPatchField<Type>::increasingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    startValue_(p.size()),
    endValue_(p.size()),
    dT_(0.0),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{}


template<class Type>
increasingFixedValueFvPatchField<Type>::increasingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    startValue_("startValue", dict, p.size()),
    endValue_("endValue",dict, p.size()),
    dT_(readScalar(dict.lookup("dT"))),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(startValue_ + (endValue_ - startValue_)*currentScale());
    }
}


template<class Type>
increasingFixedValueFvPatchField<Type>::increasingFixedValueFvPatchField
(
    const increasingFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    startValue_(ptf.startValue_, mapper),
    endValue_(ptf.endValue_, mapper),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{}


template<class Type>
increasingFixedValueFvPatchField<Type>::increasingFixedValueFvPatchField
(
    const increasingFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    startValue_(ptf.startValue_),
    endValue_(ptf.endValue_),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{}


template<class Type>
increasingFixedValueFvPatchField<Type>::increasingFixedValueFvPatchField
(
    const increasingFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    startValue_(ptf.startValue_),
    endValue_(ptf.endValue_),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void increasingFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    startValue_.autoMap(m);
}


template<class Type>
void increasingFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const increasingFixedValueFvPatchField<Type>& tiptf =
        refCast<const increasingFixedValueFvPatchField<Type> >(ptf);

    startValue_.rmap(tiptf.startValue_, addr);
}


template<class Type>
void increasingFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        patchField = startValue_ + (endValue_ - startValue_)*currentScale();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void increasingFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    startValue_.writeEntry("startValue", os);
    endValue_.writeEntry("endValue", os);
    os.writeKeyword("dT")
        << dT_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
