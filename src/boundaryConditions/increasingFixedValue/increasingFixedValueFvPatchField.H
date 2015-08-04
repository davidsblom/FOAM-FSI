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

Class
    Foam::increasingFixedValueFvPatchField

Description
    Foam::increasingFixedValueFvPatchField

SourceFiles
    increasingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef increasingFixedValueFvPatchField_H
#define increasingFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class increasingFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class increasingFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Reference value
        Field<Type> startValue_;

        //- Amplitude
        Field<Type> endValue_;

        //- time to reach endValue
        scalar dT_;

        //- Current time index
        label curTimeIndex_;

        //- Time when constructed: increase will hapen from time0_ to time0_+dT_
        scalar time0_;        


    // Private member functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("increasingFixedValue");


    // Constructors

        //- Construct from patch and internal field
        increasingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        increasingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given increasingFixedValueFvPatchField
        //  onto a new patch
        increasingFixedValueFvPatchField
        (
            const increasingFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        increasingFixedValueFvPatchField
        (
            const increasingFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new increasingFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        increasingFixedValueFvPatchField
        (
            const increasingFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new increasingFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& startValue() const
            {
                return startValue_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& startValue()
            {
                return startValue_;
            }

            //- Return amplitude
            scalar endValue() const
            {
                return endValue_;
            }

            scalar& endValue()
            {
                return endValue_;
            }

            //- Return frequency
            scalar dT() const
            {
                return dT_;
            }

            scalar& dT()
            {
                return dT_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "increasingFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
