
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef setInletVelocityFvPatchVectorField_H
#define setInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
    *             Class parabolicVelocityFvPatchField Declaration
    \*---------------------------------------------------------------------------*/

    class setInletVelocityFvPatchVectorField
        :
        public fixedValueFvPatchVectorField
    {
        // Private data

        // - Peak velocity magnitude
        scalar maxValue_;

public:

        // - Runtime type information
        TypeName( "setInletVelocity" );


        // Constructors

        // - Construct from patch and internal field
        setInletVelocityFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &
            );

        // - Construct from patch, internal field and dictionary
        setInletVelocityFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &
            );

        // - Construct by mapping given
        // setInletVelocityFvPatchVectorField
        // onto a new patch
        setInletVelocityFvPatchVectorField(
            const setInletVelocityFvPatchVectorField &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &
            );

        // - Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
                   (
                new setInletVelocityFvPatchVectorField( *this )
                   );
        }

        // - Construct as copy setting internal field reference
        setInletVelocityFvPatchVectorField(
            const setInletVelocityFvPatchVectorField &,
            const DimensionedField<vector, volMesh> &
            );

        // - Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
            ( const DimensionedField<vector, volMesh> & iF ) const
        {
            return tmp<fvPatchVectorField>
                   (
                new setInletVelocityFvPatchVectorField( *this, iF )
                   );
        }

        // Member functions

        // - Return max value
        scalar & maxValue()
        {
            return maxValue_;
        }

        // - Update coefficients
        virtual void updateCoeffs();

        // - Write
        virtual void write( Ostream & ) const;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
