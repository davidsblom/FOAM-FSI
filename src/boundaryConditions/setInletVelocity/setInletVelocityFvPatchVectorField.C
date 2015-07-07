
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "setInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    setInletVelocityFvPatchVectorField::
    setInletVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        maxValue_( 0 )
    {}


    setInletVelocityFvPatchVectorField::
    setInletVelocityFvPatchVectorField(
        const setInletVelocityFvPatchVectorField & ptf,
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const fvPatchFieldMapper & mapper
        )
        :
        fixedValueFvPatchVectorField( ptf, p, iF, mapper ),
        maxValue_( ptf.maxValue_ )
    {}


    setInletVelocityFvPatchVectorField::
    setInletVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const dictionary & dict
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        maxValue_( readScalar( dict.lookup( "maxValue" ) ) )
    {
        evaluate();
    }

    setInletVelocityFvPatchVectorField::
    setInletVelocityFvPatchVectorField(
        const setInletVelocityFvPatchVectorField & fcvpvf,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( fcvpvf, iF ),
        maxValue_( fcvpvf.maxValue_ )
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void setInletVelocityFvPatchVectorField::updateCoeffs()
    {
        if ( updated() )
        {
            return;
        }

        scalar Umax = maxValue_;

        scalar t = this->db().time().value();

        scalar T = 2; // 0.2;

        if ( t < T )
        {
            Umax = maxValue_ * ( 1.0 - cos( M_PI * t / T ) ) / 2.0;
        }

        const vectorField & Cf = patch().Cf();

        scalarField y( Cf.size() );
        scalarField z( Cf.size() );
        forAll( Cf, i )
        {
            y[i] = Cf[i][1];
            z[i] = Cf[i][2];
        }

        vectorField::operator=( Umax * y * (0.4 - y) * ( sqr( 0.4 ) - sqr( z ) ) * vector( 1, 0, 0 ) / ( sqr( 0.2 ) * sqr( 0.4 ) ) );
    }

    // Write
    void setInletVelocityFvPatchVectorField::write( Ostream & os ) const
    {
        fvPatchVectorField::write( os );
        os.writeKeyword( "maxValue" )
        << maxValue_ << token::END_STATEMENT << nl;
        writeEntry( "value", os );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
        fvPatchVectorField,
        setInletVelocityFvPatchVectorField
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
