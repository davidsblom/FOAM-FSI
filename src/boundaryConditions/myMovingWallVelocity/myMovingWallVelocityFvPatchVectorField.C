
#include "myMovingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    myMovingWallVelocityFvPatchVectorField::myMovingWallVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        myTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
        Fc_( p.patch().size(), vector::zero ),
        oldFc_( p.patch().size(), vector::zero ),
        oldoldFc_( p.patch().size(), vector::zero )
    {}


    myMovingWallVelocityFvPatchVectorField::myMovingWallVelocityFvPatchVectorField(
        const myMovingWallVelocityFvPatchVectorField & ptf,
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const fvPatchFieldMapper & mapper
        )
        :
        fixedValueFvPatchVectorField( ptf, p, iF, mapper ),
        myTimeIndex_( ptf.myTimeIndex_ ),
        Fc_( p.patch().size(), vector::zero ),
        oldFc_( p.patch().size(), vector::zero ),
        oldoldFc_( p.patch().size(), vector::zero )
    {}


    myMovingWallVelocityFvPatchVectorField::myMovingWallVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const dictionary & dict
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        myTimeIndex_( dimensionedInternalField().mesh().time().timeIndex() ),
        Fc_( p.patch().size(), vector::zero ),
        oldFc_( p.patch().size(), vector::zero ),
        oldoldFc_( p.patch().size(), vector::zero )
    {
        fvPatchVectorField::operator=( vectorField( "value", dict, p.size() ) );

        Fc_ = p.patch().faceCentres();
        oldFc_ = p.patch().faceCentres();
        oldoldFc_ = p.patch().faceCentres();
    }

    myMovingWallVelocityFvPatchVectorField::myMovingWallVelocityFvPatchVectorField( const myMovingWallVelocityFvPatchVectorField & pivpvf )
        :
        fixedValueFvPatchVectorField( pivpvf ),
        myTimeIndex_( pivpvf.myTimeIndex_ ),
        Fc_( pivpvf.Fc_ ),
        oldFc_( pivpvf.oldFc_ ),
        oldoldFc_( pivpvf.oldoldFc_ )
    {}


    myMovingWallVelocityFvPatchVectorField::myMovingWallVelocityFvPatchVectorField(
        const myMovingWallVelocityFvPatchVectorField & pivpvf,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( pivpvf, iF ),
        myTimeIndex_( pivpvf.myTimeIndex_ ),
        Fc_( pivpvf.oldFc_ ),
        oldFc_( pivpvf.oldFc_ ),
        oldoldFc_( pivpvf.oldoldFc_ )
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void myMovingWallVelocityFvPatchVectorField::updateCoeffs()
    {
        if ( updated() )
        {
            return;
        }

        // Info << "void myMovingWallVelocityFvPatchVectorField::updateCoeffs()" <<endl;
        const fvMesh & mesh = dimensionedInternalField().mesh();
        const fvPatch & p = patch();
        const polyPatch & pp = p.patch();

        if ( myTimeIndex_ < mesh.time().timeIndex() )
        {
            oldoldFc_ = oldFc_;
            oldFc_ = Fc_;
            Fc_ = pp.faceCentres();
            myTimeIndex_ = mesh.time().timeIndex();
        }

        // const pointField& oldPoints = mesh.oldPoints();
        const volVectorField & U = mesh.lookupObject<volVectorField>( dimensionedInternalField().name() );

        scalar deltaT = mesh.time().deltaT().value();
        scalar deltaT0 = mesh.time().deltaT0().value();

        if ( U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex() || U.oldTime().oldTime().timeIndex() < 0 )
        {
            deltaT0 = GREAT;
        }

        // Set coefficients based on deltaT and deltaT0
        scalar coefft = 1 + deltaT / (deltaT + deltaT0);
        scalar coefft00 = deltaT * deltaT / ( deltaT0 * (deltaT + deltaT0) );
        scalar coefft0 = coefft + coefft00;

        vectorField Up = (coefft * Fc_ - coefft0 * oldFc_ + coefft00 * oldoldFc_) / mesh.time().deltaT().value();

        scalarField phip = p.patchField<surfaceScalarField, scalar>( fvc::meshPhi( U ) );

        vectorField n = p.nf();
        const scalarField & magSf = p.magSf();
        scalarField Un = phip / (magSf + VSMALL);

        vectorField::operator=( Up + n *( Un - (n & Up) ) );

        fixedValueFvPatchVectorField::updateCoeffs();
    }

    void myMovingWallVelocityFvPatchVectorField::write( Ostream & os ) const
    {
        fvPatchVectorField::write( os );
        writeEntry( "value", os );
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
        fvPatchVectorField,
        myMovingWallVelocityFvPatchVectorField
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
