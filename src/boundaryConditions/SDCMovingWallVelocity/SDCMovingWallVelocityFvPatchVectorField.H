
#ifndef SDCMovingWallVelocityFvPatchVectorField_H
#define SDCMovingWallVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

#include <Eigen/Dense>
#include "SDC.H"
#include "ESDIRK.H"
#include "TimeIntegrationScheme.H"
#include <yaml-cpp/yaml.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
    *                Class SDCMovingWallVelocityFvPatch Declaration
    \*---------------------------------------------------------------------------*/

    class SDCMovingWallVelocityFvPatchVectorField
        :
        public fixedValueFvPatchVectorField
    {
public:

        // - Runtime type information
        TypeName( "SDCMovingWallVelocity" );


        // Constructors

        // - Construct from patch and internal field
        SDCMovingWallVelocityFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &
            );

        // - Construct from patch, internal field and dictionary
        SDCMovingWallVelocityFvPatchVectorField(
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const dictionary &
            );

        // - Construct by mapping given SDCMovingWallVelocityFvPatchVectorField
        // onto a new patch
        SDCMovingWallVelocityFvPatchVectorField(
            const SDCMovingWallVelocityFvPatchVectorField &,
            const fvPatch &,
            const DimensionedField<vector, volMesh> &,
            const fvPatchFieldMapper &
            );

        // - Construct as copy
        SDCMovingWallVelocityFvPatchVectorField( const SDCMovingWallVelocityFvPatchVectorField & );

        // - Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
                   (
                new SDCMovingWallVelocityFvPatchVectorField( *this )
                   );
        }

        // - Construct as copy setting internal field reference
        SDCMovingWallVelocityFvPatchVectorField(
            const SDCMovingWallVelocityFvPatchVectorField &,
            const DimensionedField<vector, volMesh> &
            );

        // - Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
            ( const DimensionedField<vector, volMesh> & iF ) const
        {
            return tmp<fvPatchVectorField>
                   (
                new SDCMovingWallVelocityFvPatchVectorField( *this, iF )
                   );
        }

        // Member functions

        // - Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // - Write
        virtual void write( Ostream & ) const;

        void setSDCInfo(
            bool corrector,
            int k,
            int sweep
            );

        std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme;

private:

        void init();

        int k;
        int sweep;
        bool corrector;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
