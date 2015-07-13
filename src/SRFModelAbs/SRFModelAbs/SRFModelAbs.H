
/*
 *
 * Namespace
 * Foam::SRF
 *
 * Description
 * Namespace for single rotating frame (SRF) models
 *
 * Class
 * Foam::SRF::SRFModelAbs
 *
 * Description
 * Top level model for single rotating frame
 * - Steady state only - no time derivatives included
 *
 * SourceFiles
 * SRFModelAbs
 *
 */

#ifndef SRFModelAbs_H
#define SRFModelAbs_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "volFields.H"
#include "vectorField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace SRF
    {
        /*---------------------------------------------------------------------------*\
        *                          Class SRFModelAbs Declaration
        \*---------------------------------------------------------------------------*/

        class SRFModelAbs
            :
            public IOdictionary
        {
protected:

            // Protected data

            // - Reference to the relative velocity field
            const volVectorField & U_;

            // - Reference to the mesh
            const fvMesh & mesh_;

            // - Axis of rotation
            vector axis_;

            // - Origin of rotation
            vector origin_;

            // - SRF model coeficients dictionary
            dictionary SRFModelAbsCoeffs_;

            // - Angular velocity of the frame (rad/s)
            dimensionedVector omega_;

            // - List of excluded patches (nonRotatingPatches)
            const wordList excludedPatchNames_;
            labelList excludedPatchLabels_;

            // - Patches that do move with SRF
            labelList includedPatches_;

            // - Excluded patches (nonRotatingPatches)
            labelList excludedPatches_;

private:

            // Private Member Functions

            void setPatchProperties();

            // - Disallow default bitwise copy construct
            SRFModelAbs( const SRFModelAbs & );

            // - Disallow default bitwise assignment
            void operator=( const SRFModelAbs & );

public:

            // - Runtime type information
            TypeName( "SRFModelAbs" );


            // Declare runtime constructor selection table

            declareRunTimeSelectionTable
            (
                autoPtr,
                SRFModelAbs,
                dictionary,
                (
                    const volVectorField &U
                ),
                (U)
            );


            // Constructors

            // - Construct from components
            SRFModelAbs(
                const word & type,
                const volVectorField & U
                );


            // Selectors

            // - Return a reference to the selected SRF model
            static autoPtr<SRFModelAbs> New( const volVectorField & U );


            // Destructor

            virtual ~SRFModelAbs();


            // Member Functions

            // Edit

            // - Read radiationProperties dictionary
            virtual bool read();


            // Access

            // - Return the axis of rotation
            const vector & axis() const;

            // - Return the angular velocity field [rad/s]
            const dimensionedVector & omega() const;

            // - Correct boundary velocities on walls
            void correctBoundaryVelocity( volVectorField & U ) const;

            // - Return the coriolis force
            tmp<DimensionedField<vector, volMesh> > Fcoriolis() const;

            // - Add coriolis term to matrix
            void addCoriolis( fvVectorMatrix & UEqn ) const;

            // - Make the given absolute flux relative
            void relativeFlux( surfaceScalarField & phi ) const;

            // - Make the given relative flux absolute
            void absoluteFlux( surfaceScalarField & phi ) const;
        };


        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace SRF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
