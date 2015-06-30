
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFMeshRigidMotionSolver_H
#define RBFMeshRigidMotionSolver_H

#include "motionSolver.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "RBFInterpolation.H"
#include "RBFCoarsening.H"
#include "TPSFunction.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"
#include "twoDPointCorrectorRBF.H"
#include <memory>
#include <assert.h>

namespace Foam
{
    class RBFMeshRigidMotionSolver : public motionSolver
    {
private:

        // Disallow default bitwise copy construct
        RBFMeshRigidMotionSolver( const RBFMeshRigidMotionSolver & );

        // Disallow default bitwise assignment
        void operator=( const RBFMeshRigidMotionSolver & );

        Foam::vector calcTransformation();

        Field<vectorField> motionCenters;

        wordList staticPatches;
        labelList staticPatchIDs;
        wordList movingPatches;
        labelList movingPatchIDs;
        wordList fixedPatches;
        labelList fixedPatchIDs;

        pointField newPoints;

        std::shared_ptr<rbf::RBFCoarsening> rbf;

        labelList nbGlobalFaceCenters;
        labelList nbGlobalMovingFaceCenters;
        labelList nbGlobalStaticFaceCenters;
        labelList nbGlobalFixedFaceCenters;

        List<labelList> globalMovingPointsLabelList;

        // 2d corrector
        twoDPointCorrectorRBF twoDCorrector;
        int nbPoints;
        bool faceCellCenters;
        bool cpu;

        // Rigid body motion
        double translationAmplitude;
        double translationFrequency;
        Foam::vector translationDirection;
        double rotationAmplitude;
        double rotationFrequency;
        Foam::vector rotationOrigin;

public:

        // Runtime type information
        TypeName( "RBFMeshRigidMotionSolver" );

        // Constructors

        // Construct from polyMesh
        RBFMeshRigidMotionSolver(
            const polyMesh & mesh,
            Istream & msData
            );

        // Destructor
        virtual ~RBFMeshRigidMotionSolver();

        // Set the motion of every boundary patch, where m is equal to number of patches and with non-empty vectorFields for moving patches.
        // The motion is defined at the face centers of the boundary patch.
        void setMotion( const Field<vectorField> & motion );

        // Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        // Solve for motion
        virtual void solve();

        // Update the mesh corresponding to given map
        virtual void updateMesh( const mapPolyMesh & );
    };
}

#endif
