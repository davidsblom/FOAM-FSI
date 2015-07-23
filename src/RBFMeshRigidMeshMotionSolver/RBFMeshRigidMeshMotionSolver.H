
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFMeshRigidMeshMotionSolver_H
#define RBFMeshRigidMeshMotionSolver_H

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

#include "RBFMeshMotionSolver.H"
#include "SDC.H"
#include "ESDIRK.H"
#include "TimeIntegrationScheme.H"
#include <yaml-cpp/yaml.h>

namespace Foam
{
    class RBFMeshRigidMeshMotionSolver : public RBFMeshMotionSolver
    {
private:

        // Disallow default bitwise copy construct
        RBFMeshRigidMeshMotionSolver( const RBFMeshRigidMeshMotionSolver & );

        // Disallow default bitwise assignment
        void operator=( const RBFMeshRigidMeshMotionSolver & );

        Foam::vector calcTransformation( scalar t );
        Foam::vector calcVelocity();

        // Rigid body motion
        scalar translationAmplitude;
        scalar translationFrequency;
        Foam::vector translationDirection;
        scalar rotationAmplitude;
        scalar rotationFrequency;
        Foam::vector rotationOrigin;
        Foam::vector oldTransformation;

        scalar told;

public:

        // Runtime type information
        TypeName( "RBFMeshRigidMeshMotionSolver" );

        // Constructors

        // Construct from polyMesh
        RBFMeshRigidMeshMotionSolver(
            const polyMesh & mesh,
            Istream & msData
            );

        // Destructor
        virtual ~RBFMeshRigidMeshMotionSolver();

        // Solve for motion
        virtual void solve();
    };
}

#endif
