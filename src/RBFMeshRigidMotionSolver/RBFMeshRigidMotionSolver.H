
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

#include "RBFMeshMotionSolver.H"
#include "SDC.H"
#include "ESDIRK.H"
#include "TimeIntegrationScheme.H"
#include <yaml-cpp/yaml.h>

namespace Foam
{
    class RBFMeshRigidMotionSolver : public RBFMeshMotionSolver
    {
private:

        // Disallow default bitwise copy construct
        RBFMeshRigidMotionSolver( const RBFMeshRigidMotionSolver & );

        // Disallow default bitwise assignment
        void operator=( const RBFMeshRigidMotionSolver & );

        Foam::vector calcTransformation( double t );
        Foam::vector calcVelocity();

        // Rigid body motion
        double translationAmplitude;
        double translationFrequency;
        Foam::vector translationDirection;
        double rotationAmplitude;
        double rotationFrequency;
        Foam::vector rotationOrigin;
        Foam::vector oldTransformation;

        scalar told;

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

        // Solve for motion
        virtual void solve();
    };
}

#endif
