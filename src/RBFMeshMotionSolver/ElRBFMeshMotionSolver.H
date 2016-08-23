
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "ElRBFInterpolation.H"
#include "motionSolver.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "twoDPointCorrectorRBF.H"
#include <unordered_map>

struct Vertex
{
    bool owner;
    bool moving;
    int pointId;
    int globalPointId;
    int id;
    std::vector<double> coord;
};

namespace Foam
{
    class ElRBFMeshMotionSolver : public motionSolver
    {
        protected:
            // Disallow default bitwise copy construct
            ElRBFMeshMotionSolver( const ElRBFMeshMotionSolver & );

            // Disallow default bitwise assignment
            void operator=( const ElRBFMeshMotionSolver & );

            Field<vectorField> motionCenters;

            pointField newPoints;

        public:
            // Runtime type information
            TypeName( "ElRBFMeshMotionSolver" );

            // Constructors

            // Construct from polyMesh
            ElRBFMeshMotionSolver(
                const polyMesh & mesh,
                Istream & msData
                );

            // Destructor
            virtual ~ElRBFMeshMotionSolver();

            // Set the motion of every boundary patch, where m is equal to number of patches and with non-empty vectorFields for moving patches.
            // The motion is defined at the face centers of the boundary patch.
            void setMotion( const Field<vectorField> & motion );

            // Return point location obtained from the current motion field
            virtual tmp<pointField> curPoints() const;

            // Solve for motion
            virtual void solve();

            // Update the mesh corresponding to given map
            virtual void updateMesh( const mapPolyMesh & );

        private:
            labelList movingPatchIDs;
            labelList staticPatchIDs;
            std::unique_ptr<rbf::ElRBFInterpolation> rbf;
            twoDPointCorrectorRBF twoDCorrector;
            std::unordered_map<int, Vertex> boundaryPoints;
    };
}
