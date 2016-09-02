
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
#include "Coarsener.H"

struct Vertex
{
    bool moving;
    int id;
    std::vector<double> data;
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
            El::Grid grid;
            std::unique_ptr<rbf::Coarsener> rbf;
            std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
            twoDPointCorrectorRBF twoDCorrector;
            std::map<int, Vertex> boundaryPoints;
            std::unique_ptr<rbf::ElDistVector> data;
    };
}
