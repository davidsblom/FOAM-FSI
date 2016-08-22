
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFMeshMotionSolver.H"
#include <cassert>

using namespace Foam;

defineTypeNameAndDebug( ElRBFMeshMotionSolver, 0 );

addToRunTimeSelectionTable
(
    motionSolver,
    ElRBFMeshMotionSolver,
    dictionary
);

ElRBFMeshMotionSolver::ElRBFMeshMotionSolver(
    const polyMesh & mesh,
    Istream & msData
    )
    :
    motionSolver( mesh ),
    motionCenters( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    newPoints( mesh.points().size(), vector::zero ),
    movingPatchIDs( wordList( lookup( "movingPatches" ) ).size(), 0 ),
    staticPatchIDs( wordList( lookup( "staticPatches" ) ).size(), 0 )
{
    wordList staticPatches( lookup( "staticPatches" ) );
    wordList movingPatches( lookup( "movingPatches" ) );
    assert( staticPatches.size() == staticPatchIDs.size() );
    assert( movingPatches.size() == movingPatchIDs.size() );

    // Find IDs of staticPatches
    forAll( staticPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( staticPatches[patchI] );

        assert( patchIndex >= 0 );

        staticPatchIDs[patchI] = patchIndex;
    }

    // Find IDs of movingPatches
    forAll( movingPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( movingPatches[patchI] );

        assert( patchIndex >= 0 );

        movingPatchIDs[patchI] = patchIndex;
    }

    // Verify that a patch is not defined as a static and a moving patch

    forAll( staticPatchIDs, staticPatchI )
    {
        // Search the moving patches for static patchI
        forAll( movingPatchIDs, movingPatchI )
        {
            assert( movingPatchIDs[movingPatchI] != staticPatchIDs[staticPatchI] );
        }
    }

    // Initialize zero motion

    forAll( movingPatchIDs, patchId )
    {
        motionCenters[movingPatchIDs[patchId]] = vectorField( mesh.boundaryMesh()[movingPatchIDs[patchId]].meshPoints().size(), Foam::vector::zero );
    }
}

ElRBFMeshMotionSolver::~ElRBFMeshMotionSolver()
{}

tmp<pointField> ElRBFMeshMotionSolver::curPoints() const
{
    Pout << "ElRBFMeshMotionSolver::curPoints start" << endl;

    // Prepare new points: same as old point
    tmp<pointField> tnewPoints
    (
        new vectorField( mesh().nPoints(), vector::zero )
    );

    pointField & newPoints = tnewPoints();

    newPoints = this->newPoints;

    // Add old point positions
    newPoints += mesh().points();

    Pout << "ElRBFMeshMotionSolver::curPoints end" << endl;

    return tnewPoints;
}

// As a first step, the motion is defined in the
void ElRBFMeshMotionSolver::setMotion( const Field<vectorField> & motion )
{
    Pout << "ElRBFMeshMotionSolver::setMotion start" << endl;

    // Input checking

    assert( motion.size() == mesh().boundaryMesh().size() );

    forAll( motion, ipatch )
    {
        const vectorField & mpatch = motion[ipatch];

        if ( mpatch.size() > 0 )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size() );

        // Check whether the size of a moving patch is equal to the number of face centers in the patch
        // First check if patchid is a moving patch
        bool movingPatch = false;
        forAll( movingPatchIDs, movingPatchI )
        {
            if ( movingPatchIDs[movingPatchI] == ipatch )
                movingPatch = true;
        }

        if ( movingPatch )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size() );
    }

    motionCenters = motion;
    Pout << "ElRBFMeshMotionSolver::setMotion end" << endl;
}

void ElRBFMeshMotionSolver::updateMesh( const mapPolyMesh & )
{
    assert( false );
}

void ElRBFMeshMotionSolver::solve()
{
    Pout << "ElRBFMeshMotionSolver::solve start" << endl;

    Pout << "ElRBFMeshMotionSolver::solve end" << endl;
}
