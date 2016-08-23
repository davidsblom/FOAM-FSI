
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ElRBFMeshMotionSolver.H"
#include "TPSFunction.H"
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
    const Foam::polyMesh & mesh,
    Foam::Istream & msData
    )
    :
    motionSolver( mesh ),
    motionCenters( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    newPoints( mesh.points().size(), vector::zero ),
    movingPatchIDs( wordList( lookup( "movingPatches" ) ).size(), 0 ),
    staticPatchIDs( wordList( lookup( "staticPatches" ) ).size(), 0 ),
    rbf( new rbf::ElRBFInterpolation() ),
    twoDCorrector( mesh )
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
    // Prepare new points: same as old point
    tmp<pointField> tnewPoints
    (
        new vectorField( mesh().nPoints(), vector::zero )
    );

    pointField & newPoints = tnewPoints();

    newPoints = this->newPoints;

    // Add old point positions
    newPoints += mesh().points();

    return tnewPoints;
}

// As a first step, the motion is defined in the
void ElRBFMeshMotionSolver::setMotion( const Field<vectorField> & motion )
{
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
}

void ElRBFMeshMotionSolver::updateMesh( const mapPolyMesh & )
{
    assert( false );
}

void ElRBFMeshMotionSolver::solve()
{
    if ( not rbf->initialized() )
    {
        forAll( staticPatchIDs, patchId )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[staticPatchIDs[patchId]].meshPoints();

            forAll( meshPoints, i )
            {
                if ( twoDCorrector.marker()[meshPoints[i]] != 0 )
                    continue;

                if ( staticPoints.find( meshPoints[i] ) != staticPoints.end() )
                    continue;

                Vertex vertex;
                vertex.owner = true;
                vertex.patchId = patchId;
                vertex.meshPointId = i;
                vertex.pointId = meshPoints[i];
                vertex.globalPointId = 0;
                vertex.id = staticPoints.size();

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                    vertex.coord.push_back( mesh().points()[vertex.pointId][j] );

                staticPoints[vertex.pointId] = vertex;
            }
        }

        assert( movingPoints.size() == 0 );

        forAll( movingPatchIDs, patchId )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchId]].meshPoints();

            forAll( meshPoints, i )
            {
                if ( twoDCorrector.marker()[meshPoints[i]] != 0 )
                    continue;

                if ( staticPoints.find( meshPoints[i] ) != staticPoints.end() )
                    continue;

                if ( movingPoints.find( meshPoints[i] ) != movingPoints.end() )
                    continue;

                Vertex vertex;
                vertex.owner = true;
                vertex.patchId = patchId;
                vertex.meshPointId = i;
                vertex.pointId = meshPoints[i];
                vertex.globalPointId = 0;
                vertex.id = movingPoints.size();

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                    vertex.coord.push_back( mesh().points()[vertex.pointId][j] );

                movingPoints[vertex.pointId] = vertex;
            }
        }

        std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>() );
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double> () );
        El::Zeros( *positions, staticPoints.size() + movingPoints.size(), mesh().nGeometricD() );

        int nbInterpolationPoints = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
                nbInterpolationPoints++;
        }

        El::Zeros( *positionsInterpolation, nbInterpolationPoints, mesh().nGeometricD() );

        positions->Reserve( positions->LocalHeight() * positions->LocalWidth() );

        for ( const auto & vertex : movingPoints )
        {
            const int globalRow = positions->GlobalRow( vertex.second.id );

            for ( int j = 0; j < mesh().nGeometricD(); j++ )
            {
                const int globalCol = positions->GlobalCol( j );

                positions->QueueUpdate( globalRow, globalCol, vertex.second.coord[j] );
            }
        }

        for ( const auto & vertex : staticPoints )
        {
            const int globalRow = positions->GlobalRow( movingPoints.size() + vertex.second.id );

            for ( int j = 0; j < mesh().nGeometricD(); j++ )
            {
                const int globalCol = positions->GlobalCol( j );

                positions->QueueUpdate( globalRow, globalCol, vertex.second.coord[j] );
            }
        }

        positionsInterpolation->Reserve( positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth() );

        int index = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] != 0 )
                continue;

            const int globalRow = positionsInterpolation->GlobalRow( index );

            for ( int j = 0; j < mesh().nGeometricD(); j++ )
            {
                const int globalCol = positionsInterpolation->GlobalCol( j );

                positionsInterpolation->QueueUpdate( globalRow, globalCol, mesh().points()[i][j] );
            }

            index++;
        }

        std::unique_ptr<rbf::RBFFunctionInterface> rbfFunction( new rbf::TPSFunction() );
        rbf->compute( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );
    }

    std::unique_ptr<El::DistMatrix<double> > data( new El::DistMatrix<double>() );
    El::Zeros( *data, staticPoints.size() + movingPoints.size(), mesh().nGeometricD() );

    data->Reserve( movingPoints.size() * data->LocalWidth() );

    forAll( movingPatchIDs, patchId )
    {
        const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchId]].meshPoints();

        forAll( meshPoints, i )
        {
            int pointId = meshPoints[i];

            if ( movingPoints.find( pointId ) != movingPoints.end() )
            {
                const auto & vertex = movingPoints[pointId];
                const int globalRow = data->GlobalRow( vertex.id );

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                {
                    const int globalCol = data->GlobalCol( j );
                    data->QueueUpdate( globalRow, globalCol, motionCenters[movingPatchIDs[patchId]][i][j] );
                }
            }
        }
    }

    std::unique_ptr<El::DistMatrix<double> > result = rbf->interpolate( data );

    vectorField valuesInterpolationField( mesh().points().size(), Foam::vector::zero );

    std::vector<double> buffer;
    result->ReservePulls( result->LocalHeight() * result->LocalWidth() );

    int index = 0;

    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] != 0 )
            continue;

        const int globalRow = result->GlobalRow( index );

        for ( int j = 0; j < mesh().nGeometricD(); j++ )
        {
            const int globalCol = result->GlobalCol( j );
            result->QueuePull( globalRow, globalCol );
        }

        index++;
    }
    result->ProcessPullQueue( buffer );

    index = 0;

    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] != 0 )
            continue;

        for ( int j = 0; j < mesh().nGeometricD(); j++ )
        {
            valuesInterpolationField[i][j] = buffer[index];
            index++;
        }
    }

    twoDCorrector.setShadowSide( valuesInterpolationField );

    assert( newPoints.size() == valuesInterpolationField.size() );

    newPoints = valuesInterpolationField;
}
