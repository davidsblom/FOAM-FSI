
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
        assert( boundaryPoints.empty() );

        forAll( staticPatchIDs, patchId )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[staticPatchIDs[patchId]].meshPoints();

            forAll( meshPoints, i )
            {
                if ( twoDCorrector.marker()[meshPoints[i]] != 0 )
                    continue;

                if ( boundaryPoints.find( meshPoints[i] ) != boundaryPoints.end() )
                    continue;

                Vertex vertex;
                vertex.moving = false;
                vertex.id = meshPoints[i];

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                    vertex.data.push_back( mesh().points()[vertex.id][j] );

                boundaryPoints[vertex.id] = vertex;
            }
        }

        forAll( movingPatchIDs, patchId )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchId]].meshPoints();

            forAll( meshPoints, i )
            {
                if ( twoDCorrector.marker()[meshPoints[i]] != 0 )
                    continue;

                if ( boundaryPoints.find( meshPoints[i] ) != boundaryPoints.end() )
                    continue;

                Vertex vertex;
                vertex.moving = true;
                vertex.id = meshPoints[i];

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                    vertex.data.push_back( mesh().points()[vertex.id][j] );

                boundaryPoints[vertex.id] = vertex;
            }
        }

        std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>() );
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double> () );
        El::Zeros( *positions, boundaryPoints.size(), mesh().nGeometricD() );

        int nbInterpolationPoints = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
                nbInterpolationPoints++;
        }

        El::Zeros( *positionsInterpolation, nbInterpolationPoints, mesh().nGeometricD() );

        positions->Reserve( boundaryPoints.size() * mesh().nGeometricD() );

        {
            int i = 0;

            for ( const auto & vertex : boundaryPoints )
            {
                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                {
                    positions->QueueUpdate( i, j, vertex.second.data[j] );
                }

                i++;
            }
        }

        positionsInterpolation->Reserve( nbInterpolationPoints * mesh().nGeometricD() );

        int index = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] != 0 )
                continue;

            for ( int j = 0; j < mesh().nGeometricD(); j++ )
            {
                positionsInterpolation->QueueUpdate( index, j, mesh().points()[i][j] );
            }

            index++;
        }

        std::unique_ptr<rbf::RBFFunctionInterface> rbfFunction( new rbf::TPSFunction() );
        rbf->compute( std::move( rbfFunction ), std::move( positions ), std::move( positionsInterpolation ) );
    }

    std::unique_ptr<El::DistMatrix<double> > data( new El::DistMatrix<double>() );
    El::Zeros( *data, boundaryPoints.size(), mesh().nGeometricD() );

    for ( auto & vertex : boundaryPoints )
    {
        if ( not vertex.second.moving )
        {
            for ( auto & el: vertex.second.data )
                el = 0.0;
        }
    }

    forAll( movingPatchIDs, patchId )
    {
        const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchId]].meshPoints();

        forAll( meshPoints, i )
        {
            int pointId = meshPoints[i];

            if ( boundaryPoints.find( pointId ) != boundaryPoints.end() )
            {
                auto & vertex = boundaryPoints[pointId];

                for ( int j = 0; j < mesh().nGeometricD(); j++ )
                {
                    vertex.data[j] = motionCenters[movingPatchIDs[patchId]][i][j];
                }
            }
        }
    }

    data->Reserve( boundaryPoints.size() * mesh().nGeometricD() );

    {
        int i = 0;

        for ( const auto & vertex : boundaryPoints )
        {
            for ( int j = 0; j < mesh().nGeometricD(); j++ )
                data->QueueUpdate( i, j, vertex.second.data[j] );

            i++;
        }
    }

    std::unique_ptr<El::DistMatrix<double> > result = rbf->interpolate( data );

    vectorField valuesInterpolationField( mesh().points().size(), Foam::vector::zero );

    int nbInterpolationPoints = 0;
    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] != 0 )
            continue;

        nbInterpolationPoints++;
    }

    std::vector<double> buffer;
    result->ReservePulls( nbInterpolationPoints * mesh().nGeometricD() );

    int index = 0;

    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] != 0 )
            continue;

        for ( int j = 0; j < mesh().nGeometricD(); j++ )
        {
            result->QueuePull( index, j );
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
