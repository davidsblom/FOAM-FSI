
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <mxx/sort.hpp>
#include "ElRBFMeshMotionSolver.H"
#include "TPSFunction.H"
#include <cassert>
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"
#include "LinearFunction.H"
#include "NoCoarsening.H"
#include "AdaptiveCoarsening.H"
#include "UnitCoarsening.H"

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
    Foam::Istream &
    )
    :
    motionSolver( mesh ),
    motionCenters( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    newPoints( mesh.points().size(), vector::zero ),
    movingPatchIDs( wordList( lookup( "movingPatches" ) ).size(), 0 ),
    staticPatchIDs( wordList( lookup( "staticPatches" ) ).size(), 0 ),
    rbf( nullptr ),
    rbfFunction( nullptr ),
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

    // Initialize RBF interpolator

    dictionary & dict = subDict( "interpolation" );

    word function = dict.lookup( "function" );

    assert( function == "TPS" || function == "WendlandC0" || function == "WendlandC2" || function == "WendlandC4" || function == "WendlandC6" || function == "linear" );

    Info << "Radial Basis Function interpolation: Selecting RBF function: " << function << endl;

    if ( function == "TPS" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::TPSFunction() );

    if ( function == "WendlandC0" )
    {
        scalar radius = readScalar( dict.lookup( "radius" ) );
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC0Function( radius ) );
    }

    if ( function == "WendlandC2" )
    {
        scalar radius = readScalar( dict.lookup( "radius" ) );
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC2Function( radius ) );
    }

    if ( function == "WendlandC4" )
    {
        scalar radius = readScalar( dict.lookup( "radius" ) );
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC4Function( radius ) );
    }

    if ( function == "WendlandC6" )
    {
        scalar radius = readScalar( dict.lookup( "radius" ) );
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC6Function( radius ) );
    }

    if ( function == "linear" )
    {
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::LinearFunction() );
    }

    assert( rbfFunction );

    // Setup coarsening strategy

    word coarseningStrategy = lookup( "coarseningStrategy" );

    assert( coarseningStrategy == "NoCoarsening" || coarseningStrategy == "UnitCoarsening" || coarseningStrategy == "AdaptiveCoarsening" );

    if ( coarseningStrategy == "NoCoarsening" )
    {
        rbf = std::unique_ptr<rbf::Coarsener>( new rbf::NoCoarsening() );
    }

    if ( coarseningStrategy == "UnitCoarsening" )
    {
        double tol = double( readScalar( subDict( "UnitCoarsening" ).lookup( "tol" ) ) );
        int minPoints = readLabel( subDict( "UnitCoarsening" ).lookup( "minPoints" ) );
        int maxPoints = readLabel( subDict( "UnitCoarsening" ).lookup( "maxPoints" ) );
        rbf = std::unique_ptr<rbf::Coarsener>( new rbf::UnitCoarsening( tol, minPoints, maxPoints ) );
    }

    if ( coarseningStrategy == "AdaptiveCoarsening" )
    {
        double tol = double( readScalar( subDict( "AdaptiveCoarsening" ).lookup( "tol" ) ) );
        double reselectionTol = double( readScalar( subDict( "AdaptiveCoarsening" ).lookup( "reselectionTol" ) ) );
        int minPoints = readLabel( subDict( "AdaptiveCoarsening" ).lookup( "minPoints" ) );
        int maxPoints = readLabel( subDict( "AdaptiveCoarsening" ).lookup( "maxPoints" ) );
        rbf = std::unique_ptr<rbf::Coarsener>( new rbf::AdaptiveCoarsening( tol, reselectionTol, minPoints, maxPoints ) );
    }

    assert( rbf );
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
    typedef std::tuple<int, double, double, double> tuple_type;

    // comparator for sorting
    auto cmp = []( const tuple_type & x, const tuple_type & y ) {
            return std::get<0>( x ) < std::get<0>( y );
        };

    auto cmp_unique = []( const tuple_type & x, const tuple_type & y ) {
            return std::get<0>( x ) == std::get<0>( y );
        };

    // If the points are already selected, the rbf.compute() is already
    // called.
    if ( boundaryPoints.empty() )
    {
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

        if ( Pstream::parRun() )
        {
            IOobject addrHeader
            (
                "pointProcAddressing",
                mesh().facesInstance(),
                mesh().meshSubDir,
                mesh(),
                IOobject::MUST_READ
            );

            assert( addrHeader.headerOk() );
            labelIOList pointProcAddressing( addrHeader );

            assert( pointProcAddressing.size() == mesh().points().size() );

            for ( auto & vertex : boundaryPoints )
            {
                vertex.second.id = pointProcAddressing[vertex.second.id];
            }
        }

        std::vector<tuple_type> boundaryData;

        for ( const auto & vertex : boundaryPoints )
        {
            if ( mesh().nGeometricD() == 2 )
                boundaryData.push_back( std::make_tuple( vertex.second.id, vertex.second.data[0], vertex.second.data[1], 0 ) );
            else
                boundaryData.push_back( std::make_tuple( vertex.second.id, vertex.second.data[0], vertex.second.data[1], vertex.second.data[2] ) );
        }

        mxx::sort( boundaryData.begin(), boundaryData.end(), cmp );

        auto it = mxx::unique( boundaryData.begin(), boundaryData.end(), cmp_unique );

        boundaryData.resize( std::distance( boundaryData.begin(), it ) );

        std::unique_ptr<El::DistMatrix<double> > positions( new El::DistMatrix<double>() );
        std::unique_ptr<El::DistMatrix<double> > positionsInterpolation( new El::DistMatrix<double> () );

        std::vector<size_t> allBoundaryPoints = mxx::allgather( boundaryData.size() );
        size_t totalNbBoundaryPoints = 0;

        for ( size_t n : allBoundaryPoints )
            totalNbBoundaryPoints += n;

        El::Zeros( *positions, totalNbBoundaryPoints, mesh().nGeometricD() );

        int nbInterpolationPoints = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
                nbInterpolationPoints++;
        }

        std::vector<int> allInterpolationPoints = mxx::allgather( nbInterpolationPoints );
        int totalNbInterpolationPoints = 0;

        for ( int n : allInterpolationPoints )
            totalNbInterpolationPoints += n;

        El::Zeros( *positionsInterpolation, totalNbInterpolationPoints, mesh().nGeometricD() );

        positions->Reserve( boundaryPoints.size() * mesh().nGeometricD() );

        size_t localBoundaryOffset = 0;

        for ( int i = 0; i < Pstream::myProcNo(); i++ )
            localBoundaryOffset += allBoundaryPoints[i];

        {
            int i = 0;

            for ( const auto & vertex : boundaryData )
            {
                positions->QueueUpdate( localBoundaryOffset + i, 0, std::get<1>( vertex ) );
                positions->QueueUpdate( localBoundaryOffset + i, 1, std::get<2>( vertex ) );

                if ( mesh().nGeometricD() == 3 )
                    positions->QueueUpdate( localBoundaryOffset + i, 2, std::get<3>( vertex ) );

                i++;
            }
        }

        positionsInterpolation->Reserve( nbInterpolationPoints * mesh().nGeometricD() );

        size_t localInterpolationOffset = 0;

        for ( int i = 0; i < Pstream::myProcNo(); i++ )
            localInterpolationOffset += allInterpolationPoints[i];

        int index = 0;

        forAll( mesh().points(), i )
        {
            if ( twoDCorrector.marker()[i] != 0 )
                continue;

            for ( int j = 0; j < mesh().nGeometricD(); j++ )
            {
                positionsInterpolation->QueueUpdate( localInterpolationOffset + index, j, mesh().points()[i][j] );
            }

            index++;
        }

        rbf->compute( rbfFunction, std::move( positions ), std::move( positionsInterpolation ) );
    }

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

    std::vector<tuple_type> boundaryData;

    for ( const auto & vertex : boundaryPoints )
    {
        if ( mesh().nGeometricD() == 2 )
            boundaryData.push_back( std::make_tuple( vertex.second.id, vertex.second.data[0], vertex.second.data[1], 0 ) );
        else
            boundaryData.push_back( std::make_tuple( vertex.second.id, vertex.second.data[0], vertex.second.data[1], vertex.second.data[2] ) );
    }

    mxx::sort( boundaryData.begin(), boundaryData.end(), cmp );

    auto it = mxx::unique( boundaryData.begin(), boundaryData.end(), cmp_unique );

    boundaryData.resize( std::distance( boundaryData.begin(), it ) );

    std::vector<size_t> allBoundaryPoints = mxx::allgather( boundaryData.size() );
    size_t totalNbBoundaryPoints = 0;

    for ( size_t n : allBoundaryPoints )
        totalNbBoundaryPoints += n;

    std::unique_ptr<El::DistMatrix<double> > data( new El::DistMatrix<double>() );
    El::Zeros( *data, totalNbBoundaryPoints, mesh().nGeometricD() );

    data->Reserve( boundaryPoints.size() * mesh().nGeometricD() );

    size_t localBoundaryOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        localBoundaryOffset += allBoundaryPoints[i];

    {
        int i = 0;

        for ( const auto & vertex : boundaryData )
        {
            data->QueueUpdate( localBoundaryOffset + i, 0, std::get<1>( vertex ) );
            data->QueueUpdate( localBoundaryOffset + i, 1, std::get<2>( vertex ) );

            if ( mesh().nGeometricD() == 3 )
                data->QueueUpdate( localBoundaryOffset + i, 2, std::get<3>( vertex ) );

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

    std::vector<int> allInterpolationPoints = mxx::allgather( nbInterpolationPoints );

    std::vector<double> buffer;
    result->ReservePulls( nbInterpolationPoints * mesh().nGeometricD() );

    size_t localInterpolationOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        localInterpolationOffset += allInterpolationPoints[i];

    int index = 0;

    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] != 0 )
            continue;

        for ( int j = 0; j < mesh().nGeometricD(); j++ )
        {
            result->QueuePull( localInterpolationOffset + index, j );
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
