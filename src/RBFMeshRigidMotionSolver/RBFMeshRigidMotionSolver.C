
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFMeshRigidMotionSolver.H"

using namespace Foam;

defineTypeNameAndDebug( RBFMeshRigidMotionSolver, 0 );

addToRunTimeSelectionTable
(
    motionSolver,
    RBFMeshRigidMotionSolver,
    dictionary
);

RBFMeshRigidMotionSolver::RBFMeshRigidMotionSolver(
    const polyMesh & mesh,
    Istream & msData
    )
    :
    motionSolver( mesh ),
    motionCenters( mesh.boundaryMesh().size(), vectorField( 0 ) ),
    staticPatches( lookup( "staticPatches" ) ),
    staticPatchIDs( staticPatches.size() ),
    movingPatches( lookup( "movingPatches" ) ),
    movingPatchIDs( movingPatches.size() ),
    fixedPatches( lookup( "fixedPatches" ) ),
    fixedPatchIDs( fixedPatches.size() ),
    newPoints( mesh.points().size(), vector::zero ),
    nbGlobalFaceCenters( Pstream::nProcs(), 0 ),
    nbGlobalMovingFaceCenters( Pstream::nProcs(), 0 ),
    nbGlobalStaticFaceCenters( Pstream::nProcs(), 0 ),
    nbGlobalFixedFaceCenters( Pstream::nProcs(), 0 ),
    globalMovingPointsLabelList( mesh.boundaryMesh().size(), labelList( 0 ) ),
    twoDCorrector( mesh ),
    nbPoints( 0 ),
    faceCellCenters( true ),
    cpu( false ),
    translationAmplitude( 0 ),
    translationFrequency( 0 ),
    translationDirection( Foam::vector::zero ),
    rotationAmplitude( 0 ),
    rotationFrequency( 0 ),
    rotationOrigin( Foam::vector::zero ),
    told(0)
{
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

    assert( movingPatches.size() == 1 );

    // Find IDs of fixedPatches
    forAll( fixedPatches, patchI )
    {
        label patchIndex = mesh.boundaryMesh().findPatchID( fixedPatches[patchI] );

        assert( patchIndex >= 0 );

        fixedPatchIDs[patchI] = patchIndex;
    }

    // Verify that a patch is not defined as a static and a moving patch

    forAll( staticPatchIDs, staticPatchI )
    {
        // Search the moving patches for static patchI
        forAll( movingPatchIDs, movingPatchI )
        {
            assert( movingPatchIDs[movingPatchI] != staticPatchIDs[staticPatchI] );
        }

        // Search the fixed patches for static patchI
        forAll( fixedPatchIDs, fixedPatchI )
        {
            assert( fixedPatchIDs[fixedPatchI] != staticPatchIDs[staticPatchI] );
        }
    }

    forAll( fixedPatchIDs, fixedPatchI )
    {
        // Search the moving patches for fixed patchI
        forAll( movingPatchIDs, movingPatchI )
        {
            assert( movingPatchIDs[movingPatchI] != fixedPatchIDs[fixedPatchI] );
        }
    }

    // Initialize RBF interpolator

    dictionary & dict = subDict( "interpolation" );

    word function = dict.lookup( "function" );

    assert( function == "TPS" || function == "WendlandC0" || function == "WendlandC2" || function == "WendlandC4" || function == "WendlandC6" );

    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;

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

    assert( rbfFunction );

    bool polynomialTerm = dict.lookupOrDefault( "polynomial", false );
    bool cpu = dict.lookupOrDefault( "cpu", false );
    this->cpu = dict.lookupOrDefault( "fullCPU", false );
    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator( new rbf::RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );

    if ( this->cpu == true )
        assert( cpu == true );

    bool coarsening = readBool( subDict( "coarsening" ).lookup( "enabled" ) );
    double tol = 0.1;
    double tolLivePointSelection = 0.1;
    bool livePointSelection = false;
    bool exportSelectedPoints = false;
    int coarseningMinPoints = 1;
    int coarseningMaxPoints = 2;
    bool twoPointSelection = false;
    bool surfaceCorrection = false;
    double ratioRadiusError = 10.0;

    if ( coarsening )
    {
        tol = readScalar( subDict( "coarsening" ).lookup( "tol" ) );
        coarseningMinPoints = readLabel( subDict( "coarsening" ).lookup( "minPoints" ) );
        coarseningMaxPoints = readLabel( subDict( "coarsening" ).lookup( "maxPoints" ) );
        livePointSelection = readBool( subDict( "coarsening" ).lookup( "livePointSelection" ) );
        exportSelectedPoints = readBool( subDict( "coarsening" ).lookup( "exportSelectedPoints" ) );
        twoPointSelection = subDict( "coarsening" ).lookupOrDefault( "twoPointSelection", false );
    }

    if ( livePointSelection )
    {
        tolLivePointSelection = readScalar( subDict( "coarsening" ).lookup( "tolLivePointSelection" ) );
        surfaceCorrection = subDict( "coarsening" ).lookupOrDefault( "surfaceCorrection", false );

        if ( surfaceCorrection )
        {
            ratioRadiusError = subDict( "coarsening" ).lookupOrDefault( "ratioRadiusError", 10.0 );
        }
    }

    rbf = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, true, tol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, twoPointSelection, surfaceCorrection, ratioRadiusError, exportSelectedPoints ) );

    faceCellCenters = lookupOrDefault( "faceCellCenters", true );

    Info << "RBF mesh deformation settings:" << endl;
    Info << "    interpolation function = " << function << endl;
    Info << "    interpolation polynomial term = " << polynomialTerm << endl;
    Info << "    interpolation cpu formulation = " << cpu << endl;
    Info << "    coarsening = " << coarsening << endl;
    Info << "        coarsening tolerance = " << tol << endl;
    Info << "        coarsening reselection tolerance = " << tolLivePointSelection << endl;
    Info << "        coarsening two-point selection = " << twoPointSelection << endl;

    // Rigid body motion parameters
    assert( faceCellCenters == false );

    dictionary & rigidBodyMotionDict = subDict( "rigidBodyMotion" );
    translationAmplitude = readScalar( rigidBodyMotionDict.lookup( "translationAmplitude" ) );
    translationFrequency = readScalar( rigidBodyMotionDict.lookup( "translationFrequency" ) );
    rotationAmplitude = readScalar( rigidBodyMotionDict.lookup( "rotationAmplitude" ) );
    rotationFrequency = readScalar( rigidBodyMotionDict.lookup( "rotationFrequency" ) );
    translationDirection = rigidBodyMotionDict.lookup( "translationDirection" );
    rotationOrigin = rigidBodyMotionDict.lookup( "rotationOrigin" );

    string filename = static_cast<std::string>( mesh.time().rootPath() ) + "/" + static_cast<std::string>( mesh.time().caseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );

    assert( config["sdc"] || config["esdirk"] );

    if ( config["sdc"] )
    {
        YAML::Node sdcConfig( config["sdc"] );
        assert( sdcConfig["convergence-tolerance"] );
        assert( sdcConfig["number-of-points"] );
        assert( sdcConfig["quadrature-rule"] );

        int n = sdcConfig["number-of-points"].as<int>();
        scalar tol = sdcConfig["convergence-tolerance"].as<scalar>();
        std::string quadratureRule = sdcConfig["quadrature-rule"].as<std::string>();

        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( quadratureRule, n, tol ) );
    }

    if ( config["esdirk"] )
    {
        YAML::Node esdirkConfig( config["esdirk"] );

        assert( esdirkConfig["method"] );

        std::string method = esdirkConfig["method"].as<std::string>();

        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::ESDIRK( method ) );
    }
}

RBFMeshRigidMotionSolver::~RBFMeshRigidMotionSolver()
{}

tmp<pointField> RBFMeshRigidMotionSolver::curPoints() const
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
void RBFMeshRigidMotionSolver::setMotion( const Field<vectorField> & motion )
{
    // Input checking

    assert( motion.size() == mesh().boundaryMesh().size() );

    forAll( motion, ipatch )
    {
        const vectorField & mpatch = motion[ipatch];

        // Check whether the size of patch motion is equal to number of face centers in patch
        if ( faceCellCenters && mpatch.size() > 0 )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size() );

        if ( not faceCellCenters && mpatch.size() > 0 )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size() );

        // Check whether the size of a moving patch is equal to the number of face centers in the patch
        // First check if patchid is a moving patch
        bool movingPatch = false;
        forAll( movingPatchIDs, movingPatchI )
        {
            if ( movingPatchIDs[movingPatchI] == ipatch )
                movingPatch = true;
        }

        if ( faceCellCenters && movingPatch )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size() );

        if ( not faceCellCenters && movingPatch )
            assert( mpatch.size() == mesh().boundaryMesh()[ipatch].meshPoints().size() );
    }

    motionCenters = motion;
}

void RBFMeshRigidMotionSolver::updateMesh( const mapPolyMesh & )
{
    assert( false );
}

Foam::vector RBFMeshRigidMotionSolver::calcTransformation( double t )
{
    scalar smoothStartup = 1;
    bool smoothStart = false;

    if ( smoothStart )
    {
        if ( t < 1.0 / (translationFrequency + SMALL) )
        {
            smoothStartup = 0.5 - 0.5 * Foam::cos( M_PI * translationFrequency * t );
        }
    }

    Foam::vector transformation = smoothStartup * translationAmplitude * translationDirection * ( 0.5 - 0.5 * Foam::cos( M_PI * translationFrequency * t ) );

    return transformation;
}

void RBFMeshRigidMotionSolver::setTimeIntegrationInfo( const bool corrector, const int k )
{
    this->corrector = corrector;
    this->k = k;
}

Foam::vector RBFMeshRigidMotionSolver::calcVelocity()
{
    assert( timeIntegrationScheme );

    scalar t = mesh().time().value();
    scalar dt = mesh().time().deltaT().value();

    Foam::vector transformation = calcTransformation( t ) - calcTransformation( told );

    Foam::vector transformationT = calcTransformation( t );
    Foam::vector transformationTold = calcTransformation( told );
    fsi::vector disp ( 3 ), dudt( 3 ), dispold( 3 );
    disp << transformationT(0), transformationT(1), transformationT(2);
    dispold << transformationTold(0), transformationTold(1), transformationTold(2);

    fsi::vector rhs( 3 ), result( 3 ), f( 3 ), qold;

    if ( not corrector && k == 0 )
        timeIntegrationScheme->setOldSolution( mesh().time().timeIndex(), dispold );

    timeIntegrationScheme->getSourceTerm( corrector, k, dt, rhs, qold );

    scalar smoothStartup = 1;

    scalar t_offset = 0.01;

    if ( t < t_offset )
        smoothStartup = 0;

    f( 0 ) = 0;
    f( 1 ) = 0.5 * smoothStartup * translationAmplitude * std::sin( M_PI * translationFrequency * (t-t_offset) ) * M_PI * translationFrequency;
    f( 2 ) = 0;
    result = dt * f + qold + rhs;
    dudt = ( result - qold - rhs ) / dt;

    timeIntegrationScheme->setFunction( k, f, result );

    for ( int i = 0 ; i < 3; i++ )
        transformation[i] = result(i) - oldTransformation[i];

    told = t;
    //oldTransformation = calcTransformation( told );
    oldTransformation(1) = result(1);

    return transformation;
}

void RBFMeshRigidMotionSolver::solve()
{
    Field<vectorField> motion( mesh().boundaryMesh().size(), vectorField( 0 ) );
    const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[0]].meshPoints();
    motion[movingPatchIDs[0]] = vectorField( meshPoints.size(), Foam::vector::zero );

    Foam::vector transformation = calcVelocity();
    Info << "transformation = " << transformation << endl;
    //transformation = Foam::vector::zero;
    Info << "transformation = " << transformation << endl;

    forAll( motion[movingPatchIDs[0]], i )
    {
        motion[movingPatchIDs[0]][i] = transformation;
    }

    motionCenters = motion;

    assert( motionCenters.size() == mesh().boundaryMesh().size() );

    /*
     * RBF interpolator from face centers to local complete mesh vertices
     * The interpolation consists of the following steps:
     * 1. Build a matrix with the face center positions of the static patches and the moving patches
     * 2. Build a matrix with the positions of every vertex in the local mesh
     * 3. Build a matrix with the displacement/motion of the face center positions of the static patches and the moving patches
     * 4. Perform the interpolation from the face centers to the complete mesh
     * 5. Correct the mesh vertices of the static patches. Set these displacement to zero.
     * 6. Set the motion of the mesh vertices
     */

    /*
     * Step 1: Build a matrix with the face center positions of the static patches and the moving patches
     * The order of the matrix is defined as first a list of the moving patch face centers,
     * thereafter the static patch face centers. These are the control points used by the
     * radial basis function interpolation.
     * The control points should be exactly the same at each processor, and are therefore communicated
     * to each process. As only an absolute RBF interpolation is implemented, this communication step is only
     * performed once per simulation.
     * The global ordering of the data is first the information of the moving patches,
     * thereafter the static patches.
     */

    int nbFaceCenters = 0;
    int nbMovingFaceCenters = 0;
    int nbStaticFaceCenters = 0;
    int nbFixedFaceCenters = 0;

    std::vector<int> staticControlPointLabels;
    std::vector<int> staticControlGlobalPointLabels;
    std::vector<int> fixedControlPointLabels;
    std::vector<int> fixedControlGlobalPointLabels;
    std::vector<int> movingControlPointLabels;
    std::vector<int> movingControlPointPatchIds;
    std::vector<int> movingControlPointIndices;
    std::vector<int> movingControlGlobalPointLabels;

    labelList globalStaticPointsListEnabled( nbStaticFaceCenters, 0 );
    labelList globalFixedPointsListEnabled( nbFixedFaceCenters, 0 );
    labelList globalMovingPointsListEnabled( nbMovingFaceCenters, 0 );
    int globalStaticOffsetNonUnique = 0;
    int globalFixedOffsetNonUnique = 0;
    int globalMovingOffsetNonUnique = 0;

    if ( sum( nbGlobalFaceCenters ) == 0 )
    {
        // Determine the number of face centers
        // The total number of face centers is simply the sum of the face centers
        // on each processor.
        nbFaceCenters = 0;

        // First add the static patches, thereafter the fixed patches, and
        // the moving patches as last.

        forAll( staticPatchIDs, i )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[staticPatchIDs[i]].meshPoints();

            forAll( meshPoints, j )
            {
                if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                    continue;

                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                {
                    staticControlPointLabels.push_back( meshPoints[j] );
                }
            }
        }

        nbStaticFaceCenters = staticControlPointLabels.size();

        forAll( fixedPatchIDs, i )
        {
            const labelList & meshPoints = mesh().boundaryMesh()[fixedPatchIDs[i]].meshPoints();

            forAll( meshPoints, j )
            {
                if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                    continue;

                // Only add the static vertex point if it's not already added to the list
                // and if it's not a static point
                std::vector<int>::iterator it;

                it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                {
                    it = std::find( fixedControlPointLabels.begin(), fixedControlPointLabels.end(), meshPoints[j] );

                    if ( it == fixedControlPointLabels.end() || fixedControlPointLabels.size() == 0 )
                    {
                        fixedControlPointLabels.push_back( meshPoints[j] );
                    }
                }
            }
        }

        nbFixedFaceCenters = fixedControlPointLabels.size();

        if ( faceCellCenters )
        {
            forAll( movingPatchIDs, i )
            {
                nbMovingFaceCenters += mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres().size();
            }
        }

        if ( not faceCellCenters )
        {
            forAll( movingPatchIDs, patchI )
            {
                const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList( meshPoints.size(), 0 );

                forAll( meshPoints, j )
                {
                    if ( twoDCorrector.marker()[meshPoints[j]] != 0 )
                        continue;

                    // Only add the static vertex point if it's not already added to the list
                    // and if it's not a static point
                    std::vector<int>::iterator it;

                    it = std::find( staticControlPointLabels.begin(), staticControlPointLabels.end(), meshPoints[j] );

                    if ( it == staticControlPointLabels.end() || staticControlPointLabels.size() == 0 )
                    {
                        it = std::find( fixedControlPointLabels.begin(), fixedControlPointLabels.end(), meshPoints[j] );

                        if ( it == fixedControlPointLabels.end() || fixedControlPointLabels.size() == 0 )
                        {
                            it = std::find( movingControlPointLabels.begin(), movingControlPointLabels.end(), meshPoints[j] );

                            if ( it == movingControlPointLabels.end() || movingControlPointLabels.size() == 0 )
                            {
                                movingControlPointLabels.push_back( meshPoints[j] );
                                movingControlPointPatchIds.push_back( movingPatchIDs[patchI] );
                                movingControlPointIndices.push_back( j );
                                globalMovingPointsLabelList[movingPatchIDs[patchI]][j] = 1;
                            }
                        }
                    }
                }
            }

            nbMovingFaceCenters = movingControlPointLabels.size();
        }

        if ( Pstream::nProcs() == 1 )
        {
            globalStaticPointsListEnabled.resize( nbStaticFaceCenters );
            globalStaticPointsListEnabled = 1;
            globalFixedPointsListEnabled.resize( nbFixedFaceCenters );
            globalFixedPointsListEnabled = 1;
            globalMovingPointsListEnabled.resize( nbMovingFaceCenters );
            globalMovingPointsListEnabled = 1;
        }

        if ( Pstream::nProcs() > 1 )
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

            // Count the number of global static points including double points
            nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
            nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
            nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
            reduce( nbGlobalStaticFaceCenters, sumOp<labelList>() );
            reduce( nbGlobalFixedFaceCenters, sumOp<labelList>() );
            reduce( nbGlobalMovingFaceCenters, sumOp<labelList>() );
            nbStaticFaceCenters = sum( nbGlobalStaticFaceCenters );
            nbFixedFaceCenters = sum( nbGlobalFixedFaceCenters );

            if ( not faceCellCenters )
                nbMovingFaceCenters = sum( nbGlobalMovingFaceCenters );

            // Construct a list with all the global point labels, thus including
            // also double points. Thereafter, construct a list of static control
            // list which indicates whether the point is already included or not.
            // Use this later to build a list of the unique static control points.
            labelList globalStaticPointsList( nbStaticFaceCenters, 0 );
            labelList globalFixedPointsList( nbFixedFaceCenters, 0 );
            labelList globalMovingPointsList( nbMovingFaceCenters, 0 );
            labelList globalMovingPointsPatchIds( nbMovingFaceCenters, 0 );
            labelList globalMovingPointsIndices( nbMovingFaceCenters, 0 );

            globalStaticOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalStaticOffsetNonUnique += nbGlobalStaticFaceCenters[i];

            globalFixedOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalFixedOffsetNonUnique += nbGlobalFixedFaceCenters[i];

            globalMovingOffsetNonUnique = 0;

            for ( int i = 0; i < Pstream::myProcNo(); i++ )
                globalMovingOffsetNonUnique += nbGlobalMovingFaceCenters[i];

            for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
            {
                globalStaticPointsList[i + globalStaticOffsetNonUnique] = pointProcAddressing[staticControlPointLabels[i]];
            }

            for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
            {
                globalFixedPointsList[i + globalFixedOffsetNonUnique] = pointProcAddressing[fixedControlPointLabels[i]];
            }

            for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
            {
                globalMovingPointsList[i + globalMovingOffsetNonUnique] = pointProcAddressing[movingControlPointLabels[i]];
                globalMovingPointsPatchIds[i + globalMovingOffsetNonUnique] = movingControlPointPatchIds[i];
                globalMovingPointsIndices[i + globalMovingOffsetNonUnique] = movingControlPointIndices[i];
            }

            reduce( globalStaticPointsList, sumOp<labelList>() );
            reduce( globalFixedPointsList, sumOp<labelList>() );

            if ( not faceCellCenters )
            {
                reduce( globalMovingPointsList, sumOp<labelList>() );
                reduce( globalMovingPointsPatchIds, sumOp<labelList>() );
                reduce( globalMovingPointsIndices, sumOp<labelList>() );
            }

            // Construct a list of static control points which indicate whether
            // should be included or not.

            globalStaticPointsListEnabled.resize( nbStaticFaceCenters );
            globalStaticPointsListEnabled = 0;
            globalFixedPointsListEnabled.resize( nbFixedFaceCenters );
            globalFixedPointsListEnabled = 0;
            globalMovingPointsListEnabled.resize( nbMovingFaceCenters );
            globalMovingPointsListEnabled = 0;
            forAll( globalStaticPointsList, i )
            {
                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalStaticPointsList[i] );

                if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                {
                    staticControlGlobalPointLabels.push_back( globalStaticPointsList[i] );
                    globalStaticPointsListEnabled[i] = 1;
                }
            }

            forAll( globalFixedPointsList, i )
            {
                // Only add the static vertex point if it's not already added to the list
                std::vector<int>::iterator it;

                it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalFixedPointsList[i] );

                if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                {
                    it = std::find( fixedControlGlobalPointLabels.begin(), fixedControlGlobalPointLabels.end(), globalFixedPointsList[i] );

                    if ( it == fixedControlGlobalPointLabels.end() || fixedControlGlobalPointLabels.size() == 0 )
                    {
                        fixedControlGlobalPointLabels.push_back( globalFixedPointsList[i] );
                        globalFixedPointsListEnabled[i] = 1;
                    }
                }
            }

            if ( not faceCellCenters )
            {
                forAll( movingPatchIDs, patchI )
                {
                    const labelList & meshPoints = mesh().boundaryMesh()[movingPatchIDs[patchI]].meshPoints();
                    globalMovingPointsLabelList[movingPatchIDs[patchI]] = labelList( meshPoints.size(), 0 );
                }

                forAll( globalMovingPointsList, i )
                {
                    // Only add the static vertex point if it's not already added to the list
                    std::vector<int>::iterator it;

                    it = std::find( staticControlGlobalPointLabels.begin(), staticControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                    if ( it == staticControlGlobalPointLabels.end() || staticControlGlobalPointLabels.size() == 0 )
                    {
                        it = std::find( fixedControlGlobalPointLabels.begin(), fixedControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                        if ( it == fixedControlGlobalPointLabels.end() || fixedControlGlobalPointLabels.size() == 0 )
                        {
                            it = std::find( movingControlGlobalPointLabels.begin(), movingControlGlobalPointLabels.end(), globalMovingPointsList[i] );

                            if ( it == movingControlGlobalPointLabels.end() || movingControlGlobalPointLabels.size() == 0 )
                            {
                                movingControlGlobalPointLabels.push_back( globalMovingPointsList[i] );
                                globalMovingPointsListEnabled[i] = 1;

                                if ( i < static_cast<int>( movingControlPointLabels.size() ) + globalMovingOffsetNonUnique
                                    && i >= globalMovingOffsetNonUnique )
                                {
                                    label patchId = globalMovingPointsPatchIds[i];
                                    label index = globalMovingPointsIndices[i];
                                    globalMovingPointsLabelList[patchId][index] = 1;
                                }
                            }
                        }
                    }
                }
            }

            // Count the number of local unique static points
            nbStaticFaceCenters = 0;

            for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
            {
                if ( globalStaticPointsListEnabled[i + globalStaticOffsetNonUnique] == 1 )
                    nbStaticFaceCenters++;
            }

            nbFixedFaceCenters = 0;

            for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
            {
                if ( globalFixedPointsListEnabled[i + globalFixedOffsetNonUnique] == 1 )
                    nbFixedFaceCenters++;
            }

            if ( not faceCellCenters )
            {
                nbMovingFaceCenters = 0;

                for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
                {
                    if ( globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1 )
                        nbMovingFaceCenters++;
                }
            }
        }

        // Calculate sum of all faces on each processor
        nbGlobalStaticFaceCenters = 0;
        nbGlobalFixedFaceCenters = 0;
        nbGlobalMovingFaceCenters = 0;
        nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
        nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
        nbGlobalFixedFaceCenters[Pstream::myProcNo()] = nbFixedFaceCenters;
        nbGlobalFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters + nbStaticFaceCenters + nbFixedFaceCenters;

        reduce( nbGlobalMovingFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalStaticFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalFixedFaceCenters, sumOp<labelList>() );
        reduce( nbGlobalFaceCenters, sumOp<labelList>() );
    }

    nbMovingFaceCenters = sum( nbGlobalMovingFaceCenters );
    nbStaticFaceCenters = sum( nbGlobalStaticFaceCenters );
    nbFixedFaceCenters = sum( nbGlobalFixedFaceCenters );
    nbFaceCenters = sum( nbGlobalFaceCenters );

    // Determine the offset taking into account multiple processors

    int globalMovingOffset = 0;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalMovingOffset += nbGlobalMovingFaceCenters[i];

    int globalStaticOffset = nbMovingFaceCenters;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalStaticOffset += nbGlobalStaticFaceCenters[i];

    int globalFixedOffset = nbMovingFaceCenters + nbStaticFaceCenters;

    for ( int i = 0; i < Pstream::myProcNo(); i++ )
        globalFixedOffset += nbGlobalFixedFaceCenters[i];

    if ( !rbf->rbf->computed )
    {
        rbf::matrix positions( nbFaceCenters, mesh().nGeometricD() );
        positions.setZero();
        int index = 0;

        const Foam::pointField & points = mesh().points();

        vectorField positionsField( positions.rows(), vector::zero );

        if ( faceCellCenters )
        {
            int offset = 0;

            forAll( movingPatchIDs, i )
            {
                const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

                // Set the positions for patch i
                forAll( faceCentres, j )
                {
                    positionsField[j + offset + globalMovingOffset] = faceCentres[j];
                }

                offset += faceCentres.size();
            }
        }

        if ( not faceCellCenters )
        {
            for ( unsigned int i = 0; i < movingControlPointLabels.size(); i++ )
            {
                if ( globalMovingPointsListEnabled[i + globalMovingOffsetNonUnique] == 1 )
                {
                    positionsField[index + globalMovingOffset] = points[movingControlPointLabels[i]];
                    index++;
                }
            }

            assert( index == nbGlobalMovingFaceCenters[Pstream::myProcNo()] );
        }

        index = 0;

        for ( unsigned int i = 0; i < staticControlPointLabels.size(); i++ )
        {
            if ( globalStaticPointsListEnabled[i + globalStaticOffsetNonUnique] == 1 )
            {
                positionsField[index + globalStaticOffset] = points[staticControlPointLabels[i]];
                index++;
            }
        }

        assert( index == nbGlobalStaticFaceCenters[Pstream::myProcNo()] );

        index = 0;

        for ( unsigned int i = 0; i < fixedControlPointLabels.size(); i++ )
        {
            if ( globalFixedPointsListEnabled[i + globalFixedOffsetNonUnique] == 1 )
            {
                positionsField[index + globalFixedOffset] = points[fixedControlPointLabels[i]];
                index++;
            }
        }

        assert( index == nbGlobalFixedFaceCenters[Pstream::myProcNo()] );

        reduce( positionsField, sumOp<vectorField>() );

        // Copy the FOAM vector field to an Eigen matrix
        for ( int i = 0; i < positions.rows(); i++ )
            for ( int j = 0; j < positions.cols(); j++ )
                positions( i, j ) = positionsField[i][j];

        /*
         * Step 2: Build a matrix with the positions of every vertex in the local mesh.
         * This is only local information and does not need to be communicated to other
         * processors.
         */

        // Determine the number of points by using the 2d corrector
        nbPoints = 0;
        forAll( points, i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
                nbPoints++;
        }

        rbf::matrix positionsInterpolation( nbPoints, positions.cols() );

        index = 0;
        forAll( points, i )
        {
            if ( twoDCorrector.marker()[i] == 0 )
            {
                for ( int j = 0; j < positionsInterpolation.cols(); j++ )
                    positionsInterpolation( index, j ) = points[i][j];

                index++;
            }
        }

        rbf->compute( positions, positionsInterpolation );

        rbf->setNbMovingAndStaticFaceCenters( nbMovingFaceCenters, nbStaticFaceCenters + nbFixedFaceCenters );
    }

    /*
     * Step 3: Build a matrix with the displacement/motion of the face center
     * positions of the static patches and the moving patches.
     * The motion needs to be communicated to every process at every mesh deformation.
     * This is considered to be the most expensive step with regards to parallel
     * scalability of the overall algorithm.
     */

    rbf::matrix values( nbFaceCenters, mesh().nGeometricD() );
    values.setZero();

    vectorField valuesField( values.rows(), vector::zero );

    if ( faceCellCenters )
    {
        int offset = 0;

        forAll( movingPatchIDs, i )
        {
            const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

            forAll( motionCenters[movingPatchIDs[i]], j )
            {
                valuesField[j + offset + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
            }

            offset += faceCentres.size();
        }
    }

    if ( not faceCellCenters )
    {
        int index = 0;

        forAll( movingPatchIDs, i )
        {
            forAll( globalMovingPointsLabelList[movingPatchIDs[i]], j )
            {
                if ( globalMovingPointsLabelList[movingPatchIDs[i]][j] == 1 )
                {
                    valuesField[index + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
                    index++;
                }
            }
        }

        assert( index == nbGlobalMovingFaceCenters[Pstream::myProcNo()] );
    }

    reduce( valuesField, sumOp<vectorField>() );

    // Copy the FOAM vector field to an Eigen matrix
    for ( int i = 0; i < values.rows(); i++ )
        for ( int j = 0; j < values.cols(); j++ )
            values( i, j ) = valuesField[i][j];

    /*
     * Step 4: Perform the interpolation from the face centers to the complete mesh
     */

    rbf::matrix valuesInterpolation( nbPoints, values.cols() );
    valuesInterpolation.setZero();

    if ( cpu )
        rbf->rbf->computed = false;

    rbf->interpolate( values, valuesInterpolation );

    // Apply the 2d correction

    vectorField valuesInterpolationField( mesh().points().size(), Foam::vector::zero );
    int index = 0;
    forAll( valuesInterpolationField, i )
    {
        if ( twoDCorrector.marker()[i] == 0 )
        {
            for ( int j = 0; j < valuesInterpolation.cols(); j++ )
                valuesInterpolationField[i][j] = valuesInterpolation( index, j );

            index++;
        }
    }

    twoDCorrector.setShadowSide( valuesInterpolationField );

    /*
     * Step 5: Correct the mesh vertices of the fixed patches. Set these displacements to zero.
     */

    // Loop over all the patches, and set the fixed patches to zero.

    forAll( mesh().boundaryMesh(), i )
    {
        const labelList & meshPoints = mesh().boundaryMesh()[i].meshPoints();

        bool isFixedPatch = false;
        forAll( fixedPatchIDs, j )
        {
            if ( i == fixedPatchIDs[j] )
                isFixedPatch = true;
        }

        if ( isFixedPatch )
        {
            for ( int j = 0; j < meshPoints.size(); j++ )
                valuesInterpolationField[meshPoints[j]] = Foam::vector::zero;
        }
    }

    /*
     * Step 6: Set the motion of the mesh vertices
     */

    assert( newPoints.size() == valuesInterpolationField.size() );

    newPoints = valuesInterpolationField;
}
