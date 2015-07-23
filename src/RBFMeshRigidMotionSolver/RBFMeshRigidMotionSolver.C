
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
    RBFMeshMotionSolver( mesh, msData ),
    translationAmplitude(0),
    translationFrequency(0),
    translationDirection(Foam::vector::zero),
    rotationAmplitude(0),
    rotationFrequency(0),
    rotationOrigin(Foam::vector::zero),
    oldTransformation(Foam::vector::zero),
    told(0)
{
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

    forAll( motion[movingPatchIDs[0]], i )
    {
        motion[movingPatchIDs[0]][i] = transformation;
    }

    motionCenters = motion;

    RBFMeshMotionSolver::solve();
}
