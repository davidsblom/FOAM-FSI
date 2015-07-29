
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

    int nbTimeSteps = mesh.time().endTime().value() / mesh.time().deltaT().value();
    scalar dt = mesh.time().deltaT().value();
    fsi::vector q0 ( 2 );
    q0 << 1, 0;
    scalar amplitude = 0;
    scalar frequency = 1;
    scalar m = 1;
    scalar k = 100;
    oscillator = std::shared_ptr<Oscillator>( new Oscillator( nbTimeSteps, dt, q0, amplitude, frequency, m, k ) );
}

RBFMeshRigidMotionSolver::~RBFMeshRigidMotionSolver()
{}

Foam::vector RBFMeshRigidMotionSolver::calcVelocity()
{
    assert( timeIntegrationScheme );

    scalar t = mesh().time().value();
    scalar dt = mesh().time().deltaT().value();

    fsi::vector rhs( 2 ), result( 2 ), f( 2 ), qold;

    if ( not corrector && k == 0 )
        timeIntegrationScheme->setOldSolution( mesh().time().timeIndex(), oscillator->sol );

    if ( not corrector && k == 0 )
    {
        oscillator->evaluateFunction( 0, oscillator->sol, t - dt, f );
        timeIntegrationScheme->setFunction( k-1, f, oscillator->sol );
    }

    timeIntegrationScheme->getSourceTerm( corrector, k, dt, rhs, qold );

    oscillator->implicitSolve( corrector, k, k, t, dt, qold, rhs, f, result );

    timeIntegrationScheme->setFunction( k, f, result );

    scalar disp = oscillator->sol(0);
    // change reference frame
    disp *= -1;
    disp += 1;
    disp *= translationAmplitude;

    Foam::vector transformation = Foam::vector::zero;
    transformation[1] = disp - oldTransformation[1];

    Info << "dudt(1) = " << oscillator->sol(0) << endl;
    Info << "dudt(0) = " << oscillator->sol(1) << endl;

    told = t;
    oldTransformation[1] = disp;

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
