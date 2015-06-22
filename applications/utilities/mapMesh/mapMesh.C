
#include <iostream>
#include <memory>
#include <yaml-cpp/yaml.h>
#include "fvCFD.H"
#include "SolidSolver.H"
#include "CoupledFluidSolver.H"
#include "FluidSolver.H"
#include "RBFInterpolation.H"
#include "RBFCoarsening.H"
#include "TPSFunction.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"

class MapMeshSolidSolver : public SolidSolver
{
public:

    MapMeshSolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime
        )
        :
        SolidSolver( name, args, runTime )
    {}

    void getDisplacementLocal( matrix & displacement )
    {
        displacement.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

        int offset = 0;

        forAll( movingPatchIDs, patchI )
        {
            int size = U.boundaryField()[movingPatchIDs[patchI]].size();

            for ( int i = 0; i < size; i++ )
                for ( int j = 0; j < displacement.cols(); j++ )
                    displacement( i + offset, j ) = U.boundaryField()[movingPatchIDs[patchI]][i][j];

            offset += size;
        }
    }

    void getWritePositionsLocal( matrix & writePositions )
    {
        writePositions.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

        int offset = 0;

        forAll( movingPatchIDs, patchI )
        {
            const vectorField faceCentres( mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres() );

            for ( int i = 0; i < faceCentres.size(); i++ )
                for ( int j = 0; j < writePositions.cols(); j++ )
                    writePositions( i + offset, j ) = faceCentres[i][j];

            offset += faceCentres.size();
        }
    }
};

std::shared_ptr<rbf::RBFInterpolation> createRBFInterpolator(
    std::string interpolationFunction,
    double radius,
    bool cpu
    )
{
    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;

    if ( interpolationFunction == "thin-plate-spline" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );

    if ( interpolationFunction == "wendland-c0" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC0Function( radius ) );

    if ( interpolationFunction == "wendland-c2" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function( radius ) );

    if ( interpolationFunction == "wendland-c4" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC4Function( radius ) );

    if ( interpolationFunction == "wendland-c6" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC6Function( radius ) );

    assert( rbfFunction );

    bool polynomialTerm = true;
    return std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );
}

int main(
    int argc,
    char * argv[]
    )
{
    std::shared_ptr<argList> args( new argList( argc, argv ) );

    if ( !args->checkRootCase() )
    {
        FatalError.exit();
    }

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    assert( Pstream::nProcs() == 1 );

    std::shared_ptr<MapMeshSolidSolver> solid( new MapMeshSolidSolver( "solid", args, runTime ) );

    std::shared_ptr<foamFluidSolver> fluid;
    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";
    YAML::Node config = YAML::LoadFile( filename );
    YAML::Node configInterpolation( config["interpolation"] );

    std::string fluidSolver = config["fluid-solver"].as<std::string>();
    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" );

    assert( configInterpolation["coarsening"] );
    assert( configInterpolation["coarsening"]["enabled"] );
    assert( configInterpolation["radial-basis-function"] );
    assert( configInterpolation["radial-basis-function"]["function"] );
    bool coarsening = configInterpolation["coarsening"]["enabled"].as<bool>();
    std::string interpolationFunction = configInterpolation["radial-basis-function"]["function"].as<std::string>();
    double coarseningTol = 1.0e-5;
    int coarseningMinPoints = 1;
    int coarseningMaxPoints = 10;
    double radius = 1;
    bool livePointSelection = false;
    double tolLivePointSelection = 1.0e-5;
    bool cpu = false;

    assert( interpolationFunction == "thin-plate-spline" || interpolationFunction == "wendland-c0" || interpolationFunction == "wendland-c2" || interpolationFunction == "wendland-c4" || interpolationFunction == "wendland-c6" );

    if ( configInterpolation["radial-basis-function"]["cpu"] )
    {
        cpu = configInterpolation["radial-basis-function"]["cpu"].as<bool>();
    }

    if ( interpolationFunction != "thin-plate-spline" )
    {
        assert( configInterpolation["radial-basis-function"]["radius"] );
        radius = configInterpolation["radial-basis-function"]["radius"].as<double>();
    }

    if ( coarsening )
    {
        assert( configInterpolation["coarsening"]["tol"] );
        assert( configInterpolation["coarsening"]["min-points"] );
        assert( configInterpolation["coarsening"]["max-points"] );
        assert( configInterpolation["coarsening"]["live-point-selection"] );

        coarseningTol = configInterpolation["coarsening"]["tol"].as<double>();
        coarseningMinPoints = configInterpolation["coarsening"]["min-points"].as<int>();
        coarseningMaxPoints = configInterpolation["coarsening"]["max-points"].as<int>();
        livePointSelection = configInterpolation["coarsening"]["live-point-selection"].as<bool>();

        if ( livePointSelection )
        {
            assert( configInterpolation["coarsening"]["tol-live-point-selection"] );
            tolLivePointSelection = configInterpolation["coarsening"]["tol-live-point-selection"].as<double>();
        }
    }

    if ( fluidSolver == "coupled-pressure-velocity-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "pimple-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
    /*
     * 1. Interpolate displacement from solid interface to fluid interface
     * 2. Mesh deformation fluid mesh
     */

    // Step 1: Interpolate displacement from solid interface to fluid interface

    // Settings coarsening;
    bool livePointSelectionSumValues = false;
    bool exportTxt = false;

    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
    std::shared_ptr<rbf::RBFCoarsening> rbfSolidToFluidInterface;

    rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
    rbfSolidToFluidInterface = std::shared_ptr<rbf::RBFCoarsening>( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, livePointSelectionSumValues, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, exportTxt ) );

    Eigen::MatrixXd positions, positionsInterpolation, values, valuesInterpolation;

    solid->getWritePositionsLocal( positions );
    fluid->getReadPositions( positionsInterpolation );

    rbfSolidToFluidInterface->compute( positions, positionsInterpolation );

    // valuesInterpolation = displacement of fluid interface
    // values = displacement of solid interface

    solid->getDisplacementLocal( values );

    rbfSolidToFluidInterface->interpolate( values, valuesInterpolation );

    // Step 2
    // Fluid mesh deformation

    fluid->setDisplacementLocal( valuesInterpolation );
    fluid->moveMesh();

    runTime->loop();
    runTime->write();
}
