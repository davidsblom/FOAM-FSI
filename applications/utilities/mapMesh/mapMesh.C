
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

    //std::shared_ptr<foamFluidSolver> fluid( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
// Capability to use different fluid solvers (copuled-pressure-velocity-solver or pimple-solver)
    std::shared_ptr<foamFluidSolver> fluid;
    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";
    YAML::Node config = YAML::LoadFile( filename );
    std::string fluidSolver = config["fluid-solver"].as<std::string>();
    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" );


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
    bool enabled = true;
    bool livePointSelection = true;
    bool livePointSelectionSumValues = false;
    double tol = 1.0e-4;
    double tolLivePointSelection = 1.0e-4;
    int coarseningMinPoints = 5;
    int coarseningMaxPoints = 500;
    bool exportTxt = false;
    bool polynomialTerm = false;
    bool cpu = true;

    std::shared_ptr<TPSFunction> tpsFunction( new TPSFunction() );
    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator( new rbf::RBFInterpolation( tpsFunction, polynomialTerm, cpu ) );
    std::shared_ptr<rbf::RBFCoarsening> rbfSolidToFluidInterface( new rbf::RBFCoarsening( rbfInterpolator, enabled, livePointSelection, livePointSelectionSumValues, tol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, exportTxt ) );

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
