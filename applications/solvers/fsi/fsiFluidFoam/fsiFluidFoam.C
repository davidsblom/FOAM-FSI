
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <yaml-cpp/yaml.h>

#include "version.H"
#include "CoupledFluidSolver.H"
#include "PreciceFluidSolver.H"
#include "FluidSolver.H"
#include "CompressibleFluidSolver.H"
#include "SteadyStateFluidSolver.H"

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

    Info << nl << "FOAM-FSI build: " << word( FOAMFSIbuild ) << nl << endl;

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );

    assert( config["fluid-solver"] );

    std::string fluidSolver = config["fluid-solver"].as<std::string>();

    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "compressible-solver" || fluidSolver == "steady-state-pimple-solver" );

    std::shared_ptr<foamFluidSolver> fluid;

    if ( fluidSolver == "coupled-pressure-velocity-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "pimple-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "steady-state-pimple-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new SteadyStateFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "compressible-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    assert( fluid );

    PreciceFluidSolver solver( fluid );

    solver.run();

    Info << "End\n" << endl;

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    return (0);
}

// ************************************************************************* //
