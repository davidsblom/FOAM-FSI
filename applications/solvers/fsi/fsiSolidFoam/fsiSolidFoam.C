
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <yaml-cpp/yaml.h>

#include "version.H"
#include "PreciceSolidSolver.H"
#include "SolidSolver.H"
#include "ElasticSolidSolver.H"
#include "SteadyStateSolidSolver.H"

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

    assert( config["solid-solver"] );

    std::string solidSolver = config["solid-solver"].as<std::string>();

    assert( solidSolver == "nonlinear-elastic-solver" || solidSolver == "linear-elastic-solver" || solidSolver == "steady-state-nonlinear-elastic-solver" );

    std::shared_ptr<foamSolidSolver> solid;

    if ( solidSolver == "linear-elastic-solver" )
        solid = std::shared_ptr<foamSolidSolver> ( new ElasticSolidSolver( fvMesh::defaultRegion, args, runTime ) );

    if ( solidSolver == "nonlinear-elastic-solver" )
        solid = std::shared_ptr<foamSolidSolver>( new SolidSolver( fvMesh::defaultRegion, args, runTime ) );

    if ( solidSolver == "steady-state-nonlinear-elastic-solver" )
        solid = std::shared_ptr<foamSolidSolver>( new SteadyStateSolidSolver( fvMesh::defaultRegion, args, runTime ) );

    assert( solid );

    PreciceSolidSolver solver( solid );

    solver.run();

    Info << "End\n" << endl;

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    return (0);
}

// ************************************************************************* //
