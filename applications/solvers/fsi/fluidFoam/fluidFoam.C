
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <yaml-cpp/yaml.h>

#include "FluidSolver.H"
#include "CoupledFluidSolver.H"
#include "CompressibleFluidSolver.H"
#include "SDCFluidSolver.H"
#include "SDC.H"

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

  string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";

  YAML::Node config = YAML::LoadFile( filename );

  assert( config["fluid-solver"] );

  std::string fluidSolver = config["fluid-solver"].as<std::string>();

  assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "compressible-solver" || fluidSolver == "sdc-pimple-solver" );

  std::shared_ptr<foamFluidSolver> fluid;
  std::shared_ptr<sdc::SDC> sdc;

  if ( fluidSolver == "coupled-pressure-velocity-solver" )
    fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

  if ( fluidSolver == "pimple-solver" )
    fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

  if ( fluidSolver == "compressible-solver" )
    fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

  if ( fluidSolver == "sdc-pimple-solver" )
  {
    YAML::Node sdcConfig( config["sdc"] );
    assert( sdcConfig["convergence-tolerance"] );
    assert( sdcConfig["number-of-points"] );

    int n = sdcConfig["number-of-points"].as<int>();
    double tol = sdcConfig["convergence-tolerance"].as<double>();

    std::shared_ptr<sdc::SDCSolver> solver( new SDCFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
    sdc = std::shared_ptr<sdc::SDC> ( new sdc::SDC( solver, n, tol ) );
  }

  assert( fluid || sdc );

  if ( fluid )
    fluid->run();

  if ( sdc )
    sdc->run();

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
