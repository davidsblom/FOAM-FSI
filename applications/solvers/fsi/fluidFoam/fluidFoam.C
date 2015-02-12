
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include "CoupledFluidSolver.H"

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

  std::shared_ptr<CoupledFluidSolver> fluid;

  fluid = std::shared_ptr<CoupledFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

  fluid->run();

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
