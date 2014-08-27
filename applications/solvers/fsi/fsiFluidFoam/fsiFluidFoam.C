
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>

#include "CoupledFluidSolver.H"
#include "PreciceFluidSolver.H"

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

  std::shared_ptr<foamFluidSolver> fluid( new CoupledFluidSolver( dynamicFvMesh::defaultRegion, args, runTime ) );

  PreciceFluidSolver solver( fluid );

  solver.run();

  Info << "End\n" << endl;

  return (0);
}

// ************************************************************************* //
