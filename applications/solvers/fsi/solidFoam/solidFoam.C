
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include "SolidSolver.H"

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

  std::shared_ptr<foamSolidSolver> solid( new SolidSolver( fvMesh::defaultRegion, args, runTime ) );

  solid->run();

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
