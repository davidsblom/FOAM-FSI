
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include "SolidSolver.H"
#include "TPSFunction.H"

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

  std::shared_ptr<rbf::RBFFunctionInterface>rbfFunction( new rbf::TPSFunction() );
  std::shared_ptr<rbf::RBFInterpolation>rbfInterpolator( new rbf::RBFInterpolation( rbfFunction ) );
  std::shared_ptr<rbf::RBFCoarsening> interpolator( new rbf::RBFCoarsening( rbfInterpolator, false, 1.0e-3, 10, 20 ) );

  std::shared_ptr<foamSolidSolver> solid( new SolidSolver( fvMesh::defaultRegion, args, runTime, interpolator ) );

  solid->run();

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
