
#include "fvCFD.H"
#include "Random.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(
  int argc,
  char * argv[]
  )
{
  argList::validArgs.append( "randomPercentage" );

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  scalar randomPercentage( readScalar( IStringStream( args.additionalArgs()[0] ) () ) );

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "Reading field U\n" << endl;
  volVectorField U
  (
    IOobject
    (
      "U",
      runTime.timeName(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    ),
    mesh
  );

  vector meanU = sum( U.internalField() ) / U.internalField().size();

  Random random( 1234567 );
  forAll( U.internalField(), iCell ){
    U.internalField()[iCell] += (randomPercentage / 100) * meanU * random.scalar01();
  }

  U.write();

  Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
       << "  ClockTime = " << runTime.elapsedClockTime() << " s"
       << nl << endl;

  Info << "End\n" << endl;

  return 0;
}

// ************************************************************************* //
