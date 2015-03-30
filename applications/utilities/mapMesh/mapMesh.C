
#include <iostream>
#include "fvCFD.H"
#include <memory>
#include "SolidSolver.H"
#include "CoupledFluidSolver.H"
#include "RBFInterpolation.H"
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

    std::shared_ptr<foamSolidSolver> solid( new SolidSolver( "solid", args, runTime ) );

    std::shared_ptr<foamFluidSolver> fluid ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    /*
    1. Interpolate displacement from solid interface to fluid interface
    2. Mesh deformation fluid mesh
    */

    // Step 1: Interpolate displacement from solid interface to fluid interface

    std::shared_ptr<TPSFunction> tpsFunction( new TPSFunction() );
    rbf::RBFInterpolation rbfSolidToFluidInterface ( tpsFunction );

    Eigen::MatrixXd positions, positionsInterpolation, values, valuesInterpolation;

    solid->getWritePositions( positions );
    fluid->getReadPositions( positionsInterpolation );

    rbfSolidToFluidInterface.compute(positions, positionsInterpolation);

    // valuesInterpolation = displacement of fluid interface
    // values = displacement of solid interface

    solid->getDisplacementLocal( values );

    rbfSolidToFluidInterface.interpolate(values, valuesInterpolation);

    std::cout << "step 1 done!" << std::endl;

    // Step 2
    // Fluid mesh deformation

    fluid->setDisplacementLocal( valuesInterpolation );
    fluid->moveMesh( fluid->movingPatchesDispl - fluid->movingPatchesDisplOld );

  }
