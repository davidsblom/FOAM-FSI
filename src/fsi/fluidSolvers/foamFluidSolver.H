
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */


#ifndef foamFluidSolver_H
#define foamFluidSolver_H

#include <memory>

#include "BaseMultiLevelSolver.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"

using namespace fsi;
using std::shared_ptr;

class foamFluidSolver : public BaseMultiLevelSolver
{
public:

  foamFluidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    );

  ~foamFluidSolver();

  virtual void finalizeTimeStep();

  virtual void getReadPositions( matrix & readPositions );

  virtual void getWritePositions( matrix & writePositions );

  virtual void getTractionLocal( matrix & traction ) = 0;

  virtual void solve() = 0;

  virtual void solve(
    const matrix & input,
    matrix & output
    );

  int getInterfaceCentersSize();

  int getInterfacePointsSize();

  int getInterfacePointsSizeLocal();

  void setDisplacementLocal( const matrix & displacement );

  void getReadPositionsLocal( matrix & readPositions );

  void getWritePositionsLocal( matrix & writePositions );

  string name;
  shared_ptr<argList> args;
  shared_ptr<Time> runTime;

protected:

  virtual void moveMesh( vectorField motion ) = 0;

  // Mesh
  autoPtr<dynamicFvMesh> meshPtr;
  dynamicFvMesh & mesh;

  // Dictionaries
  IOdictionary couplingProperties;

  // FSI interface displacement
  word fluidPatchName;
  label fluidPatchID;
  vectorField fluidPatchPointsDispl;
  vectorField fluidPatchPointsDisplOld;

  // Parallel FSI coupling variables
  labelList nPoints;
  labelList nGlobalPoints;
  labelList nGlobalCenters;
  labelList globalPointsUnique;
  labelList globalPointsNonUnique;
};

#endif
