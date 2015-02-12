
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */


#ifndef foamSolidSolver_H
#define foamSolidSolver_H

#include <memory>

#include "BaseMultiLevelSolver.H"
#include "fvCFD.H"
#include "solidTractionFvPatchVectorField.H"
#include "leastSquaresVolPointInterpolation.H"
#include <time.h>

using namespace fsi;
using std::shared_ptr;

class foamSolidSolver : public BaseMultiLevelSolver
{
public:

  foamSolidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    );

  ~foamSolidSolver();

  virtual void finalizeTimeStep();

  virtual void getReadPositions( matrix & readPositions );

  virtual void getWritePositions( matrix & writePositions );

  virtual void solve() = 0;

  virtual void solve(
    const matrix & input,
    matrix & output
    );

  int getInterfacePointsSize();

  int getInterfaceCentersSize();

  int getInterfaceCentersSizeLocal();

  void getDisplacementLocal( matrix & displacement );

  void getReadPositionsLocal( matrix & readPositions );

  void getWritePositionsLocal( matrix & writePositions );

  void setTractionLocal( const matrix & traction );

  string name;
  std::shared_ptr<argList> args;
  std::shared_ptr<Time> runTime;

  Foam::fvMesh mesh;

protected:

  void readCouplingProperties();

  // Fields
  volVectorField U;

  // Dictionaries
  IOdictionary couplingProperties;

  // Solid patch data
  word solidPatchName;
  label solidPatchID;

  solidTractionFvPatchVectorField * tForce;

  labelList nPoints;
  labelList nGlobalPoints;
  labelList nGlobalCenters;
  labelList globalPointsUnique;
  labelList globalPointsNonUnique;

  float totalRunTime;
  int totalNbIterations;
};

#endif
