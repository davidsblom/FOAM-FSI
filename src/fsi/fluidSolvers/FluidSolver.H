
/*
 * Class
 *   FluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef FluidSolver_H
#define FluidSolver_H

#include <memory>

#include "foamFluidSolver.H"
#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

#include "RBFMotionSolverExt.H"
#include "faceTetPolyPatch.H"
#include "fixedValuePointPatchFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetFemMatrices.H"
#include "tetPointFields.H"
#include "tetPolyPatchInterpolation.H"

using namespace fsi;
using std::shared_ptr;

class FluidSolver : public foamFluidSolver
{
public:

  FluidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    );

  virtual ~FluidSolver();

  virtual void getTractionLocal( matrix & traction );

  virtual void initTimeStep();

  virtual bool isRunning();

  virtual void solve();

private:

  virtual void moveMesh( vectorField motion );

  void checkTimeDiscretisationScheme();

  void checkTotalVolume();

  void continuityErrs();

  void correctPhiField();

  void courantNo();

  void createFields();

  void createSf();

  void initContinuityErrs();

  void initialize();

  void initTotalVolume();

  void postProcessing();

  void readControls();

  void readCouplingProperties();

  void readPIMPLEControls();

  void readTimeControls();

  void setDeltaT();

  void updateSf();

  // Dictionaries
  IOdictionary transportProperties;
  dictionary pimple;

  // Dimensioned scalars
  dimensionedScalar nu;
  dimensionedScalar rhoFluid;

  // Fields
  volScalarField p;
  volVectorField U;
  surfaceScalarField phi;
  volScalarField AU;
  volVectorField HU;

  // PIMPLE controls
  int nOuterCorr;
  int nCorr;
  int nNonOrthCorr;
  bool momentumPredictor;
  bool transonic;

  // Continuity errors
  scalar sumLocalContErr;
  scalar globalContErr;
  scalar cumulativeContErr;

  label pRefCell;
  scalar pRefValue;
  singlePhaseTransportModel laminarTransport;
  autoPtr<incompressible::turbulenceModel> turbulence;

  IOobject SfHeader;
  surfaceVectorField Sf;

  // Controls
  bool correctPhi;
  scalar convergenceTolerance;

  // Time discretization scheme
  word ddtScheme;

  // Time controls
  Switch adjustTimeStep;
  scalar maxCo;
  scalar maxDeltaT;

  // Total volume
  scalar totalVolume;

  // Courant no
  scalar CoNum;
  scalar meanCoNum;
  scalar velMag;
};

#endif
