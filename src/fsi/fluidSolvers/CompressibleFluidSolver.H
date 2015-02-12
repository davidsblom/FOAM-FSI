
/*
 * Class
 *   CompressibleFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef CompressibleFluidSolver_H
#define CompressibleFluidSolver_H

#include "foamFluidSolver.H"
#include "basicPsiThermo.H"
#include "compressible/RAS/lnInclude/RASModel.H"
#include "fvBlockMatrix.H"

class CompressibleFluidSolver : public foamFluidSolver
{
public:

  CompressibleFluidSolver(
    string name,
    shared_ptr<argList> args,
    shared_ptr<Time> runTime
    );

  virtual ~CompressibleFluidSolver();

  virtual void getAcousticsDensityLocal( matrix & data );

  virtual void getAcousticsVelocityLocal( matrix & data );

  virtual void getAcousticsPressureLocal( matrix & data );

  virtual void getTractionLocal( matrix & traction );

  virtual void getWritePositionsLocalAcoustics( matrix & writePositions );

  virtual void initTimeStep();

  virtual bool isRunning();

  virtual void resetSolution();

  virtual void solve();

protected:

  void continuityErrs();

  // Fields
  autoPtr<basicPsiThermo> pThermo;
  basicPsiThermo & thermo;

  volScalarField & p;
  volScalarField & h;
  const volScalarField & T;
  const volScalarField & psi;
  volScalarField rho;
  volVectorField U;
  surfaceScalarField phi;
  autoPtr<compressible::RASModel> turbulence;
  volVector4Field Up;
  volScalarField DpDt;

  // Continuity errors
  scalar cumulativeContErr;

  // Convergence tolerance fluid domain
  scalar convergenceTolerance;
  label nOuterCorr;

  // Acoustics coupling patch

  word acousticsPatchName;
  label acousticsPatchID;
};

#endif
