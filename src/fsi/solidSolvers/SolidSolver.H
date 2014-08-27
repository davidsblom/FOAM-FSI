
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SolidSolver_H
#define SolidSolver_H

#include <memory>

#include "foamSolidSolver.H"
#include "constitutiveModel.H"

using namespace fsi;

class SolidSolver : public foamSolidSolver
{
public:

  SolidSolver(
    string name,
    std::shared_ptr<argList> args,
    std::shared_ptr<Time> runTime
    );

  virtual ~SolidSolver();

  virtual void initTimeStep();

  virtual bool isRunning();

  virtual void solve();

private:

  void calculateEpsilonSigma();

  void initialize();

  void readSolidMechanicsControls();

  // Fields
  volTensorField gradU;
  volSymmTensorField epsilon;
  volSymmTensorField sigma;

  // Constitutive model
  constitutiveModel rheology;
  volScalarField rho;
  volScalarField mu;
  volScalarField lambda;
  surfaceScalarField muf;
  surfaceScalarField lambdaf;
  surfaceVectorField n;

  // Solid mechanics controls
  int nCorr;
  scalar convergenceTolerance;

  // Dictionaries
  IOdictionary couplingProperties;
};

#endif
