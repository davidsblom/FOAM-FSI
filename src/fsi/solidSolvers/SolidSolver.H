
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SolidSolver_H
#define SolidSolver_H

#include <memory>

#include "foamSolidSolver.H"
#include "constitutiveModel.H"
#include "RBFCoarsening.H"
#include "meshToMesh.H"

using namespace fsi;
using namespace rbf;

class SolidSolver : public foamSolidSolver
{
public:

    SolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime
        );

    SolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime,
        std::shared_ptr<rbf::RBFCoarsening> interpolator
        );

    virtual ~SolidSolver();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void resetSolution();

    virtual void solve();

    virtual bool interpolateVolField( std::shared_ptr<BaseMultiLevelSolver> fineModel );

protected:

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
    int minIter;
    int maxIter;
    scalar absoluteTolerance;
    scalar relativeTolerance;

    std::shared_ptr<rbf::RBFCoarsening> interpolator;
};

#endif
