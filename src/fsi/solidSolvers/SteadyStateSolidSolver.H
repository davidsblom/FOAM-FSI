
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SteadyStateSolidSolver_H
#define SteadyStateSolidSolver_H

#include <memory>

#include "foamSolidSolver.H"
#include "constitutiveModel.H"
#include "RBFCoarsening.H"
#include "meshToMesh.H"

using namespace fsi;
using namespace rbf;

class SteadyStateSolidSolver : public foamSolidSolver {
    public:
        SteadyStateSolidSolver(const std::string & name,
            const std::shared_ptr<argList> & args,
            const std::shared_ptr<Time> & runTime
            );

        SteadyStateSolidSolver(const std::string & name,
            const std::shared_ptr<argList> & args,
            const std::shared_ptr<Time> & runTime,
            const std::shared_ptr<rbf::RBFCoarsening> & interpolator
            );

        virtual ~SteadyStateSolidSolver();

        virtual void initTimeStep();

        virtual bool isRunning();

        virtual void resetSolution();

        virtual void solve();

        virtual bool interpolateVolField(std::shared_ptr<BaseMultiLevelSolver> fineModel);

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
