
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include <memory>

#include "foamSolidSolver.H"
#include "constitutiveModel.H"

using namespace fsi;

class ElasticSolidSolver : public foamSolidSolver
{
    public:
        ElasticSolidSolver(
            string name,
            std::shared_ptr<argList> args,
            std::shared_ptr<Time> runTime
            );

        virtual ~ElasticSolidSolver();

        virtual void finalizeTimeStep();

        virtual void initTimeStep();

        virtual bool isRunning();

        virtual void resetSolution();

        virtual void solve();

    private:
        // Fields
        volTensorField gradU;
        surfaceVectorField snGradU;
        volVectorField V;
        volTensorField gradV;
        surfaceVectorField snGradV;
        volSymmTensorField epsilon;
        volSymmTensorField sigma;
        volVectorField divSigmaExp;
        constitutiveModel rheology;
        volScalarField rho;
        volScalarField mu;
        volScalarField lambda;
        surfaceScalarField muf;
        surfaceScalarField lambdaf;
        surfaceVectorField n;

        // Controls
        int nCorr;
        scalar convergenceTolerance;
        word divSigmaExpMethod;
};
