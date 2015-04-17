
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef CoupledSolidSolver_H
#define CoupledSolidSolver_H

#include <memory>

#include "foamSolidSolver.H"
#include "solidPolyMesh.H"
#include "constitutiveModel.H"
#include "blockLduSolvers.H"
#include "volVectorNFields.H"
#include "nonLinearGeometry.H"
#include "processorFvPatchFields.H"
#include "emptyFvPatchFields.H"

using namespace fsi;

class CoupledSolidSolver : public foamSolidSolver
{
public:

    CoupledSolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime
        );

    virtual ~CoupledSolidSolver();

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void solve();

private:

    void createNonLinearGeometry();

    // Mesh
    solidPolyMesh solidMesh;
    pointMesh pMesh;
    leastSquaresVolPointInterpolation volToPointInterp;

    // Fields
    volTensorField gradU;
    surfaceTensorField gradUf;
    volSymmTensorField epsilon;
    volSymmTensorField sigma;
    surfaceVectorField n;

    // Constitutive model
    constitutiveModel rheology;
    const volScalarField rho;
    const volScalarField mu;
    const volScalarField lambda;
    surfaceScalarField muf;
    surfaceScalarField lambdaf;

    // Solid mechanics controls
    const dictionary & stressControl;
    int nCorr;
};

#endif
