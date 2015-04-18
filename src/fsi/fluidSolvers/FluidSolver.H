
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

    virtual void getAcousticsDensityLocal( matrix & data );

    virtual void getAcousticsVelocityLocal( matrix & data );

    virtual void getAcousticsPressureLocal( matrix & data );

    virtual void getTractionLocal( matrix & traction );

    virtual void getWritePositionsLocalAcoustics( matrix & writePositions );

    virtual void initTimeStep();

    virtual bool isRunning();

    virtual void resetSolution();

    virtual void solve();

private:

    void checkTimeDiscretisationScheme();

    void continuityErrs();

    void courantNo();

    double evaluateMomentumResidual();

    // Dictionaries
    IOdictionary transportProperties;

    // Dimensioned scalars
    dimensionedScalar nu;
    dimensionedScalar rho;

    // Fields
    volScalarField p;
    volVectorField U;
    surfaceScalarField phi;
    IOobject UfHeader;
    surfaceVectorField Uf;
    volScalarField AU;
    volVectorField HU;

    // PIMPLE controls
    int nOuterCorr;
    int nCorr;
    int nNonOrthCorr;
    scalar absoluteTolerance;
    scalar relativeTolerance;

    // Continuity errors
    scalar sumLocalContErr;
    scalar globalContErr;
    scalar cumulativeContErr;

    label pRefCell;
    scalar pRefValue;
    singlePhaseTransportModel laminarTransport;
    autoPtr<incompressible::turbulenceModel> turbulence;

    // Courant no
    scalar CoNum;
    scalar meanCoNum;
    scalar velMag;
};

#endif
