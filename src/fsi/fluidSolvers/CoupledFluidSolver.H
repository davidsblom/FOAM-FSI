
/*
 * Class
 *   CoupledFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef CoupledFluidSolver_H
#define CoupledFluidSolver_H

#include <memory>

#include "foamFluidSolver.H"
#include "RASModel.H"
#include "singlePhaseTransportModel.H"
#include "fvBlockMatrix.H"

using namespace fsi;
using std::shared_ptr;

class CoupledFluidSolver : public foamFluidSolver
{
    public:
        CoupledFluidSolver(
            string name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        virtual ~CoupledFluidSolver();

        virtual void getAcousticsDensityLocal( matrix & data );

        virtual void getAcousticsVelocityLocal( matrix & data );

        virtual void getAcousticsPressureLocal( matrix & data );

        virtual void getAcousticsDensityGradientLocal( matrix & data );

        virtual void getAcousticsVelocityGradientLocal( matrix & data );

        virtual void getAcousticsPressureGradientLocal( matrix & data );

        virtual void setAcousticsData(
            const matrix & pressure,
            const matrix & temperature,
            const matrix & velocity
            );

        virtual void getTractionLocal( matrix & traction );

        virtual void getWritePositionsLocalAcoustics( matrix & writePositions );

        virtual void initTimeStep();

        virtual bool isRunning();

        virtual void resetSolution();

        virtual void solve();

    protected:
        void checkTotalVolume();

        void continuityErrs();

        void courantNo();

        void readBlockSolverControls();

        // Fields
        volScalarField p;
        volVectorField U;
        surfaceScalarField phi;
        volVector4Field Up;

        // Turbulence
        singlePhaseTransportModel laminarTransport;
        autoPtr<incompressible::turbulenceModel> turbulence;

        // Continuity errors
        scalar sumLocalContErr;
        scalar globalContErr;
        scalar cumulativeContErr;

        // Pressure reference cell
        label pRefCell;
        scalar pRefValue;

        // Dictionaries
        IOdictionary transportProperties;

        // Dimensioned scalars
        dimensionedScalar nu;
        dimensionedScalar rho;

        // Convergence tolerance fluid domain
        scalar convergenceTolerance;
        label nOuterCorr;

        // Total volume
        scalar totalVolume;

        // Courant no
        scalar CoNum;
        scalar meanCoNum;
        scalar velMag;
};

#endif
