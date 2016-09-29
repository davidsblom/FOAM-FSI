
/*
 * Class
 *   SteadyStateFluidSolver
 *
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef SteadyStateFluidSolver_H
#define SteadyStateFluidSolver_H

#include <memory>

#include "foamFluidSolver.H"
#include "dynamicFvMesh.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

using namespace fsi;
using std::shared_ptr;

class SteadyStateFluidSolver : public foamFluidSolver
{
    public:
        SteadyStateFluidSolver(
            const std::string & name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        virtual ~SteadyStateFluidSolver();

        void getAcousticsDensityLocal( matrix & data ) override;

        void getAcousticsTemperatureLocal( matrix & data ) override;

        void getAcousticsVelocityLocal( matrix & data )  override;

        void getAcousticsPressureLocal( matrix & data ) override;

        void getAcousticsTemperatureGradientLocal( matrix & data ) override;

        void getAcousticsVelocityGradientLocal( matrix & data )  override;

        void getAcousticsPressureGradientLocal( matrix & data ) override;

        void getTractionLocal( matrix & traction ) override;

        void getWritePositionsLocalAcoustics( matrix & writePositions ) override;

        void initTimeStep() override;

        bool isRunning() override;

        void resetSolution() override;

        void setAcousticsTemperatureGradient( const matrix & temperatureGradient ) override;

        void setAcousticsPressureGradient( const matrix & pressureGradient ) override;

        void setAcousticsVelocityGradient( const matrix & velocityGradient ) override;

        void setAcousticsTemperature( const matrix & temperatureGradient ) override;

        void setAcousticsPressure( const matrix & pressureGradient ) override;

        void setAcousticsVelocity( const matrix & velocityGradient ) override;

        void solve() override;

    private:
        void checkTimeDiscretisationScheme();

        void continuityErrs();

        void courantNo();

        scalar evaluateMomentumResidual();

        // Dictionaries
        IOdictionary transportProperties;

        // Dimensioned scalars
        dimensionedScalar nu;
        dimensionedScalar rho;

        // Fields
        volScalarField p;
        volVectorField U;
        surfaceScalarField phi;
        volScalarField AU;
        volVectorField HU;

        // PIMPLE controls
        int nCorr;
        int nNonOrthCorr;
        int minIter;
        int maxIter;
        scalar absoluteTolerance;
        scalar relativeTolerance;
        scalar pisoTolerance;

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

        bool turbulenceSwitch;
};

#endif
