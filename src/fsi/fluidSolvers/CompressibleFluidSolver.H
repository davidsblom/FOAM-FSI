
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
#include "scalarIOList.H"

namespace Foam
{
    typedef IOList<vector> vectorIOList;
}

class CompressibleFluidSolver : public foamFluidSolver
{
    public:
        CompressibleFluidSolver(
            const string & name,
            shared_ptr<argList> args,
            shared_ptr<Time> runTime
            );

        virtual ~CompressibleFluidSolver();

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

        void setAcousticsTemperature( const matrix & temperature ) override;

        void setAcousticsPressure( const matrix & pressure ) override;

        void setAcousticsVelocity( const matrix & velocity ) override;

        void solve() override;

        void finalizeTimeStep() override;

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
        volScalarField ddtp;
        volScalarField ddtrho;

        // Continuity errors
        scalar cumulativeContErr;

        // Convergence tolerance fluid domain
        scalar convergenceTolerance;
        label nOuterCorr;

        // Acoustics coupling patch

        word acousticsPatchName;
        label acousticsPatchID;

        scalarIOList pressureAcousticsBC;
        scalarIOList temperatureAcousticsBC;
        vectorIOList velocityAcousticsBC;
};

#endif
