
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowSolidSolver.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RBFCoarsening.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "AitkenPostProcessing.H"
#include <iomanip>

class SDCFsiSolverTest : public ::testing::Test
{
    protected:
        virtual void SetUp()
        {
            scalar r0 = 0.2;
            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.1;
            scalar p0 = 0;
            scalar dt = 1;
            int N = 5;
            scalar L = 1;
            scalar T = 1;
            scalar dx = L / N;
            scalar rho = 1.225;
            scalar E = 490;
            scalar h = 1.0e-3;
            scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
            scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
            scalar kappa = c0 / u0;

            bool parallel = false;
            int extrapolation = 0;
            scalar tol = 1.0e-2;
            int maxIter = 20;
            scalar initialRelaxation = 1.0e-3;
            int maxUsedIterations = 50;
            int nbReuse = 0;

            scalar singularityLimit = 1.0e-13;
            int reuseInformationStartingFromTimeIndex = 0;
            bool scaling = false;
            bool updateJacobian = false;
            scalar beta = 0.5;
            int minIter = 5;

            ASSERT_NEAR( kappa, 10, 1.0e-13 );
            ASSERT_TRUE( dx > 0 );

            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
            std::shared_ptr<tubeflow::SDCTubeFlowSolidSolver> solid( new tubeflow::SDCTubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

            shared_ptr<RBFFunctionInterface> rbfFunction;
            shared_ptr<RBFInterpolation> rbfInterpolator;
            shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
            shared_ptr<RBFCoarsening> rbfInterpToMesh;

            rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
            rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
            rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

            rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
            rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
            rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

            shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

            rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
            rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
            rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

            rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
            rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
            rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

            shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

            std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
            convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, false, tol ) ) );
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );

            shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

            shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

            assert( sdcFluidSolver );
            assert( sdcSolidSolver );

            std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing ) );

            int nbNodes = 3;

            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
            quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussLobatto<scalar>( nbNodes ) );

            sdc = std::shared_ptr<sdc::SDC> ( new sdc::SDC( fsiSolver, quadrature, 1.0e-10, 1, 20 ) );
        }

        virtual void TearDown()
        {
            sdc.reset();
        }

        std::shared_ptr<sdc::SDC> sdc;
};

TEST_F( SDCFsiSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_F( SDCFsiSolverTest, timeStep )
{
    std::cout << std::setprecision( 9 );
    sdc->solveTimeStep( 0 );
    ASSERT_TRUE( sdc->isConverged() );
}

TEST_F( SDCFsiSolverTest, run )
{
    sdc->run();
    ASSERT_TRUE( sdc->isConverged() );
}

TEST( SDCFsiTest, order )
{
    int N = 5;
    scalar r0 = 0.2;
    scalar a0 = M_PI * r0 * r0;
    scalar u0 = 0.1;
    scalar p0 = 0;
    scalar L = 1;
    scalar T = 1;
    scalar dx = L / N;
    scalar rho = 1.225;
    scalar E = 490;
    scalar h = 1.0e-3;
    scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
    scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
    scalar kappa = c0 / u0;

    bool parallel = false;
    int extrapolation = 0;
    scalar tol = 1.0e-5;
    int maxIter = 20;
    scalar initialRelaxation = 1.0e-3;
    int maxUsedIterations = 50;
    int nbReuse = 0;
    scalar singularityLimit = 1.0e-13;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;
    scalar beta = 0.01;
    bool updateJacobian = false;
    int minIter = 5;

    ASSERT_NEAR( kappa, 10, 1.0e-13 );
    ASSERT_TRUE( dx > 0 );

    int nbComputations = 6;

    std::deque<std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> > fluidSolvers;
    std::deque<int> nbTimeStepsList;

    for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
    {
        int nbTimeSteps = 4 * std::pow( 2, iComputation );
        std::cout << "nbTimeSteps = " << nbTimeSteps << std::endl;
        scalar dt = T / nbTimeSteps;

        std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        std::shared_ptr<tubeflow::SDCTubeFlowSolidSolver> solid( new tubeflow::SDCTubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        shared_ptr<RBFFunctionInterface> rbfFunction;
        shared_ptr<RBFInterpolation> rbfInterpolator;
        shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
        shared_ptr<RBFCoarsening> rbfInterpToMesh;

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
        rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

        shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, false, tol ) ) );
        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );

        shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

        shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
        std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

        assert( sdcFluidSolver );
        assert( sdcSolidSolver );

        std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing ) );

        int nbNodes = 3;

        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
        quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussLobatto<scalar>( nbNodes ) );

        std::shared_ptr<sdc::SDC> sdc( new sdc::SDC( fsiSolver, quadrature, 1.0e-13, 10, 20 ) );

        sdc->run();
        ASSERT_TRUE( sdc->isConverged() );

        fluidSolvers.push_back( fluid );
        nbTimeStepsList.push_back( nbTimeSteps );
    }

    for ( int i = 0; i < 2; i++ )
    {
        fsi::vector ref;

        if ( i == 0 )
            ref = fluidSolvers.back()->u;
        else
            ref = fluidSolvers.back()->a;

        std::deque<scalar> errors;

        for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
        {
            fsi::vector data;

            if ( i == 0 )
                data = fluidSolvers.at( iComputation )->u;
            else
                data = fluidSolvers.at( iComputation )->a;

            scalar error = (ref - data).norm() / data.norm();
            std::cout << "error = " << error << std::endl;
            errors.push_back( error );
        }

        for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
        {
            scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
            std::cout << "order = " << order << std::endl;
            ASSERT_GE( order, 3.8 );
        }
    }
}
