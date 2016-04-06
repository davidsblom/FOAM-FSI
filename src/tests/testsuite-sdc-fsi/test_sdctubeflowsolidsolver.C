
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowLinearSolidSolver.H"
#include "SDCTubeFlowLinearizedSolidSolver.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RBFCoarsening.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "AitkenPostProcessing.H"
#include "AbsoluteConvergenceMeasure.H"

TEST( SDCSolidTest, order )
{
    scalar r0 = 0.2;
    scalar h = 1.0e-3;
    scalar L = 1;
    scalar rho_s = 1.225;
    scalar E0 = 490;
    scalar G = 490;
    scalar nu = 0.5;
    scalar p0 = 0.001;
    scalar T = 1;

    int N = 10;
    int nbComputations = 3;

    std::deque<std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver> > solidSolvers;
    std::deque<int> nbTimeStepsList;

    for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
    {
        int nbTimeSteps = 200 * std::pow( 2, iComputation );
        scalar dt = T / nbTimeSteps;

        std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver> solid( new tubeflow::SDCTubeFlowLinearSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T ) );
        solid->p.fill( p0 );

        int nbNodes = 2;

        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
        quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( nbNodes ) );

        std::shared_ptr<sdc::SDC> sdc( new sdc::SDC( solid, quadrature, 1.0e-13, 10, 50 ) );

        sdc->run();
        ASSERT_TRUE( sdc->isConverged() );

        solidSolvers.push_back( solid );
        nbTimeStepsList.push_back( nbTimeSteps );
    }

    fsi::vector ref = solidSolvers.back()->r;

    std::deque<scalar> errors;

    for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
    {
        fsi::vector data;

        data = solidSolvers.at( iComputation )->r;

        scalar error = (ref - data).norm() / data.norm();
        std::cout << "error = " << error << std::endl;
        errors.push_back( error );
    }

    for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
    {
        scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
        std::cout << "order = " << order << std::endl;
        ASSERT_GE( order, 2 );
    }
}


class SDCFsiSolidSolverTest : public ::testing::TestWithParam<std::string>
{
    protected:
        virtual void SetUp()
        {
            scalar r0 = 3.0e-3;
            scalar h = 3.0e-4;
            scalar L = 0.126;
            scalar rho_s = 1000;
            scalar E0 = 4.0e5;
            scalar G = 4.0e5;
            scalar nu = 0.5;

            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.26;
            scalar p0 = 0;
            scalar dt = 1;
            int N = 5;
            scalar T = 1;
            scalar rho_f = 1060;
            scalar E = 490;
            scalar cmk = std::sqrt( E * h / (2 * rho_f * r0) );

            bool parallel = false;
            int extrapolation = 0;
            scalar tol = 1.0e-3;
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

            std::string solidSolverSetting = GetParam();

            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );

            std::shared_ptr<fsi::BaseMultiLevelSolver> solid;

            if ( solidSolverSetting == "linear" )
                solid = std::shared_ptr<fsi::BaseMultiLevelSolver>( new tubeflow::SDCTubeFlowLinearSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T ) );

            if ( solidSolverSetting == "linearized" )
                solid = std::shared_ptr<fsi::BaseMultiLevelSolver>( new tubeflow::SDCTubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0 ) );

            assert( solid );

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

            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, true, tol ) ) );
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

            sdc = std::shared_ptr<sdc::SDC> ( new sdc::SDC( fsiSolver, quadrature, 1.0e-15, 1, 250 ) );
        }

        virtual void TearDown()
        {
            sdc.reset();
        }

        std::shared_ptr<sdc::SDC> sdc;
};

INSTANTIATE_TEST_CASE_P( tests, SDCFsiSolidSolverTest, ::testing::Values( "linearized", "linear" ) );

TEST_P( SDCFsiSolidSolverTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( SDCFsiSolidSolverTest, timeStep )
{
    sdc->solveTimeStep( 0 );
    ASSERT_TRUE( sdc->isConverged() );
}

TEST_P( SDCFsiSolidSolverTest, run )
{
    sdc->run();
    ASSERT_TRUE( sdc->isConverged() );
}

TEST( SDCFsiSolidTest, order )
{
    scalar r0 = 3.0e-3;
    scalar h = 3.0e-4;
    scalar L = 0.126;
    scalar rho_s = 1000;
    scalar E0 = 4.0e5;
    scalar G = 4.0e5;
    scalar nu = 0.5;

    scalar a0 = M_PI * r0 * r0;
    scalar u0 = 0.26;
    scalar p0 = 0;
    int N = 5;
    scalar T = 1;
    scalar rho_f = 1060;
    scalar cmk = std::sqrt( E0 * h / (2 * rho_f * r0) );

    bool parallel = false;
    int extrapolation = 0;
    scalar tol = 1.0e-3;
    int maxIter = 50;
    scalar initialRelaxation = 1.0e-3;
    int maxUsedIterations = 50;
    int nbReuse = 0;
    scalar singularityLimit = 1.0e-13;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;
    scalar beta = 0.01;
    bool updateJacobian = false;
    int minIter = 5;

    int nbComputations = 3;

    std::deque<std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> > fluidSolvers;
    std::deque<std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver> > solidSolvers;
    std::deque<int> nbTimeStepsList;

    for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
    {
        int nbTimeSteps = 80 * std::pow( 2, iComputation );
        std::cout << "nbTimeSteps = " << nbTimeSteps << std::endl;
        scalar dt = T / nbTimeSteps;

        std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
        std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver> solid( new tubeflow::SDCTubeFlowLinearSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T ) );

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

        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, true, tol ) ) );
        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new AbsoluteConvergenceMeasure( 0, true, 1.0e-14 ) ) );
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

        std::shared_ptr<sdc::SDC> sdc( new sdc::SDC( fsiSolver, quadrature, 1.0e-12, 10, 50 ) );

        sdc->run();
        ASSERT_TRUE( sdc->isConverged() );

        fluidSolvers.push_back( fluid );
        solidSolvers.push_back( solid );
        nbTimeStepsList.push_back( nbTimeSteps );
    }

    for ( int i = 0; i < 2; i++ )
    {
        fsi::vector ref;

        if ( i == 0 )
            ref = solidSolvers.back()->r;
        else
            ref = solidSolvers.back()->u;

        std::deque<scalar> errors;

        for ( int iComputation = 0; iComputation < nbComputations - 1; iComputation++ )
        {
            fsi::vector data;

            if ( i == 0 )
                data = solidSolvers.at( iComputation )->r;
            else
                data = solidSolvers.at( iComputation )->u;

            scalar error = (ref - data).norm() / data.norm();
            std::cout << "error = " << error << std::endl;
            errors.push_back( error );
        }

        for ( int iComputation = 0; iComputation < nbComputations - 2; iComputation++ )
        {
            scalar order = ( std::log10( errors.at( iComputation ) ) - std::log10( errors.at( iComputation + 1 ) ) ) / ( std::log10( nbTimeStepsList.at( iComputation + 1 ) ) - std::log10( nbTimeStepsList.at( iComputation ) ) );
            std::cout << "order = " << order << std::endl;
        }
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
            ASSERT_NEAR( order, 4, 0.3 );
        }
    }
}
