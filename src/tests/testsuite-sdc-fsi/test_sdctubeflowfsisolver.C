
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowSolidSolver.H"
#include "GaussRadau.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RBFCoarsening.H"
#include "RelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"

class SDCFsiSolverTest : public::testing::Test
 {
 protected:

     virtual void SetUp()
     {
         scalar r0 = 0.2;
         scalar a0 = M_PI * r0 * r0;
         scalar u0 = 0.1;
         scalar p0 = 0;
         scalar dt = 0.1;
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
         scalar tau = u0 * dt / L;

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
         scalar beta = 0.5;
         bool updateJacobian = false;
         int minIter = 5;

         ASSERT_NEAR( tau, 0.01, 1.0e-13 );
         ASSERT_NEAR( kappa, 10, 1.0e-13 );
         ASSERT_TRUE( dx > 0 );

         std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid ( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
         std::shared_ptr<tubeflow::SDCTubeFlowSolidSolver> solid ( new tubeflow::SDCTubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

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

         int nbNodes = 4;

         std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
         quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

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
     sdc->solveTimeStep( 0 );
     ASSERT_TRUE( sdc->isConverged() );
 }

 TEST_F( SDCFsiSolverTest, run )
 {
     sdc->run();
     ASSERT_TRUE( sdc->isConverged() );
 }
