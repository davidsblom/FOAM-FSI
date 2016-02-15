
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "GaussRadau.H"
#include "SDC.H"

class SDCFluidSolverTest : public::testing::Test
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

         ASSERT_NEAR( tau, 0.01, 1.0e-13 );
         ASSERT_NEAR( kappa, 10, 1.0e-13 );
         ASSERT_TRUE( dx > 0 );

         std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid ( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );

         int nbNodes = 4;
         scalar tol = 1.0e-10;

         std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
         quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

         sdc = std::shared_ptr<sdc::SDC> ( new sdc::SDC( fluid, quadrature, tol, 3, 10 ) );
     }

     virtual void TearDown()
     {
         sdc.reset();
     }

     std::shared_ptr<sdc::SDC> sdc;
 };

 TEST_F( SDCFluidSolverTest, object )
 {
     ASSERT_TRUE( true );
 }

 TEST_F( SDCFluidSolverTest, timeStep )
 {
     sdc->solveTimeStep( 0 );
 }

 TEST_F( SDCFluidSolverTest, run )
 {
     sdc->run();
 }
