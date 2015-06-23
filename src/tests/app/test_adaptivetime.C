
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "ESDIRK.H"
#include "Piston.H"
#include "AdaptiveTimeStepper.H"
#include "gtest/gtest.h"

using namespace sdc;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

/*
 * Implements the piston problem for which a exact solution exists.
 * With the exact solution, the convergence orders can be confirmed
 * exactly for this test case.
 *
 * Reference: A. H. van Zuijlen, Efficient Higher Order Time
 * Integration of Partitioned Systems, PhD thesis, page 22.
 *
 * Notes Gaussian quadrature:
 * http://undergraduate.csse.uwa.edu.au/units/GENG2140/Mech2012/Gaussian_Quadrature.pdf
 * http://en.wikipedia.org/wiki/Gaussian_quadrature
 * http://en.wikipedia.org/wiki/Legendre_polynomials
 */

class AdaptiveTimeSteppingTest : public TestWithParam< std::tr1::tuple<int, std::string> >
{
protected:

    virtual void SetUp()
    {
        double dt, q0, qdot0, As, Ac, omega, endTime;

        int nbTimeSteps = std::tr1::get<0>( GetParam() );
        std::string method = std::tr1::get<1>( GetParam() );

        endTime = 100;
        dt = endTime / nbTimeSteps;
        As = 100;
        Ac = As;
        omega = 1;
        q0 = -As;
        qdot0 = -As;

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( true, "h211b", 1.0e-3, 0.5, 2 ) );

        piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        esdirk = std::shared_ptr<ESDIRK> ( new ESDIRK( piston, method, adaptiveTimeStepper ) );
    }

    virtual void TearDown()
    {
        piston.reset();
        esdirk.reset();
        adaptiveTimeStepper.reset();
    }

    std::shared_ptr<Piston> piston;
    std::shared_ptr<ESDIRK> esdirk;
    std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;
};

INSTANTIATE_TEST_CASE_P( testParameters, AdaptiveTimeSteppingTest, ::testing::Combine( Values( 2, 100, 200, 400, 800, 1600, 3200 ), Values( "SDIRK2", "SDIRK3", "SDIRK4" ) ) );

TEST_P( AdaptiveTimeSteppingTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( AdaptiveTimeSteppingTest, solveTimeStep )
{
    esdirk->solveTimeStep( 0 );
}

TEST_P( AdaptiveTimeSteppingTest, run )
{
    esdirk->run();

    ASSERT_NEAR( piston->t, 100, 1.0e-5 );
}
