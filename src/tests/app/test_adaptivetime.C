
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "ESDIRK.H"
#include "SDC.H"
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

class AdaptiveTimeSteppingESDIRKTest : public TestWithParam< std::tr1::tuple<int, std::string, scalar> >
{
protected:

    virtual void SetUp()
    {
        scalar dt, q0, qdot0, As, Ac, omega, endTime;

        int nbTimeSteps = std::tr1::get<0>( GetParam() );
        std::string method = std::tr1::get<1>( GetParam() );
        scalar tol = std::tr1::get<2>( GetParam() );

        endTime = 100;
        dt = endTime / nbTimeSteps;
        As = 100;
        Ac = As;
        omega = 1;
        q0 = -As;
        qdot0 = -As;

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( true, "h211b", tol, 0.5 ) );

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

INSTANTIATE_TEST_CASE_P( testParameters, AdaptiveTimeSteppingESDIRKTest, ::testing::Combine( Values( 2, 100, 200, 400, 800, 1600, 3200 ), Values( "SDIRK2", "SDIRK3", "SDIRK4", "ESDIRK3", "ESDIRK4", "ESDIRK5" ), Values( 1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8 ) ) );

TEST_P( AdaptiveTimeSteppingESDIRKTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( AdaptiveTimeSteppingESDIRKTest, solveTimeStep )
{
    esdirk->solveTimeStep( 0 );
}

TEST_P( AdaptiveTimeSteppingESDIRKTest, run )
{
    esdirk->run();

    fsi::vector solution( 2 );
    piston->getSolution( solution );

    scalar result = solution( 1 );
    scalar ref = piston->referenceSolution( 100 );
    scalar error = std::abs( result - ref ) / std::abs( ref );

    scalar tol = std::tr1::get<2>( GetParam() );

    std::string method = std::tr1::get<1>( GetParam() );

    if ( method != "ESDIRK5" )
        ASSERT_LT( error, 1000 * tol );
}

class AdaptiveTimeSteppingSDCTest : public TestWithParam< std::tr1::tuple<scalar, int> >
{
protected:

    virtual void SetUp()
    {
        scalar dt, q0, qdot0, As, Ac, omega, endTime;

        int nbTimeSteps = 2;
        std::string rule = "clenshaw-curtis";

        scalar tol = std::tr1::get<0>( GetParam() );
        int nbNodes = std::tr1::get<1>( GetParam() );

        endTime = 100;
        dt = endTime / nbTimeSteps;
        As = 100;
        Ac = As;
        omega = 1;
        q0 = -As;
        qdot0 = -As;

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( true, "h211b", tol, 0.9 ) );

        piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        sdc = std::shared_ptr<SDC> ( new SDC( piston, adaptiveTimeStepper, rule, nbNodes, tol * 1.0e-2 ) );
    }

    virtual void TearDown()
    {
        piston.reset();
        sdc.reset();
        adaptiveTimeStepper.reset();
    }

    std::shared_ptr<Piston> piston;
    std::shared_ptr<SDC> sdc;
    std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;
};

INSTANTIATE_TEST_CASE_P( testParameters, AdaptiveTimeSteppingSDCTest, ::testing::Combine( Values( 1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8 ), Values( 5, 7, 9 ) ) );

TEST_P( AdaptiveTimeSteppingSDCTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( AdaptiveTimeSteppingSDCTest, solveTimeStep )
{
    sdc->solveTimeStep( 0 );
}

TEST_P( AdaptiveTimeSteppingSDCTest, run )
{
    sdc->run();

    ASSERT_NEAR( piston->t, 100, 1.0e-5 );

    fsi::vector solution( 2 );
    piston->getSolution( solution );

    scalar result = solution( 1 );
    scalar ref = piston->referenceSolution( 100 );
    scalar error = std::abs( result - ref ) / std::abs( ref );

    scalar tol = std::tr1::get<0>( GetParam() );

    ASSERT_LT( error, tol );
}
