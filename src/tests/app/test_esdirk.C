
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "ESDIRK.H"
#include "Piston.H"
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

class ESDIRKTest : public TestWithParam< std::tr1::tuple<int, std::string> >
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

        piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        esdirk = std::shared_ptr<ESDIRK> ( new ESDIRK( piston, method ) );
    }

    virtual void TearDown()
    {
        piston.reset();
        esdirk.reset();
    }

    std::shared_ptr<Piston> piston;
    std::shared_ptr<ESDIRK> esdirk;
};

INSTANTIATE_TEST_CASE_P( testParameters, ESDIRKTest, ::testing::Combine( Values( 2, 100, 200, 400, 800, 1600, 3200 ), Values( "SDIRK2", "SDIRK3", "ESDIRK3", "ESDIRK4", "ESDIRK5" ) ) );

TEST_P( ESDIRKTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( ESDIRKTest, solveTimeStep )
{
    esdirk->solveTimeStep( 0 );
}

TEST_P( ESDIRKTest, run )
{
    esdirk->run();

    Eigen::VectorXd solution( 2 );
    piston->getSolution( solution );

    double result = solution( 1 );
    double ref = piston->referenceSolution( 100 );
    double error = std::abs( result - ref ) / std::abs( ref );

    int nbTimeSteps = std::tr1::get<0>( GetParam() );
    std::string method = std::tr1::get<1>( GetParam() );

    ASSERT_NEAR( piston->t, 100, 1.0e-5 );

    // ---------------------------------------------

    if ( nbTimeSteps == 2 && method == "SDIRK2" )
        ASSERT_NEAR( error, 6983.04523052, 1.0e-8 );

    if ( nbTimeSteps == 100 && method == "SDIRK2" )
        ASSERT_NEAR( error, 1.51043910e+00, 1.0e-8 );

    if ( nbTimeSteps == 200 && method == "SDIRK2" )
        ASSERT_NEAR( error, 5.56619516e-01, 1.0e-7 );

    if ( nbTimeSteps == 400 && method == "SDIRK2" )
        ASSERT_NEAR( error, 1.61165601e-01, 1.0e-7 );

    if ( nbTimeSteps == 800 && method == "SDIRK2" )
        ASSERT_NEAR( error, 4.30446278e-02, 1.0e-7 );

    if ( nbTimeSteps == 1600 && method == "SDIRK2" )
        ASSERT_NEAR( error, 1.11061831e-02, 1.0e-7 );

    // ---------------------------------------------

    if ( nbTimeSteps == 100 && method == "SDIRK3" )
        ASSERT_NEAR( error, 2.44907705e+00, 1.0e-8 );

    if ( nbTimeSteps == 200 && method == "SDIRK3" )
        ASSERT_NEAR( error, 2.89782371e-01, 1.0e-7 );

    if ( nbTimeSteps == 400 && method == "SDIRK3" )
        ASSERT_NEAR( error, 3.54034132e-02, 1.0e-7 );

    if ( nbTimeSteps == 800 && method == "SDIRK3" )
        ASSERT_NEAR( error, 4.36890963e-03, 1.0e-7 );

    if ( nbTimeSteps == 1600 && method == "SDIRK3" )
        ASSERT_NEAR( error, 5.31256827e-04, 1.0e-7 );

    // ---------------------------------------------

    if ( nbTimeSteps == 100 && method == "ESDIRK3" )
        ASSERT_NEAR( error, 2.16401700e-01, 1.0e-8 );

    if ( nbTimeSteps == 200 && method == "ESDIRK3" )
        ASSERT_NEAR( error, 1.59883384e-02, 1.0e-7 );

    if ( nbTimeSteps == 400 && method == "ESDIRK3" )
        ASSERT_NEAR( error, 1.31386639e-03, 1.0e-7 );

    if ( nbTimeSteps == 800 && method == "ESDIRK3" )
        ASSERT_NEAR( error, 1.21257715e-04, 1.0e-7 );

    if ( nbTimeSteps == 1600 && method == "ESDIRK3" )
        ASSERT_NEAR( error, 1.24594213e-05, 1.0e-7 );

    // ---------------------------------------------

    if ( nbTimeSteps == 100 && method == "ESDIRK4" )
        ASSERT_NEAR( error, 1.06552842e-01, 1.0e-8 );

    if ( nbTimeSteps == 200 && method == "ESDIRK4" )
        ASSERT_NEAR( error, 6.62505388e-03, 1.0e-7 );

    if ( nbTimeSteps == 400 && method == "ESDIRK4" )
        ASSERT_NEAR( error, 4.15160887e-04, 1.0e-7 );

    if ( nbTimeSteps == 800 && method == "ESDIRK4" )
        ASSERT_NEAR( error, 2.60153951e-05, 1.0e-7 );

    if ( nbTimeSteps == 1600 && method == "ESDIRK4" )
        ASSERT_NEAR( error, 1.62861015e-06, 1.0e-7 );

    // ---------------------------------------------

    if ( nbTimeSteps == 100 && method == "ESDIRK5" )
        ASSERT_NEAR( error, 9.93536421e-03, 1.0e-8 );

    if ( nbTimeSteps == 200 && method == "ESDIRK5" )
        ASSERT_NEAR( error, 3.04059916e-04, 1.0e-7 );

    if ( nbTimeSteps == 400 && method == "ESDIRK5" )
        ASSERT_NEAR( error, 9.46071494e-06, 1.0e-7 );

    if ( nbTimeSteps == 800 && method == "ESDIRK5" )
        ASSERT_NEAR( error, 2.95451539e-07, 1.0e-7 );

    if ( nbTimeSteps == 1600 && method == "ESDIRK5" )
        ASSERT_NEAR( error, 9.23400561e-09, 1.0e-7 );
}
