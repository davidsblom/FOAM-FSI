
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "SDC.H"
#include "Piston.H"
#include "gtest/gtest.h"

using namespace sdc;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class SDCTest : public TestWithParam< std::tr1::tuple<int, int, std::string> >
{
protected:

    virtual void SetUp()
    {
        double dt, q0, qdot0, As, Ac, omega, endTime;

        int nbNodes = std::tr1::get<0>( GetParam() );
        int nbTimeSteps = std::tr1::get<1>( GetParam() );
        std::string rule = std::tr1::get<2>( GetParam() );

        endTime = 100;
        dt = endTime / nbTimeSteps;
        As = 100;
        Ac = As;
        omega = 1;
        q0 = -As;
        qdot0 = -As;

        piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        sdc = std::shared_ptr<SDC> ( new SDC( piston, rule, nbNodes, 1.0e-10 ) );
    }

    virtual void TearDown()
    {
        piston.reset();
        sdc.reset();
    }

    std::shared_ptr<Piston> piston;
    std::shared_ptr<SDC> sdc;
};

INSTANTIATE_TEST_CASE_P( testParameters, SDCTest, ::testing::Combine( Values( 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ), Values( 10, 20, 40, 80, 100, 160 ), Values( "gauss-radau", "gauss-lobatto" ) ) );

TEST_P( SDCTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( SDCTest, referenceSolution )
{
    double result = piston->referenceSolution( 100 );

    ASSERT_NEAR( result, -35.5953231178, 1.0e-11 );
}

TEST_P( SDCTest, solve )
{
    piston->run();

    Eigen::VectorXd solution( 2 );
    piston->getSolution( solution );

    double result = solution( 1 );

    int nbTimeSteps = std::tr1::get<1>( GetParam() );

    if ( nbTimeSteps == 10 )
        ASSERT_NEAR( result, -75814.5607609, 1.0e-8 );

    ASSERT_NEAR( piston->t, 100, 1.0e-10 );
}

TEST_P( SDCTest, evaluateFunction )
{
    Eigen::VectorXd q( 2 ), f( 2 );
    double t;

    q << -100, -100;
    t = 1;

    piston->evaluateFunction( 0, q, t, f );

    ASSERT_NEAR( f( 0 ), 138.17732907, 1.0e-8 );
    ASSERT_NEAR( f( 1 ), -100, 1.0e-9 );
}

TEST_P( SDCTest, solveTimeStep )
{
    sdc->solveTimeStep( 0 );
}

TEST_P( SDCTest, run )
{
    sdc->run();

    Eigen::VectorXd solution( 2 );
    piston->getSolution( solution );

    double result = solution( 1 );
    double ref = piston->referenceSolution( 100 );
    double error = std::abs( result - ref ) / std::abs( ref );

    int nbTimeSteps = std::tr1::get<1>( GetParam() );
    int nbNodes = std::tr1::get<0>( GetParam() );
    std::string rule = std::tr1::get<2>( GetParam() );

    ASSERT_NEAR( piston->t, 100, 1.0e-10 );

    if ( nbTimeSteps == 100 && nbNodes > 4 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 1.0e-5, 1.0e-4 );

    if ( nbTimeSteps == 10 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 1.72134678e+03, 1.0e-5 );

    if ( nbTimeSteps == 20 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 5.38712906e+02, 1.0e-6 );

    if ( nbTimeSteps == 40 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 2.62943671e+01, 1.0e-7 );

    if ( nbTimeSteps == 80 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 2.74294224e+00, 1.0e-8 );

    if ( nbTimeSteps == 160 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 3.25178345e-01, 1.0e-9 );

    if ( nbTimeSteps == 10 && nbNodes == 3 && rule == "gauss-lobatto" )
        ASSERT_NEAR( error, 1.37211827e+03, 1.0e-5 );

    if ( nbTimeSteps == 20 && nbNodes == 3 && rule == "gauss-lobatto" )
        ASSERT_NEAR( error, 1.79890465e+02, 1.0e-5 );

    if ( nbTimeSteps == 40 && nbNodes == 3 && rule == "gauss-lobatto" )
        ASSERT_NEAR( error, 4.55455671e+00, 1.0e-5 );

    if ( nbTimeSteps == 80 && nbNodes == 3 && rule == "gauss-lobatto" )
        ASSERT_NEAR( error, 2.44720474e-01, 1.0e-5 );

    if ( nbTimeSteps == 160 && nbNodes == 3 && rule == "gauss-lobatto" )
        ASSERT_NEAR( error, 1.47661061e-02, 1.0e-5 );
}
