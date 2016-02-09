
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "SDC.H"
#include "Piston.H"
#include "Cos.H"
#include "Oscillator.H"
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
        scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

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
        tol = 1.0e-10;

        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( false ) );

        piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        sdc = std::shared_ptr<SDC> ( new SDC( piston, adaptiveTimeStepper, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );

        std::shared_ptr<sdc::SDC> sdc( new SDC( rule, nbNodes, tol ) );
        piston_sdc = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega, sdc, sdc->nodes.rows() ) );
    }

    virtual void TearDown()
    {
        piston.reset();
        sdc.reset();
        piston_sdc.reset();
    }

    std::shared_ptr<Piston> piston;
    std::shared_ptr<SDC> sdc;
    std::shared_ptr<Piston> piston_sdc;
};

INSTANTIATE_TEST_CASE_P( testParameters, SDCTest, ::testing::Combine( Values( 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ), Values( 1, 10, 20, 40, 80, 100, 160 ), Values( "gauss-radau", "gauss-lobatto", "clenshaw-curtis", "uniform", "uniform-right-sided" ) ) );

class SDCEstimateOrderTest : public TestWithParam< std::tr1::tuple<int, std::string> >
{
protected:

    virtual void SetUp()
    {
        scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

        int nbNodes = std::tr1::get<0>( GetParam() );
        int nbTimeSteps = 100;
        std::string rule = std::tr1::get<1>( GetParam() );

        endTime = 100;
        dt = endTime / nbTimeSteps;
        As = 100;
        Ac = As;
        omega = 1;
        q0 = -As;
        qdot0 = -As;
        tol = 1.0e-15;

        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper1;
        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper2;

        adaptiveTimeStepper1 = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( false ) );
        adaptiveTimeStepper2 = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( false ) );

        piston1 = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        sdc1 = std::shared_ptr<SDC> ( new SDC( piston1, adaptiveTimeStepper1, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );

        nbTimeSteps *= 2;
        dt = endTime / nbTimeSteps;

        piston2 = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
        sdc2 = std::shared_ptr<SDC> ( new SDC( piston2, adaptiveTimeStepper2, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );
    }

    virtual void TearDown()
    {
        piston1.reset();
        piston2.reset();
        sdc1.reset();
        sdc2.reset();
    }

    std::shared_ptr<Piston> piston1;
    std::shared_ptr<Piston> piston2;
    std::shared_ptr<SDC> sdc1;
    std::shared_ptr<SDC> sdc2;
};

INSTANTIATE_TEST_CASE_P( testParameters, SDCEstimateOrderTest, ::testing::Combine( Values( 2, 4 ), Values( "gauss-radau", "gauss-lobatto", "clenshaw-curtis", "uniform", "uniform-right-sided" ) ) );

TEST_P( SDCTest, object )
{
    ASSERT_TRUE( true );
}

TEST_P( SDCTest, referenceSolution )
{
    scalar result = piston->referenceSolution( 100 );

    ASSERT_NEAR( result, -35.5953231178, 1.0e-11 );
}

TEST_P( SDCTest, solve )
{
    piston->run();

    fsi::vector solution( 2 ), f;
    piston->getSolution( solution, f );

    scalar result = solution( 1 );

    int nbTimeSteps = std::tr1::get<1>( GetParam() );

    if ( nbTimeSteps == 10 )
        ASSERT_NEAR( result, -75814.5607609, 1.0e-8 );

    ASSERT_NEAR( piston->t, 100, 1.0e-10 );
}

TEST_P( SDCTest, evaluateFunction )
{
    fsi::vector q( 2 ), f( 2 );
    scalar t;

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

    fsi::vector solution( 2 ), f;
    piston->getSolution( solution, f );

    scalar result = solution( 1 );
    scalar ref = piston->referenceSolution( 100 );
    scalar error = std::abs( result - ref ) / std::abs( ref );

    int nbTimeSteps = std::tr1::get<1>( GetParam() );
    int nbNodes = std::tr1::get<0>( GetParam() );
    std::string rule = std::tr1::get<2>( GetParam() );

    ASSERT_NEAR( piston->t, 100, 1.0e-10 );

    if ( nbTimeSteps == 100 && nbNodes > 4 && rule == "gauss-radau" )
        ASSERT_NEAR( error, 1.0e-5, 1.0e-4 );

    if ( nbTimeSteps == 10 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_LE( error, 1.72134678e+03 );

    if ( nbTimeSteps == 20 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_LE( error, 5.38712906e+02 );

    if ( nbTimeSteps == 40 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_LE( error, 2.62943671e+01 );

    if ( nbTimeSteps == 80 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_LE( error, 2.74294224e+00 );

    if ( nbTimeSteps == 160 && nbNodes == 3 && rule == "gauss-radau" )
        ASSERT_LE( error, 3.25178345e-01 );

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

TEST_P( SDCTest, runCompareSDC )
{
    piston_sdc->run();
    sdc->run();

    fsi::vector solution_piston_sdc( 2 ), solution_piston( 2 ), f;
    piston->getSolution( solution_piston, f );
    piston_sdc->getSolution( solution_piston_sdc, f );

    ASSERT_NEAR( solution_piston_sdc( 0 ), solution_piston( 0 ), 1.0e-8 );
    ASSERT_NEAR( solution_piston_sdc( 1 ), solution_piston( 1 ), 1.0e-7 );
}

TEST_P( SDCEstimateOrderTest, order )
{
    sdc1->run();
    sdc2->run();

    fsi::vector solution1( 2 ), solution2( 2 ), f;
    piston1->getSolution( solution1, f );
    piston2->getSolution( solution2, f );

    scalar result1 = solution1( 1 );
    scalar result2 = solution2( 1 );

    scalar ref = piston1->referenceSolution( 100 );
    scalar error1 = std::abs( result1 - ref ) / std::abs( ref );
    scalar error2 = std::abs( result2 - ref ) / std::abs( ref );

    Info << "error1 = " << error1 << ", error2 = " << error2 << endl;

    int nbTimeSteps = 10;
    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    int nbNodes = std::tr1::get<0>( GetParam() );
    std::string rule = std::tr1::get<1>( GetParam() );

    scalar tol = 3.0e-1;

    if ( rule == "gauss-radau" )
        ASSERT_NEAR( order, nbNodes * 2 - 1, tol );

    if ( rule == "gauss-lobatto" )
        ASSERT_NEAR( order, nbNodes * 2 - 2, tol );

    if ( rule == "clenshaw-curtis" )
        ASSERT_NEAR( order, nbNodes, tol );

    if ( rule == "uniform" )
        ASSERT_NEAR( order, nbNodes, tol );

    if ( rule == "uniform-right-sided" )
        ASSERT_NEAR( order, nbNodes - 1, tol );
}

TEST( CosTest, SDC )
{
    std::string rule = "gauss-lobatto";

    // rule = "uniform";
    // rule = "clenshaw-curtis";
    int nbNodes = 3;
    scalar tol = 1.0e-25;

    int nbTimeSteps = 5;
    scalar endTime = 0.05;
    scalar dt = endTime / nbTimeSteps;
    scalar amplitude = 0.2;
    scalar frequency = 5;

    std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
    std::shared_ptr<Cos> cos1( new Cos( nbTimeSteps, dt, endTime, amplitude, frequency ) );
    std::shared_ptr<SDC> sdc1( new SDC( cos1, adaptiveTimeStepper, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );
    std::shared_ptr<Cos> cos2( new Cos( nbTimeSteps * 2, dt / 2, endTime, amplitude, frequency ) );
    std::shared_ptr<SDC> sdc2( new SDC( cos2, adaptiveTimeStepper, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );

    sdc1->run();
    sdc2->run();

    scalar ref = 0.5 * amplitude * std::sin( M_PI * frequency * endTime ) * M_PI * frequency;
    scalar error1 = std::abs( cos1->f - ref ) / std::abs( ref );
    scalar error2 = std::abs( cos2->f - ref ) / std::abs( ref );
    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    std::cout << "error1 = " << error1 << std::endl;
    std::cout << "error2 = " << error2 << std::endl;
    std::cout << "order = " << order << std::endl;

    error1 = std::abs( cos1->sol - ref ) / std::abs( ref );
    error2 = std::abs( cos2->sol - ref ) / std::abs( ref );
    order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    std::cout << "error1 = " << error1 << std::endl;
    std::cout << "error2 = " << error2 << std::endl;
    std::cout << "order = " << order << std::endl;
}

TEST( OscillatorTest, SDC )
{
    std::string rule = "gauss-lobatto";
    int nbNodes = 3;
    scalar tol = 1.0e-25;

    int nbTimeSteps = 100;
    scalar endTime = 10;
    scalar dt = endTime / nbTimeSteps;
    scalar amplitude = 0;
    scalar frequency = 1;
    scalar m = 1;
    scalar k = 1;
    fsi::vector q0( 2 );
    q0 << 1, 0;

    std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
    std::shared_ptr<Oscillator> oscillator1( new Oscillator( nbTimeSteps, dt, q0, amplitude, frequency, m, k ) );
    std::shared_ptr<SDC> sdc1( new SDC( oscillator1, adaptiveTimeStepper, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );
    std::shared_ptr<Oscillator> oscillator2( new Oscillator( nbTimeSteps * 2, dt / 2, q0, amplitude, frequency, m, k ) );
    std::shared_ptr<SDC> sdc2( new SDC( oscillator2, adaptiveTimeStepper, rule, nbNodes, tol, nbNodes, 10 * nbNodes ) );

    sdc1->run();
    sdc2->run();

    scalar ref = std::cos( endTime );
    scalar error1 = std::abs( oscillator1->sol( 0 ) - ref ) / std::abs( ref );
    scalar error2 = std::abs( oscillator2->sol( 0 ) - ref ) / std::abs( ref );
    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    std::cout << "error1 = " << error1 << std::endl;
    std::cout << "error2 = " << error2 << std::endl;
    std::cout << "order = " << order << std::endl;

    ASSERT_NEAR( order, 4, 1.0e-3 );
}
