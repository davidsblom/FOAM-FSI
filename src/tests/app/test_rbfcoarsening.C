
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"
#include "RBFCoarsening.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"
#include "TPSFunction.H"
#include "gtest/gtest.h"

using namespace rbf;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class RBFCoarseningParametrizedTest : public TestWithParam < std::tr1::tuple<bool, int, int, int, double, bool> >
{
protected:

    virtual void SetUp()
    {
        int rbfFunctionId = std::tr1::get<1>( GetParam() );

        std::shared_ptr<RBFFunctionInterface> rbfFunction;

        if ( rbfFunctionId == 0 )
            rbfFunction = std::shared_ptr<RBFFunctionInterface>( new WendlandC0Function( 5 ) );

        else
        if ( rbfFunctionId == 1 )
            rbfFunction = std::shared_ptr<RBFFunctionInterface>( new WendlandC2Function( 5 ) );

        else
        if ( rbfFunctionId == 2 )
            rbfFunction = std::shared_ptr<RBFFunctionInterface>( new WendlandC4Function( 5 ) );

        else
        if ( rbfFunctionId == 3 )
            rbfFunction = std::shared_ptr<RBFFunctionInterface>( new WendlandC6Function( 5 ) );

        else
        if ( rbfFunctionId == 4 )
            rbfFunction = std::shared_ptr<RBFFunctionInterface>( new TPSFunction() );

        else
            assert( false );

        bool polynomialTerm = false;
        bool cpu = std::tr1::get<5>( GetParam() );
        std::shared_ptr<RBFInterpolation> rbfInterpolator( new RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );

        bool enabled = std::tr1::get<0>( GetParam() );
        int coarseningMinPoints = std::tr1::get<2>( GetParam() );
        int coarseningMaxPoints = std::tr1::get<3>( GetParam() );
        double tol = std::tr1::get<4>( GetParam() );

        rbf = std::shared_ptr<RBFCoarsening>( new RBFCoarsening( rbfInterpolator, enabled, false, false, tol, 0.1, coarseningMinPoints, coarseningMaxPoints, false ) );
    }

    virtual void TearDown()
    {
        rbf.reset();
    }

    std::shared_ptr<RBFCoarsening> rbf;
};

INSTANTIATE_TEST_CASE_P( RBFTest, RBFCoarseningParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 2, 3, 4 ), Values( 10, 25 ), Values( 25, 100 ), Values( 1.0e-6, 1.0e-30 ), Bool() ) );

TEST_P( RBFCoarseningParametrizedTest, rbf1d_regularity )
{
    // Check that the RBF function approximates a smooth function well
    // away from the nodes

    matrix x( 30, 1 ), y( 30, 1 ), xnew( 100, 1 ), ynew;
    x.col( 0 ).setLinSpaced( 30, 0, 10 );
    y.col( 0 ) = x.array().sin();
    xnew.col( 0 ).setLinSpaced( 100, 0, 10 );

    rbf->compute( x, xnew );
    rbf->interpolate( y, ynew );

    for ( int i = 0; i < ynew.rows(); i++ )
        ASSERT_NEAR( ynew( i, 0 ), std::sin( xnew( i, 0 ) ), 1.0e-1 );
}

TEST_P( RBFCoarseningParametrizedTest, rbf1d_regularity2 )
{
    matrix x( 30, 1 ), y( 30, 1 ), xnew( 30, 1 ), ynew;
    x.col( 0 ).setLinSpaced( 30, 0, 10 );
    y.col( 0 ) = x.array().sin();
    xnew.col( 0 ).setLinSpaced( 30, 0, 10 );

    rbf->compute( x, xnew );
    rbf->interpolate( y, ynew );

    for ( int i = 0; i < ynew.rows(); i++ )
        ASSERT_NEAR( ynew( i, 0 ), std::sin( xnew( i, 0 ) ), 1.0e-1 );

    double tol = std::tr1::get<4>( GetParam() );
    int coarseningMaxPoints = std::tr1::get<3>( GetParam() );

    if ( tol < 1.0e-20 && coarseningMaxPoints == 100 )
    {
        for ( int i = 0; i < y.rows(); i++ )
            ASSERT_NEAR( ynew( i, 0 ), std::sin( xnew( i, 0 ) ), 1.0e-10 );
    }
}

TEST_P( RBFCoarseningParametrizedTest, rbf_interpolation_comparison )
{
    matrix x( 30, 1 ), y( 30, 1 ), xnew( 100, 1 ), ynew;
    x.col( 0 ).setLinSpaced( 30, 0, 10 );
    y.col( 0 ) = x.array().sin();
    xnew.col( 0 ).setLinSpaced( 100, 0, 10 );

    rbf->compute( x, xnew );
    rbf->interpolate( y, ynew );

    double tol = std::tr1::get<4>( GetParam() );
    int coarseningMaxPoints = std::tr1::get<3>( GetParam() );

    if ( tol < 1.0e-20 && coarseningMaxPoints == 100 )
    {
        bool cpu = std::tr1::get<5>( GetParam() );
        bool polynomialTerm = false;

        std::shared_ptr<RBFFunctionInterface> rbfFunction = rbf->rbf->rbfFunction;
        std::shared_ptr<RBFInterpolation> rbfInterpolator( new RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );

        matrix ynew2( ynew.rows(), ynew.cols() );

        rbfInterpolator->compute( x, xnew );
        rbfInterpolator->interpolate( y, ynew2 );

        for ( int i = 0; i < y.rows(); i++ )
            ASSERT_NEAR( ynew( i, 0 ), ynew2( i, 0 ), 1.0e-10 );
    }
}
