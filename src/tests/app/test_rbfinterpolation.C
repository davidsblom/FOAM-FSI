
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"
#include "TPSFunction.H"
#include "DataValues.H"
#include "gtest/gtest.h"

using namespace rbf;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class RBFInterpolationParametrizedTest : public TestWithParam < std::tr1::tuple<int, bool> >
{
protected:

    virtual void SetUp()
    {
        int rbfFunctionId = std::tr1::get<0>( GetParam() );

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

        bool polynomialTerm = std::tr1::get<1>( GetParam() );
        bool cpu = false;
        rbf = std::shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );

        cpu = true;
        rbfCPU = std::shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );
    }

    virtual void TearDown()
    {
        rbf.reset();
        rbfCPU.reset();
    }

    std::shared_ptr<RBFInterpolation> rbf;
    std::shared_ptr<RBFInterpolation> rbfCPU;
};

INSTANTIATE_TEST_CASE_P( RBFTest, RBFInterpolationParametrizedTest, ::testing::Combine( Values( 0, 1, 2, 3, 4 ), Bool() ) );

TEST_P( RBFInterpolationParametrizedTest, rbf3d_directly_interpolate )
{
    // Check that the RBF function interpolates through the nodes (2D)

    matrix x, y, ynew, ynew2;

    x = fsi::matrix::Random( 50, 3 ).array() * 4 - 2;
    y = fsi::matrix::Random( 50, 3 ).array() * 4 - 2;

    rbf->compute( x, x );
    rbf->interpolate( y, ynew );
    rbfCPU->compute( x, x );
    rbfCPU->interpolate( y, ynew2 );

    // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
    // interpolation.
    rbf::vector rowsum = rbf->Hhat.rowwise().sum();

    for ( int i = 0; i < rowsum.rows(); i++ )
        ASSERT_NEAR( rowsum( i ), 1, 1.0e-12 );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-6 );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( ynew( i, j ), ynew2( i, j ), 1.0e-11 );
}

TEST_P( RBFInterpolationParametrizedTest, rbf1d_regularity )
{
    // Check that the RBF function approximates a smooth function well
    // away from the nodes

    matrix x( 30, 1 ), y( 30, 1 ), xnew( 100, 1 ), ynew, ynew2;
    x.col( 0 ).setLinSpaced( 30, 0, 10 );
    y.col( 0 ) = x.array().sin();
    xnew.col( 0 ).setLinSpaced( 100, 0, 10 );

    rbf->compute( x, xnew );
    rbf->interpolate( y, ynew );
    rbfCPU->compute( x, xnew );
    rbfCPU->interpolate( y, ynew2 );

    for ( int i = 0; i < ynew.rows(); i++ )
        ASSERT_NEAR( ynew( i, 0 ), std::sin( xnew( i, 0 ) ), 1.0e-1 );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( ynew( i, j ), ynew2( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, oneD_small )
{
    // Initialize variables

    int n_A = 7;
    int n_B = 6;

    matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

    // Location of original values
    x << 0, 1, 2, 3, 4, 5, 6;
    xnew << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;

    // Values
    y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    rbf.compute( x, xnew );
    rbf.interpolate( y, ynew );

    ASSERT_NEAR( 0.44626242465329269, ynew( 0 ), 1.0e-14 );
    ASSERT_NEAR( 0.98177874060900239, ynew( 1 ), 1.0e-14 );
    ASSERT_NEAR( 0.58634017875503497, ynew( 2 ), 1.0e-14 );
    ASSERT_NEAR( -0.34305440473301307, ynew( 3 ), 1.0e-14 );
    ASSERT_NEAR( -0.96344050977208218, ynew( 4 ), 1.0e-14 );
    ASSERT_NEAR( -0.65619975149611487, ynew( 5 ), 1.0e-14 );
    ASSERT_EQ( 1, ynew.cols() );
    ASSERT_EQ( 6, ynew.rows() );
    ASSERT_EQ( xnew.cols(), ynew.cols() );
    ASSERT_EQ( xnew.rows(), ynew.rows() );

    // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
    // interpolation.
    rbf::vector rowsum = rbf.Hhat.rowwise().sum();

    assert( rowsum.rows() > 0 );

    for ( int i = 0; i < rowsum.rows(); i++ )
        ASSERT_NEAR( rowsum( i ), 1, 1.0e-14 );
}

TEST( RBFInterpolationTest, oneD_small_cpu )
{
    // Initialize variables

    int n_A = 7;
    int n_B = 6;

    matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

    // Location of original values
    x << 0, 1, 2, 3, 4, 5, 6;
    xnew << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;

    // Values
    y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

    bool polynomialTerm = true;
    bool cpu = true;
    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction, polynomialTerm, cpu );

    rbf.compute( x, xnew );
    rbf.interpolate( y, ynew );

    ASSERT_NEAR( 0.44626242465329269, ynew( 0 ), 1.0e-14 );
    ASSERT_NEAR( 0.98177874060900239, ynew( 1 ), 1.0e-14 );
    ASSERT_NEAR( 0.58634017875503497, ynew( 2 ), 1.0e-14 );
    ASSERT_NEAR( -0.34305440473301307, ynew( 3 ), 1.0e-14 );
    ASSERT_NEAR( -0.96344050977208218, ynew( 4 ), 1.0e-14 );
    ASSERT_NEAR( -0.65619975149611487, ynew( 5 ), 1.0e-14 );
    ASSERT_EQ( 1, ynew.cols() );
    ASSERT_EQ( 6, ynew.rows() );
    ASSERT_EQ( xnew.cols(), ynew.cols() );
    ASSERT_EQ( xnew.rows(), ynew.rows() );

    // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
    // interpolation.
    rbf::vector rowsum = rbf.Hhat.rowwise().sum();

    ASSERT_EQ( rowsum.rows(), 0 );
}

TEST( RBFInterpolationTest, oneD_small_directly_interpolate )
{
    // Initialize variables

    int n_A = 7;
    int n_B = 6;

    matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

    // Location of original values
    x << 0, 1, 2, 3, 4, 5, 6;
    xnew << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;

    // Values
    y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    // rbf.interpolate( x, xnew, y, ynew );
    rbf.compute( x, xnew );
    rbf.interpolate( y, ynew );

    ASSERT_NEAR( 0.44626242465329269, ynew( 0 ), 1.0e-14 );
    ASSERT_NEAR( 0.98177874060900239, ynew( 1 ), 1.0e-14 );
    ASSERT_NEAR( 0.58634017875503497, ynew( 2 ), 1.0e-14 );
    ASSERT_NEAR( -0.34305440473301307, ynew( 3 ), 1.0e-14 );
    ASSERT_NEAR( -0.96344050977208218, ynew( 4 ), 1.0e-14 );
    ASSERT_NEAR( -0.65619975149611487, ynew( 5 ), 1.0e-14 );
    ASSERT_EQ( 1, ynew.cols() );
    ASSERT_EQ( 6, ynew.rows() );
    ASSERT_EQ( xnew.cols(), ynew.cols() );
    ASSERT_EQ( xnew.rows(), ynew.rows() );
}

TEST( RBFInterpolationTest, oneD_smallSamePosition )
{
    // Initialize variables

    int n_A = 7;
    int n_B = n_A;

    matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

    x << 0, 1, 2, 3, 4, 5, 6;
    xnew = x;

    y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    rbf.compute( x, xnew );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < ynew.rows(); i++ )
        ASSERT_NEAR( y( i ), ynew( i ), 1.0e-14 );

    ASSERT_EQ( 1, ynew.cols() );
    ASSERT_EQ( 7, ynew.rows() );
    ASSERT_EQ( xnew.cols(), ynew.cols() );
    ASSERT_EQ( xnew.rows(), ynew.rows() );

    // Verify consistent interpolation. Rowsum of Hhat is one for a consistent
    // interpolation.
    rbf::vector rowsum = rbf.Hhat.rowwise().sum();

    for ( int i = 0; i < rowsum.rows(); i++ )
        ASSERT_NEAR( rowsum( i ), 1, 1.0e-14 );
}

TEST( RBFInterpolationTest, oneD_smallSamePosition_directlyInterpolate )
{
    // Initialize variables

    int n_A = 7;
    int n_B = n_A;

    matrix x( n_A, 1 ), y( n_A, 1 ), xnew( n_B, 1 ), ynew( n_B, 1 );

    x << 0, 1, 2, 3, 4, 5, 6;
    xnew = x;

    y << 0, 8.414709848078965e-01, 9.092974268256817e-01, 1.411200080598672e-01, -7.568024953079282e-01, -9.589242746631385e-01, -2.794154981989259e-01;

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    rbf.compute( x, xnew );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < ynew.rows(); i++ )
        ASSERT_NEAR( y( i ), ynew( i ), 1.0e-14 );

    ASSERT_EQ( 1, ynew.cols() );
    ASSERT_EQ( 7, ynew.rows() );
    ASSERT_EQ( xnew.cols(), ynew.cols() );
    ASSERT_EQ( xnew.rows(), ynew.rows() );
}

TEST( RBFInterpolationTest, rbf1d )
{
    // Check that the RBF function interpolates through the nodes (1D)

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    matrix x( 9, 1 ), y( 9, 1 ), ynew;
    x.col( 0 ).setLinSpaced( 9, 0, 10 );

    y.col( 0 ) = x.array().sin();

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, rbf1d_directly_interpolate )
{
    // Check that the RBF function interpolates through the nodes (1D)

    matrix x( 9, 1 ), y( 9, 1 ), ynew;
    x.col( 0 ).setLinSpaced( 9, 0, 10 );

    y.col( 0 ) = x.array().sin();

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, rbf2d )
{
    // Check that the RBF function interpolates through the nodes (2D)

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    matrix x, y, ynew;

    x = fsi::matrix::Random( 50, 2 ).array() * 4 - 2;
    y = fsi::matrix::Random( 50, 2 ).array() * 4 - 2;

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-10 );
}

TEST( RBFInterpolationTest, rbf2d_directly_interpolate )
{
    // Check that the RBF function interpolates through the nodes (2D)

    matrix x, y, ynew;

    x = fsi::matrix::Random( 50, 2 ).array() * 4 - 2;
    y = fsi::matrix::Random( 50, 2 ).array() * 4 - 2;

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-10 );
}

TEST( RBFInterpolationTest, rbf3d )
{
    // Check that the RBF function interpolates through the nodes (2D)

    std::shared_ptr<RBFFunctionInterface> rbfFunction( new TPSFunction() );
    RBFInterpolation rbf( rbfFunction );

    matrix x, y, ynew;

    x = fsi::matrix::Random( 50, 3 ).array() * 4 - 2;
    y = fsi::matrix::Random( 50, 3 ).array() * 4 - 2;

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-11 );
}

TEST( RBFInterpolationTest, wendlandC6 )
{
    std::shared_ptr<RBFFunctionInterface> rbfFunction( new WendlandC6Function( 5 ) );
    RBFInterpolation rbf( rbfFunction );

    matrix x, y, ynew;

    x = fsi::matrix::Random( 50, 3 ).array() * 50.0 - 2.0;
    y = fsi::matrix::Random( 50, 3 ).array() * 32.0 - 2.0;

    rbf.compute( x, x );
    rbf.interpolate( y, ynew );

    for ( int i = 0; i < y.rows(); i++ )
        for ( int j = 0; j < y.cols(); j++ )
            ASSERT_NEAR( y( i, j ), ynew( i, j ), 1.0e-12 );
}

TEST( RBFInterpolationTest, wendlandC6Unit )
{
    std::shared_ptr<RBFFunctionInterface> rbfFunction( new WendlandC6Function( 1 ) );

    ASSERT_NEAR( rbfFunction->evaluate( 2 ), 0, 1.0e-13 );
    ASSERT_NEAR( rbfFunction->evaluate( 2.0e5 ), 0, 1.0e-13 );
    ASSERT_NEAR( rbfFunction->evaluate( 0.5 ), 0.0595703125, 1.0e-9 );
    ASSERT_NEAR( rbfFunction->evaluate( 0.69 ), 0.00246782213555, 1.0e-9 );
}
