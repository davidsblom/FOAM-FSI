
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "fvCFD.H"
#include "gtest/gtest.h"
#include "Uniform.H"

TEST( UniformTest, object )
{
    const size_t num_nodes = 2;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );
}

TEST( UniformTest, nodes )
{
    const size_t num_nodes = 5;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );
    const std::vector<double> nodes = Uniform.get_nodes();

    double tol = 1.0e-13;
    ASSERT_NEAR( nodes[0], 0, tol );
    ASSERT_NEAR( nodes[1], 0.2, tol );
    ASSERT_NEAR( nodes[2], 0.4, tol );
    ASSERT_NEAR( nodes[3], 0.6, tol );
    ASSERT_NEAR( nodes[4], 0.8, tol );
    ASSERT_NEAR( nodes[5], 1, tol );
}

TEST( UniformTest, scalar )
{
    const size_t num_nodes = 5;
    fsi::quadrature::Uniform<scalar> Uniform( num_nodes );
    const std::vector<scalar> nodes = Uniform.get_nodes();

    double tol = 1.0e-13;
    ASSERT_NEAR( nodes[0], 0, tol );
    ASSERT_NEAR( nodes[1], 0.2, tol );
    ASSERT_NEAR( nodes[2], 0.4, tol );
    ASSERT_NEAR( nodes[3], 0.6, tol );
    ASSERT_NEAR( nodes[4], 0.8, tol );
    ASSERT_NEAR( nodes[5], 1, tol );
}

TEST( UniformTest, api )
{
    const size_t num_nodes = 5;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );
    Uniform.get_nodes();
    Uniform.get_q_mat();
    Uniform.get_s_mat();
    ASSERT_EQ( Uniform.get_num_nodes(), num_nodes );
}

TEST( UniformTest, smat )
{
    const size_t num_nodes = 3;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );

    const Matrix<double> smat2 = Uniform.get_s_mat();

    Eigen::MatrixXd smat( 3, 4 );

    smat( 0, 0 ) = 0;
    smat( 0, 1 ) = 0.6388888888888889235833584084250413697420726465274495076917846637083615748861992409406796930069108539;
    smat( 0, 2 ) = -0.4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444;
    smat( 0, 3 ) = 0.1388888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889;
    smat( 1, 0 ) = 0;
    smat( 1, 1 ) = 0.1388888888888887986832681380948924386706111190286312800013598743582599052958819735542327981820317799;
    smat( 1, 2 ) = 0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222;
    smat( 1, 3 ) = -0.02777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777778;
    smat( 2, 0 ) = 0;
    smat( 2, 1 ) = -0.02777777777777765981658141135485934287695300180667167384793214339157064538692258080168904377650309685;
    smat( 2, 2 ) = 0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222;
    smat( 2, 3 ) = 0.1388888888888887639887986185587711578825264346221478530000423069561328870437259517330245730385413901;

    ASSERT_EQ( smat.rows(), smat2.rows() );
    ASSERT_EQ( smat.cols(), smat2.cols() );

    for ( int i = 0; i < smat.rows(); i++ )
        for ( int j = 0; j < smat.cols(); j++ )
            ASSERT_NEAR( smat(i,j), smat2(i,j), 1.0e-13 );
}

TEST( UniformTest, qmat )
{
    const size_t num_nodes = 3;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );

    const Matrix<double> qmat2 = Uniform.get_q_mat();

    Eigen::MatrixXd qmat( 3, 4 );

    ASSERT_EQ( qmat.rows(), qmat2.rows() );
    ASSERT_EQ( qmat.cols(), qmat2.cols() );

    qmat( 0, 0 ) = 0;
    qmat( 0, 1 ) = 0.6388888888888889235833584084250413697420726465274495076917846637083615748861992409406796930069108539;
    qmat( 0, 2 ) = -0.4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444;
    qmat( 0, 3 ) = 0.1388888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889;
    qmat( 1, 0 ) = 0;
    qmat( 1, 1 ) = 0.7777777777777777222666265465199338084126837655560807876931445380666214801820812144949124911889426338;
    qmat( 1, 2 ) = -0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222;
    qmat( 1, 3 ) = 0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111;
    qmat( 2, 0 ) = 0;
    qmat( 2, 1 ) = 0.7500000000000000624500451351650744655357307637494091138452123946750508347951586336932234474124395370;
    qmat( 2, 2 ) = 0;
    qmat( 2, 3 ) = 0.2499999999999998750999097296698822689936375457332589641111534180672439981548370628441356841496525012;

    for ( int i = 0; i < qmat.rows(); i++ )
        for ( int j = 0; j < qmat.cols(); j++ )
            ASSERT_NEAR( qmat(i,j), qmat2(i,j), 1.0e-13 );
}
