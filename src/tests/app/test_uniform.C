
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
    ASSERT_NEAR( nodes[0], 0.2, tol );
    ASSERT_NEAR( nodes[1], 0.4, tol );
    ASSERT_NEAR( nodes[2], 0.6, tol );
    ASSERT_NEAR( nodes[3], 0.8, tol );
    ASSERT_NEAR( nodes[4], 1, tol );
}

TEST( UniformTest, scalar )
{
    const size_t num_nodes = 5;
    fsi::quadrature::Uniform<scalar> Uniform( num_nodes );
    const std::vector<scalar> nodes = Uniform.get_nodes();

    double tol = 1.0e-13;
    ASSERT_NEAR( nodes[0], 0.2, tol );
    ASSERT_NEAR( nodes[1], 0.4, tol );
    ASSERT_NEAR( nodes[2], 0.6, tol );
    ASSERT_NEAR( nodes[3], 0.8, tol );
    ASSERT_NEAR( nodes[4], 1, tol );
}

TEST( UniformTest, api )
{
    const size_t num_nodes = 5;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );
    Uniform.get_nodes();
    Uniform.get_q_mat();
    Uniform.get_s_mat();
    Uniform.get_b_mat();
    Uniform.get_q_vec();
    ASSERT_EQ( Uniform.get_num_nodes(), num_nodes );
}

TEST( UniformTest, expectedError )
{
    const size_t num_nodes = 7;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );

    double error = Uniform.expected_error();
    std::cout << "expected error = " << error << std::endl;
}

TEST( UniformTest, smat )
{
    const size_t num_nodes = 3;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );

    const Matrix<double> smat2 = Uniform.get_s_mat();

    Eigen::MatrixXd smat( 3, 3 );

    ASSERT_EQ( smat.rows(), smat2.rows() );
    ASSERT_EQ( smat.cols(), smat2.cols() );

    smat( 0, 0 ) = 0.638888888888888923583358408425041369742072646527449507691784663;
    smat( 0, 1 ) = -0.4444444444444444444444444444444444444444444444444444444444444;
    smat( 0, 2 ) = 0.138888888888888888888888888888888888888888888888888888888888888;
    smat( 1, 0 ) = 0.138888888888888798683268138094892438670611119028631280001359;
    smat( 1, 1 ) = 0.22222222222222222222222222222222222222222222222222222222222222222;
    smat( 1, 2 ) = -0.0277777777777777777777777777777777777777777777777777777777777777777777777777777;
    smat( 2, 0 ) = -0.0277777777777776598165814113548593428769530018066716738479321433915706453869225808;
    smat( 2, 1 ) = 0.222222222222222222222222222222222222222222222222222222222222222222222222222;
    smat( 2, 2 ) = 0.138888888888888763988798618558771157882526434622147853000042306956132887043;

    for ( int i = 0; i < smat.rows(); i++ )
        for ( int j = 0; j < smat.cols(); j++ )
            ASSERT_NEAR( smat(i,j), smat2(i,j), 1.0e-13 );
}

TEST( UniformTest, qmat )
{
    const size_t num_nodes = 3;
    fsi::quadrature::Uniform<double> Uniform( num_nodes );

    const Matrix<double> qmat2 = Uniform.get_q_mat();

    Eigen::MatrixXd qmat( 3, 3 );

    ASSERT_EQ( qmat.rows(), qmat2.rows() );
    ASSERT_EQ( qmat.cols(), qmat2.cols() );

    qmat( 0, 0 ) = 0.6388888888888889235833584084250413697420726465274495076917846637083;
    qmat( 0, 1 ) = -0.444444444444444444444444444444444444444444444444444444444444444444444444;
    qmat( 0, 2 ) = 0.13888888888888888888888888888888888888888888888888888888888888888888888;
    qmat( 1, 0 ) = 0.77777777777777772226662654651993380841268376555608078769314453806662148;
    qmat( 1, 1 ) = -0.2222222222222222222222222222222222222222222222222222222222222222222222222222;
    qmat( 1, 2 ) = 0.11111111111111111111111111111111111111111111111111111111111111111111111111;
    qmat( 2, 0 ) = 0.750000000000000062450045135165074465535730763749409113845212394675050834;
    qmat( 2, 1 ) = 0;
    qmat( 2, 2 ) = 0.24999999999999987509990972966988226899363754573325896411115341806724399815;

    for ( int i = 0; i < qmat.rows(); i++ )
        for ( int j = 0; j < qmat.cols(); j++ )
            ASSERT_NEAR( qmat(i,j), qmat2(i,j), 1.0e-13 );
}
