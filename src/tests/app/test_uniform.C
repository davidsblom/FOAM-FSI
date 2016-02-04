
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
