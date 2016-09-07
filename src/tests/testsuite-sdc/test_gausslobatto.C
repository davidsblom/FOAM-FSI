
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "fvCFD.H"
#include "gtest/gtest.h"
#include "GaussLobatto.H"

TEST(GaussLobatto, object)
{
    const size_t num_nodes = 2;
    fsi::quadrature::GaussLobatto<double> Uniform(num_nodes);
}

TEST(GaussLobatto, nodes)
{
    const size_t num_nodes = 5;
    fsi::quadrature::GaussLobatto<double> Uniform(num_nodes);
    const std::vector<double> nodes = Uniform.get_nodes();

    double tol = 1.0e-13;

    ASSERT_NEAR(nodes[0], 0, tol);
    ASSERT_NEAR(nodes[1], 0.1726731646460114281008537718765708222153959588022877212423398, tol);
    ASSERT_NEAR(nodes[2], 0.5, tol);
    ASSERT_NEAR(nodes[3], 0.8273268353539885718991462281234291777846040, tol);
    ASSERT_NEAR(nodes[4], 1, tol);
}

TEST(GaussLobatto, scalar)
{
    const size_t num_nodes = 5;
    fsi::quadrature::GaussLobatto<scalar> gaussLobatto(num_nodes);
    const std::vector<scalar> nodes = gaussLobatto.get_nodes();

    double tol = 1.0e-13;

    ASSERT_NEAR(nodes[0], 0, tol);
    ASSERT_NEAR(nodes[1], 0.1726731646460114281008537718765708222153959588022877212423398, tol);
    ASSERT_NEAR(nodes[2], 0.5, tol);
    ASSERT_NEAR(nodes[3], 0.8273268353539885718991462281234291777846040, tol);
    ASSERT_NEAR(nodes[4], 1, tol);
}

TEST(GaussLobatto, api)
{
    const size_t num_nodes = 5;
    fsi::quadrature::GaussLobatto<double> gaussLobatto(num_nodes);
    gaussLobatto.get_nodes();
    gaussLobatto.get_q_mat();
    gaussLobatto.get_s_mat();
    ASSERT_EQ(gaussLobatto.get_num_nodes(), num_nodes);
}

TEST(GaussLobatto, smat)
{
    const size_t num_nodes = 3;
    fsi::quadrature::GaussLobatto<double> gaussLobatto(num_nodes);

    const fsi::quadrature::Matrix<double> smat2 = gaussLobatto.get_s_mat();

    for (auto i : gaussLobatto.get_nodes())
        std::cout << i << std::endl;

    Eigen::MatrixXd smat(2, 3);

    smat(0, 0) = 0.2083333333333333148296162562473909929394721984863281250000000000000000000000000000000000000000000000;
    smat(0, 1) = 0.3333333333333333703407674875052180141210556030273437500000000000000000000000000000000000000000000000;
    smat(0, 2) = -0.04166666666666667129259593593815225176513195037841796875000000000000000000000000000000000000000000000;
    smat(1, 0) = -0.04166666666666668517038374375260900706052780151367187500000000000000000000000000000000000000000000000;
    smat(1, 1) = 0.3333333333333333703407674875052180141210556030273437500000000000000000000000000000000000000000000000;
    smat(1, 2) = 0.2083333333333333009518284484329342376440763473510742187500000000000000000000000000000000000000000000;

    std::cout << "smat = \n" << smat << std::endl;
    std::cout << "smat2 = \n" << smat2 << std::endl;

    ASSERT_EQ(smat.rows(), smat2.rows());
    ASSERT_EQ(smat.cols(), smat2.cols());

    for (int i = 0; i < smat.rows(); i++)
        for (int j = 0; j < smat.cols(); j++)
            ASSERT_NEAR(smat(i, j), smat2(i, j), 1.0e-13);
}

TEST(GaussLobatto, qmat)
{
    const size_t num_nodes = 3;
    fsi::quadrature::GaussLobatto<double> gaussLobatto(num_nodes);

    const fsi::quadrature::Matrix<double> qmat2 = gaussLobatto.get_q_mat();

    Eigen::MatrixXd qmat(2, 3);

    ASSERT_EQ(qmat.rows(), qmat2.rows());
    ASSERT_EQ(qmat.cols(), qmat2.cols());

    qmat(0, 0) = 0.2083333333333333148296162562473909929394721984863281250000000000000000000000000000000000000000000000;
    qmat(0, 1) = 0.3333333333333333703407674875052180141210556030273437500000000000000000000000000000000000000000000000;
    qmat(0, 2) = -0.04166666666666667129259593593815225176513195037841796875000000000000000000000000000000000000000000000;
    qmat(1, 0) = 0.1666666666666666296592325124947819858789443969726562500000000000000000000000000000000000000000000000;
    qmat(1, 1) = 0.6666666666666667406815349750104360282421112060546875000000000000000000000000000000000000000000000000;
    qmat(1, 2) = 0.1666666666666666296592325124947819858789443969726562500000000000000000000000000000000000000000000000;

    for (int i = 0; i < qmat.rows(); i++)
        for (int j = 0; j < qmat.cols(); j++)
            ASSERT_NEAR(qmat(i, j), qmat2(i, j), 1.0e-13);
}

TEST(GaussLobatto, qmat2)
{
    const size_t num_nodes = 2;
    fsi::quadrature::GaussLobatto<double> gaussLobatto(num_nodes);

    const fsi::quadrature::Matrix<double> qmat2 = gaussLobatto.get_q_mat();

    Eigen::MatrixXd qmat(1, 2);

    ASSERT_EQ(qmat.rows(), qmat2.rows());
    ASSERT_EQ(qmat.cols(), qmat2.cols());

    qmat(0, 0) = 0.5;
    qmat(0, 1) = 0.5;

    for (int i = 0; i < qmat.rows(); i++)
        for (int j = 0; j < qmat.cols(); j++)
            ASSERT_NEAR(qmat(i, j), qmat2(i, j), 1.0e-13);
}

TEST(GaussLobatto, smat2)
{
    const size_t num_nodes = 2;
    fsi::quadrature::GaussLobatto<double> gaussLobatto(num_nodes);

    const fsi::quadrature::Matrix<double> smat2 = gaussLobatto.get_s_mat();

    Eigen::MatrixXd smat(1, 2);

    smat(0, 0) = 0.5;
    smat(0, 1) = 0.5;
    ASSERT_EQ(smat.rows(), smat2.rows());
    ASSERT_EQ(smat.cols(), smat2.cols());

    for (int i = 0; i < smat.rows(); i++)
        for (int j = 0; j < smat.cols(); j++)
            ASSERT_NEAR(smat(i, j), smat2(i, j), 1.0e-13);
}
