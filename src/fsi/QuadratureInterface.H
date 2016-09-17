
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include "Polynomial.H"

namespace fsi {
namespace quadrature {
template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename scalar>
using Index = typename Matrix<scalar>::Index;

/**
 * Quadrature type descriptors.
 */
enum class QuadratureType : int {
    GaussLegendre = 0         // !< @ref pfasst::quadrature::GaussLegendre "Gauss-Legendre" quadrature
    , GaussLobatto = 1         // !< @ref pfasst::quadrature::GaussLobatto "Gauss-Lobatto" quadrature
    , GaussRadau = 2          // !< @ref pfasst::quadrature::GaussRadau "Gauss-Radau" quadrature
    , ClenshawCurtis = 3         // !< @ref pfasst::quadrature::ClenshawCurtis "Clenshaw-Curtis" quadrature
    , Uniform = 4             // !< Uniform quadrature
    , UNDEFINED = -1
};


template<typename scalar>
static Polynomial<scalar> build_polynomial(const size_t node,
    const std::vector<scalar> & nodes
    ) {
    const size_t num_nodes = nodes.size();
    Polynomial<scalar> p(num_nodes + 1), p1(num_nodes + 1);
    p[0] = 1.0;

    for (size_t m = 0; m < num_nodes; ++m) {
        if (m == node) {
            continue;
        }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = scalar(0.0);

        for (size_t j = 0; j < num_nodes; ++j) {
            p1[j + 1] = p[j];
        }

        for (size_t j = 0; j < num_nodes + 1; ++j) {
            p1[j] -= p[j] * nodes[m];
        }

        for (size_t j = 0; j < num_nodes + 1; ++j) {
            p[j] = p1[j];
        }
    }

    return p;
}

/**
 * Compute quadrature matrix \\( Q \\) between two sets of nodes.
 *
 * Computing the quadrature matrix \\( Q \\) for polynomial-based integration from one set of
 * quadrature nodes (@p from) to another set of quadrature nodes (@p to).
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] from first set of quadrature nodes
 * @param[in] to second set of quadrature nodes
 * @returns quadrature matrix \\( Q \\) with `to.size()` rows and `from.size()` colums
 *
 * @pre For correctness of the algorithm it is assumed, that both sets of nodes are in the range
 *   \\( [0, 1] \\).
 *
 * @since v0.3.0
 */
template<typename scalar>
static Matrix<scalar> compute_q_matrix(const std::vector<scalar> & from,
    const std::vector<scalar> & to
    ) {
    const size_t to_size = to.size();
    const size_t from_size = from.size();
    assert(to_size >= 1 && from_size >= 1);

    Matrix<scalar> q_mat = Matrix<scalar>::Zero(to_size, from_size);

    for (size_t m = 0; m < from_size; ++m) {
        Polynomial<scalar> p = build_polynomial(m, from);

        // evaluate integrals
        auto den = p.evaluate(from[m]);
        auto P = p.integrate();

        for (size_t j = 0; j < to_size; ++j) {
            q_mat(j, m) = (P.evaluate(to[j]) - P.evaluate(0.0)) / den;
        }
    }

    return q_mat;
}

/**
 * Compute quadrature matrix \\( Q \\) for one set of nodes.
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] nodes quadrature nodes to compute \\( Q \\) matrix for
 *
 * @overload
 */
template<typename scalar>
static Matrix<scalar> compute_q_matrix(const std::vector<scalar> & nodes) {
    return compute_q_matrix<scalar>(nodes, nodes);
}

/**
 * Compute quadrature matrix \\( Q \\) from a given node-to-node quadrature matrix \\( S \\).
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] s_mat \\( S \\) matrix to compute \\( Q \\) from
 * @see pfasst::quadrature::compute_s_matrix
 *
 * @since v0.3.0
 *
 * @overload
 */
template<typename scalar>
static Matrix<scalar> compute_q_matrix(const Matrix<scalar> & s_mat) {
    Matrix<scalar> q_mat = Matrix<scalar>::Zero(s_mat.rows(), s_mat.cols());
    q_mat.col(0) = s_mat.col(0);

    for (Index<scalar> q_mat_col = 1; q_mat_col < q_mat.cols(); ++q_mat_col) {
        q_mat.col(q_mat_col) = q_mat.col(q_mat_col - 1) + s_mat.col(q_mat_col);
    }

    return q_mat;
}

/**
 * Compute node-to-node quadrature matrix \\( S \\) from a given quadrature matrix \\( Q \\).
 *
 * The \\( S \\) matrix provides a node-to-node quadrature where the \\( i \\)-th row of
 * \\( S \\) represents a quadrature from the \\( i-1 \\)-th node to the \\( i \\)-th node.
 *
 * The procedure is simply subtracting the \\( i-1 \\)-th row of \\( Q \\) from the
 * \\( i \\)-th row of \\( Q \\).
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] q_mat \\( Q \\) matrix to compute \\( S \\) of
 * @returns \\( S \\) matrix
 *
 * @since v0.3.0
 */
template<typename scalar>
static Matrix<scalar> compute_s_matrix(const Matrix<scalar> & q_mat) {
    Matrix<scalar> s_mat = Matrix<scalar>::Zero(q_mat.rows(), q_mat.cols());
    s_mat.row(0) = q_mat.row(0);

    for (Index<scalar> row = 1; row < s_mat.rows(); ++row) {
        s_mat.row(row) = q_mat.row(row) - q_mat.row(row - 1);
    }

    return s_mat;
}

/**
 * Compute node-to-node quadrature matrix \\( S \\) from two given sets of nodes
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] from first set of quadrature nodes
 * @param[in] to second set of quadrature nodes
 *
 * @since v0.3.0
 *
 * @overload
 */
template<typename scalar>
static Matrix<scalar> compute_s_matrix(const std::vector<scalar> & from,
    const std::vector<scalar> & to
    ) {
    return compute_s_matrix(compute_q_matrix(from, to));
}

/**
 * Compute vector \\( q \\) for integration from \\( 0 \\) to \\( 1 \\) for given set of nodes.
 *
 * This equals to the last row of the quadrature matrix \\( Q \\) for the given set of nodes.
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 * @param[in] nodes quadrature nodes to compute \\( Q \\) matrix for
 * @pre For correctness of the algorithm it is assumed, that the nodes are in the range
 *   \\( [0, 1] \\).
 *
 * @since v0.3.0
 */
template<typename scalar>
static std::vector<scalar> compute_q_vec(const std::vector<scalar> & nodes) {
    const size_t num_nodes = nodes.size();
    assert(num_nodes >= 1);

    std::vector<scalar> q_vec = std::vector<scalar>(num_nodes, scalar(0.0));

    for (size_t m = 0; m < num_nodes; ++m) {
        Polynomial<scalar> p = build_polynomial(m, nodes);

        // evaluate integrals
        auto den = p.evaluate(nodes[m]);
        auto P = p.integrate();
        q_vec[m] = (P.evaluate(scalar(1.0)) - P.evaluate(scalar(0.0))) / den;
    }

    return q_vec;
}

/**
 * Interface for quadrature handlers.
 *
 * Quadrature handlers provide \\( Q \\), \\( S \\) and \\( B \\) matrices respecting the left
 * and right nodes, i.e. whether \\( 0 \\) and \\( 1 \\) are part of the nodes or not.
 *
 * Computation of the quadrature nodes and matrices (i.e. quadrature weights) is done on
 * initialization.
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 *
 * @since v0.3.0
 */
template<typename precision>
class IQuadrature {
    protected:
        // ! @{
        static const bool LEFT_IS_NODE = false;
        static const bool RIGHT_IS_NODE = false;

        size_t num_nodes;
        Matrix<precision> q_mat;
        Matrix<precision> s_mat;
        std::vector<precision> q_vec;         // XXX: black spot
        Matrix<precision> b_mat;
        std::vector<precision> nodes;

        // ! @}

    public:
        // ! @{

        /**
         * @throws invalid_argument if number of nodes is invalid for quadrature type
         */
        explicit IQuadrature(const size_t num_nodes);

        /**
         * @throws invalid_argument if number of nodes is invalid for quadrature type
         */
        IQuadrature();
        virtual ~IQuadrature() = default;

        // ! @}

        // ! @{
        virtual const Matrix<precision> & get_q_mat() const;

        virtual const Matrix<precision> & get_s_mat() const;

        virtual const std::vector<precision> & get_nodes() const;

        virtual size_t get_num_nodes() const;

        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual bool left_is_node() const;

        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual bool right_is_node() const;

    protected:
        // ! @{

        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual void compute_nodes();

        virtual void compute_weights();

        // ! @}
};
}     // namespace quadrature
} // namespace fsi

#include "QuadratureInterface.tpp"
