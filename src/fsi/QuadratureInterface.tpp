
/*
 * Copyright [2016] <David Blom>
 */

#include <vector>
#include "QuadratureInterface.H"

namespace fsi {
namespace quadrature {
template<typename precision>
IQuadrature<precision>::IQuadrature(const size_t num_nodes)
    : num_nodes(num_nodes) {
    assert(this->num_nodes > 0);
}

template<typename precision>
IQuadrature<precision>::IQuadrature()
    : num_nodes(0)
{}

template<typename precision>
const Matrix<precision> & IQuadrature<precision>::get_q_mat() const {
    return this->q_mat;
}

template<typename precision>
const Matrix<precision> & IQuadrature<precision>::get_s_mat() const {
    return this->s_mat;
}

template<typename precision>
const std::vector<precision> & IQuadrature<precision>::get_nodes() const {
    return this->nodes;
}

template<typename precision>
size_t IQuadrature<precision>::get_num_nodes() const {
    return this->nodes.size();
}

template<typename precision>
bool IQuadrature<precision>::left_is_node() const {
    assert(false);
    return LEFT_IS_NODE;
}

template<typename precision>
bool IQuadrature<precision>::right_is_node() const {
    assert(false);
    return RIGHT_IS_NODE;
}

template<typename precision>
void IQuadrature<precision>::compute_nodes() {
    assert(false);
}

/**
 * @internals
 * Computing weights means computing \\( Q \\) and \\( S \\) matrices as well as the \\( q \\)
 * vector.
 * The \\( B \\) matrix is constructed from \\( q \\) vector.
 *
 * @note As long as the weight computation of the implemented quadrature is based on polynomial
 *   interpolation, this function must not be overwritten.
 * @endinternals
 */
template<typename precision>
void IQuadrature<precision>::compute_weights() {
    this->q_mat = compute_q_matrix(this->nodes);
    this->s_mat = compute_s_matrix(this->q_mat);
    this->q_vec = compute_q_vec(this->nodes);
    this->b_mat = Matrix<precision>::Zero(1, this->num_nodes);

    for (size_t i = 0; i < this->num_nodes; i++) {
        this->b_mat(0, i) = this->q_vec[i];
    }

    if (not left_is_node()) {
        this->nodes.insert(this->nodes.begin(), 0);

        Matrix<precision> q_mat2(this->q_mat.rows(), this->q_mat.cols() + 1);
        q_mat2.setZero();
        q_mat2.rightCols(this->q_mat.cols()) = this->q_mat;
        this->q_mat = q_mat2;

        Matrix<precision> s_mat2(this->s_mat.rows(), this->s_mat.cols() + 1);
        s_mat2.setZero();
        s_mat2.rightCols(this->s_mat.cols()) = this->s_mat;
        this->s_mat = s_mat2;
    }

    if (left_is_node()) {
        Matrix<precision> q_mat2(this->q_mat.rows() - 1, this->q_mat.cols());
        q_mat2 = this->q_mat.bottomRows(q_mat2.rows());
        this->q_mat = q_mat2;

        Matrix<precision> s_mat2(this->s_mat.rows() - 1, this->s_mat.cols());
        s_mat2 = this->s_mat.bottomRows(s_mat2.rows());
        this->s_mat = s_mat2;
    }
}
}     // namespace quadrature
}  // namespace fsi
