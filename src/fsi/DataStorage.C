
/*
 * Copyright [2016] <David Blom>
 */

#include <memory>
#include <vector>
#include "DataStorage.H"

namespace sdc {
DataStorage::DataStorage(std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
    int N
    )
    :
    quadrature(quadrature),
    F(quadrature->get_num_nodes(), N),
    solStages(quadrature->get_num_nodes(), N) {
    assert(quadrature);
}

DataStorage::~DataStorage()
{}

const fsi::matrix DataStorage::integrate(const std::vector<scalar> & nodes,
    scalar dt
    ) const {
    std::vector<scalar> quadratureNodes = quadrature->get_nodes();

    if (not quadrature->left_is_node())
        quadratureNodes.erase(quadratureNodes.begin());

    fsi::quadrature::Matrix<scalar> q_matrix = fsi::quadrature::compute_q_matrix(quadratureNodes, nodes);

    if (not quadrature->left_is_node()) {
        fsi::quadrature::Matrix<scalar> q_mat2(q_matrix.rows(), q_matrix.cols() + 1);
        q_mat2.setZero();
        q_mat2.rightCols(q_matrix.cols()) = q_matrix;
        q_matrix = q_mat2;
    }

    fsi::quadrature::Matrix<scalar> data = dt * q_matrix * getFunctions();

    for (int i = 0; i < data.rows(); i++)
        data.row(i) += getSolutions().row(0);

    return data;
}

template<typename scalar>
static fsi::quadrature::Matrix<scalar> compute_interpolation_matrix(const std::vector<scalar> & from,
    const std::vector<scalar> & to
    ) {
    const size_t to_size = to.size();
    const size_t from_size = from.size();
    assert(to_size >= 1 && from_size >= 1);

    fsi::quadrature::Matrix<scalar> q_mat = fsi::quadrature::Matrix<scalar>::Zero(to_size, from_size);

    for (size_t m = 0; m < from_size; ++m) {
        fsi::quadrature::Polynomial<scalar> p = fsi::quadrature::build_polynomial(m, from);

        auto den = p.evaluate(from[m]);

        for (size_t j = 0; j < to_size; ++j) {
            q_mat(j, m) = p.evaluate(to[j]) / den;
        }
    }

    return q_mat;
}

const fsi::matrix DataStorage::interpolate(const fsi::matrix functions,
    const std::vector<scalar> & nodes
    ) const {
    std::vector<scalar> quadratureNodes = quadrature->get_nodes();
    fsi::quadrature::Matrix<scalar> q_matrix = compute_interpolation_matrix(quadratureNodes, nodes);

    return q_matrix * functions;
}

void DataStorage::copyFunctions() {
    Fold = F;
}

const fsi::matrix & DataStorage::getFunctions() const {
    assert(F.rows() > 0);
    return F;
}

const fsi::matrix & DataStorage::getOldFunctions() const {
    assert(Fold.rows() > 0);
    return Fold;
}

const fsi::matrix & DataStorage::getSolutions() const {
    assert(solStages.rows() > 0);
    return solStages;
}

const fsi::vector DataStorage::getFunction(int substep) const {
    assert(substep <= F.rows());
    assert(F.cols() > 0);
    assert(F.rows() > 0);
    return F.row(substep);
}

const fsi::vector DataStorage::getSolution(int substep) const {
    assert(substep <= solStages.rows());
    return solStages.row(substep);
}

const fsi::vector DataStorage::getLastSolution() const {
    return getSolution(quadrature->get_num_nodes() - 1);
}

void DataStorage::initialize(int k,
    int N
    ) {
    assert(N >= 0);
    assert(k >= 2);
    F.resize(k, N);
    solStages.resize(k, N);
    F.setZero();
    solStages.setZero();
}

void DataStorage::storeFunction(const fsi::vector & f,
    int substep
    ) {
    assert(f.rows() == F.cols());
    assert(substep <= F.rows());
    assert(not std::isnan(f.norm()));
    F.row(substep) = f;
}

void DataStorage::storeSolution(const fsi::vector & sol,
    int substep
    ) {
    assert(sol.rows() == solStages.cols());
    assert(substep <= solStages.rows());
    assert(not std::isnan(sol.norm()));
    solStages.row(substep) = sol;
}
} // namespace sdc
