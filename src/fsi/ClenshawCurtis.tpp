
/*
 * Copyright [2016] <David Blom>
 */

#include <boost/math/constants/constants.hpp>
#include <stdexcept>
#include <vector>
#include "Polynomial.H"
#include "ClenshawCurtis.H"

namespace fsi {
namespace quadrature {
template<typename precision>
ClenshawCurtis<precision>::ClenshawCurtis(const size_t num_nodes)
    : IQuadrature<precision>(num_nodes) {
    if (this->num_nodes < 2) {
        throw std::invalid_argument("Clenshaw-Curtis quadrature requires at least two quadrature nodes.");
    }

    this->compute_nodes();
    this->compute_weights();
}

template<typename precision>
bool ClenshawCurtis<precision>::left_is_node() const {
    return LEFT_IS_NODE;
}

template<typename precision>
bool ClenshawCurtis<precision>::right_is_node() const {
    return RIGHT_IS_NODE;
}

template<typename precision>
void ClenshawCurtis<precision>::compute_nodes() {
    this->nodes = std::vector<precision>(this->num_nodes, precision(0.0));
    auto roots = Polynomial<precision>::legendre(this->num_nodes).roots();

    for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = 0.5 * (1.0 - std::cos(j * boost::math::constants::pi<precision>() / (this->num_nodes - 1)));
    }
}
}     // namespace quadrature
} // namespace fsi
