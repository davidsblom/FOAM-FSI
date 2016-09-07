
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <cassert>
#include <vector>
#include "QuadratureInterface.H"

namespace fsi {
namespace quadrature {
/**
 * Quadrature handler for Clenshaw-Curtis quadrature.
 *
 * @tparam scalar precision of quadrature (i.e. `double`)
 *
 * @since v0.3.0
 */
template<typename precision>
class ClenshawCurtis
      : public IQuadrature<precision>{
    protected:
        // ! @{
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;

        // ! @}

    public:
        // ! @{

        /**
         * @throws invalid_argument if less than two nodes are requested
         */
        explicit ClenshawCurtis(const size_t num_nodes);
        ClenshawCurtis() = default;
        virtual ~ClenshawCurtis() = default;

        // ! @}

        // ! @{
        bool left_is_node() const override;

        bool right_is_node() const override;

        // ! @}

    protected:
        // ! @{
        void compute_nodes() override;

        // ! @}
};
}     // namespace quadrature
} // namespace fsi

#include "ClenshawCurtis.tpp"
