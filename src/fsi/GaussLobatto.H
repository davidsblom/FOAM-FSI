
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include "QuadratureInterface.H"

namespace fsi
{
    namespace quadrature
    {
        /**
         * Quadrature handler for Gauss-Lobatto quadrature.
         *
         * @tparam scalar precision of quadrature (i.e. `double`)
         */
        template<typename precision>
        class GaussLobatto
              : public IQuadrature<precision>
        {
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
                explicit GaussLobatto( const size_t num_nodes );
                GaussLobatto() = default;
                virtual ~GaussLobatto() = default;

                // ! @}

                // ! @{
                virtual bool left_is_node() const override;

                virtual bool right_is_node() const override;

                // ! @}

            protected:
                // ! @{
                virtual void compute_nodes() override;

                // ! @}
        };

} // namespace quadrature

} // namespace fsi

#include "GaussLobatto.tpp"
