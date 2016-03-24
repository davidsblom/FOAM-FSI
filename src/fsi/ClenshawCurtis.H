
#pragma once

#include <cassert>
#include <vector>
#include "QuadratureInterface.H"

namespace fsi
{
    namespace quadrature
    {
        using namespace std;

        /**
         * Quadrature handler for Clenshaw-Curtis quadrature.
         *
         * @tparam scalar precision of quadrature (i.e. `double`)
         *
         * @since v0.3.0
         */
        template<typename precision>
        class ClenshawCurtis
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
                explicit ClenshawCurtis( const size_t num_nodes );
                ClenshawCurtis() = default;
                virtual ~ClenshawCurtis() = default;

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
    }
}

#include "ClenshawCurtis.tpp"
