
#pragma once

#include "QuadratureInterface.H"


namespace fsi
{
    namespace quadrature
    {
        /**
         * Quadrature handler for uniform distributed nodes.
         *
         * @tparam scalar precision of quadrature (i.e. `double`)
         */
        template<typename precision>
        class Uniform
              : public IQuadrature<precision>
        {
            protected:
                // ! @{
                static const bool LEFT_IS_NODE = false;
                static const bool RIGHT_IS_NODE = true;

                // ! @}

            public:
                // ! @{

                /**
                 * @throws invalid_argument if less than two nodes are requested
                 */
                explicit Uniform( const size_t num_nodes );
                Uniform() = default;
                virtual ~Uniform() = default;

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

#include "Uniform.tpp"
