
#pragma once

#include "QuadratureInterface.H"


namespace fsi
{
    namespace quadrature
    {

        template<typename precision>
        class UserDefinedNodes
              : public IQuadrature<precision>
        {
            protected:
                static const bool LEFT_IS_NODE = true;
                static const bool RIGHT_IS_NODE = true;

            public:

                explicit UserDefinedNodes( const std::vector<precision> nodes );
                UserDefinedNodes() = default;
                virtual ~UserDefinedNodes() = default;

                virtual bool left_is_node() const override;

                virtual bool right_is_node() const override;
        };
    }
}

#include "UserDefinedNodes.tpp"
