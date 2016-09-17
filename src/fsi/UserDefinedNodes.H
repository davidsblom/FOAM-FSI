
/*
 * Copyright [2016] <David Blom>
 */

#pragma once

#include <vector>
#include "QuadratureInterface.H"

namespace fsi {
namespace quadrature {
template<typename precision>
class UserDefinedNodes
      : public IQuadrature<precision>{
    protected:
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;

    public:
        explicit UserDefinedNodes(const std::vector<precision> & nodes);
        UserDefinedNodes() = default;
        virtual ~UserDefinedNodes() = default;

        bool left_is_node() const override;

        bool right_is_node() const override;
};
}     // namespace quadrature
} // namespace fsi

#include "UserDefinedNodes.tpp"
