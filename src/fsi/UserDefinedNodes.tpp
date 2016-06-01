
#include "UserDefinedNodes.H"

#include <stdexcept>
#include <vector>
using namespace std;

#include "Polynomial.H"

namespace fsi
{
    namespace quadrature
    {
        template<typename precision>
        UserDefinedNodes<precision>::UserDefinedNodes( const std::vector<precision> & nodes )
            : IQuadrature<precision>( nodes.size() )
        {
            this->nodes = nodes;
            this->compute_weights();
        }

        template<typename precision>
        bool UserDefinedNodes<precision>::left_is_node() const
        {
            return LEFT_IS_NODE;
        }

        template<typename precision>
        bool UserDefinedNodes<precision>::right_is_node() const
        {
            return RIGHT_IS_NODE;
        }
    }
}
