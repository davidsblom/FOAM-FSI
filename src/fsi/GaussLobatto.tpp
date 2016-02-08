
#include <stdexcept>
#include <vector>

#include "Polynomial.H"
#include "GaussLobatto.H"

namespace fsi
{
    namespace quadrature
    {
        using namespace std;

        template<typename precision>
        GaussLobatto<precision>::GaussLobatto( const size_t num_nodes )
            : IQuadrature<precision>(num_nodes)
        {
            if ( this->num_nodes < 2 )
            {
                throw invalid_argument( "Gauss-Lobatto quadrature requires at least two quadrature nodes." );
            }

            this->compute_nodes();
            this->compute_weights();
        }

        template<typename precision>
        bool GaussLobatto<precision>::left_is_node() const
        {
            return LEFT_IS_NODE;
        }

        template<typename precision>
        bool GaussLobatto<precision>::right_is_node() const
        {
            return RIGHT_IS_NODE;
        }

        template<typename precision>
        void GaussLobatto<precision>::compute_nodes()
        {
            this->nodes = std::vector<precision>( this->num_nodes, precision( 0.0 ) );
            auto roots = Polynomial<precision>::legendre( this->num_nodes - 1 ).differentiate().roots();

            for ( size_t j = 0; j < this->num_nodes - 2; j++ )
            {
                this->nodes[j + 1] = 0.5 * (1.0 + roots[j]);
            }

            this->nodes.front() = 0.0;
            this->nodes.back() = 1.0;
        }
    }
}
