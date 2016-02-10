
#include "GaussRadau.H"

#include <stdexcept>
#include <vector>
using namespace std;

#include "Polynomial.H"

namespace fsi
{
    namespace quadrature
    {
        template<typename precision>
        GaussRadau<precision>::GaussRadau( const size_t num_nodes )
            : IQuadrature<precision>( num_nodes )
        {
            if ( this->num_nodes < 2 )
            {
                throw invalid_argument( "Gauss-Radau quadrature requires at least two quadrature nodes." );
            }

            this->compute_nodes();
            this->compute_weights();
        }

        template<typename precision>
        bool GaussRadau<precision>::left_is_node() const
        {
            return LEFT_IS_NODE;
        }

        template<typename precision>
        bool GaussRadau<precision>::right_is_node() const
        {
            return RIGHT_IS_NODE;
        }

        template<typename precision>
        void GaussRadau<precision>::compute_nodes()
        {
            this->nodes = std::vector<precision>( this->num_nodes, precision( 0.0 ) );
            auto l = Polynomial<precision>::legendre( this->num_nodes );
            auto lm1 = Polynomial<precision>::legendre( this->num_nodes - 1 );

            for ( size_t i = 0; i < this->num_nodes; i++ )
            {
                l[i] += lm1[i];
            }

            auto roots = l.roots();

            for ( size_t j = 1; j < this->num_nodes; j++ )
            {
                this->nodes[j - 1] = 0.5 * (1.0 - roots[this->num_nodes - j]);
            }

            this->nodes.back() = 1.0;
        }
    }
}
