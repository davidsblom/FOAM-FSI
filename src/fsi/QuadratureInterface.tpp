
#include "QuadratureInterface.H"

namespace fsi
{
    namespace quadrature
    {
        template<typename precision>
        IQuadrature<precision>::IQuadrature( const size_t num_nodes )
            : num_nodes( num_nodes )
        {
            assert( this->num_nodes > 0 );
        }

        template<typename precision>
        IQuadrature<precision>::IQuadrature()
            : num_nodes( 0 )
        {}

        template<typename precision>
        const Matrix<precision> & IQuadrature<precision>::get_q_mat() const
        {
            return this->q_mat;
        }

        template<typename precision>
        const Matrix<precision> & IQuadrature<precision>::get_s_mat() const
        {
            return this->s_mat;
        }

        template<typename precision>
        const Matrix<precision> & IQuadrature<precision>::get_b_mat() const
        {
            return this->b_mat;
        }

        template<typename precision>
        const vector<precision> & IQuadrature<precision>::get_q_vec() const
        {
            return this->q_vec;
        }

        template<typename precision>
        const vector<precision> & IQuadrature<precision>::get_nodes() const
        {
            return this->nodes;
        }

        template<typename precision>
        size_t IQuadrature<precision>::get_num_nodes() const
        {
            return this->num_nodes;
        }

        template<typename precision>
        bool IQuadrature<precision>::left_is_node() const
        {
            assert( false );
            return LEFT_IS_NODE;
        }

        template<typename precision>
        bool IQuadrature<precision>::right_is_node() const
        {
            assert( false );
            return RIGHT_IS_NODE;
        }

        template<typename precision>
        void IQuadrature<precision>::compute_nodes()
        {
            assert( false );
        }

        template<typename precision>
        precision IQuadrature<precision>::expected_error() const
        {
            using cvec = Eigen::Array<precision, Eigen::Dynamic, 1>;
            const cvec row_sums = this->get_q_mat().rowwise().sum();
            Eigen::Map<const cvec> nodes( this->get_nodes().data(), this->get_nodes().size() );
            return (row_sums - nodes).maxCoeff();
        }

        /**
         * @internals
         * Computing weights means computing \\( Q \\) and \\( S \\) matrices as well as the \\( q \\)
         * vector.
         * The \\( B \\) matrix is constructed from \\( q \\) vector.
         *
         * @note As long as the weight computation of the implemented quadrature is based on polynomial
         *   interpolation, this function must not be overwritten.
         * @endinternals
         */
        template<typename precision>
        void IQuadrature<precision>::compute_weights()
        {
            this->q_mat = compute_q_matrix( this->nodes );
            this->s_mat = compute_s_matrix( this->q_mat );
            this->q_vec = compute_q_vec( this->nodes );
            this->b_mat = Matrix<precision>::Zero( 1, this->num_nodes );

            for ( size_t i = 0; i < this->num_nodes; i++ )
            {
                this->b_mat( 0, i ) = this->q_vec[i];
            }
        }
    } // ::pfasst::quadrature
}  // ::pfasst
