
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

 #include "DataStorage.H"

 namespace sdc
 {
     DataStorage::DataStorage(
         std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
         int N )
     :
     quadrature( quadrature ),
     F( quadrature->get_num_nodes(), N ),
     solStages( quadrature->get_num_nodes(), N )
     {
         assert( quadrature );
     }

     DataStorage::~DataStorage()
     {

     }

     const fsi::matrix DataStorage::integrate( const std::vector<scalar> & nodes, scalar dt ) const
     {
         std::vector<scalar> quadratureNodes = quadrature->get_nodes();

         if ( not quadrature->left_is_node() )
         quadratureNodes.erase( quadratureNodes.begin() );

         fsi::quadrature::Matrix<scalar> q_matrix = fsi::quadrature::compute_q_matrix( quadratureNodes, nodes );

         if ( not quadrature->left_is_node() )
         {
             fsi::quadrature::Matrix<scalar> q_mat2( q_matrix.rows(), q_matrix.cols() + 1 );
             q_mat2.setZero();
             q_mat2.rightCols( q_matrix.cols() ) = q_matrix;
             q_matrix = q_mat2;
         }

         fsi::quadrature::Matrix<scalar> data = dt * q_matrix * getFunctions();

         for ( int i = 0; i < data.rows(); i++ )
             data.row( i ) += getSolutions().row( 0 );

         return data;
     }

     const fsi::matrix & DataStorage::getFunctions() const
     {
         assert( F.rows() > 0 );
         assert( F.cols() > 0 );
         return F;
     }

     const fsi::matrix & DataStorage::getSolutions() const
     {
         assert( solStages.rows() > 0 );
         return solStages;
     }

     const fsi::vector DataStorage::getFunction( int substep ) const
     {
         assert( substep <= F.rows() );
         assert( F.cols() > 0 );
         assert( F.rows() > 0 );
         return F.row( substep );
     }

     const fsi::vector DataStorage::getSolution( int substep ) const
     {
         assert( substep <= solStages.rows() );
         assert( solStages.cols() > 0 );
         assert( solStages.rows() > 0 );
         return solStages.row( substep );
     }

     const fsi::vector DataStorage::getLastSolution() const
     {
         return getSolution( quadrature->get_num_nodes() - 1 );
     }

     void DataStorage::initialize( int k, int N )
     {
         assert( N > 0 );
         assert( k >= 2 );
         F.resize( k, N );
         solStages.resize( k, N );
         F.setZero();
         solStages.setZero();
     }

     void DataStorage::storeFunction( const fsi::vector & f, int substep )
     {
         assert( f.rows() == F.cols() );
         assert( substep <= F.rows() );
         assert( F.cols() > 0 );
         assert( F.rows() > 0 );
         F.row( substep ) = f;
     }

     void DataStorage::storeSolution( const fsi::vector & sol, int substep )
     {
         assert( sol.rows() == solStages.cols() );
         assert( substep <= solStages.rows() );
         assert( solStages.cols() > 0 );
         assert( solStages.rows() > 0 );
         solStages.row( substep ) = sol;
     }
 }
