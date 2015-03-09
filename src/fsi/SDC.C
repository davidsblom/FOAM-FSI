
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"
#include "gauss_radau.H"

namespace sdc
{
  SDC::SDC(
    std::shared_ptr<SDCSolver> solver,
    int nbNodes
    )
    :
    solver( solver ),
    nbNodes( nbNodes ),
    nodes(),
    smat(),
    qmat(),
    dsdc(),
    dt( solver->getTimeStep() ),
    N( solver->getDOF() )
  {
    assert( solver );
    assert( nbNodes > 1 );
    assert( nbNodes < 15 );
    assert( dt > 0 );

    quadrature::rules( nbNodes, nodes, smat, qmat );

    dsdc.resize( nodes.rows() - 1 );

    for ( int i = 0; i < dsdc.rows(); i++ )
      dsdc( i ) = nodes( i + 1 ) - nodes( i );
  }

  SDC::~SDC()
  {}

  void SDC::run()
  {
    int i = 0;

    while ( solver->isRunning() )
    {
      solveTimeStep( dt * i );
      i++;
    }
  }

  void SDC::solveTimeStep( const double t0 )
  {
    Eigen::VectorXd dtsdc = this->dt * dsdc;
    Eigen::MatrixXd solStages( nbNodes, N ), F( nbNodes, N );

    Eigen::VectorXd sol( N );
    solver->getSolution( sol );
    solStages.row( 0 ) = sol;

    double t = t0;

    for ( int j = 0; j < nbNodes - 1; j++ )
    {
      double dt = dtsdc( j );
      t += dt;

      Eigen::VectorXd f( N ), rhs( N ), result( N ), qold( N );
      f.setZero();
      rhs.setZero();
      result.setZero();
      qold = solStages.row( j );

      solver->initTimeStep();
      solver->implicitSolve( t, dt, qold, rhs, f, result );
      solver->finalizeTimeStep();

      solStages.row( j + 1 ) = result;
      F.row( j + 1 ) = f;
    }

    // Compute successive corrections

    for ( int j = 0; j < nbNodes; j++ )
    {
      t = t0;
      Eigen::MatrixXd Fold = F;

      Eigen::VectorXd f( N ), rhs( N ), result( N ), qold( N ), q( N );
      f.setZero();
      rhs.setZero();
      result.setZero();

      q = solStages.row( 0 );
      solver->evaluateFunction( q, t, f );
      F.row( 0 ) = f;

      Eigen::MatrixXd Sj = this->dt * (smat * F);

      // SDC sweep
      for ( int p = 0; p < nbNodes - 1; p++ )
      {
        double dt = dtsdc( p );
        t += dt;
        f.setZero();
        result.setZero();

        qold = solStages.row( p );

        // Form right hand side
        rhs = -dt * Fold.row( p + 1 ) + Sj.row( p );

        // rhs /= dt;

        solver->initTimeStep();
        solver->implicitSolve( t, dt, qold, rhs, f, result );
        solver->finalizeTimeStep();

        solStages.row( p + 1 ) = result;
        F.row( p + 1 ) = f;
      }
    }

    // Compute the SDC residual

    Eigen::MatrixXd Qj = dt * (qmat * F);
    Eigen::MatrixXd residual = solStages.row( 0 ) + Qj.row( nbNodes - 2 ) - solStages.row( nbNodes - 1 );
  }
}
