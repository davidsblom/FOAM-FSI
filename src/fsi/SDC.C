
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
    int nbNodes,
    double tol
    )
    :
    solver( solver ),
    nodes(),
    smat(),
    qmat(),
    dsdc(),
    dt( solver->getTimeStep() ),
    N( solver->getDOF() ),
    k( 0 ),
    tol( tol )
  {
    assert( solver );
    assert( nbNodes > 1 );
    assert( nbNodes < 15 );
    assert( dt > 0 );
    assert( tol > 0 );
    assert( tol < 1 );

    quadrature::rules( nbNodes, nodes, smat, qmat );

    k = nodes.rows();

    dsdc.resize( nodes.rows() - 1 );

    for ( int i = 0; i < dsdc.rows(); i++ )
      dsdc( i ) = nodes( i + 1 ) - nodes( i );

    solver->setNumberOfStages( k );
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
    solver->nextTimeStep();

    Eigen::VectorXd dtsdc = this->dt * dsdc;
    Eigen::MatrixXd solStages( k, N ), F( k, N );

    Eigen::VectorXd sol( N );
    solver->getSolution( sol );
    solStages.row( 0 ) = sol;

    double t = t0;

    for ( int j = 0; j < k - 1; j++ )
    {
      double dt = dtsdc( j );
      t += dt;

      Eigen::VectorXd f( N ), rhs( N ), result( N ), qold( N );
      f.setZero();
      rhs.setZero();
      result.setZero();
      qold = solStages.row( j );

      Info << "\nTime = " << t << ", SDC sweep = 0, SDC stage = " << j + 1 << nl << endl;

      solver->initTimeStep();
      solver->implicitSolve( j, t, dt, qold, rhs, f, result );
      solver->finalizeTimeStep();

      solStages.row( j + 1 ) = result;
      F.row( j + 1 ) = f;
    }

    // Compute successive corrections

    for ( int j = 0; j < 5 * k; j++ )
    {
      t = t0;
      Eigen::MatrixXd Fold = F;

      Eigen::VectorXd f( N ), rhs( N ), result( N ), qold( N ), q( N );
      f.setZero();
      rhs.setZero();
      result.setZero();

      q = solStages.row( 0 );
      solver->evaluateFunction( 0, q, t, f );
      F.row( 0 ) = f;

      Eigen::MatrixXd Sj = this->dt * (smat * F);

      // SDC sweep
      for ( int p = 0; p < k - 1; p++ )
      {
        double dt = dtsdc( p );
        t += dt;
        f.setZero();
        result.setZero();

        Info << "\nTime = " << t << ", SDC sweep = " << j + 1 << ", SDC stage = " << p + 1 << nl << endl;

        qold = solStages.row( p );

        // Form right hand side
        rhs = -dt * Fold.row( p + 1 ) + Sj.row( p );

        solver->initTimeStep();
        solver->implicitSolve( p, t, dt, qold, rhs, f, result );
        solver->finalizeTimeStep();

        solStages.row( p + 1 ) = result;
        F.row( p + 1 ) = f;
      }

      // Compute the SDC residual

      Eigen::MatrixXd Qj = dt * (qmat * F);
      Eigen::MatrixXd residual = solStages.row( 0 ) + Qj.row( k - 2 ) - solStages.row( k - 1 );

      double error = std::sqrt( residual.squaredNorm() / ( residual.rows() * residual.cols() ) );
      bool convergence = error < tol;

      labelList convergenceList( Pstream::nProcs(), 0 );
      convergenceList[Pstream::myProcNo()] = convergence;
      reduce( convergenceList, sumOp<labelList>() );
      convergence = min( convergenceList );

      Info << "SDC residual = " << error << ", tol = " << tol << ", convergence = ";

      if ( convergence )
        Info << "true";
      else
        Info << "false";

      Info << endl;

      if ( convergence )
        break;
    }
  }
}
