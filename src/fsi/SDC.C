
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"
#include "QuadratureRules.H"

namespace sdc
{
    SDC::SDC(
        std::shared_ptr<SDCSolver> solver,
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
        std::string rule,
        int nbNodes,
        double tol
        )
        :
        solver( solver ),
        adaptiveTimeStepper( adaptiveTimeStepper ),
        nodes(),
        smat(),
        qmat(),
        dsdc(),
        dt( solver->getTimeStep() ),
        N( solver->getDOF() ),
        k( 0 ),
        tol( tol )
    {
        assert( adaptiveTimeStepper );
        assert( solver );
        assert( nbNodes > 1 );
        assert( nbNodes < 15 );
        assert( dt > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( rule == "gauss-radau" || rule == "gauss-lobatto" || rule == "clenshaw-curtis" || rule == "uniform" );

        int refine = 1;
        quadrature::rules( rule, nbNodes, refine, nodes, smat, qmat );

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );

        solver->setNumberOfStages( k );

        if ( adaptiveTimeStepper->isEnabled() )
        {
            refine = 2;
            quadrature::rules( rule, nbNodes, refine, nodesEmbedded, smatEmbedded, qmatEmbedded );

            int orderEmbedded = 0;

            if ( rule == "gauss-radau" )
                orderEmbedded = (nodesEmbedded.rows() - 1) * 2 - 1;

            if ( rule == "gauss-lobatto" )
                orderEmbedded = nodesEmbedded.rows() * 2 - 3;

            if ( rule == "clenshaw-curtis" )
                orderEmbedded = nodesEmbedded.rows() - 1;

            if ( rule == "uniform" )
                orderEmbedded = nodesEmbedded.rows();

            adaptiveTimeStepper->setOrderEmbeddedMethod( orderEmbedded );
            adaptiveTimeStepper->setEndTime( solver->getEndTime() );
        }
    }

    SDC::~SDC()
    {}

    void SDC::run()
    {
        double t = 0;

        while ( t < solver->getEndTime() )
        {
            double computedTimeStep = dt;

            solveTimeStep( t );

            if ( adaptiveTimeStepper->isAccepted() )
                t += computedTimeStep;
        }
    }

    void SDC::solveTimeStep( const double t0 )
    {
        if ( adaptiveTimeStepper->isPreviousStepAccepted() )
            solver->nextTimeStep();

        Eigen::VectorXd dtsdc = this->dt * dsdc;
        Eigen::MatrixXd solStages( k, N ), F( k, N );

        Eigen::VectorXd sol( N ), f( N );
        solver->getSolution( sol );
        solStages.row( 0 ) = sol;

        double t = t0;

        solver->evaluateFunction( 0, sol, t, f );
        F.row( 0 ) = f;

        Eigen::VectorXd rhs( N ), result( N ), qold( N );
        rhs.setZero();

        solver->initTimeStep();

        for ( int j = 0; j < k - 1; j++ )
        {
            double dt = dtsdc( j );
            t += dt;

            qold = solStages.row( j );

            Info << "\nTime = " << t << ", SDC sweep = 0, SDC substep = " << j + 1 << nl << endl;

            solver->implicitSolve( false, j, t, dt, qold, rhs, f, result );

            solStages.row( j + 1 ) = result;
            F.row( j + 1 ) = f;
        }

        // Compute successive corrections

        bool convergence = false;

        for ( int j = 0; j < 10 * k; j++ )
        {
            t = t0;

            Eigen::MatrixXd Sj = this->dt * (smat * F);

            // SDC sweep
            for ( int p = 0; p < k - 1; p++ )
            {
                double dt = dtsdc( p );
                t += dt;

                Info << "\nTime = " << t << ", SDC sweep = " << j + 1 << ", SDC substep = " << p + 1 << nl << endl;

                qold = solStages.row( p );

                // Form right hand side
                rhs = -dt * F.row( p + 1 ) + Sj.row( p );

                solver->implicitSolve( true, p, t, dt, qold, rhs, f, result );

                solStages.row( p + 1 ) = result;
                F.row( p + 1 ) = f;
            }

            // Compute the SDC residual

            // Eigen::MatrixXd Qj = dt * (qmat * F);
            // Eigen::MatrixXd residual = solStages.row( 0 ) + Qj.row( k - 2 ) - solStages.row( k - 1 );

            // Only compute row k-2 of matrix Qj for efficiency
            Eigen::MatrixXd qj( 1, solStages.cols() );

            int ii = k - 2, jj, kk;

            for ( jj = 0; jj < F.cols(); ++jj )
            {
                qj( 0, jj ) = 0;

                for ( kk = 0; kk < F.rows(); ++kk )
                    qj( 0, jj ) += qmat( ii, kk ) * F( kk, jj );

                qj( 0, jj ) *= dt;
            }

            Eigen::MatrixXd residual = solStages.row( 0 ) + qj.row( 0 ) - solStages.row( k - 1 );

            scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
            squaredNorm[Pstream::myProcNo()] = residual.squaredNorm();
            reduce( squaredNorm, sumOp<scalarList>() );
            double error = std::sqrt( sum( squaredNorm ) / N );
            error /= solver->getScalingFactor();
            convergence = error < tol && j >= k - 1;

            Info << "SDC residual = " << error;
            Info << ", tol = " << tol;
            Info << ", time = " << t;
            Info << ", sweep = " << j + 1;
            Info << ", convergence = ";

            if ( convergence )
                Info << "true";
            else
                Info << "false";

            Info << endl;

            if ( convergence )
                break;
        }

        if ( adaptiveTimeStepper->isEnabled() )
        {
            double newTimeStep = 0;
            Eigen::VectorXd errorEstimate( N );

            Eigen::MatrixXd Fembedded( nodesEmbedded.rows(), N );

            for ( int i = 0; i < Fembedded.rows(); i++ )
                Fembedded.row( i ) = F.row( i * 2 );

            Eigen::MatrixXd Qj = dt * (qmat * F);
            Eigen::MatrixXd QjEmbedded = dt * (qmatEmbedded * Fembedded);

            errorEstimate = Qj.row( k - 2 ) - QjEmbedded.row( nodesEmbedded.rows() - 2 );

            bool accepted = adaptiveTimeStepper->determineNewTimeStep( errorEstimate, result, dt, newTimeStep );

            dt = newTimeStep;

            if ( not accepted )
                solver->setSolution( solStages.row( 0 ), F.row( 0 ) );
        }

        if ( adaptiveTimeStepper->isAccepted() )
            solver->finalizeTimeStep();
    }
}
