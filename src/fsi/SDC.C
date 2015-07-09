
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
        N( solver->getDOF() ),
        k( 0 ),
        dt( solver->getTimeStep() ),
        tol( tol ),
        nodes(),
        smat(),
        qmat(),
        dsdc()
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
        double t = solver->getStartTime();

        while ( std::abs( t - solver->getEndTime() ) > 1.0e-13 && t < solver->getEndTime() )
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
        Eigen::MatrixXd solStages( k, N ), F( k, N ), Fembedded( nodesEmbedded.rows(), N ), residual;
        Eigen::MatrixXd qj( 1, solStages.cols() ), qjEmbedded( 1, solStages.cols() );
        Eigen::VectorXd errorEstimate( N );

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

        computeResidual( qmat, F, dt, qj );

        residual.noalias() = solStages.row( 0 ) + qj.row( 0 ) - solStages.row( k - 1 );

        scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
        squaredNorm[Pstream::myProcNo()] = residual.squaredNorm();
        reduce( squaredNorm, sumOp<scalarList>() );
        double error = std::sqrt( sum( squaredNorm ) / N );
        error /= solver->getScalingFactor();
        bool convergence = error < tol;

        Info << "SDC residual = " << error;
        Info << ", tol = " << tol;
        Info << ", time = " << t;
        Info << ", sweep = 0";
        Info << ", convergence = ";

        if ( convergence )
            Info << "true";
        else
            Info << "false";

        Info << endl;

        Eigen::VectorXd variables;
        solver->getVariableDOF( variables );
        assert( variables.sum() == N );

        if ( variables.rows() > 1 )
        {
            int index = 0;

            for ( int i = 0; i < variables.rows(); i++ )
            {
                scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );

                for ( int j = 0; j < variables( i ); j++ )
                {
                    squaredNorm[Pstream::myProcNo()] += residual( 0, index ) * residual( 0, index );
                    index++;
                }

                reduce( squaredNorm, sumOp<scalarList>() );
                double error = std::sqrt( sum( squaredNorm ) / variables( i ) );
                Info << "SDC residual variable " << i << " = " << error << endl;
            }

            assert( index == N );
        }

        // Compute successive corrections

        convergence = false;

        for ( int j = 0; j < 4 * k; j++ )
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
                rhs.noalias() = -dt * F.row( p + 1 ) + Sj.row( p );

                solver->implicitSolve( true, p, t, dt, qold, rhs, f, result );

                solStages.row( p + 1 ) = result;
                F.row( p + 1 ) = f;
            }

            // Compute the SDC residual

            // Eigen::MatrixXd Qj = dt * (qmat * F);
            // Eigen::MatrixXd residual = solStages.row( 0 ) + Qj.row( k - 2 ) - solStages.row( k - 1 );

            // Only compute row k-2 of matrix Qj for efficiency
            computeResidual( qmat, F, dt, qj );

            residual.noalias() = solStages.row( 0 ) + qj.row( 0 ) - solStages.row( k - 1 );

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

            Eigen::VectorXd variables;
            solver->getVariableDOF( variables );

            if ( variables.rows() > 1 )
            {
                int index = 0;

                for ( int i = 0; i < variables.rows(); i++ )
                {
                    scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );

                    for ( int j = 0; j < variables( i ); j++ )
                    {
                        squaredNorm[Pstream::myProcNo()] += residual( 0, index ) * residual( 0, index );
                        index++;
                    }

                    reduce( squaredNorm, sumOp<scalarList>() );
                    double error = std::sqrt( sum( squaredNorm ) / variables( i ) );
                    Info << "SDC residual variable " << i << " = " << error << endl;
                }

                assert( index == N );
            }

            if ( convergence )
                break;
        }

        if ( adaptiveTimeStepper->isEnabled() )
        {
            double newTimeStep = 0;

            for ( int i = 0; i < Fembedded.rows(); i++ )
                Fembedded.row( i ) = F.row( i * 2 );

            computeResidual( qmat, F, dt, qj );
            computeResidual( qmatEmbedded, Fembedded, dt, qjEmbedded );

            errorEstimate.noalias() = qj.row( 0 ) - qjEmbedded.row( 0 );

            bool accepted = adaptiveTimeStepper->determineNewTimeStep( errorEstimate, result, dt, newTimeStep );

            dt = newTimeStep;

            if ( not accepted )
                solver->setSolution( solStages.row( 0 ), F.row( 0 ) );
        }

        if ( adaptiveTimeStepper->isAccepted() )
            solver->finalizeTimeStep();
    }

    void SDC::computeResidual(
        const Eigen::MatrixXd & qmat,
        const Eigen::MatrixXd & F,
        const double dt,
        Eigen::MatrixXd & qj
        )
    {
        // Eigen::MatrixXd Qj = dt * (qmat * F);

        // Only compute row k-2 of matrix Qj for efficiency
        int k = F.rows();
        int ii = k - 2, jj, kk;

        for ( jj = 0; jj < F.cols(); ++jj )
        {
            qj( 0, jj ) = 0;

            for ( kk = 0; kk < F.rows(); ++kk )
                qj( 0, jj ) += qmat( ii, kk ) * F( kk, jj );

            qj( 0, jj ) *= dt;
        }
    }
}
