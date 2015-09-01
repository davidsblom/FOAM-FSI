
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"
#include "QuadratureRules.H"

namespace sdc
{
    SDC::SDC(
        std::string rule,
        int nbNodes,
        scalar tol
        )
        :
        solver( false ),
        adaptiveTimeStepper( false ),
        nbNodes( nbNodes ),
        N( 0 ),
        k( 0 ),
        dt( -1 ),
        tol( tol ),
        nodes(),
        smat(),
        qmat(),
        nodesEmbedded(),
        smatEmbedded(),
        qmatEmbedded(),
        dsdc(),
        corrector( false ),
        stageIndex( 0 ),
        F(),
        Fold(),
        Sj(),
        solStages(),
        convergence( false ),
        timeIndex( 0 ),
        minSweeps( 0 ),
        maxSweeps( 0 )
    {
        assert( tol > 0 );
        assert( tol < 1 );
        assert( rule == "gauss-radau" || rule == "gauss-lobatto" || rule == "clenshaw-curtis" || rule == "uniform" || rule == "uniform-right-sided" );

        int refine = 1;
        quadrature::rules( rule, nbNodes, refine, nodes, smat, qmat );

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );
    }

    SDC::SDC(
        std::shared_ptr<SDCSolver> solver,
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper,
        std::string rule,
        int nbNodes,
        scalar tol,
        int minSweeps,
        int maxSweeps
        )
        :
        solver( solver ),
        adaptiveTimeStepper( adaptiveTimeStepper ),
        nbNodes( nbNodes ),
        N( solver->getDOF() ),
        k( 0 ),
        dt( solver->getTimeStep() ),
        tol( tol ),
        nodes(),
        smat(),
        qmat(),
        nodesEmbedded(),
        smatEmbedded(),
        qmatEmbedded(),
        dsdc(),
        corrector( false ),
        stageIndex( 0 ),
        F(),
        Fold(),
        Sj(),
        solStages(),
        convergence( false ),
        timeIndex( 0 ),
        minSweeps( minSweeps ),
        maxSweeps( maxSweeps )
    {
        assert( adaptiveTimeStepper );
        assert( solver );
        assert( nbNodes > 1 );
        assert( nbNodes < 15 );
        assert( dt > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( rule == "gauss-radau" || rule == "gauss-lobatto" || rule == "clenshaw-curtis" || rule == "uniform" || rule == "uniform-right-sided" );
        assert( maxSweeps > 0 );
        assert( maxSweeps >= minSweeps );
        assert( minSweeps > 0 );
        assert( N > 0 );

        int refine = 1;
        quadrature::rules( rule, nbNodes, refine, nodes, smat, qmat );

        k = nodes.rows();

        dsdc.resize( nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = nodes( i + 1 ) - nodes( i );

        solver->setNumberOfImplicitStages( k - 1 );

        assert( std::abs( nodes.tail( 1 ) ( 0 ) - 1 ) < 1.0e-13 );

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

            if ( rule == "uniform-right-sided" )
                orderEmbedded = nodesEmbedded.rows() - 1;

            adaptiveTimeStepper->setOrderEmbeddedMethod( orderEmbedded );
            adaptiveTimeStepper->setEndTime( solver->getEndTime() );

            assert( std::abs( nodesEmbedded.tail( 1 ) ( 0 ) - 1 ) < 1.0e-13 );
        }
    }

    SDC::~SDC()
    {}

    void SDC::run()
    {
        assert( solver );

        scalar t = solver->getStartTime();

        while ( std::abs( t - solver->getEndTime() ) > 1.0e-13 && t < solver->getEndTime() )
        {
            scalar computedTimeStep = dt;

            solveTimeStep( t );

            if ( adaptiveTimeStepper->isAccepted() )
                t += computedTimeStep;
        }
    }

    void SDC::solveTimeStep( const scalar t0 )
    {
        assert( qmat.rows() == smat.rows() );
        assert( qmat.cols() == smat.cols() );
        assert( qmat.rows() == k - 1 );
        assert( qmat.cols() == k );
        assert( solver );

        if ( adaptiveTimeStepper->isPreviousStepAccepted() )
            solver->nextTimeStep();

        fsi::vector dtsdc = this->dt * dsdc;
        fsi::matrix solStages( k, N ), F( k, N ), Fembedded( nodesEmbedded.rows(), N ), residual;
        fsi::matrix qj( 1, solStages.cols() ), qjEmbedded( 1, solStages.cols() );
        fsi::vector errorEstimate( N );

        fsi::vector sol( N ), f( N );
        solver->getSolution( sol, f );
        solStages.row( 0 ) = sol;

        scalar t = t0;

        solver->evaluateFunction( 0, sol, t, f );
        F.row( 0 ) = f;

        fsi::vector rhs( N ), result( N ), qold( N );
        rhs.setZero();

        solver->initTimeStep();

        for ( int j = 0; j < k - 1; j++ )
        {
            scalar dt = dtsdc( j );
            t += dt;

            qold = solStages.row( j );

            Info << "\nTime = " << t << ", SDC sweep = 0, SDC substep = " << j + 1 << nl << endl;

            solver->implicitSolve( false, j, j, t, dt, qold, rhs, f, result );

            solStages.row( j + 1 ) = result;
            F.row( j + 1 ) = f;
        }

        // Compute successive corrections

        for ( int j = 0; j < maxSweeps; j++ )
        {
            t = t0;

            fsi::matrix Sj = this->dt * (smat * F);

            // SDC sweep
            for ( int p = 0; p < k - 1; p++ )
            {
                scalar dt = dtsdc( p );
                t += dt;

                Info << "\nTime = " << t << ", SDC sweep = " << j + 1 << ", SDC substep = " << p + 1 << nl << endl;

                qold = solStages.row( p );

                // Form right hand side
                rhs.noalias() = -dt * F.row( p + 1 ) + Sj.row( p );

                solver->implicitSolve( true, p, p, t, dt, qold, rhs, f, result );

                solStages.row( p + 1 ) = result;
                F.row( p + 1 ) = f;
            }

            // Compute the SDC residual

            residual = dt * (qmat * F);

            for ( int i = 0; i < residual.rows(); i++ )
                residual.row( i ) += solStages.row( 0 ) - solStages.row( i + 1 );

            scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
            squaredNorm[Pstream::myProcNo()] = residual.squaredNorm();
            reduce( squaredNorm, sumOp<scalarList>() );
            scalar error = std::sqrt( sum( squaredNorm ) / N );
            error /= solver->getScalingFactor();
            bool convergence = error < tol && j >= minSweeps - 1;

            std::deque<int> dofVariables;
            std::deque<bool> enabledVariables;
            std::deque<std::string> namesVariables;
            solver->getVariablesInfo( dofVariables, enabledVariables, namesVariables );

            assert( enabledVariables.size() == dofVariables.size() );
            assert( enabledVariables.size() == namesVariables.size() );

            std::deque<bool> convergenceVariables;

            for ( unsigned int i = 0; i < enabledVariables.size(); i++ )
                convergenceVariables.push_back( not enabledVariables.at( i ) );

            bool solverConverged = solver->isConverged();

            if ( dofVariables.size() == 1 )
            {
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

                if ( convergence && solverConverged )
                    break;
            }

            if ( dofVariables.size() > 1 )
            {
                assert( std::accumulate( dofVariables.begin(), dofVariables.end(), 0 ) == N );

                bool convergence = solverConverged;

                for ( unsigned int substep = 0; substep < residual.rows(); substep++ )
                {
                    int index = 0;

                    for ( unsigned int i = 0; i < dofVariables.size(); i++ )
                    {
                        assert( dofVariables.at( i ) > 0 );

                        scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
                        labelList dofVariablesGlobal( Pstream::nProcs(), 0 );
                        dofVariablesGlobal[Pstream::myProcNo()] = dofVariables.at( i );
                        reduce( dofVariablesGlobal, sumOp<labelList>() );

                        for ( int j = 0; j < dofVariables.at( i ); j++ )
                        {
                            squaredNorm[Pstream::myProcNo()] += residual( substep, index ) * residual( substep, index );
                            index++;
                        }

                        reduce( squaredNorm, sumOp<scalarList>() );
                        scalar error = std::sqrt( sum( squaredNorm ) / sum( dofVariablesGlobal ) );

                        if ( error > tol || j < minSweeps - 1 )
                            convergence = false;

                        if ( enabledVariables.at( i ) )
                        {
                            Info << "SDC " << namesVariables.at( i ).c_str();
                            Info << " residual = " << error;
                            Info << " substep = " << substep + 1;
                            Info << ", time = " << t;
                            Info << ", sweep = " << j + 1;
                            Info << ", convergence = ";

                            if ( error < tol && j >= minSweeps - 1 )
                                Info << "true";
                            else
                                Info << "false";

                            Info << endl;
                        }

                        if ( enabledVariables.at( i ) )
                            convergenceVariables.at( i ) = convergence;
                    }

                    assert( index == N );
                }

                convergence = solverConverged;

                for ( unsigned int i = 0; i < convergenceVariables.size(); i++ )
                    if ( not convergenceVariables.at( i ) )
                        convergence = false;

                if ( convergence )
                    break;
            }
        }

        if ( adaptiveTimeStepper->isEnabled() )
        {
            scalar newTimeStep = 0;

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
        const fsi::matrix & qmat,
        const fsi::matrix & F,
        const scalar dt,
        fsi::matrix & qj
        )
    {
        // fsi::matrix Qj = dt * (qmat * F);

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

    void SDC::getSourceTerm(
        const bool corrector,
        const int k,
        const scalar deltaT,
        fsi::vector & rhs,
        fsi::vector & qold
        )
    {
        assert( k <= this->k - 1 );
        assert( solStages.rows() > 0 );
        assert( solStages.cols() > 0 );

        qold = solStages.row( k );

        // Compute the time step from the stage deltaT
        if ( dt < 0 )
        {
            // first time step, first prediction step
            dt = deltaT / dsdc( 0 );
        }

        assert( dt > 0 );

        if ( not corrector )
            rhs.setZero();

        if ( corrector )
        {
            if ( (k == 0 && stageIndex != 0) || (k == 0 && nbNodes == 2) )
            {
                Sj = dt * (smat * F);
                Fold = F;
            }

            rhs.noalias() = -dt * dsdc( k ) * Fold.row( k + 1 ) + Sj.row( k );
        }

        this->stageIndex = k;
        this->corrector = corrector;

        assert( rhs.rows() == qold.rows() );
        assert( rhs.rows() > 0 );
    }

    void SDC::setFunction(
        const int k,
        const fsi::vector & f,
        const fsi::vector & result
        )
    {
        assert( f.rows() == result.rows() );
        assert( k <= this->k - 1 );

        if ( F.cols() == 0 )
        {
            F.resize( this->k, f.rows() );
            F.setZero();
        }

        if ( solStages.cols() == 0 )
        {
            solStages.resize( this->k, f.rows() );
            solStages.setZero();
        }

        F.row( k + 1 ) = f;
        solStages.row( k + 1 ) = result;
    }

    void SDC::setOldSolution(
        int timeIndex,
        const fsi::vector & result
        )
    {
        assert( timeIndex >= this->timeIndex );

        if ( solStages.cols() == 0 )
        {
            solStages.resize( this->k, result.rows() );
            solStages.setZero();
            solStages.row( 0 ) = result;
        }
        else
        if ( timeIndex > this->timeIndex )
            solStages.row( 0 ) = solStages.bottomRows( 1 );

        this->timeIndex = timeIndex;
    }

    void SDC::outputResidual( std::string name )
    {
        fsi::matrix Qj = dt * (qmat * F);
        fsi::matrix residual = solStages.row( 0 ) + Qj.row( k - 2 ) - solStages.row( k - 1 );

        scalarList squaredNorm( Pstream::nProcs(), scalar( 0 ) );
        squaredNorm[Pstream::myProcNo()] = residual.squaredNorm();
        reduce( squaredNorm, sumOp<scalarList>() );
        scalar error = std::sqrt( sum( squaredNorm ) / F.cols() );
        convergence = error < tol;

        Info << "SDC " << name.c_str();
        Info << " residual = " << error;
        Info << ", convergence = ";

        if ( convergence )
            Info << "true";
        else
            Info << "false";

        Info << endl;
    }

    bool SDC::isConverged()
    {
        return convergence;
    }

    int SDC::getNbImplicitStages()
    {
        return this->k - 1;
    }
}
