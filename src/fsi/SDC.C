
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"
#include "GaussRadau.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "ClenshawCurtis.H"

namespace sdc
{
    SDC::SDC(
        std::shared_ptr<SDCSolver> solver,
        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
        scalar tol,
        int minSweeps,
        int maxSweeps
        )
        :
        solver( solver ),
        nbNodes( quadrature->get_num_nodes() ),
        N( solver->getDOF() ),
        k( quadrature->get_num_nodes() ),
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
        maxSweeps( maxSweeps ),
        quadrature( quadrature )
    {
        assert( solver );
        assert( dt > 0 );
        assert( tol > 0 );
        assert( tol < 1 );
        assert( maxSweeps >= minSweeps );
        assert( minSweeps > 0 );
        assert( N > 0 );
        assert( quadrature );
        assert( quadrature->right_is_node() );

        init();

        solver->setNumberOfImplicitStages( k - 1 );
    }

    SDC::SDC(
        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
        scalar tol
        )
        :
        solver( nullptr ),
        nbNodes( quadrature->get_num_nodes() ),
        N( 0 ),
        k( quadrature->get_num_nodes() ),
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
        maxSweeps( 0 ),
        quadrature( quadrature )
    {
        assert( tol > 0 );
        assert( tol < 1 );
        assert( quadrature );
        assert( quadrature->right_is_node() );

        init();
    }

    SDC::~SDC()
    {}

    void SDC::init()
    {
        smat = quadrature->get_s_mat();
        qmat = quadrature->get_q_mat();
        const std::vector<scalar> nodes = quadrature->get_nodes();

        this->nodes.resize( nodes.size() );

        for ( unsigned int i = 0; i < nodes.size(); i++ )
            this->nodes( i ) = nodes[i];

        dsdc.resize( this->nodes.rows() - 1 );

        for ( int i = 0; i < dsdc.rows(); i++ )
            dsdc( i ) = this->nodes( i + 1 ) - this->nodes( i );
    }

    void SDC::run()
    {
        assert( solver );

        scalar t = solver->getStartTime();

        while ( std::abs( t - solver->getEndTime() ) > 1.0e-13 && t < solver->getEndTime() )
        {
            scalar computedTimeStep = dt;

            solveTimeStep( t );

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
        assert( N > 0 );
        assert( k > 0 );

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

            Info << "\nTime = " << t << ", SDC sweep = 1, SDC substep = " << j + 1 << nl << endl;

            solver->implicitSolve( false, j, j, t, dt, qold, rhs, f, result );

            solStages.row( j + 1 ) = result;
            F.row( j + 1 ) = f;
        }

        // Compute successive corrections

        for ( int j = 0; j < maxSweeps - 1; j++ )
        {
            t = t0;

            fsi::matrix Sj = this->dt * (smat * F);

            // SDC sweep
            for ( int p = 0; p < k - 1; p++ )
            {
                scalar dt = dtsdc( p );
                t += dt;

                Info << "\nTime = " << t << ", SDC sweep = " << j + 2 << ", SDC substep = " << p + 1 << nl << endl;

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

            scalarList squaredNormResidual( Pstream::nProcs(), scalar( 0 ) );
            labelList dof( Pstream::nProcs(), label( 0 ) );
            squaredNormResidual[Pstream::myProcNo()] = residual.squaredNorm();
            dof[Pstream::myProcNo()] = residual.rows() * residual.cols();
            reduce( squaredNormResidual, sumOp<scalarList>() );
            reduce( dof, sumOp<labelList>() );
            scalar error = std::sqrt( sum( squaredNormResidual ) / sum( dof ) );
            bool convergence = error < tol && j >= minSweeps - 2;

            std::deque<int> dofVariables;
            std::deque<bool> enabledVariables;
            std::deque<std::string> namesVariables;
            solver->getVariablesInfo( dofVariables, enabledVariables, namesVariables );

            assert( enabledVariables.size() == dofVariables.size() );
            assert( enabledVariables.size() == namesVariables.size() );

            std::deque<bool> convergenceVariables;

            for ( unsigned int i = 0; i < enabledVariables.size(); i++ )
                convergenceVariables.push_back( true );

            bool solverConverged = solver->isConverged();

            if ( dofVariables.size() == 1 )
            {
                Info << "SDC residual = " << error;
                Info << ", tol = " << tol;
                Info << ", time = " << t;
                Info << ", sweep = " << j + 2;
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

                for ( unsigned int substep = 0; substep < residual.rows(); substep++ )
                {
                    int index = 0;

                    for ( unsigned int i = 0; i < dofVariables.size(); i++ )
                    {
                        assert( dofVariables.at( i ) > 0 );

                        scalarList squaredNormResidual( Pstream::nProcs(), scalar( 0 ) );
                        labelList dof( Pstream::nProcs(), label( 0 ) );

                        for ( int j = 0; j < dofVariables.at( i ); j++ )
                        {
                            squaredNormResidual[Pstream::myProcNo()] += residual( substep, index ) * residual( substep, index );
                            dof[Pstream::myProcNo()] += 1;
                            index++;
                        }

                        reduce( squaredNormResidual, sumOp<scalarList>() );
                        reduce( dof, sumOp<labelList>() );
                        scalar error = std::sqrt( sum( squaredNormResidual ) / sum( dof ) );

                        bool convergence = convergenceVariables.at( i );

                        if ( error > tol || j < minSweeps - 2 )
                            convergence = false;

                        if ( enabledVariables.at( i ) )
                        {
                            Info << "SDC " << namesVariables.at( i ).c_str();
                            Info << " substep = " << substep + 1;
                            Info << ", residual = " << error;
                            Info << ", time = " << t;
                            Info << ", sweep = " << j + 2;
                            Info << ", convergence = ";

                            if ( error < tol && j >= minSweeps - 2 )
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

        scalarList squaredNormResidual( Pstream::nProcs(), scalar( 0 ) );
        labelList dof( Pstream::nProcs(), label( 0 ) );
        squaredNormResidual[Pstream::myProcNo()] = residual.squaredNorm();
        reduce( squaredNormResidual, sumOp<scalarList>() );
        dof[Pstream::myProcNo()] = residual.rows() * residual.cols();
        scalar error = std::sqrt( sum( squaredNormResidual ) / sum( dof ) );
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
