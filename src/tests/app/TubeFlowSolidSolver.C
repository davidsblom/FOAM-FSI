
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TubeFlowSolidSolver.H"

namespace tubeflow
{
    TubeFlowSolidSolver::TubeFlowSolidSolver(
        double a0,
        double cmk,
        double p0,
        double rho,
        double L,
        int N
        )
        :
        BaseMultiLevelSolver( N, 1, a0 ),
        a0( a0 ),
        cmk( cmk ),
        p0( p0 ),
        rho( rho ),
        dx( L / N ),
        grid()
    {
        assert( cmk > 0 );
        assert( N > 0 );
        assert( rho > 0 );
        assert( a0 > 0 );
        assert( dx > 0 );
        data.fill( a0 );
    }

    TubeFlowSolidSolver::~TubeFlowSolidSolver()
    {}

    void TubeFlowSolidSolver::calcGrid()
    {
        if ( grid.rows() != N )
        {
            grid.resize( N, 1 );
            grid.setZero();

            for ( int i = 0; i < N; i++ )
                grid( i, 0 ) = dx * i + 0.5 * dx;
        }
    }

    void TubeFlowSolidSolver::finalizeTimeStep()
    {
        assert( init );

        init = false;
    }

    void TubeFlowSolidSolver::getWritePositions( matrix & writePositions )
    {
        calcGrid();

        writePositions = grid;
    }

    void TubeFlowSolidSolver::getReadPositions( matrix & readPositions )
    {
        calcGrid();

        readPositions = grid;
    }

    void TubeFlowSolidSolver::initTimeStep()
    {
        assert( !init );

        init = true;
    }

    bool TubeFlowSolidSolver::isRunning()
    {
        return true;
    }

    void TubeFlowSolidSolver::resetSolution()
    {}

    void TubeFlowSolidSolver::solve(
        const matrix & input,
        matrix & output
        )
    {
        assert( input.rows() == N );
        assert( input.cols() == 1 );

        // Map the matrices to fsi::vectors
        fsi::vector p = Eigen::Map< const fsi::vector> ( input.data(), input.rows() );
        fsi::vector a;

        // Call solve with the correct arguments
        solve( p, a );

        output = a;

        assert( output.rows() == N );
        assert( output.cols() == 1 );
    }

    void TubeFlowSolidSolver::solve(
        const fsi::vector & p,
        fsi::vector & a
        )
    {
        std::cout << "Solve solid domain with size " << N << std::endl;

        // Verify input variables
        assert( init );
        assert( p.rows() == N );
        assert( data.rows() == N );
        a.resize( N );
        a.setZero();

        for ( int i = 0; i < a.rows(); i++ )
        {
            a( i ) = p0 / (2 * rho) - cmk * cmk;
            a( i ) /= p( i ) / (2 * rho) - cmk * cmk;
            a( i ) = a0 * a( i ) * a( i );
        }

        data.col( 0 ) = a;
    }
}
