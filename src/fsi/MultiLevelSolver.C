
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MultiLevelSolver.H"

namespace fsi
{
    MultiLevelSolver::MultiLevelSolver(
        shared_ptr<BaseMultiLevelSolver> solver,
        shared_ptr<BaseMultiLevelSolver> couplingGridSolver,
        int participantId,
        int level
        )
        :
        solver( solver ),
        couplingGridSolver( couplingGridSolver ),
        rbfInterpToCouplingMesh( shared_ptr<RBFCoarsening> ( new RBFCoarsening() ) ),
        rbfInterpToMesh( shared_ptr<RBFCoarsening> ( new RBFCoarsening() ) ),
        participantId( participantId ),
        level( level ),
        couplingGridSize( 0 )
    {
        assert( solver );
        assert( couplingGridSolver );
        assert( participantId == 0 || participantId == 1 );
        assert( level >= 0 );
        assert( rbfInterpToCouplingMesh );
        assert( rbfInterpToMesh );

        matrix couplingGridPositions;

        if ( participantId == 0 )
            couplingGridSolver->getWritePositions( couplingGridPositions );

        if ( participantId == 1 )
            couplingGridSolver->getReadPositions( couplingGridPositions );

        couplingGridSize = couplingGridPositions.rows();

        assert( couplingGridSize > 0 );
    }

    MultiLevelSolver::MultiLevelSolver(
        shared_ptr<BaseMultiLevelSolver> solver,
        shared_ptr<BaseMultiLevelSolver> couplingGridSolver,
        shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh,
        shared_ptr<RBFCoarsening> rbfInterpToMesh,
        int participantId,
        int level
        )
        :
        solver( solver ),
        couplingGridSolver( couplingGridSolver ),
        rbfInterpToCouplingMesh( rbfInterpToCouplingMesh ),
        rbfInterpToMesh( rbfInterpToMesh ),
        participantId( participantId ),
        level( level ),
        couplingGridSize( 0 )
    {
        assert( solver );
        assert( couplingGridSolver );
        assert( participantId == 0 || participantId == 1 );
        assert( level >= 0 );
        assert( rbfInterpToCouplingMesh );
        assert( rbfInterpToMesh );

        matrix couplingGridPositions;

        if ( participantId == 0 )
            couplingGridSolver->getWritePositions( couplingGridPositions );

        if ( participantId == 1 )
            couplingGridSolver->getReadPositions( couplingGridPositions );

        couplingGridSize = couplingGridPositions.rows();

        assert( couplingGridSize > 0 );
    }

    void MultiLevelSolver::interpToCouplingMesh(
        matrix & data,
        matrix & dataInterpolated
        )
    {
        assert( data.cols() > 0 );
        assert( data.rows() > 0 );

        rbfInterpToCouplingMesh->interpolate( data, dataInterpolated );

        assert( dataInterpolated.cols() > 0 );
        assert( dataInterpolated.rows() > 0 );
    }

    void MultiLevelSolver::interpToMesh(
        matrix & data,
        matrix & dataInterpolated
        )
    {
        assert( data.cols() > 0 );
        assert( data.rows() > 0 );

        rbfInterpToMesh->interpolate( data, dataInterpolated );

        assert( dataInterpolated.cols() > 0 );
        assert( dataInterpolated.rows() > 0 );
    }

    void MultiLevelSolver::solve(
        matrix & input,
        matrix & outputInterpolated
        )
    {
        // Initialize variables
        matrix inputInterpolated, output;

        // Interpolate data from coupling mesh to own mesh

        interpToMesh( input, inputInterpolated );

        // Call solver

        solver->solve( inputInterpolated, output );

        // Interpolate result from own mesh to coupling mesh

        interpToCouplingMesh( output, outputInterpolated );
    }

    void MultiLevelSolver::updateInterpolationMatrices()
    {
        matrix writePositions, readPositions, couplingGridPositions;

        solver->getWritePositions( writePositions );
        solver->getReadPositions( readPositions );

        if ( participantId == 0 )
            couplingGridSolver->getWritePositions( couplingGridPositions );

        if ( participantId == 1 )
            couplingGridSolver->getReadPositions( couplingGridPositions );

        rbfInterpToCouplingMesh->compute( writePositions, couplingGridPositions );

        if ( participantId == 0 )
            couplingGridSolver->getReadPositions( couplingGridPositions );

        if ( participantId == 1 )
            couplingGridSolver->getWritePositions( couplingGridPositions );

        rbfInterpToMesh->compute( couplingGridPositions, readPositions );
    }
}
