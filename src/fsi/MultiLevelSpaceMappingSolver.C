
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MultiLevelSpaceMappingSolver.H"

using namespace fsi;

MultiLevelSpaceMappingSolver::MultiLevelSpaceMappingSolver(
    shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > solvers,
    shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models,
    bool synchronization
    )
    :
    solvers( solvers ),
    models( models ),
    fineModel( models->back() ),
    init( false ),
    synchronization( synchronization )
{
    assert( solvers->size() > 0 );
    assert( models->size() == solvers->size() + 1 );

    int level = 0;

    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = models->begin(); it != models->end(); ++it )
    {
        shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;

        assert( model->fsi->fluidSolver->level == model->fsi->solidSolver->level );
        assert( model->fsi->fluidSolver->level == level );

        level++;
    }

    level = 0;

    for ( std::deque<shared_ptr<SpaceMappingSolver> >::iterator it = solvers->begin(); it != solvers->end(); ++it )
    {
        shared_ptr<SpaceMappingSolver> solver = *it;

        assert( solver->fineModel->fsi->fluidSolver->level == level + 1 );
        assert( solver->coarseModel->fsi->fluidSolver->level == level );

        level++;
    }
}

void MultiLevelSpaceMappingSolver::initTimeStep()
{
    assert( !init );

    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = models->begin(); it != models->end(); ++it )
    {
        shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;
        model->initTimeStep();
    }

    for ( std::deque<shared_ptr<SpaceMappingSolver> >::iterator it = solvers->begin(); it != solvers->end(); ++it )
    {
        shared_ptr<SpaceMappingSolver> solver = *it;
        solver->init = true;
    }

    init = true;
}

void MultiLevelSpaceMappingSolver::finalizeTimeStep()
{
    assert( init );

    Info << endl << "Synchronize solvers" << endl;

    int level = 1;

    for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = models->begin(); it != models->end(); ++it )
    {
        shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;

        if ( level < static_cast<int>( models->size() ) )
        {
            assert( model->fsi->x.rows() == fineModel->fsi->x.rows() );

            if ( synchronization )
            {
                matrix input = Eigen::Map<const matrix> ( fineModel->fsi->x.head( model->fsi->solidSolver->couplingGridSize * model->fsi->solid->dim ).data(), model->fsi->solidSolver->couplingGridSize, model->fsi->solid->dim );

                matrix output( fineModel->fsi->solid->data.rows(), fineModel->fsi->solid->data.cols() );

                if ( std::abs( model->fsi->x.norm() - fineModel->fsi->x.norm() ) > 1.0e-14 )
                    model->fsi->fluidSolver->solve( input, output );

                if ( model->fsi->parallel )
                    input = Eigen::Map<const matrix> ( fineModel->fsi->x.tail( model->fsi->fluidSolver->couplingGridSize * model->fsi->fluid->dim ).data(), model->fsi->fluidSolver->couplingGridSize, model->fsi->fluid->dim );

                if ( !model->fsi->parallel )
                    input = fineModel->fsi->fluid->data;

                output.resize( fineModel->fsi->fluid->data.rows(), fineModel->fsi->fluid->data.cols() );

                bool interpolated = model->fsi->solid->interpolateVolField( fineModel->fsi->solid );

                if ( !interpolated )
                    if ( std::abs( model->fsi->x.norm() - fineModel->fsi->x.norm() ) > 1.0e-14 || !model->fsi->parallel )
                        model->fsi->solidSolver->solve( input, output );

                model->fsi->x = fineModel->fsi->x;
            }

            if ( !synchronization )
            {
                model->solve();
            }
        }

        model->finalizeTimeStep();

        level++;
    }

    for ( std::deque<shared_ptr<SpaceMappingSolver> >::iterator it = solvers->begin(); it != solvers->end(); ++it )
    {
        shared_ptr<SpaceMappingSolver> solver = *it;

        solver->spaceMapping->finalizeTimeStep();
    }

    init = false;
}

bool MultiLevelSpaceMappingSolver::isRunning()
{
    return models->back()->fsi->fluid->isRunning();
}

void MultiLevelSpaceMappingSolver::run()
{
    assert( !init );

    time = std::clock();

    while ( isRunning() )
        solveTimeStep();
}

void MultiLevelSpaceMappingSolver::solve()
{
    assert( init );

    solvers->back()->solve();
}

void MultiLevelSpaceMappingSolver::solveTimeStep()
{
    assert( !init );

    initTimeStep();
    solve();
    finalizeTimeStep();
}
