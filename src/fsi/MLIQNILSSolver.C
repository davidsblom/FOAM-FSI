
/*
 * Copyright [2016] <David Blom>
 */

#include <memory>
#include <deque>
#include "MLIQNILSSolver.H"

namespace fsi {
MLIQNILSSolver::MLIQNILSSolver(shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models,
    bool synchronization
    )
    :
    models(models),
    fineModel(models->back()),
    init(false),
    synchronization(synchronization) {
    assert(models);
    assert(models->size() >= 2);
}

void MLIQNILSSolver::finalizeTimeStep() {
    assert(init);

    Info << endl << "Synchronize solvers" << endl;

    int level = 1;

    for (auto && model : *models) {
        if (synchronization) {
            if (level < static_cast<int>(models->size())) {
                assert(model->fsi->x.rows() == fineModel->fsi->x.rows());

                matrix input = Eigen::Map<const matrix> (fineModel->fsi->x.head(model->fsi->solidSolver->couplingGridSize * model->fsi->solid->dim).data(), model->fsi->solidSolver->couplingGridSize, model->fsi->solid->dim);

                matrix output(fineModel->fsi->solid->data.rows(), fineModel->fsi->solid->data.cols());
                model->fsi->fluidSolver->solve(input, output);

                if (model->fsi->parallel)
                    input = Eigen::Map<const matrix> (fineModel->fsi->x.tail(model->fsi->fluidSolver->couplingGridSize * model->fsi->fluid->dim).data(), model->fsi->fluidSolver->couplingGridSize, model->fsi->fluid->dim);

                if (!model->fsi->parallel)
                    input = fineModel->fsi->fluid->data;

                output.resize(fineModel->fsi->fluid->data.rows(), fineModel->fsi->fluid->data.cols());

                bool interpolated = model->fsi->solid->interpolateVolField(fineModel->fsi->solid);

                if (!interpolated)
                    model->fsi->solidSolver->solve(input, output);

                model->fsi->x = fineModel->fsi->x;
            }
        }

        if (!synchronization) {
            model->solve();
        }

        model->finalizeTimeStep();

        level++;
    }

    init = false;
}

void MLIQNILSSolver::initTimeStep() {
    assert(!init);

    for (auto && model : *models)
        model->initTimeStep();

    init = true;
}

bool MLIQNILSSolver::isRunning() {
    return fineModel->fsi->isRunning();
}

void MLIQNILSSolver::run() {
    assert(!init);

    time = std::clock();

    while (isRunning())
        solveTimeStep();
}

void MLIQNILSSolver::solve() {
    assert(init);

    // Initialize variables

    fsi::vector x0 = fineModel->fsi->x;
    fsi::vector xk = x0;
    int level = 1;

    for (auto && model : *models) {
        // Initialize variables

        x0 = xk;

        // Optimize the model

        model->optimize(x0, xk);

        model->fsi->x = xk;

        // Copy Jacobian information from the current level to the next level

        if (level < static_cast<int>(models->size())) {
            // Copy residuals of current time step

            for (auto && residuals : model->postProcessing->residualsList)
                models->at(level)->postProcessing->residualsList.push_back(residuals);

            for (auto && sols : model->postProcessing->solsList)
                models->at(level)->postProcessing->solsList.push_back(sols);
        }

        level++;
    }
}

void MLIQNILSSolver::solveTimeStep() {
    assert(!init);

    initTimeStep();
    solve();
    finalizeTimeStep();
}
} // namespace fsi
