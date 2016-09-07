
/*
 * Copyright [2016] <David Blom>
 */

#include <numeric>
#include <memory>
#include <string>
#include <vector>
#include <deque>
#include "PstreamReduceOps.H"
#include "SDC.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "ClenshawCurtis.H"

namespace sdc {
SDC::SDC(std::shared_ptr<SDCSolver> solver,
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
    std::shared_ptr<sdc::DataStorage> data,
    scalar tol,
    int minSweeps,
    int maxSweeps
    )
    :
    solver(solver),
    nbNodes(quadrature->get_num_nodes()),
    N(solver->getDOF()),
    k(quadrature->get_num_nodes()),
    dt(solver->getTimeStep()),
    tol(tol),
    nodes(),
    smat(),
    qmat(),
    nodesEmbedded(),
    smatEmbedded(),
    qmatEmbedded(),
    dsdc(),
    corrector(false),
    stageIndex(0),
    Sj(),
    convergence(false),
    timeIndex(0),
    minSweeps(minSweeps),
    maxSweeps(maxSweeps),
    quadrature(quadrature),
    data(data) {
    assert(solver);
    assert(dt > 0);
    assert(tol > 0);
    assert(tol < 1);
    assert(maxSweeps >= minSweeps);
    assert(minSweeps > 0);
    assert(N > 0);
    assert(quadrature);
    assert(quadrature->right_is_node());
    assert(data);

    init();

    solver->setNumberOfImplicitStages(k - 1);
}

SDC::SDC(std::shared_ptr<SDCSolver> solver,
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
    scalar tol,
    int minSweeps,
    int maxSweeps
    )
    :
    solver(solver),
    nbNodes(quadrature->get_num_nodes()),
    N(solver->getDOF()),
    k(quadrature->get_num_nodes()),
    dt(solver->getTimeStep()),
    tol(tol),
    nodes(),
    smat(),
    qmat(),
    nodesEmbedded(),
    smatEmbedded(),
    qmatEmbedded(),
    dsdc(),
    corrector(false),
    stageIndex(0),
    Sj(),
    convergence(false),
    timeIndex(0),
    minSweeps(minSweeps),
    maxSweeps(maxSweeps),
    quadrature(quadrature),
    data(new DataStorage(quadrature, N)) {
    assert(solver);
    assert(dt > 0);
    assert(tol > 0);
    assert(tol < 1);
    assert(maxSweeps >= minSweeps);
    assert(minSweeps > 0);
    assert(N > 0);
    assert(quadrature);
    assert(quadrature->right_is_node());

    init();

    solver->setNumberOfImplicitStages(k - 1);
}

SDC::SDC(std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature,
    scalar tol
    )
    :
    solver(nullptr),
    nbNodes(quadrature->get_num_nodes()),
    N(0),
    k(quadrature->get_num_nodes()),
    dt(-1),
    tol(tol),
    nodes(),
    smat(),
    qmat(),
    nodesEmbedded(),
    smatEmbedded(),
    qmatEmbedded(),
    dsdc(),
    corrector(false),
    stageIndex(0),
    Sj(),
    convergence(false),
    timeIndex(0),
    minSweeps(0),
    maxSweeps(0),
    quadrature(quadrature),
    data(new DataStorage(quadrature, 0)) {
    assert(tol > 0);
    assert(tol < 1);
    assert(quadrature);
    assert(quadrature->right_is_node());

    init();
}

SDC::~SDC()
{}

void SDC::init() {
    smat = quadrature->get_s_mat();
    qmat = quadrature->get_q_mat();
    const std::vector<scalar> nodes = quadrature->get_nodes();

    this->nodes.resize(nodes.size());

    for (unsigned int i = 0; i < nodes.size(); i++)
        this->nodes(i) = nodes[i];

    dsdc.resize(this->nodes.rows() - 1);

    for (int i = 0; i < dsdc.rows(); i++)
        dsdc(i) = this->nodes(i + 1) - this->nodes(i);
}

void SDC::run() {
    assert(solver);

    scalar t = solver->getStartTime();

    while (std::abs(t - solver->getEndTime()) > 1.0e-13 && t < solver->getEndTime()) {
        scalar computedTimeStep = dt;

        solveTimeStep(t);

        t += computedTimeStep;
    }
}

void SDC::solveTimeStep(const scalar t0) {
    assert(qmat.rows() == smat.rows());
    assert(qmat.cols() == smat.cols());
    assert(qmat.rows() == k - 1);
    assert(qmat.cols() == k);
    assert(solver);
    assert(N > 0);
    assert(k > 0);

    fsi::vector dtsdc = this->dt * dsdc;
    fsi::matrix residual;
    fsi::vector errorEstimate(N);

    scalar t = t0;

    fsi::vector sol(N), f(N);
    solver->getSolution(sol, f);
    solver->evaluateFunction(0, sol, t, f);

    data->storeFunction(f, 0);
    data->storeSolution(sol, 0);

    fsi::vector rhs(N), result(N);
    rhs.setZero();

    solver->nextTimeStep();
    solver->initTimeStep();

    for (int j = 0; j < k - 1; j++) {
        scalar dt = dtsdc(j);
        t += dt;

        Info << "\nTime = " << t << ", SDC sweep = 1, SDC substep = " << j + 1 << nl << endl;

        solver->implicitSolve(false, j, j, t, dt, data->getSolution(j), rhs, f, result);

        assert((1.0 / dt * (result - data->getSolution(j) - rhs) - f).array().abs().maxCoeff() < 1.0e-8);

        data->storeFunction(f, j + 1);
        data->storeSolution(result, j + 1);
    }

    // Compute successive corrections

    for (int j = 0; j < maxSweeps - 1; j++) {
        t = t0;

        fsi::matrix Sj = this->dt * (smat * data->getFunctions());

        // SDC sweep
        for (int p = 0; p < k - 1; p++) {
            scalar dt = dtsdc(p);
            t += dt;

            Info << "\nTime = " << t << ", SDC sweep = " << j + 2 << ", SDC substep = " << p + 1 << nl << endl;

            // Form right hand side
            rhs.noalias() = -dt * data->getFunctions().row(p + 1) + Sj.row(p);

            solver->implicitSolve(true, p, p, t, dt, data->getSolution(p), rhs, f, result);

            assert((1.0 / dt * (result - data->getSolution(p) - rhs) - f).array().abs().maxCoeff() < 1.0e-8);

            data->storeFunction(f, p + 1);
            data->storeSolution(result, p + 1);
        }

        // Compute the SDC residual

        residual = dt * (qmat * data->getFunctions());

        for (int i = 0; i < residual.rows(); i++)
            residual.row(i) += data->getSolutions().row(0) - data->getSolutions().row(i + 1);

        assert(not std::isnan(residual.norm()));

        scalarList squaredNormResidual(Pstream::nProcs(), scalar(0));
        labelList dof(Pstream::nProcs(), label(0));
        squaredNormResidual[Pstream::myProcNo()] = residual.squaredNorm();
        dof[Pstream::myProcNo()] = residual.rows() * residual.cols();
        reduce(squaredNormResidual, sumOp<scalarList>());
        reduce(dof, sumOp<labelList>());
        scalar error = std::sqrt(sum(squaredNormResidual) / sum(dof));
        convergence = error < tol && j >= minSweeps - 2;

        std::deque<int> dofVariables;
        std::deque<bool> enabledVariables;
        std::deque<std::string> namesVariables;
        solver->getVariablesInfo(dofVariables, enabledVariables, namesVariables);

        assert(enabledVariables.size() == dofVariables.size());
        assert(enabledVariables.size() == namesVariables.size());

        std::deque<bool> convergenceVariables;

        for (unsigned int i = 0; i < enabledVariables.size(); i++)
            convergenceVariables.push_back(true);

        bool solverConverged = solver->isConverged();

        if (dofVariables.size() == 1) {
            Info << "SDC residual = " << error;
            Info << ", tol = " << tol;
            Info << ", time = " << t;
            Info << ", sweep = " << j + 2;
            Info << ", convergence = ";

            if (convergence)
                Info << "true";
            else
                Info << "false";

            Info << endl;

            if (convergence && solverConverged)
                break;
        }

        if (dofVariables.size() > 1) {
            assert(std::accumulate(dofVariables.begin(), dofVariables.end(), 0) == N);

            for (unsigned int substep = 0; substep < residual.rows(); substep++) {
                int index = 0;

                for (unsigned int i = 0; i < dofVariables.size(); i++) {
                    scalarList squaredNormResidual(Pstream::nProcs(), scalar(0));
                    labelList dof(Pstream::nProcs(), label(0));

                    for (int j = 0; j < dofVariables.at(i); j++) {
                        squaredNormResidual[Pstream::myProcNo()] += residual(substep, index) * residual(substep, index);
                        dof[Pstream::myProcNo()] += 1;
                        index++;
                    }

                    reduce(squaredNormResidual, sumOp<scalarList>());
                    reduce(dof, sumOp<labelList>());
                    scalar error = std::sqrt(sum(squaredNormResidual) / sum(dof));

                    bool convergence = convergenceVariables.at(i);

                    if (error > tol || j < minSweeps - 2)
                        convergence = false;

                    if (enabledVariables.at(i)) {
                        Info << "SDC " << namesVariables.at(i).c_str();
                        Info << " substep = " << substep + 1;
                        Info << ", residual = " << error;
                        Info << ", time = " << t;
                        Info << ", sweep = " << j + 2;
                        Info << ", convergence = ";

                        if (error < tol && j >= minSweeps - 2)
                            Info << "true";
                        else
                            Info << "false";

                        Info << endl;
                    }

                    if (enabledVariables.at(i))
                        convergenceVariables.at(i) = convergence;
                }

                assert(index == N);
            }

            convergence = solverConverged;

            for (unsigned int i = 0; i < convergenceVariables.size(); i++)
                if (not convergenceVariables.at(i))
                    convergence = false;

            if (convergence)
                break;
        }
    }

    solver->finalizeTimeStep();
}

void SDC::computeResidual(const fsi::matrix & qmat,
    const scalar dt,
    fsi::matrix & qj
    ) {
    // fsi::matrix Qj = dt * (qmat * F);

    // Only compute row k-2 of matrix Qj for efficiency
    int k = data->getFunctions().rows();
    int ii = k - 2, jj, kk;

    for (jj = 0; jj < data->getFunctions().cols(); ++jj) {
        qj(0, jj) = 0;

        for (kk = 0; kk < data->getFunctions().rows(); ++kk)
            qj(0, jj) += qmat(ii, kk) * data->getFunctions() (kk, jj);

        qj(0, jj) *= dt;
    }
}

void SDC::getSourceTerm(const bool corrector,
    const int k,
    const int sweep,
    const scalar deltaT,
    fsi::vector & rhs,
    fsi::vector & qold
    ) {
    assert(k <= this->k - 1);

    qold = data->getSolution(k);

    // Compute the time step from the stage deltaT
    if (dt < 0) {
        // first time step, first prediction step
        dt = deltaT / dsdc(0);
    }

    assert(dt > 0);

    if (not corrector)
        rhs.setZero();

    if (corrector) {
        if ((this->stageIndex != k || this->sweep != sweep) && k == 0) {
            data->copyFunctions();
            Sj = dt * (smat * data->getFunctions());
        }

        rhs.noalias() = -dt * dsdc(k) * data->getOldFunctions().row(k + 1) + Sj.row(k);
    }

    this->stageIndex = k;
    this->sweep = sweep;
    this->corrector = corrector;

    assert(rhs.rows() == qold.rows());
}

void SDC::setFunction(const int k,
    const fsi::vector & f,
    const fsi::vector & result
    ) {
    assert(f.rows() == result.rows());
    assert(k <= this->k - 1);
    assert(not std::isnan(f.norm()));
    assert(not std::isnan(result.norm()));

    if (data->getFunctions().cols() == 0)
        data->initialize(this->k, f.rows());

    data->storeFunction(f, k + 1);
    data->storeSolution(result, k + 1);
}

void SDC::setOldSolution(int timeIndex,
    const fsi::vector & result
    ) {
    assert(timeIndex >= this->timeIndex);
    assert(not std::isnan(result.norm()));

    if (data->getSolutions().cols() == 0) {
        data->initialize(this->k, result.rows());
        data->storeSolution(result, 0);
    } else if (timeIndex > this->timeIndex) {
        data->storeSolution(data->getLastSolution(), 0);
    }

    this->timeIndex = timeIndex;
}

void SDC::outputResidual(const std::string & name) {
    fsi::matrix residual = dt * (qmat * data->getFunctions());

    for (int i = 0; i < residual.rows(); i++)
        residual.row(i) += data->getSolutions().row(0) - data->getSolutions().row(i + 1);

    scalarList squaredNormResidual(Pstream::nProcs(), scalar(0));
    labelList dof(Pstream::nProcs(), label(0));
    squaredNormResidual[Pstream::myProcNo()] = residual.squaredNorm();
    reduce(squaredNormResidual, sumOp<scalarList>());
    dof[Pstream::myProcNo()] = residual.rows() * residual.cols();
    reduce(dof, sumOp<labelList>());
    scalar error = std::sqrt(sum(squaredNormResidual) / sum(dof));
    convergence = error < tol;

    Info << "SDC " << name.c_str();
    Info << " residual = " << error;
    Info << ", convergence = ";

    if (convergence)
        Info << "true";
    else
        Info << "false";

    Info << endl;
}

bool SDC::isConverged() {
    return convergence;
}

int SDC::getNbImplicitStages() {
    return this->k - 1;
}
} // namespace sdc
