
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> EigenMatrix;

template <int dimension>
dealiiSolidSolver<dimension>::dealiiSolidSolver(dealiifsi::DataStorage & data)
    :
    dealiifsi::LinearElasticity<dimension>(data),
    BaseMultiLevelSolver(5, dimension),
    k(0),
    kindex(0),
    UStages(),
    VStages(),
    ndofs(Pstream::nProcs(), 0),
    globalOffset(0) {
    matrix tmp;
    getWritePositions(tmp);

    BaseMultiLevelSolver::N = tmp.rows();
    BaseMultiLevelSolver::data.resize(N, dimension);
    BaseMultiLevelSolver::data.setZero();
}

template <int dimension>
dealiiSolidSolver<dimension>::dealiiSolidSolver(double time_step,
    double final_time,
    double theta,
    double degree,
    double gravity,
    double distributed_load,
    double rho,
    double E,
    double nu,
    unsigned int n_global_refines
    )
    :
    dealiifsi::LinearElasticity<dimension>(time_step, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines),
    BaseMultiLevelSolver(5, dimension, 0),
    k(0),
    kindex(0),
    UStages(),
    VStages(),
    ndofs(Pstream::nProcs(), 0),
    globalOffset(0) {
    matrix tmp;
    getWritePositions(tmp);

    BaseMultiLevelSolver::N = tmp.rows();
    BaseMultiLevelSolver::data.resize(N, dimension);
    BaseMultiLevelSolver::data.setZero();
}

template <int dimension>
dealiiSolidSolver<dimension>::~dealiiSolidSolver()
{}

template <int dimension>
void dealiiSolidSolver<dimension>::finalizeTimeStep() {
    assert(BaseMultiLevelSolver::init);

    dealiifsi::LinearElasticity<dimension>::finalizeTimeStep();

    BaseMultiLevelSolver::init = false;
}

template <int dimension>
void dealiiSolidSolver<dimension>::getReadPositions(matrix & readPositions) {
    getWritePositions(readPositions);
}

template <int dimension>
void dealiiSolidSolver<dimension>::getWritePositions(matrix & writePositions) {
    EigenMatrix local_writePositionsEigen;
    dealiifsi::LinearElasticity<dimension>::getWritePositions(local_writePositionsEigen);
    matrix local_writePositions = local_writePositionsEigen.cast<scalar>();

    ndofs = 0;
    ndofs[Pstream::myProcNo()] = local_writePositions.rows();
    reduce(ndofs, sumOp<labelList>());

    globalOffset = 0;

    for (int i = 0; i < Pstream::myProcNo(); i++)
        globalOffset += ndofs[i];

    vectorField writePositionsField(sum(ndofs), Foam::vector::zero);

    for (int i = 0; i < local_writePositions.rows(); i++)
        for (int j = 0; j < local_writePositions.cols(); j++)
            writePositionsField[i + globalOffset][j] = local_writePositions(i, j);

    reduce(writePositionsField, sumOp<vectorField>());

    writePositions.resize(writePositionsField.size(), local_writePositions.cols());

    for (int i = 0; i < writePositions.rows(); i++)
        for (int j = 0; j < writePositions.cols(); j++)
            writePositions(i, j) = writePositionsField[i][j];
}

template <int dimension>
void dealiiSolidSolver<dimension>::initTimeStep() {
    assert(!BaseMultiLevelSolver::init);

    BaseMultiLevelSolver::t += dealiifsi::LinearElasticity<dimension>::time_step;

    dealiifsi::LinearElasticity<dimension>::initTimeStep();

    BaseMultiLevelSolver::init = true;
}

template <int dimension>
bool dealiiSolidSolver<dimension>::isRunning() {
    return dealiifsi::LinearElasticity<dimension>::isRunning();
}

template <int dimension>
void dealiiSolidSolver<dimension>::resetSolution()
{}

template <int dimension>
void dealiiSolidSolver<dimension>::solve(const matrix & input,
    matrix & output
    ) {
    Info << "Solve solid domain with deal.II" << endl;

    assert(ndofs[0] > 0);

    matrix local_input(ndofs[Pstream::myProcNo()], input.cols());

    for (int i = 0; i < local_input.rows(); i++)
        for (int j = 0; j < local_input.cols(); j++)
            local_input(i, j) = input(i + globalOffset, j);

    EigenMatrix inputEigen = local_input.cast<double>();

    dealiifsi::LinearElasticity<dimension>::setTraction(inputEigen);

    dealiifsi::LinearElasticity<dimension>::solve();

    EigenMatrix local_outputEigen(ndofs[Pstream::myProcNo()], local_input.cols());
    dealiifsi::LinearElasticity<dimension>::getDisplacement(local_outputEigen);
    matrix local_output = local_outputEigen.cast<scalar>();

    assert(local_outputEigen.rows() == ndofs[Pstream::myProcNo()]);

    vectorField displacementField(sum(ndofs), Foam::vector::zero);

    for (int i = 0; i < local_output.rows(); i++)
        for (int j = 0; j < local_output.cols(); j++)
            displacementField[i + globalOffset][j] = local_output(i, j);

    reduce(displacementField, sumOp<vectorField>());

    output.resize(displacementField.size(), local_output.cols());

    for (int i = 0; i < output.rows(); i++)
        for (int j = 0; j < output.cols(); j++)
            output(i, j) = displacementField[i][j];

    BaseMultiLevelSolver::data = output;
}

template <int dimension>
void dealiiSolidSolver<dimension>::evaluateFunction(const int,
    const fsi::vector &,
    const scalar,
    fsi::vector & f
    ) {
    copy(dealiifsi::LinearElasticity<dimension>::u_f, f, 0);
    copy(dealiifsi::LinearElasticity<dimension>::v_f, f, dealiifsi::LinearElasticity<dimension>::u_f.size());
}

template <int dimension>
int dealiiSolidSolver<dimension>::getDOF() {
    return dealiifsi::LinearElasticity<dimension>::solution_u.size() + dealiifsi::LinearElasticity<dimension>::solution_v.size();
}

template <int dimension>
void dealiiSolidSolver<dimension>::getSolution(fsi::vector & solution,
    fsi::vector & f
    ) {
    copy(dealiifsi::LinearElasticity<dimension>::solution_u, solution, 0);
    copy(dealiifsi::LinearElasticity<dimension>::solution_v, solution, dealiifsi::LinearElasticity<dimension>::solution_u.size());
    copy(dealiifsi::LinearElasticity<dimension>::u_f, f, 0);
    copy(dealiifsi::LinearElasticity<dimension>::v_f, f, dealiifsi::LinearElasticity<dimension>::u_f.size());
}

template <int dimension>
void dealiiSolidSolver<dimension>::setSolution(const fsi::vector & solution,
    const fsi::vector & f
    ) {
    copy(solution, dealiifsi::LinearElasticity<dimension>::solution_u, 0);
    copy(solution, dealiifsi::LinearElasticity<dimension>::solution_v, dealiifsi::LinearElasticity<dimension>::solution_u.size());
    copy(f, dealiifsi::LinearElasticity<dimension>::u_f, 0);
    copy(f, dealiifsi::LinearElasticity<dimension>::v_f, dealiifsi::LinearElasticity<dimension>::u_f.size());
}

template <int dimension>
scalar dealiiSolidSolver<dimension>::getEndTime() {
    return dealiifsi::LinearElasticity<dimension>::final_time;
}

template <int dimension>
scalar dealiiSolidSolver<dimension>::getTimeStep() {
    return dealiifsi::LinearElasticity<dimension>::time_step;
}

template <int dimension>
void dealiiSolidSolver<dimension>::nextTimeStep() {
    for (int i = 0; i < k; i++) {
        UStages.at(i) = dealiifsi::LinearElasticity<dimension>::solution_u;
        VStages.at(i) = dealiifsi::LinearElasticity<dimension>::solution_v;
    }
}

template <int dimension>
void dealiiSolidSolver<dimension>::setNumberOfImplicitStages(int k) {
    this->k = k + 1;

    UStages.clear();
    VStages.clear();

    for (int i = 0; i < k + 1; i++) {
        UStages.push_back(dealiifsi::LinearElasticity<dimension>::solution_u);
        VStages.push_back(dealiifsi::LinearElasticity<dimension>::solution_v);
    }
}

template <int dimension>
void dealiiSolidSolver<dimension>::prepareImplicitSolve(bool corrector,
    const int k,
    const int,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs
    ) {
    dealiifsi::LinearElasticity<dimension>::time_step = dt;
    dealiifsi::LinearElasticity<dimension>::time = t;
    BaseMultiLevelSolver::t = t;

    if (corrector) {
        dealiifsi::LinearElasticity<dimension>::solution_u = UStages.at(k + 1);
        dealiifsi::LinearElasticity<dimension>::solution_v = VStages.at(k + 1);
    }

    copy(qold, dealiifsi::LinearElasticity<dimension>::old_solution_u, 0);
    copy(qold, dealiifsi::LinearElasticity<dimension>::old_solution_v, dealiifsi::LinearElasticity<dimension>::old_solution_u.size());
    copy(rhs, dealiifsi::LinearElasticity<dimension>::u_rhs, 0);
    copy(rhs, dealiifsi::LinearElasticity<dimension>::v_rhs, dealiifsi::LinearElasticity<dimension>::u_rhs.size());
}

template <int dimension>
void dealiiSolidSolver<dimension>::implicitSolve(bool corrector,
    const int k,
    const int kold,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    ) {
    kindex = k;
    prepareImplicitSolve(corrector, k, kold, t, dt, qold, rhs);

    dealiifsi::LinearElasticity<dimension>::solve();

    finalizeImplicitSolve(k);

    getSolution(result, f);
}

template <int dimension>
void dealiiSolidSolver<dimension>::finalizeImplicitSolve(int k) {
    UStages.at(k + 1) = dealiifsi::LinearElasticity<dimension>::solution_u;
    VStages.at(k + 1) = dealiifsi::LinearElasticity<dimension>::solution_v;
}

template <int dimension>
void dealiiSolidSolver<dimension>::getVariablesInfo(std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    ) {
    dof.push_back(dealiifsi::LinearElasticity<dimension>::solution_u.size());
    dof.push_back(dealiifsi::LinearElasticity<dimension>::solution_v.size());
    enabled.push_back(true);
    enabled.push_back(true);
    names.push_back("solid U");
    names.push_back("solid V");
}

void copy(const dealii::PETScWrappers::MPI::Vector & source,
    fsi::vector & target,
    unsigned int targetOffset
    ) {
    dealii::Vector<double> localized_vector(source);

    for (unsigned int i = 0; i < localized_vector.size(); ++i)
        target(i + targetOffset) = localized_vector[i];
}

void copy(const fsi::vector & source,
    dealii::PETScWrappers::MPI::Vector & target,
    unsigned int sourceOffset
    ) {
    for (unsigned int i = 0; i < target.size(); ++i)
        target[i] = source(i + sourceOffset);

    target.compress(dealii::VectorOperation::insert);
}
